#!/usr/bin/env python3
"""
MSPP Data Plotter - Flask Backend API
Provides REST endpoints for MS proteomics data visualization
"""

# Importing the necessary libraries
import base64
import contextlib
import io
import logging
import re
import tempfile
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from flask import Flask, jsonify, request, send_from_directory
from flask_cors import CORS

matplotlib.use('Agg')  # Non-interactive backend

app = Flask(__name__, static_folder='../frontend/dist', static_url_path='')
CORS(app)  # Enable CORS for React dev server

# Set dark mode for matplotlib
plt.style.use('dark_background')


# Shared utility functions
def fig_to_base64(fig):
    """Convert matplotlib figure to base64 encoded PNG."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=100, bbox_inches='tight')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return img_base64


class DataProcessor:
    """Handles all data loading, processing, and calculation logic."""

    # Class-level constants
    ORGANISMS = ["HeLa", "E.coli", "Yeast"]
    ORGANISM_PATTERNS = {
        "HeLa": ["_HUMAN", "HOMO_SAPIENS"],
        "E.coli": [
            "_ECOLI", "_ECOL", "_ECO2", "_ECO5", "_ECO7",
            "_SHIF", "_SHIB", "_SHIS", "ESCHERICHIA",
        ],
        "Yeast": ["_YEAST", "SACCHAROMYCES", "CEREVISIAE"],
    }

    def __init__(self):
        self.file_to_raw_column = {}
        self.cached_data = None
        self.cached_file_list = []

    def identify_organism_vectorized(self, series):
        """Vectorized organism identification - much faster than row-by-row apply."""
        upper = series.fillna("").astype(str).str.upper()
        result = pd.Series("Unknown", index=series.index)

        for organism, patterns in self.ORGANISM_PATTERNS.items():
            mask = upper.str.contains("|".join(patterns), regex=True)
            result = result.where(~mask, organism)

        return pd.Categorical(result, categories=self.ORGANISMS + ["Unknown"])

    def load_data(self, file_paths):
        """Load data from selected files with caching."""
        if not file_paths:
            raise ValueError("No files provided")

        # Return cached data if file list unchanged
        if self.cached_data is not None and self.cached_file_list == file_paths:
            return self.cached_data

        all_data = []
        self.file_to_raw_column = {}

        for filepath in file_paths:
            df = pd.read_csv(filepath, sep="\t", low_memory=False)
            source_name = Path(filepath).stem
            df["Source_File"] = source_name

            # Find and map the .raw column
            raw_cols = [col for col in df.columns if ".raw" in col.lower()]
            if raw_cols:
                self.file_to_raw_column[source_name] = raw_cols[0]

            # Identify organism from protein name column
            protein_col = next(
                (col for col in ["Protein.Names", "Protein.Group"] if col in df.columns), None
            ) or next((col for col in df.columns if "protein" in col.lower()), None)

            df["Organism"] = (
                self.identify_organism_vectorized(df[protein_col]) if protein_col else "Unknown"
            )
            all_data.append(df)

        # Cache the result
        self.cached_data = pd.concat(all_data, ignore_index=True)
        self.cached_file_list = file_paths.copy()
        return self.cached_data

    def clear_cache(self):
        """Clear cached data."""
        self.cached_data = None
        self.cached_file_list = []

    def calculate_intensity_ratios(self, data, e25_file, e100_file, organism):
        """Calculate log2 intensity ratios (E25/E100) for consensus proteins.
        For each consensus protein, calculate log2(E25_intensity / E100_intensity).
        This shows the fold change between replicates without HeLa normalization.
        Args:
            data: Full DataFrame containing both samples
            e25_file: E25 sample file name
            e100_file: E100 sample file name
            organism: 'HeLa', 'E.coli', or 'Yeast'
        Returns:
            Array of log2 ratios for consensus proteins, or None if insufficient data
        """
        if e25_file not in self.file_to_raw_column or e100_file not in self.file_to_raw_column:
            return None

        e25_intensity_col = self.file_to_raw_column[e25_file]
        e100_intensity_col = self.file_to_raw_column[e100_file]

        # Get data for each sample
        e25_data = data[data["Source_File"] == e25_file].copy()
        e100_data = data[data["Source_File"] == e100_file].copy()

        # Filter by organism
        e25_org = e25_data[e25_data["Organism"] == organism]
        e100_org = e100_data[e100_data["Organism"] == organism]

        if len(e25_org) == 0 or len(e100_org) == 0:
            return None

        # Find protein identifier column
        protein_col = next(
            (col for col in ["Protein.Group", "Protein.Ids", "Protein.Names"] if col in e25_org.columns),
            None
        )

        if protein_col is None:
            return None

        # Get valid proteins from each sample (non-zero, non-NaN intensities)
        e25_valid = e25_org[
            (e25_org[e25_intensity_col].notna()) &
            (e25_org[e25_intensity_col] > 0)
        ]
        e100_valid = e100_org[
            (e100_org[e100_intensity_col].notna()) &
            (e100_org[e100_intensity_col] > 0)
        ]

        # Find consensus proteins present in BOTH samples
        e25_proteins = set(e25_valid[protein_col])
        e100_proteins = set(e100_valid[protein_col])
        consensus_proteins = e25_proteins & e100_proteins

        if len(consensus_proteins) == 0:
            return None

        # Filter to consensus proteins only
        e25_consensus = e25_valid[e25_valid[protein_col].isin(consensus_proteins)].set_index(protein_col)
        e100_consensus = e100_valid[e100_valid[protein_col].isin(consensus_proteins)].set_index(protein_col)

        # Align by protein ID to ensure paired comparison
        common_proteins = e25_consensus.index.intersection(e100_consensus.index)
        e25_aligned = e25_consensus.loc[common_proteins, e25_intensity_col]
        e100_aligned = e100_consensus.loc[common_proteins, e100_intensity_col]

        # Calculate log2(E25 / E100) for each protein
        ratios = np.log2(e25_aligned.values / e100_aligned.values)

        # Filter out invalid values (inf, nan)
        valid_ratios = ratios[np.isfinite(ratios)]

        return valid_ratios if len(valid_ratios) > 0 else None

    def calculate_protein_id_counts(self, data):
        """Calculate protein ID counts grouped by organism and source file."""
        counts = data.groupby(["Source_File", "Organism"]).size().unstack(fill_value=0)
        org_order = self.ORGANISMS + ["Unknown"]
        counts = counts.reindex(
            columns=[col for col in org_order if col in counts.columns], fill_value=0
        )
        return counts

    def calculate_sample_comparison_data(self, data):
        """Calculate log2 intensity ratios between E25/E100 sample groups.

        Strategy:
        1. First, try to detect E25/E100 from explicit filename patterns
        2. If not found, automatically categorize using E.coli median intensities
        3. Remember: E25 ↔ Y150 (higher yeast), E100 ↔ Y75 (lower yeast)

        Returns:
            Dict with 'hela_results', 'ecoli_results', 'yeast_results', 'sample_pairs'
            Each result is a list of (ratio_array, sample_pair_label) tuples
        """
        # Get all sample files
        sample_files = sorted(data["Source_File"].unique())

        if len(sample_files) == 0:
            raise ValueError("No sample files found in data")

        if len(sample_files) < 2:
            raise ValueError("Need at least 2 samples to create E25/E100 comparisons")

        # Strategy 1: Try to detect E25/E100 from filenames (explicit naming)
        e25_samples_explicit = []
        e100_samples_explicit = []
        unclassified_samples = []

        for source_file in sample_files:
            upper = source_file.upper()
            # Check for E25 or Y150 patterns (same samples)
            if re.search(r'E[-_]?25|Y[-_]?150', upper):
                e25_samples_explicit.append(source_file)
            # Check for E100 or Y75 patterns (same samples)
            elif re.search(r'E[-_]?100|Y[-_]?75', upper):
                e100_samples_explicit.append(source_file)
            else:
                unclassified_samples.append(source_file)

        # If we found explicit naming for all samples, use it
        if len(e25_samples_explicit) > 0 and len(e100_samples_explicit) > 0 and len(unclassified_samples) == 0:
            e25_samples = sorted(e25_samples_explicit)
            e100_samples = sorted(e100_samples_explicit)
            categorization_method = "explicit naming (E25/E100/Y150/Y75)"
            logging.info(f"Using explicit sample naming: E25/Y150={[Path(s).stem for s in e25_samples]}, "
                        f"E100/Y75={[Path(s).stem for s in e100_samples]}")
        else:
            # Strategy 2: Fallback to automatic categorization using E.coli median
            logging.info("No explicit E25/E100 naming found. Using E.coli median-based categorization.")
            categorization_method = "E.coli median auto-detection"

            # Calculate median E.coli intensity for each sample
            ecoli_medians = {}
            for source_file in sample_files:
                sample_data = data[data["Source_File"] == source_file]

                # Filter for E.coli proteins only
                ecoli_data = sample_data[sample_data["Organism"] == "E.coli"]

                # Get the intensity column for this sample
                intensity_col = self.file_to_raw_column.get(source_file)
                if intensity_col and intensity_col in ecoli_data.columns:
                    # Calculate median of non-zero E.coli intensities
                    intensities = ecoli_data[intensity_col]
                    valid_intensities = intensities[(intensities > 0) & (intensities.notna())]
                    if len(valid_intensities) > 0:
                        ecoli_medians[source_file] = valid_intensities.median()
                    else:
                        ecoli_medians[source_file] = 0
                else:
                    ecoli_medians[source_file] = 0

            # Sort samples by E.coli median intensity
            sorted_samples = sorted(ecoli_medians.keys(), key=lambda x: ecoli_medians[x])

            # Split into E25 (lower half) and E100 (upper half) groups
            n_samples = len(sorted_samples)
            split_point = n_samples // 2

            e25_samples = sorted_samples[:split_point] if split_point > 0 else []
            e100_samples = sorted_samples[split_point:]

            if len(e25_samples) == 0 or len(e100_samples) == 0:
                raise ValueError("Could not split samples into E25/E100 groups. Need at least 2 samples.")

            logging.info(f"Auto-categorized samples: E25/Y150={[Path(s).stem for s in e25_samples]}, "
                        f"E100/Y75={[Path(s).stem for s in e100_samples]}")

        # Pair E25 with E100 samples
        # For uneven groups, pair as many as possible
        sample_pairs = []
        for i in range(min(len(e25_samples), len(e100_samples))):
            e25_sample = e25_samples[i]
            e100_sample = e100_samples[i]
            sample_pairs.append((e25_sample, e100_sample))

        # Calculate intensity ratios for each organism and sample pair
        hela_results = []
        ecoli_results = []
        yeast_results = []

        for e25_sample, e100_sample in sample_pairs:
            # Create label based on categorization method
            e25_name = Path(e25_sample).stem
            e100_name = Path(e100_sample).stem

            # Only create labels if using explicit naming
            if categorization_method == "explicit naming (E25/E100/Y150/Y75)":
                # Use condition names from filenames for labels
                if len(sample_pairs) == 1:
                    pair_label = f"{e25_name} vs {e100_name}"
                else:
                    pair_label = f"{e25_name} vs {e100_name}"
            else:
                # No label for auto-detection (no legend will be shown)
                pair_label = None

            # Calculate log2 ratios for each organism (E25/E100)
            # All ratios calculated using CONSENSUS PROTEINS only
            # HeLa: Expected ~0 (constant concentration, sanity check)
            # E.coli: Expected ~-2 (log2(25/100) ≈ -2.0)
            # Yeast: Expected ~1 (log2(150/75) = 1.0)
            hela_ratios = self.calculate_intensity_ratios(
                data, e25_sample, e100_sample, "HeLa"
            )
            ecoli_ratios = self.calculate_intensity_ratios(
                data, e25_sample, e100_sample, "E.coli"
            )
            yeast_ratios = self.calculate_intensity_ratios(
                data, e25_sample, e100_sample, "Yeast"
            )

            # Add results (one box per sample pair)
            if hela_ratios is not None:
                hela_results.append((hela_ratios, pair_label))
            if ecoli_ratios is not None:
                ecoli_results.append((ecoli_ratios, pair_label))
            if yeast_ratios is not None:
                yeast_results.append((yeast_ratios, pair_label))

        if not hela_results and not ecoli_results and not yeast_results:
            raise ValueError("No valid E25/E100 sample pairs with sufficient data found")

        logging.info(f"Sample comparison using {categorization_method}: {len(sample_pairs)} pairs created")

        return {
            'hela_results': hela_results,
            'ecoli_results': ecoli_results,
            'yeast_results': yeast_results,
            'sample_pairs': sample_pairs
        }



class PlotGenerator:
    """Handles all matplotlib plotting and visualization logic for web backend."""

    # Plot color scheme
    COLORS = {"HeLa": "#9b59b6", "E.coli": "#e67e22", "Yeast": "#16a085", "Unknown": "#95a5a6"}

    def __init__(self, processor):
        """Initialize with a DataProcessor instance."""
        self.processor = processor

    def create_bar_chart_figure(self, data, figsize=(12, 7)):
        """Create bar chart matplotlib figure (reusable for display and export).
        Args:
            data: DataFrame with protein data
            figsize: Tuple of (width, height) for figure size
        Returns:
            Matplotlib figure object
        """
        counts = self.processor.calculate_protein_id_counts(data)
        org_order = self.processor.ORGANISMS + ["Unknown"]
        counts = counts.reindex(
            columns=[col for col in org_order if col in counts.columns], fill_value=0
        )

        # Sort samples by numeric value extracted from .raw filename
        def get_numeric_value(sample_name):
            """Extract numeric value from a single .raw filename for sorting."""
            raw_col = self.processor.file_to_raw_column.get(sample_name, sample_name)
            raw_filename = Path(raw_col).stem if raw_col else sample_name
            # Extract first number found in filename (e.g., "57367" from "57367.raw")
            match = re.search(r'(\d+)', raw_filename)
            return int(match.group(1)) if match else 0

        # Create sorted list of sample names, then reindex
        sorted_samples = sorted(counts.index, key=get_numeric_value)
        counts = counts.reindex(sorted_samples)

        fig, ax = plt.subplots(figsize=figsize)

        # Create stacked bar chart
        bottom = np.zeros(len(counts))
        for organism in counts.columns:
            ax.bar(
                range(len(counts)), counts[organism], bottom=bottom,
                label=organism, color=self.COLORS.get(organism, "#95a5a6"), alpha=0.8
            )
            bottom += counts[organism].values

        # Add count labels on each bar segment
        for i, sample in enumerate(counts.index):
            y_offset = 0
            for organism in counts.columns:
                count = counts.loc[sample, organism]
                if count > 0:
                    bar_height = count
                    y_pos = y_offset + bar_height / 2
                    ax.text(
                        i, y_pos, str(int(count)),
                        ha='center', va='center',
                        fontsize=9, fontweight='bold',
                        color='white'
                    )
                    y_offset += bar_height

        # Get the .raw column names for x-axis labels
        x_labels = []
        for sample in counts.index:
            # Get the corresponding .raw column name, or use the sample name if not found
            raw_col = self.processor.file_to_raw_column.get(sample, sample)
            # Extract just the filename from the full path (e.g., "57367.raw" from "D:\QC_DIA\57367.raw")
            raw_filename = Path(raw_col).name if raw_col else sample
            x_labels.append(raw_filename)

        # Configure axes and styling
        ax.set_xlabel("Sample", fontsize=12, fontweight="bold")
        ax.set_ylabel("Number of Protein IDs", fontsize=12, fontweight="bold")
        ax.set_title("Protein ID Counts by Organism", fontsize=14, fontweight="bold")
        ax.set_xticks(range(len(counts)))
        ax.set_xticklabels(x_labels, rotation=45, ha='right')
        ax.legend(title="Organism", fontsize=10, loc="upper right")
        ax.grid(axis="y", alpha=0.3)
        plt.tight_layout()

        return fig

    def create_bar_chart(self, data):
        """Create and return bar chart as base64 image."""
        fig = self.create_bar_chart_figure(data)
        return fig_to_base64(fig)

    def create_comparison_figure(self, data, figsize=(18, 16)):
        """Create sample comparison matplotlib figure (reusable for display and export).
        Args:
            data: DataFrame with protein data
            figsize: Tuple of (width, height) for figure size
        Returns:
            Matplotlib figure object
        """
        # Get prepared data from processor
        comparison_data = self.processor.calculate_sample_comparison_data(data)

        hela_results = comparison_data['hela_results']
        ecoli_results = comparison_data['ecoli_results']
        yeast_results = comparison_data['yeast_results']

        # Create figure with 3 vertical subplots
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=figsize)

        # Plot HeLa comparison (top)
        if hela_results:
            self.plot_ratio_comparison(
                ax1, hela_results,
                title="HeLa Log2 Intensity Ratio (E25/E100 = constant)",
                color="#9b59b6",  # Purple
                reference_line=0  # Expected: 0 (constant concentration - sanity check)
            )
        else:
            ax1.text(0.5, 0.5, "No HeLa data", ha="center", va="center",
                    transform=ax1.transAxes, fontsize=14)

        # Plot E.coli comparison (middle)
        if ecoli_results:
            self.plot_ratio_comparison(
                ax2, ecoli_results,
                title="E.coli Log2 Intensity Ratio (E25/E100)",
                color="#e67e22",  # Orange
                reference_line=-2  # Expected: -2 (log2(25/100) ≈ -2.0)
            )
        else:
            ax2.text(0.5, 0.5, "No E.coli data", ha="center", va="center",
                    transform=ax2.transAxes, fontsize=14)

        # Plot Yeast comparison (bottom)
        if yeast_results:
            self.plot_ratio_comparison(
                ax3, yeast_results,
                title="Yeast Log2 Intensity Ratio (Y150/Y75)",
                color="#16a085",  # Teal
                reference_line=1  # Expected: 1 (log2(150/75) = 1.0)
            )
        else:
            ax3.text(0.5, 0.5, "No Yeast data", ha="center", va="center",
                    transform=ax3.transAxes, fontsize=14)

        plt.suptitle(
            "Intensity Ratio Comparison by Run",
            fontsize=16,
            fontweight="bold",
            y=0.995
        )
        plt.tight_layout()

        return fig

    def create_sample_comparison_plot(self, data):
        """Create and return stacked 3-panel comparison plot as base64 image."""
        fig = self.create_comparison_figure(data)
        return fig_to_base64(fig)

    def plot_ratio_comparison(self, ax, results, title, color, reference_line):
        """Helper method to plot box plot for intensity ratios of one organism.
        Args:
            ax: Matplotlib axis
            results: List of (ratio_array, sample_pair_label) tuples
            title: Plot title
            color: Box plot color
            reference_line: Y-value for expected ratio reference line
        """
        data_arrays = [r[0] for r in results]
        mix_labels = [r[1] if r[1] is not None else "" for r in results]
        positions = np.arange(1, len(data_arrays) + 1)

        # Create box plot
        bp = ax.boxplot(
            data_arrays,
            positions=positions,
            widths=0.6,
            patch_artist=True,
            showfliers=True,
            showmeans=True,
            flierprops={
                'marker': "o", 'markerfacecolor': color, 'markersize': 3,
                'alpha': 0.4, 'markeredgecolor': "none"
            },
            meanprops={
                'marker': "s", 'markerfacecolor': "white", 'markeredgecolor': "white", 'markersize': 5
            },
        )

        # Color all boxes the same
        for patch in bp["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
            patch.set_edgecolor("white")
            patch.set_linewidth(1.5)

        plt.setp(bp["whiskers"], color="white", linewidth=1.5)
        plt.setp(bp["caps"], color="white", linewidth=1.5)
        plt.setp(bp["medians"], color="#2c3e50", linewidth=2.5)

        # Add median value annotations
        for i, data_arr in enumerate(data_arrays):
            median_val = np.median(data_arr)
            ax.text(
                i + 1, median_val, f"{median_val:.2f}",
                fontsize=9, va="bottom", ha="center", color="white",
                fontweight="bold", bbox={'boxstyle': "round,pad=0.3",
                'facecolor': "black", 'alpha': 0.5, 'edgecolor': "none"}
            )

        # Add reference line for expected ratio
        ax.axhline(y=reference_line, color="#f39c12", linestyle="--",
                   linewidth=2, alpha=0.9,
                   label=f"Expected: {reference_line}")

        # Configure axes
        ax.set_ylabel("Log2 Intensity Ratio", fontsize=11, fontweight="bold")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xticks(positions)
        ax.set_xticklabels(mix_labels, rotation=45, ha="right", fontsize=9)
        ax.grid(axis="y", alpha=0.3)

        # Always show expected value legend
        ax.legend(fontsize=9, loc="upper right")


'''
The main Python script for the heavy lifting ends here. The rest is how the Flask app is set up and routes defined for the webapp.
If you are interested in tweaking or extending the functionality of the webapp interfaces, check out the frontend code in programs/mspp_web/frontend.
'''

# Global instances
processor = DataProcessor()
plotter = PlotGenerator(processor)
uploaded_files = {}  # Store uploaded files temporarily


# Flask routes
@app.route('/')
def serve_react_app():
    """Serve React frontend."""
    return send_from_directory(app.static_folder, 'index.html') # type: ignore


@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({'status': 'ok', 'message': 'MSPP Backend API running'})


@app.route('/api/upload', methods=['POST'])
def upload_files():
    """Handle TSV file uploads."""
    if 'files' not in request.files:
        return jsonify({'error': 'No files provided'}), 400

    files = request.files.getlist('files')
    temp_paths = []

    for file in files:
        if file.filename.endswith(('.tsv', '.txt')):
            # Save to temp directory
            temp_path = Path(tempfile.gettempdir()) / file.filename
            file.save(temp_path)
            temp_paths.append(str(temp_path))
            uploaded_files[file.filename] = str(temp_path)

    return jsonify({
        'message': f'{len(temp_paths)} files uploaded successfully',
        'files': [Path(p).name for p in temp_paths]
    })


@app.route('/api/files', methods=['GET'])
def list_files():
    """List uploaded files."""
    return jsonify({'files': list(uploaded_files.keys())})


@app.route('/api/files', methods=['DELETE'])
def clear_files():
    """Clear all uploaded files."""
    for filepath in uploaded_files.values():
        with contextlib.suppress(Exception):
            Path(filepath).unlink(missing_ok=True)
    uploaded_files.clear()
    return jsonify({'message': 'All files cleared'})


@app.route('/api/plot/bar-chart', methods=['POST'])
def generate_bar_chart():
    """Generate protein ID bar chart."""
    if not uploaded_files:
        return jsonify({'error': 'No files uploaded'}), 400

    try:
        data = processor.load_data(list(uploaded_files.values()))
        img_base64 = plotter.create_bar_chart(data)
        return jsonify({'image': img_base64})
    except Exception:
        logging.exception("Unexpected error while generating bar chart")
        return jsonify({
            'error': 'An internal error occurred while generating the bar chart.'
        }), 500


@app.route('/api/plot/sample-comparison', methods=['POST'])
def generate_sample_comparison():
    """Generate E25 vs E100 intensity comparison plot."""
    if not uploaded_files:
        return jsonify({'error': 'No files uploaded'}), 400

    try:
        data = processor.load_data(list(uploaded_files.values()))
        img_base64 = plotter.create_sample_comparison_plot(data)
        return jsonify({'image': img_base64})
    except ValueError as e:
        logging.warning(
            "Validation error while generating sample comparison plot: %s",
            e,
        )
        return jsonify({
            'error': 'Invalid input or data for sample comparison plot. Please verify your files and parameters.'
        }), 400
    except Exception:
        logging.exception("Unexpected error while generating sample comparison plot")
        return jsonify({
            'error': 'An internal error occurred while generating the sample comparison plot.'
        }), 500


@app.route('/api/export/bar-chart', methods=['POST'])
def export_bar_chart():
    """Export protein ID bar chart as PNG file."""
    if not uploaded_files:
        return jsonify({'error': 'No files uploaded'}), 400

    try:
        from flask import send_file
        data = processor.load_data(list(uploaded_files.values()))

        # Generate high-DPI plot using reusable method
        fig = plotter.create_bar_chart_figure(data, figsize=(10, 6))

        # Save to bytes buffer
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
        buf.seek(0)
        plt.close(fig)

        return send_file(buf, mimetype='image/png', as_attachment=True, download_name='protein_id_bar_chart.png')
    except Exception:
        logging.exception("Unexpected error while exporting bar chart")
        return jsonify({
            'error': 'An internal error occurred while exporting the bar chart.'
        }), 500


@app.route('/api/export/sample-comparison', methods=['POST'])
def export_sample_comparison():
    """Export sample comparison plot as PNG file."""
    if not uploaded_files:
        return jsonify({'error': 'No files uploaded'}), 400

    try:
        from flask import send_file
        data = processor.load_data(list(uploaded_files.values()))

        # Generate high-DPI plot using reusable method
        fig = plotter.create_comparison_figure(data, figsize=(18, 16))

        # Save to bytes buffer with high DPI
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
        buf.seek(0)
        plt.close(fig)

        return send_file(buf, mimetype='image/png', as_attachment=True, download_name='intensity_ratio_comparison.png')
    except Exception:
        logging.exception("Unexpected error while exporting sample comparison")
        return jsonify({
            'error': 'An internal error occurred while exporting the sample comparison plot.'
        }), 500


@app.route('/api/export/all', methods=['POST'])
def export_all_plots():
    """Export all plots as a ZIP file."""
    if not uploaded_files:
        return jsonify({'error': 'No files uploaded'}), 400

    try:
        import zipfile

        from flask import send_file

        data = processor.load_data(list(uploaded_files.values()))

        # Create ZIP file in memory
        zip_buffer = io.BytesIO()

        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            # Generate and add bar chart using reusable method
            fig = plotter.create_bar_chart_figure(data, figsize=(10, 6))
            bar_buf = io.BytesIO()
            fig.savefig(bar_buf, format='png', dpi=300, bbox_inches='tight')
            bar_buf.seek(0)
            plt.close(fig)
            zip_file.writestr('protein_id_bar_chart.png', bar_buf.getvalue())

            # Generate and add sample comparison plot using reusable method
            fig = plotter.create_comparison_figure(data, figsize=(18, 16))
            comp_buf = io.BytesIO()
            fig.savefig(comp_buf, format='png', dpi=300, bbox_inches='tight')
            comp_buf.seek(0)
            plt.close(fig)
            zip_file.writestr('intensity_ratio_comparison.png', comp_buf.getvalue())

        zip_buffer.seek(0)
        return send_file(
            zip_buffer,
            mimetype='application/zip',
            as_attachment=True,
            download_name='mspp_plots.zip'
        )
    except Exception:
        logging.exception("Unexpected error while exporting all plots")
        return jsonify({
            'error': 'An internal error occurred while exporting plots.'
        }), 500
