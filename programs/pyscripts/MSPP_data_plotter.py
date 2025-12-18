#!/usr/bin/env python3
"""
MS Protein & Peptide Data Plotter
Quick visualization of protein IDs and relative intensitys from TSV files
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Set dark mode for matplotlib
plt.style.use("dark_background")

class DataProcessor:
    """Handles all data loading, processing, and calculation logic."""
    
    # Class-level constants
    ORGANISMS = ["HeLa", "E.coli", "Yeast"]
    ORGANISM_PATTERNS = {
        "HeLa": ["_HUMAN", "HOMO_SAPIENS"],
        "E.coli": [
            "_ECOLI",
            "_ECOL",
            "_ECO2",
            "_ECO5",
            "_ECO7",
            "_SHIF",
            "_SHIB",
            "_SHIS",
            "ESCHERICHIA",
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
    
    def _get_organism_data(self, file_data, intensity_col, organism):
        """Extract positive numeric intensity data for an organism."""
        data = pd.to_numeric(
            file_data[file_data["Organism"] == organism][intensity_col], errors="coerce"
        )
        return data[data > 0].dropna()

    def _get_hela_median(self, file_data, intensity_col):
        """Get HeLa median for normalization, with fallback."""
        hela_data = self._get_organism_data(file_data, intensity_col, "HeLa")
        return hela_data.median() if len(hela_data) > 0 else 1.0
    
    def _extract_mix_identifier(self, filename):
        """Extract mix identifier from filename, excluding E25/E100 prefix.
        
        For 'report.pg_matrix_E25_30_4_440960_800.tsv', returns '30_4_440960_800'.
        This ensures E25 and E100 files from the same mix are grouped together.
        """
        # Remove E25/E100 prefix first to avoid capturing it in the mix ID
        cleaned = re.sub(r'E[-_]?\d+[-_]?', '', filename, flags=re.IGNORECASE)
        
        # Extract the numeric pattern (e.g., 30_4_440960_800)
        match = re.search(r'(\d+_\d+_\d+_\d+)', cleaned)
        if match:
            return match.group(1)
        
        # Fallback: return cleaned filename
        cleaned = cleaned.replace('report.pg_matrix_', '').replace('.tsv', '')
        return cleaned if cleaned else filename
    
    def _calculate_consensus_fold_changes(self, data, e25_file, e100_file, organism):
        """Calculate log2 fold changes for consensus proteins present in BOTH samples.
        
        Args:
            data: Full DataFrame containing both samples
            e25_file: E25 sample file name
            e100_file: E100 sample file name
            organism: 'E.coli' or 'Yeast'
        
        Returns:
            Tuple of (e25_log2_values, e100_log2_values) for consensus proteins only
        """
        if e25_file not in self.file_to_raw_column or e100_file not in self.file_to_raw_column:
            return None, None
        
        e25_intensity_col = self.file_to_raw_column[e25_file]
        e100_intensity_col = self.file_to_raw_column[e100_file]
        
        # Get data for each sample
        e25_data = data[data["Source_File"] == e25_file].copy()
        e100_data = data[data["Source_File"] == e100_file].copy()
        
        # Filter by organism
        e25_org = e25_data[e25_data["Organism"] == organism]
        e100_org = e100_data[e100_data["Organism"] == organism]
        
        if len(e25_org) == 0 or len(e100_org) == 0:
            return None, None
        
        # Find consensus proteins: must have valid intensity in BOTH samples
        # Assuming Protein.Group or similar column exists for matching
        protein_col = next(
            (col for col in ["Protein.Group", "Protein.Ids", "Protein.Names"] if col in e25_org.columns),
            None
        )
        
        if protein_col is None:
            # Fallback to old method if no protein identifier column
            return self._calculate_sample_intensities_legacy(e25_data, e25_file, organism), \
                   self._calculate_sample_intensities_legacy(e100_data, e100_file, organism)
        
        # Get valid proteins from each sample (non-zero, non-NaN intensities)
        e25_valid = e25_org[
            (e25_org[e25_intensity_col].notna()) & 
            (e25_org[e25_intensity_col] > 0)
        ]
        e100_valid = e100_org[
            (e100_org[e100_intensity_col].notna()) & 
            (e100_org[e100_intensity_col] > 0)
        ]
        
        # Find consensus proteins present in BOTH
        e25_proteins = set(e25_valid[protein_col])
        e100_proteins = set(e100_valid[protein_col])
        consensus_proteins = e25_proteins & e100_proteins
        
        if len(consensus_proteins) == 0:
            return None, None
        
        # Filter to consensus proteins only
        e25_consensus = e25_valid[e25_valid[protein_col].isin(consensus_proteins)]
        e100_consensus = e100_valid[e100_valid[protein_col].isin(consensus_proteins)]
        
        # Calculate HeLa median for each sample
        e25_hela_median = self._get_hela_median(e25_data, e25_intensity_col)
        e100_hela_median = self._get_hela_median(e100_data, e100_intensity_col)
        
        # Calculate log2 normalized intensities for consensus proteins
        e25_intensities = e25_consensus[e25_intensity_col]
        e100_intensities = e100_consensus[e100_intensity_col]
        
        e25_normalized = np.log2(e25_intensities / e25_hela_median)
        e100_normalized = np.log2(e100_intensities / e100_hela_median)
        
        # Filter out invalid values
        e25_array = np.array(e25_normalized.values)
        e100_array = np.array(e100_normalized.values)
        
        e25_valid_vals = e25_array[np.isfinite(e25_array)]
        e100_valid_vals = e100_array[np.isfinite(e100_array)]
        
        return e25_valid_vals, e100_valid_vals
    
    def _calculate_sample_intensities_legacy(self, file_data, source_file, organism):
        """Legacy method: Calculate HeLa-normalized intensities for all proteins in a single sample.
        Used as fallback when protein identifier column is not found.
        """
        if source_file not in self.file_to_raw_column:
            return None
        
        intensity_col = self.file_to_raw_column[source_file]
        
        org_data = self._get_organism_data(file_data, intensity_col, organism)
        if len(org_data) == 0:
            return None
        
        hela_median = self._get_hela_median(file_data, intensity_col)
        normalized = np.log2(org_data / hela_median)
        
        intensities = np.array(normalized.values)
        return intensities[np.isfinite(intensities)]
    
    def calculate_protein_id_counts(self, data):
        """Calculate protein ID counts grouped by organism and source file."""
        counts = data.groupby(["Source_File", "Organism"]).size().unstack(fill_value=0)
        org_order = self.ORGANISMS + ["Unknown"]
        counts = counts.reindex(
            columns=[col for col in org_order if col in counts.columns], fill_value=0
        )
        return counts
    
    def calculate_sample_comparison_data(self, data):
        """Calculate E25 vs E100 intensity comparison data for all mixes.
        
        Returns:
            Dict with 'ecoli_results', 'yeast_results', 'mix_boundaries', 'sorted_mixes'
        """
        # Get all sample files
        sample_files = sorted(data["Source_File"].unique())
        
        if len(sample_files) == 0:
            raise ValueError("No sample files found in data")
        
        # Group samples by mix identifier
        mix_groups = {}
        for source_file in sample_files:
            mix_id = self._extract_mix_identifier(source_file)
            if mix_id not in mix_groups:
                mix_groups[mix_id] = []
            mix_groups[mix_id].append(source_file)
        
        sorted_mixes = sorted(mix_groups.keys())
        
        # Calculate intensities for E25 and E100 separately for each mix
        ecoli_results = []
        yeast_results = []
        mix_boundaries = []
        current_position = 0
        
        for mix_id in sorted_mixes:
            mix_samples = sorted(mix_groups[mix_id])
            
            # Find E25 and E100 files
            e25_file = None
            e100_file = None
            
            for source_file in mix_samples:
                upper = source_file.upper()
                if re.search(r'E[-_]?25', upper):
                    e25_file = source_file
                elif re.search(r'E[-_]?100', upper):
                    e100_file = source_file
            
            if e25_file and e100_file:
                # Calculate consensus fold changes for E.coli
                e25_ecoli, e100_ecoli = self._calculate_consensus_fold_changes(
                    data, e25_file, e100_file, "E.coli"
                )
                
                # Calculate consensus fold changes for Yeast
                e25_yeast, e100_yeast = self._calculate_consensus_fold_changes(
                    data, e25_file, e100_file, "Yeast"
                )
                
                # Add E.coli results (both E25 and E100 together since they're matched)
                if e25_ecoli is not None and e100_ecoli is not None:
                    ecoli_results.append(("E25", e25_ecoli, mix_id))
                    current_position += 1
                    ecoli_results.append(("E100", e100_ecoli, mix_id))
                    current_position += 1
                
                # Add Yeast results (both E25 and E100 together since they're matched)
                if e25_yeast is not None and e100_yeast is not None:
                    yeast_results.append(("E25", e25_yeast, mix_id))
                    yeast_results.append(("E100", e100_yeast, mix_id))
            else:
                if not e25_file or not e100_file:
                    print(f"Warning: Mix '{mix_id}' missing {'E25' if not e25_file else 'E100'} file. Skipping.")
            
            if current_position > 0:
                mix_boundaries.append(current_position)
        
        if not ecoli_results and not yeast_results:
            # Generate error details
            error_details = {
                'sample_files': sample_files,
                'sorted_mixes': sorted_mixes,
                'mix_groups': mix_groups
            }
            raise ValueError("No valid E25/E100 pairs found", error_details)
        
        return {
            'ecoli_results': ecoli_results,
            'yeast_results': yeast_results,
            'mix_boundaries': mix_boundaries,
            'sorted_mixes': sorted_mixes
        }


class PlotGenerator:
    """Handles all matplotlib plotting and visualization logic."""
    
    # Plot color scheme
    COLORS = {"HeLa": "#9b59b6", "E.coli": "#e67e22", "Yeast": "#16a085", "Unknown": "#95a5a6"}
    
    def __init__(self, processor):
        """Initialize with a DataProcessor instance."""
        self.processor = processor
    
    def create_bar_chart(self, data):
        """Create and display stacked bar chart of protein ID counts."""
        counts = self.processor.calculate_protein_id_counts(data)
        plot_colors = [self.COLORS[col] for col in counts.columns]

        fig, ax = plt.subplots(figsize=(12, 7))
        counts.plot(
            kind="bar",
            stacked=True,
            ax=ax,
            color=plot_colors,
            edgecolor="black",
            linewidth=0.5,
            alpha=0.8,
        )
        ax.set_xlabel("Sample", fontsize=12, fontweight="bold")
        ax.set_ylabel("Number of Protein IDs", fontsize=12, fontweight="bold")
        ax.set_title("Protein ID Counts by Organism", fontsize=14, fontweight="bold")
        ax.legend(title="Organism", fontsize=10, loc="upper right")
        ax.grid(axis="y", alpha=0.3)
        ax.tick_params(axis="x", rotation=45, labelbottom=True)
        plt.setp(ax.xaxis.get_majorticklabels(), ha="right")
        plt.tight_layout()
        plt.show()
    
    def create_sample_comparison_plot(self, data):
        """Create and display side-by-side box plots for E25 vs E100 comparison."""
        # Get prepared data from processor
        comparison_data = self.processor.calculate_sample_comparison_data(data)
        
        ecoli_results = comparison_data['ecoli_results']
        yeast_results = comparison_data['yeast_results']
        mix_boundaries = comparison_data['mix_boundaries']
        sorted_mixes = comparison_data['sorted_mixes']
        
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))
        
        # Plot E.coli comparison
        if ecoli_results:
            self._plot_organism_comparison(
                ax1, ecoli_results, mix_boundaries, sorted_mixes,
                title="E.coli Intensity Comparison",
                colors=("#e67e22", "#3498db"),  # Orange, Blue
                label_map={"E25": "E25", "E100": "E100"}
            )
        else:
            ax1.text(0.5, 0.5, "No E.coli data", ha="center", va="center",
                    transform=ax1.transAxes, fontsize=14)
        
        # Plot Yeast comparison
        if yeast_results:
            self._plot_organism_comparison(
                ax2, yeast_results, mix_boundaries, sorted_mixes,
                title="Yeast Intensity Comparison",
                colors=("#16a085", "#9b59b6"),  # Teal, Purple
                label_map={"E25": "Y150", "E100": "Y75"}
            )
        else:
            ax2.text(0.5, 0.5, "No Yeast data", ha="center", va="center",
                    transform=ax2.transAxes, fontsize=14)
        
        plt.suptitle(
            "E25 vs E100 Intensity Comparison by Mix (HeLa-Normalized)",
            fontsize=15,
            fontweight="bold",
            y=0.98
        )
        plt.tight_layout()
        plt.show()
    
    def _plot_organism_comparison(self, ax, results, mix_boundaries, sorted_mixes, 
                                   title, colors, label_map):
        """Helper method to plot box plot comparison for one organism.
        
        Args:
            ax: Matplotlib axis
            results: List of (label, data, mix_id) tuples
            mix_boundaries: List of boundary positions between mixes
            sorted_mixes: List of mix identifiers
            title: Plot title
            colors: Tuple of (first_color, second_color) for the two conditions
            label_map: Dict mapping internal labels (E25/E100) to display labels
        """
        # Extract and convert labels
        original_labels = [r[0] for r in results]
        display_labels = [label_map.get(lbl, lbl) for lbl in original_labels]
        data_arrays = [r[1] for r in results]
        positions = np.arange(1, len(data_arrays) + 1)
        
        # Create box plot
        bp = ax.boxplot(
            data_arrays,
            positions=positions,
            widths=0.6,
            patch_artist=True,
            showfliers=True,
            showmeans=True,
            flierprops=dict(
                marker="o", markerfacecolor=colors[0], markersize=3,
                alpha=0.4, markeredgecolor="none"
            ),
            meanprops=dict(
                marker="s", markerfacecolor="white", markeredgecolor="white", markersize=5
            ),
        )
        
        # Color boxes based on original label (E25 vs E100)
        for i, (patch, orig_label) in enumerate(zip(bp["boxes"], original_labels)):
            color = colors[0] if orig_label == "E25" else colors[1]
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
            patch.set_edgecolor("white")
            patch.set_linewidth(1.5)
        
        plt.setp(bp["whiskers"], color="white", linewidth=1.5)
        plt.setp(bp["caps"], color="white", linewidth=1.5)
        plt.setp(bp["medians"], color="#2c3e50", linewidth=2.5)
        
        # Add median value annotations
        for i, data in enumerate(data_arrays):
            median_val = np.median(data)
            ax.text(
                i + 1.35, median_val, f"{median_val:.2f}",
                fontsize=9, va="center", color="white", fontweight="bold"
            )
        
        # Add visual separators and mix labels
        self._add_mix_separators(ax, mix_boundaries, sorted_mixes, len(data_arrays))
        
        # Configure axes
        ax.axhline(y=0, color="#f39c12", linestyle="--", linewidth=2, alpha=0.9, 
                   label="Reference (1:1)")
        ax.set_ylabel("Log2 Intensity (HeLa-Normalized)", fontsize=12, fontweight="bold")
        ax.set_xlabel("Sample", fontsize=12, fontweight="bold")
        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.set_xticks(positions)
        ax.set_xticklabels(display_labels, rotation=0, ha="center", fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(fontsize=9)
    
    def _add_mix_separators(self, ax, mix_boundaries, sorted_mixes, total_samples):
        """Add vertical separator lines and mix labels to plot."""
        ylim = ax.get_ylim()
        prev_boundary = 0
        
        # Add separators between mixes (but not after the last one)
        for idx, boundary in enumerate(mix_boundaries[:-1]):
            ax.axvline(x=boundary + 0.5, color="#555555", linestyle="-", 
                      linewidth=2, alpha=0.8)
            
            # Add mix label
            mid_point = (prev_boundary + boundary) / 2 + 0.5
            if idx < len(sorted_mixes):
                ax.text(
                    mid_point, ylim[1] * 0.95, f"Mix: {sorted_mixes[idx]}",
                    ha="center", va="top", fontsize=9, color="#aaaaaa",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="#2c2c2c", 
                             edgecolor="#555555", alpha=0.8)
                )
            prev_boundary = boundary
        
        # Label for last mix
        if sorted_mixes:
            mid_point = (prev_boundary + total_samples) / 2 + 0.5
            ax.text(
                mid_point, ylim[1] * 0.95, f"Mix: {sorted_mixes[-1]}",
                ha="center", va="top", fontsize=9, color="#aaaaaa",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="#2c2c2c", 
                         edgecolor="#555555", alpha=0.8)
            )


class MSPPDataPlotter:
    """GUI application for MS proteomics data visualization."""

    # UI color scheme
    DARK_BG = "#1e1e1e"
    DARK_FG = "#e0e0e0"
    DARK_ACCENT = "#3c3c3c"
    DARK_HIGHLIGHT = "#007acc"

    def __init__(self, root):
        self.root = root
        self.root.title("MSPP Data Plotter")
        self.root.geometry("400x500")

        # Initialize components
        self.processor = DataProcessor()
        self.plotter = PlotGenerator(self.processor)
        
        # UI state
        self.files = []
        
        # Apply dark mode theme and setup UI
        self._setup_dark_theme()
        self.setup_ui()

    def _setup_dark_theme(self):
        """Configure dark mode styling for ttk widgets."""
        self.root.configure(bg=self.DARK_BG)

        style = ttk.Style()
        style.theme_use("clam")

        # Configure dark colors for all ttk widgets
        style.configure(".", background=self.DARK_BG, foreground=self.DARK_FG)
        style.configure("TFrame", background=self.DARK_BG)
        style.configure("TLabel", background=self.DARK_BG, foreground=self.DARK_FG)
        style.configure(
            "TButton",
            background=self.DARK_ACCENT,
            foreground=self.DARK_FG,
            borderwidth=1,
            focuscolor=self.DARK_HIGHLIGHT,
        )
        style.map("TButton", background=[("active", self.DARK_HIGHLIGHT)])
        style.configure("TLabelframe", background=self.DARK_BG, foreground=self.DARK_FG)
        style.configure("TLabelframe.Label", background=self.DARK_BG, foreground=self.DARK_FG)

    def setup_ui(self):
        """Setup the user interface."""
        main_frame = ttk.Frame(self.root, padding="20")
        main_frame.grid(row=0, column=0, sticky="nsew")
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        # Title
        ttk.Label(
            main_frame, text="MS Protein & Peptide Data Plotter", font=("Arial", 16, "bold")
        ).grid(row=0, column=0, columnspan=2, pady=(0, 20))

        # File selection
        ttk.Label(main_frame, text="Select TSV files:").grid(row=1, column=0, sticky=tk.W, pady=5)

        self.file_listbox = tk.Listbox(
            main_frame,
            height=6,
            width=50,
            bg=self.DARK_ACCENT,
            fg=self.DARK_FG,
            selectbackground=self.DARK_HIGHLIGHT,
            selectforeground="white",
            borderwidth=0,
            highlightthickness=1,
            highlightcolor=self.DARK_HIGHLIGHT,
        )
        self.file_listbox.grid(row=2, column=0, columnspan=2, pady=5, sticky="ew")

        # Buttons
        btn_frame = ttk.Frame(main_frame)
        btn_frame.grid(row=3, column=0, columnspan=2, pady=10)

        ttk.Button(btn_frame, text="Add Files", command=self.add_files).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Clear", command=self.clear_files).pack(side=tk.LEFT, padx=5)

        # Plot buttons
        plot_frame = ttk.LabelFrame(main_frame, text="Generate Plots", padding="10")
        plot_frame.grid(row=4, column=0, columnspan=2, pady=20, sticky="ew")

        ttk.Button(
            plot_frame, text="📊 Protein ID Bar Chart", command=self.plot_protein_ids, width=35
        ).pack(pady=5, fill=tk.X, padx=10)

        ttk.Button(
            plot_frame,
            text="📊 Sample Intensity Comparison",
            command=self.plot_fold_change,
            width=35,
        ).pack(pady=5, fill=tk.X, padx=10)

        ttk.Button(plot_frame, text="📈 All Plots", command=self.plot_both, width=35).pack(
            pady=5, fill=tk.X, padx=10
        )

    def add_files(self):
        """Add TSV files."""
        filenames = filedialog.askopenfilenames(
            title="Select TSV Files",
            filetypes=[("TSV files", "*.tsv"), ("Text files", "*.txt"), ("All files", "*.*")],
        )
        for filename in filenames:
            if filename not in self.files:
                self.files.append(filename)
                self.file_listbox.insert(tk.END, Path(filename).name)

    def clear_files(self):
        """Clear file list."""
        self.files.clear()
        self.file_listbox.delete(0, tk.END)
        self.processor.clear_cache()

    def load_data(self):
        """Load data from selected files."""
        if not self.files:
            messagebox.showerror("Error", "Please add at least one TSV file")
            return None
        
        try:
            return self.processor.load_data(self.files)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load data:\n{e}")
            return None

    def plot_protein_ids(self):
        """Generate and display protein ID bar chart."""
        data = self.load_data()
        if data is not None:
            self.plotter.create_bar_chart(data)

    def plot_fold_change(self):
        """Generate and display sample comparison plot."""
        data = self.load_data()
        if data is None:
            return
        
        # Handle errors from plotting
        try:
            self.plotter.create_sample_comparison_plot(data)
        except ValueError as e:
            # Handle detailed error diagnostics
            if len(e.args) > 1 and isinstance(e.args[1], dict):
                error_details = e.args[1]
                sample_files = error_details['sample_files']
                sorted_mixes = error_details['sorted_mixes']
                mix_groups = error_details['mix_groups']
                
                error_msg = "No valid E25/E100 pairs found.\n\n"
                error_msg += f"Detected {len(sample_files)} files:\n"
                for sf in sample_files:
                    error_msg += f"  • {sf}\n"
                error_msg += f"\n{len(sorted_mixes)} mix group(s):\n"
                for mix_id, samples in mix_groups.items():
                    error_msg += f"\n{mix_id}:\n"
                    for s in samples:
                        has_e25 = 'YES' if re.search(r'E[-_]?25', s.upper()) else 'NO'
                        has_e100 = 'YES' if re.search(r'E[-_]?100', s.upper()) else 'NO'
                        error_msg += f"  - {s}\n"
                        error_msg += f"    E25: {has_e25}, E100: {has_e100}\n"
                messagebox.showerror("Error", error_msg)
            else:
                messagebox.showerror("Error", str(e))

    def plot_both(self):
        """Generate all plots simultaneously."""
        data = self.load_data()
        if data is not None:
            self.plotter.create_bar_chart(data)
            self.plot_fold_change()  # Handles its own errors

def main():
    """Launch the application."""
    root = tk.Tk()
    MSPPDataPlotter(root)
    root.mainloop()

if __name__ == "__main__":
    main()
