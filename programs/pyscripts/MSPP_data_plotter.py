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


class MSPPDataPlotter:
    """Simple GUI for plotting protein data."""

    # Class-level constants
    ORGANISMS = ["HeLa", "E.coli", "Yeast"]
    COLORS = {"HeLa": "#9b59b6", "E.coli": "#e67e22", "Yeast": "#16a085", "Unknown": "#95a5a6"}

    # Organism pattern matching for faster identification
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

    # Dark mode colors
    DARK_BG = "#1e1e1e"
    DARK_FG = "#e0e0e0"
    DARK_ACCENT = "#3c3c3c"
    DARK_HIGHLIGHT = "#007acc"

    def __init__(self, root):
        self.root = root
        self.root.title("MSPP Data Plotter")
        self.root.geometry("400x500")

        # Apply dark mode theme
        self._setup_dark_theme()

        self.files = []
        self.cached_data = None
        self.cached_file_list = []
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

    def _get_organism_data(self, file_data, intensity_col, organism) -> pd.Series:
        """Extract positive numeric intensity data for an organism."""
        data: pd.Series = pd.to_numeric(
            file_data[file_data["Organism"] == organism][intensity_col], errors="coerce"
        )  # type: ignore[assignment]
        return data[data > 0].dropna()  # type: ignore[return-value]

    def _get_hela_median(self, file_data, intensity_col):
        """Get HeLa median for normalization, with fallback."""
        hela_data = self._get_organism_data(file_data, intensity_col, "HeLa")
        return hela_data.median() if len(hela_data) > 0 else 1.0

    def load_data(self):
        """Load data from selected files with caching."""
        if not self.files:
            messagebox.showerror("Error", "Please add at least one TSV file")
            return None

        # Return cached data if file list unchanged
        if self.cached_data is not None and self.cached_file_list == self.files:
            return self.cached_data

        all_data = []
        self.file_to_raw_column = {}

        for filepath in self.files:
            try:
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

            except Exception as e:
                messagebox.showerror("Error", f"Failed to load {Path(filepath).name}:\n{e}")
                return None

        # Cache the result
        self.cached_data = pd.concat(all_data, ignore_index=True)
        self.cached_file_list = self.files.copy()
        return self.cached_data

    def _create_bar_chart(self, data, ax):
        """Helper: Create stacked bar chart."""
        counts = data.groupby(["Source_File", "Organism"]).size().unstack(fill_value=0)
        org_order = self.ORGANISMS + ["Unknown"]
        counts = counts.reindex(
            columns=[col for col in org_order if col in counts.columns], fill_value=0
        )
        plot_colors = [self.COLORS[col] for col in counts.columns]

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

    def _extract_mix_identifier(self, filename):
        """Extract mix identifier from filename, excluding E25/E100 prefix.
        
        For 'report.pg_matrix_E25_30_4_440960_600.tsv', returns '30_4_440960_600'.
        This ensures E25 and E100 files from the same mix are grouped together.
        """
        # Remove E25/E100 prefix first to avoid capturing it in the mix ID
        cleaned = re.sub(r'E[-_]?\d+[-_]?', '', filename, flags=re.IGNORECASE)
        
        # Extract the numeric pattern (e.g., 30_4_440960_600)
        match = re.search(r'(\d+_\d+_\d+_\d+)', cleaned)
        if match:
            return match.group(1)
        
        # Fallback: return cleaned filename
        cleaned = cleaned.replace('report.pg_matrix_', '').replace('.tsv', '')
        return cleaned if cleaned else filename
    
    def _calculate_sample_intensities(self, file_data, source_file, organism):
        """Calculate HeLa-normalized intensities for all proteins in a single sample.
        
        Args:
            file_data: DataFrame for the sample
            source_file: Name of the source file
            organism: 'E.coli' or 'Yeast'
        
        Returns:
            Array of normalized log2 intensities, or None if insufficient data
        """
        if source_file not in self.file_to_raw_column:
            return None
        
        intensity_col = self.file_to_raw_column[source_file]
        
        org_data = self._get_organism_data(file_data, intensity_col, organism)
        if len(org_data) == 0:
            return None
        
        hela_median = self._get_hela_median(file_data, intensity_col)
        normalized = np.log2(org_data / hela_median)
        
        # Filter out invalid values
        intensities = np.array(normalized.values)
        return intensities[np.isfinite(intensities)]

    def _create_sample_comparison_plot(self, data):
        """Create side-by-side box plots showing E25 vs E100 fold changes for each mix."""
        # Get all sample files
        sample_files = sorted(data["Source_File"].unique())
        
        if len(sample_files) == 0:
            messagebox.showerror("Error", "No sample files found in data")
            return False
        
        # Group samples by mix identifier
        mix_groups = {}
        for source_file in sample_files:
            mix_id = self._extract_mix_identifier(source_file)
            if mix_id not in mix_groups:
                mix_groups[mix_id] = []
            mix_groups[mix_id].append(source_file)
        
        # Sort mixes
        sorted_mixes = sorted(mix_groups.keys())
        
        # Calculate intensities for E25 and E100 separately for each mix
        ecoli_results = []  # Will store (label, data, mix_id) tuples
        yeast_results = []
        mix_boundaries = []  # Track where each mix ends for visual separation
        current_position = 0
        
        for mix_id in sorted_mixes:
            mix_samples = sorted(mix_groups[mix_id])
            
            # Find E25 and E100 files in this mix (case-insensitive, flexible pattern)
            e25_file = None
            e100_file = None
            
            for source_file in mix_samples:
                upper = source_file.upper()
                # Look for E25, E-25, E_25, etc.
                if re.search(r'E[-_]?25', upper):
                    e25_file = source_file
                # Look for E100, E-100, E_100, etc.
                elif re.search(r'E[-_]?100', upper):
                    e100_file = source_file
            
            if e25_file and e100_file:
                # Calculate intensities separately for E25 and E100
                # E.coli
                e25_ecoli = self._calculate_sample_intensities(
                    data[data["Source_File"] == e25_file], e25_file, "E.coli"
                )
                e100_ecoli = self._calculate_sample_intensities(
                    data[data["Source_File"] == e100_file], e100_file, "E.coli"
                )
                
                # Yeast
                e25_yeast = self._calculate_sample_intensities(
                    data[data["Source_File"] == e25_file], e25_file, "Yeast"
                )
                e100_yeast = self._calculate_sample_intensities(
                    data[data["Source_File"] == e100_file], e100_file, "Yeast"
                )
                
                # Add E25 and E100 as separate box plots
                if e25_ecoli is not None:
                    ecoli_results.append(("E25", e25_ecoli, mix_id))
                    current_position += 1
                if e100_ecoli is not None:
                    ecoli_results.append(("E100", e100_ecoli, mix_id))
                    current_position += 1
                
                if e25_yeast is not None:
                    yeast_results.append(("E25", e25_yeast, mix_id))
                if e100_yeast is not None:
                    yeast_results.append(("E100", e100_yeast, mix_id))
            else:
                # Warn about missing pairs but continue
                if not e25_file or not e100_file:
                    print(f"Warning: Mix '{mix_id}' missing {'E25' if not e25_file else 'E100'} file. Skipping.")
            
            # Mark boundary after this mix (even if no data)
            if current_position > 0:
                mix_boundaries.append(current_position)
        
        if not ecoli_results and not yeast_results:
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
            return False
        
        # Create the plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))
        
        # E.coli plot
        if ecoli_results:
            ecoli_labels = [r[0] for r in ecoli_results]
            ecoli_data = [r[1] for r in ecoli_results]
            positions = np.arange(1, len(ecoli_data) + 1)
            
            bp1 = ax1.boxplot(
                ecoli_data,
                positions=positions,
                widths=0.6,
                patch_artist=True,
                showfliers=True,
                showmeans=True,
                flierprops=dict(
                    marker="o", markerfacecolor="#e67e22", markersize=3,
                    alpha=0.4, markeredgecolor="none"
                ),
                meanprops=dict(
                    marker="s", markerfacecolor="white", markeredgecolor="white", markersize=5
                ),
            )
            # Color E25 and E100 differently
            for i, (patch, label) in enumerate(zip(bp1["boxes"], ecoli_labels)):
                if label == "E25":
                    patch.set_facecolor("#e67e22")  # Orange for E25
                else:
                    patch.set_facecolor("#3498db")  # Blue for E100
                patch.set_alpha(0.7)
                patch.set_edgecolor("white")
                patch.set_linewidth(1.5)
            plt.setp(bp1["whiskers"], color="white", linewidth=1.5)
            plt.setp(bp1["caps"], color="white", linewidth=1.5)
            plt.setp(bp1["medians"], color="#2c3e50", linewidth=2.5)
            
            # Add median value annotations
            for i, intensities in enumerate(ecoli_data):
                median_val = np.median(intensities)
                ax1.text(
                    i + 1.35, median_val, f"{median_val:.2f}",
                    fontsize=9, va="center", color="white", fontweight="bold"
                )
            
            # Add visual separators between mixes
            ylim = ax1.get_ylim()
            prev_boundary = 0
            for idx, boundary in enumerate(mix_boundaries[:-1]):  # Don't draw line after last group
                # Draw vertical separator line
                ax1.axvline(x=boundary + 0.5, color="#555555", linestyle="-", linewidth=2, alpha=0.8)
                
                # Add mix label at the top
                mid_point = (prev_boundary + boundary) / 2 + 0.5
                if idx < len(sorted_mixes):
                    ax1.text(
                        mid_point, ylim[1] * 0.95, f"Mix: {sorted_mixes[idx]}",
                        ha="center", va="top", fontsize=9, color="#aaaaaa",
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="#2c2c2c", edgecolor="#555555", alpha=0.8)
                    )
                prev_boundary = boundary
            
            # Label for last mix
            if sorted_mixes:
                mid_point = (prev_boundary + len(ecoli_data)) / 2 + 0.5
                ax1.text(
                    mid_point, ylim[1] * 0.95, f"Mix: {sorted_mixes[-1]}",
                    ha="center", va="top", fontsize=9, color="#aaaaaa",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="#2c2c2c", edgecolor="#555555", alpha=0.8)
                )
            
            ax1.axhline(y=0, color="#f39c12", linestyle="--", linewidth=2, alpha=0.9, label="Reference (1:1)")
            ax1.set_ylabel("Log2 Intensity (HeLa-Normalized)", fontsize=12, fontweight="bold")
            ax1.set_xlabel("Sample", fontsize=12, fontweight="bold")
            ax1.set_title("E.coli Intensity Comparison", fontsize=14, fontweight="bold")
            ax1.set_xticks(positions)
            ax1.set_xticklabels(ecoli_labels, rotation=0, ha="center", fontsize=10)
            ax1.grid(axis="y", alpha=0.3)
            ax1.legend(fontsize=9)
        else:
            ax1.text(0.5, 0.5, "No E.coli data", ha="center", va="center", 
                    transform=ax1.transAxes, fontsize=14)
        
        # Yeast plot
        if yeast_results:
            yeast_labels = [r[0] for r in yeast_results]
            yeast_data = [r[1] for r in yeast_results]
            positions = np.arange(1, len(yeast_data) + 1)
            
            bp2 = ax2.boxplot(
                yeast_data,
                positions=positions,
                widths=0.6,
                patch_artist=True,
                showfliers=True,
                showmeans=True,
                flierprops=dict(
                    marker="o", markerfacecolor="#16a085", markersize=3,
                    alpha=0.4, markeredgecolor="none"
                ),
                meanprops=dict(
                    marker="s", markerfacecolor="white", markeredgecolor="white", markersize=5
                ),
            )
            # Color E25 and E100 differently
            for i, (patch, label) in enumerate(zip(bp2["boxes"], yeast_labels)):
                if label == "E25":
                    patch.set_facecolor("#16a085")  # Teal for E25
                else:
                    patch.set_facecolor("#9b59b6")  # Purple for E100
                patch.set_alpha(0.7)
                patch.set_edgecolor("white")
                patch.set_linewidth(1.5)
            plt.setp(bp2["whiskers"], color="white", linewidth=1.5)
            plt.setp(bp2["caps"], color="white", linewidth=1.5)
            plt.setp(bp2["medians"], color="#2c3e50", linewidth=2.5)
            
            # Add median value annotations
            for i, intensities in enumerate(yeast_data):
                median_val = np.median(intensities)
                ax2.text(
                    i + 1.35, median_val, f"{median_val:.2f}",
                    fontsize=9, va="center", color="white", fontweight="bold"
                )
            
            # Add visual separators between mixes
            ylim = ax2.get_ylim()
            prev_boundary = 0
            for idx, boundary in enumerate(mix_boundaries[:-1]):  # Don't draw line after last group
                # Draw vertical separator line
                ax2.axvline(x=boundary + 0.5, color="#555555", linestyle="-", linewidth=2, alpha=0.8)
                
                # Add mix label at the top
                mid_point = (prev_boundary + boundary) / 2 + 0.5
                if idx < len(sorted_mixes):
                    ax2.text(
                        mid_point, ylim[1] * 0.95, f"Mix: {sorted_mixes[idx]}",
                        ha="center", va="top", fontsize=9, color="#aaaaaa",
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="#2c2c2c", edgecolor="#555555", alpha=0.8)
                    )
                prev_boundary = boundary
            
            # Label for last mix
            if sorted_mixes:
                mid_point = (prev_boundary + len(yeast_data)) / 2 + 0.5
                ax2.text(
                    mid_point, ylim[1] * 0.95, f"Mix: {sorted_mixes[-1]}",
                    ha="center", va="top", fontsize=9, color="#aaaaaa",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="#2c2c2c", edgecolor="#555555", alpha=0.8)
                )
            
            ax2.axhline(y=0, color="#f39c12", linestyle="--", linewidth=2, alpha=0.9, label="Reference (1:1)")
            ax2.set_ylabel("Log2 Intensity (HeLa-Normalized)", fontsize=12, fontweight="bold")
            ax2.set_xlabel("Sample", fontsize=12, fontweight="bold")
            ax2.set_title("Yeast Intensity Comparison", fontsize=14, fontweight="bold")
            ax2.set_xticks(positions)
            ax2.set_xticklabels(yeast_labels, rotation=0, ha="center", fontsize=10)
            ax2.grid(axis="y", alpha=0.3)
            ax2.legend(fontsize=9)
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
        
        return True

    def plot_protein_ids(self):
        """Plot protein ID counts as stacked bar chart by organism."""
        data = self.load_data()
        if data is None:
            return

        fig, ax = plt.subplots(figsize=(12, 7))
        self._create_bar_chart(data, ax)
        plt.tight_layout()
        plt.show()

    def plot_fold_change(self):
        """Plot per-sample intensity distributions."""
        data = self.load_data()
        if data is None:
            return
        
        self._create_sample_comparison_plot(data)

    def plot_both(self):
        """Generate all plots simultaneously."""
        data = self.load_data()
        if data is None:
            return

        # Plot 1: Bar chart
        fig1, ax1 = plt.subplots(figsize=(12, 7))
        self._create_bar_chart(data, ax1)
        fig1.tight_layout()

        # Plot 2: Grouped intensity distributions
        self.plot_fold_change()

        plt.show()

def main():
    """Launch the application."""
    root = tk.Tk()
    MSPPDataPlotter(root)
    root.mainloop()

if __name__ == "__main__":
    main()
