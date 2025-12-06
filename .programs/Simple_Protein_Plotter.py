#!/usr/bin/env python3
"""
MS Protein & Peptide Data Plotter
Quick visualization of protein IDs and relative abundances from TSV files
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
import re
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np

# Set dark mode for matplotlib
plt.style.use('dark_background')


class MSPPDataPlotter:
    """Simple GUI for plotting protein data."""
    
    # Class-level constants
    ORGANISMS = ['HeLa', 'E.coli', 'Yeast']
    COLORS = {'HeLa': '#9b59b6', 'E.coli': '#e67e22', 'Yeast': '#16a085', 'Unknown': '#95a5a6'}
    
    # Organism pattern matching for faster identification
    ORGANISM_PATTERNS = {
        'HeLa': ['_HUMAN', 'HOMO_SAPIENS'],
        'E.coli': ['_ECOLI', '_ECOL', '_ECO2', '_ECO5', '_ECO7', '_SHIF', '_SHIB', '_SHIS', 'ESCHERICHIA'],
        'Yeast': ['_YEAST', 'SACCHAROMYCES', 'CEREVISIAE']
    }

    # Dark mode colors
    DARK_BG = '#1e1e1e'
    DARK_FG = '#e0e0e0'
    DARK_ACCENT = '#3c3c3c'
    DARK_HIGHLIGHT = '#007acc'
    
    def __init__(self, root):
        self.root = root
        self.root.title("MSPP Data Plotter")
        self.root.geometry("400x700")
        
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
        style.theme_use('clam')
        
        # Configure dark colors for all ttk widgets
        style.configure('.', background=self.DARK_BG, foreground=self.DARK_FG)
        style.configure('TFrame', background=self.DARK_BG)
        style.configure('TLabel', background=self.DARK_BG, foreground=self.DARK_FG)
        style.configure('TButton', background=self.DARK_ACCENT, foreground=self.DARK_FG,
                        borderwidth=1, focuscolor=self.DARK_HIGHLIGHT)
        style.map('TButton', background=[('active', self.DARK_HIGHLIGHT)])
        style.configure('TLabelframe', background=self.DARK_BG, foreground=self.DARK_FG)
        style.configure('TLabelframe.Label', background=self.DARK_BG, foreground=self.DARK_FG)

    def setup_ui(self):
        """Setup the user interface."""
        main_frame = ttk.Frame(self.root, padding="20")
        main_frame.grid(row=0, column=0, sticky='nsew')
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        
        # Title
        ttk.Label(main_frame, text="MS Protein & Peptide Data Plotter", 
                 font=('Arial', 16, 'bold')).grid(row=0, column=0, columnspan=2, pady=(0, 20))
        
        # File selection
        ttk.Label(main_frame, text="Select TSV files:").grid(row=1, column=0, sticky=tk.W, pady=5)
        
        self.file_listbox = tk.Listbox(main_frame, height=6, width=50,
                                        bg=self.DARK_ACCENT, fg=self.DARK_FG,
                                        selectbackground=self.DARK_HIGHLIGHT,
                                        selectforeground='white', borderwidth=0,
                                        highlightthickness=1, highlightcolor=self.DARK_HIGHLIGHT)
        self.file_listbox.grid(row=2, column=0, columnspan=2, pady=5, sticky='ew')
        
        # Buttons
        btn_frame = ttk.Frame(main_frame)
        btn_frame.grid(row=3, column=0, columnspan=2, pady=10)
        
        ttk.Button(btn_frame, text="Add Files", command=self.add_files).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Clear", command=self.clear_files).pack(side=tk.LEFT, padx=5)
        
        # Plot buttons
        plot_frame = ttk.LabelFrame(main_frame, text="Generate Plots", padding="10")
        plot_frame.grid(row=4, column=0, columnspan=2, pady=20, sticky='ew')
        
        ttk.Button(plot_frame, text="📊 Protein ID Bar Chart", 
                  command=self.plot_protein_ids, width=35).pack(pady=5, fill=tk.X, padx=10)
        
        ttk.Button(plot_frame, text="📊 E.coli vs Yeast Fold Change", 
                  command=self.plot_fold_change, width=35).pack(pady=5, fill=tk.X, padx=10)
        
        ttk.Button(plot_frame, text="📈 Both Plots", 
                  command=self.plot_both, width=35).pack(pady=5, fill=tk.X, padx=10)
        
        # Grouping options
        group_frame = ttk.LabelFrame(main_frame, text="Group Files by Pattern", padding="10")
        group_frame.grid(row=5, column=0, columnspan=2, pady=10, sticky='ew')
        
        ttk.Label(group_frame, text="Group pattern (regex):").pack(anchor=tk.W)
        self.group_pattern = tk.StringVar(value="")
        pattern_entry = tk.Entry(group_frame, textvariable=self.group_pattern, width=40,
                                 bg=self.DARK_ACCENT, fg=self.DARK_FG,
                                 insertbackground=self.DARK_FG)
        pattern_entry.pack(fill=tk.X, pady=5)
        
        ttk.Label(group_frame, text="Examples: '(\\d+ng_\\d+ng)' or '(rep\\d+)' or '(ratio_[^_]+)'",
                 font=('Arial', 8)).pack(anchor=tk.W)
        
        ttk.Button(group_frame, text="📦 Grouped Fold Change Box Plot", 
                  command=self.plot_grouped_fold_change, width=35).pack(pady=5, fill=tk.X)

    def add_files(self):
        """Add TSV files."""
        filenames = filedialog.askopenfilenames(
            title="Select TSV Files",
            filetypes=[("TSV files", "*.tsv"), ("Text files", "*.txt"), ("All files", "*.*")]
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
        upper = series.fillna('').astype(str).str.upper()
        result = pd.Series('Unknown', index=series.index)
        
        for organism, patterns in self.ORGANISM_PATTERNS.items():
            mask = upper.str.contains('|'.join(patterns), regex=True)
            result = result.where(~mask, organism)
        
        return pd.Categorical(result, categories=self.ORGANISMS + ['Unknown'])

    def _get_organism_data(self, file_data, intensity_col, organism) -> pd.Series:
        """Extract positive numeric intensity data for an organism."""
        data: pd.Series = pd.to_numeric(
            file_data[file_data['Organism'] == organism][intensity_col],
            errors='coerce'
        )  # type: ignore[assignment]
        return data[data > 0].dropna()  # type: ignore[return-value]

    def _calculate_protein_fold_changes(self, file_data, intensity_col):
        """Calculate per-protein log2 fold changes (E.coli/Yeast) for a single file."""
        hela_data = self._get_organism_data(file_data, intensity_col, 'HeLa')
        ecoli_data = self._get_organism_data(file_data, intensity_col, 'E.coli')
        yeast_data = self._get_organism_data(file_data, intensity_col, 'Yeast')
        
        if len(ecoli_data) == 0 or len(yeast_data) == 0:
            return None
        
        hela_median = hela_data.median() if len(hela_data) > 0 else 1.0
        ecoli_norm = ecoli_data / hela_median
        yeast_median = (yeast_data / hela_median).median()
        
        protein_fcs = np.log2(ecoli_norm / yeast_median)
        return protein_fcs.replace([np.inf, -np.inf], np.nan).dropna().values

    def _style_boxplot(self, bp):
        """Apply consistent styling to boxplot elements."""
        for patch in bp['boxes']:
            patch.set_facecolor('#5dade2')
            patch.set_alpha(0.7)
            patch.set_edgecolor('white')
            patch.set_linewidth(1)
        for element in ['whiskers', 'caps']:
            plt.setp(bp[element], color='white', linewidth=1)
        plt.setp(bp['medians'], color='#2c3e50', linewidth=2)

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
                df = pd.read_csv(filepath, sep='\t', low_memory=False)
                source_name = Path(filepath).stem
                df['Source_File'] = source_name
                
                # Find and map the .raw column
                raw_cols = [col for col in df.columns if '.raw' in col.lower()]
                if raw_cols:
                    self.file_to_raw_column[source_name] = raw_cols[0]
                
                # Identify organism from protein name column
                protein_col = next((col for col in ['Protein.Names', 'Protein.Group'] 
                                   if col in df.columns), None) or \
                              next((col for col in df.columns if 'protein' in col.lower()), None)
                
                df['Organism'] = self.identify_organism_vectorized(df[protein_col]) if protein_col else 'Unknown'
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
        counts = data.groupby(['Source_File', 'Organism']).size().unstack(fill_value=0)
        org_order = self.ORGANISMS + ['Unknown']
        counts = counts.reindex(columns=[col for col in org_order if col in counts.columns], fill_value=0)
        plot_colors = [self.COLORS[col] for col in counts.columns]
        
        counts.plot(kind='bar', stacked=True, ax=ax, color=plot_colors, 
                   edgecolor='black', linewidth=0.5, alpha=0.8)
        ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Protein IDs', fontsize=12, fontweight='bold')
        ax.set_title('Protein ID Counts by Organism', fontsize=14, fontweight='bold')
        ax.legend(title='Organism', fontsize=10)
        ax.grid(axis='y', alpha=0.3)
        ax.tick_params(axis='x', rotation=45)
    
    def _calculate_fold_changes_per_sample(self, data):
        """Calculate per-protein log2 abundance ratios (E.coli/Yeast) for box plot display."""
        results = []
        
        for source_file in data['Source_File'].unique():
            if source_file not in self.file_to_raw_column:
                continue
            
            intensity_col = self.file_to_raw_column[source_file]
            file_data = data[data['Source_File'] == source_file]
            
            if intensity_col not in file_data.columns:
                continue
            
            fcs = self._calculate_protein_fold_changes(file_data, intensity_col)
            if fcs is not None and len(fcs) > 0:
                results.append((source_file, fcs, len(fcs)))
        
        if not results:
            return [], [], []
        
        names, fcs_list, counts = zip(*results)
        return list(names), list(fcs_list), list(counts)
    
    def _create_fold_change_plot(self, data, ax):
        """Helper: Create box plot showing per-protein fold change distributions."""
        sample_names, all_protein_fcs, protein_counts = self._calculate_fold_changes_per_sample(data)
        
        if not all_protein_fcs:
            messagebox.showerror("Error", "No samples with both E.coli and Yeast data found")
            return False
        
        # Calculate median FC for each sample and sort
        sample_medians = [(name, np.median(fcs), fcs, count) 
                          for name, fcs, count in zip(sample_names, all_protein_fcs, protein_counts)]
        sample_medians.sort(key=lambda x: x[1], reverse=True)
        
        sorted_names = [x[0] for x in sample_medians]
        sorted_fcs = [x[2] for x in sample_medians]
        sorted_counts = [x[3] for x in sample_medians]
        sorted_medians = [x[1] for x in sample_medians]
        
        # Create box plots
        positions = np.arange(1, len(sorted_fcs) + 1)
        bp = ax.boxplot(sorted_fcs, positions=positions, widths=0.6,
                        patch_artist=True, showfliers=True, showmeans=True,
                        flierprops=dict(marker='o', markerfacecolor='#3498db', 
                                       markersize=3, alpha=0.4, markeredgecolor='none'),
                        meanprops=dict(marker='s', markerfacecolor='white', 
                                      markeredgecolor='white', markersize=4))
        self._style_boxplot(bp)
        
        # Add expected ratio reference lines
        ax.axhline(y=0, color='#f39c12', linestyle='--', linewidth=2, alpha=0.9)
        
        # Add median annotations
        for i, (med, count) in enumerate(zip(sorted_medians, sorted_counts)):
            ax.text(i + 1.35, med, f'{med:.2f}', fontsize=8, va='center', color='#f39c12')
            ax.text(i + 1, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else max([max(fc) for fc in sorted_fcs]) + 0.5, 
                   f'{count} proteins', fontsize=7, ha='center', va='bottom', color='white', rotation=90)
        
        # Labels and formatting
        ax.set_ylabel('Log2 Abundance Ratio (E.coli / Yeast median)', fontsize=12, fontweight='bold')
        ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
        ax.set_title('Per-Protein Abundance Ratios (HeLa-Normalized)', fontsize=14, fontweight='bold')
        ax.set_xticks(positions)
        ax.set_xticklabels([s.replace('report.pg_matrix_', '') for s in sorted_names], 
                          rotation=90, fontsize=7)
        ax.grid(axis='y', alpha=0.3)
        
        # Legend
        legend_elements = [
            Line2D([0], [0], color='#f39c12', linestyle='--', linewidth=2, label='Expected 1:1 (log2=0)'),
            Line2D([0], [0], color='#2c3e50', linewidth=2, label='Median'),
            Line2D([0], [0], marker='s', color='w', markerfacecolor='white', markersize=6, label='Mean', linestyle='None')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
        
        # Overall statistics
        overall_mean = np.mean(sorted_medians)
        overall_std = np.std(sorted_medians)
        
        stats_text = f'n = {len(sorted_fcs)} samples\nMedian of medians: {np.median(sorted_medians):.3f}\nMean of medians: {overall_mean:.3f}\nSD: {overall_std:.3f}'
        ax.text(0.02, 0.02, stats_text, transform=ax.transAxes, fontsize=9, 
               verticalalignment='bottom', color='white',
               bbox=dict(boxstyle='round', facecolor='#2c2c2c', edgecolor='#555555', alpha=0.8))
        
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
        """Plot E.coli vs Yeast fold change across all samples."""
        data = self.load_data()
        if data is None:
            return
        
        fig, ax = plt.subplots(figsize=(12, 7))
        if self._create_fold_change_plot(data, ax):
            plt.tight_layout()
            plt.show()

    def plot_both(self):
        """Generate both plots simultaneously."""
        data = self.load_data()
        if data is None:
            return
        
        fig1, ax1 = plt.subplots(figsize=(12, 7))
        self._create_bar_chart(data, ax1)
        fig1.tight_layout()
        
        fig2, ax2 = plt.subplots(figsize=(12, 7))
        if self._create_fold_change_plot(data, ax2):
            fig2.tight_layout()
            plt.show()
        else:
            plt.close(fig2)
            plt.show()
    
    def plot_grouped_fold_change(self):
        """Plot fold changes grouped by filename pattern."""
        pattern = self.group_pattern.get().strip()
        if not pattern:
            messagebox.showerror("Error", "Please enter a grouping pattern (regex).\n\nExamples:\n• (\\d+ng_\\d+ng) - matches '24ng_192ng'\n• (ratio_[^_]+) - matches 'ratio_1to2'\n• ([A-Z]+\\d+) - matches 'REP01'")
            return
        
        data = self.load_data()
        if data is None:
            return
        
        try:
            regex = re.compile(pattern, re.IGNORECASE)
        except re.error as e:
            messagebox.showerror("Error", f"Invalid regex pattern: {e}")
            return
        
        # Group files by pattern and calculate fold changes
        groups = {}
        unmatched = []
        
        for source_file in data['Source_File'].unique():
            match = regex.search(source_file)
            if not match or source_file not in self.file_to_raw_column:
                if not match:
                    unmatched.append(source_file)
                continue
            
            group_name = match.group(1) if match.lastindex else match.group(0)
            file_data = data[data['Source_File'] == source_file]
            intensity_col = self.file_to_raw_column[source_file]
            
            fcs = self._calculate_protein_fold_changes(file_data, intensity_col)
            if fcs is not None and len(fcs) > 0:
                if group_name not in groups:
                    groups[group_name] = []
                groups[group_name].extend(fcs)
        
        if not groups:
            messagebox.showerror("Error", f"No files matched the pattern '{pattern}'.\n\nFiles: {', '.join(data['Source_File'].unique()[:5])}...")
            return
        
        if unmatched:
            messagebox.showwarning("Warning", f"{len(unmatched)} files did not match the pattern and will be excluded.")
        
        # Sort groups by median fold change
        sorted_groups = sorted(groups.keys(), key=lambda g: np.median(groups[g]), reverse=True)
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(10, 7))
        
        positions = np.arange(1, len(sorted_groups) + 1)
        box_data = [np.array(groups[g]) for g in sorted_groups]
        
        bp = ax.boxplot(box_data, positions=positions, widths=0.6,
                        patch_artist=True, showfliers=True, showmeans=True,
                        flierprops=dict(marker='o', markerfacecolor='#3498db', 
                                       markersize=3, alpha=0.3, markeredgecolor='none'),
                        meanprops=dict(marker='s', markerfacecolor='white', 
                                      markeredgecolor='white', markersize=5))
        self._style_boxplot(bp)
        
        # Add expected ratio reference line
        ax.axhline(y=0, color='#f39c12', linestyle='--', linewidth=2, alpha=0.9)
        
        # Add median annotations and protein counts
        for i, group in enumerate(sorted_groups):
            med = np.median(groups[group])
            ax.text(i + 1.3, med, f'{med:.2f}', fontsize=10, va='center', color='#f39c12', fontweight='bold')
            ax.text(i + 1, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] > 0 else max([max(d) for d in box_data]) + 0.3,
                   f'{len(groups[group])} proteins', fontsize=9, ha='center', va='bottom', color='white')
        
        # Labels
        ax.set_ylabel('Log2 Abundance Ratio (E.coli / Yeast)', fontsize=12, fontweight='bold')
        ax.set_xlabel('Group', fontsize=12, fontweight='bold')
        ax.set_title('Per-Protein Fold Change by Group (HeLa-Normalized)', fontsize=14, fontweight='bold')
        ax.set_xticks(positions)
        ax.set_xticklabels(sorted_groups, rotation=45, ha='right', fontsize=10)
        ax.grid(axis='y', alpha=0.3)
        
        # Legend
        legend_elements = [
            Line2D([0], [0], color='#f39c12', linestyle='--', linewidth=2, label='Expected 1:1 (log2=0)'),
            Line2D([0], [0], color='#2c3e50', linewidth=2, label='Median'),
            Line2D([0], [0], marker='s', color='w', markerfacecolor='white', markersize=6, label='Mean', linestyle='None')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
        
        # Stats
        all_medians = [np.median(groups[g]) for g in sorted_groups]
        stats_text = f'{len(sorted_groups)} groups\nMedian range: {min(all_medians):.2f} to {max(all_medians):.2f}'
        ax.text(0.02, 0.02, stats_text, transform=ax.transAxes, fontsize=9, 
               verticalalignment='bottom', color='white',
               bbox=dict(boxstyle='round', facecolor='#2c2c2c', edgecolor='#555555', alpha=0.8))
        
        plt.tight_layout()
        plt.show()


def main():
    """Launch the application."""
    root = tk.Tk()
    MSPPDataPlotter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
