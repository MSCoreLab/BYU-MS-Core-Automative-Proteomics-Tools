#!/usr/bin/env python3
"""
Astral HEY Analysis
MS Proteomics Data Analysis Script for QC of Orbitrap Astral + Vanquish Neo

HEY Mix Composition:
- HeLa: Constant concentration across all runs
- E. coli & Yeast ratios:
  * Mix 1: E100 / Y75
  * Mix 2: E25 / Y150

This script compares protein IDs and ratio fold changes between the two HEY mixtures
to verify instrument consistency and data quality.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
from pathlib import Path
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


class HEYAnalyzer:
    """Analyzer for HEY (HeLa, E. coli, Yeast) mixed sample QC data."""

    def __init__(self, progress_callback=None):
        self.mix1_data = None  # E100/Y75
        self.mix2_data = None  # E25/Y150
        self.comparison_results = None
        self.progress_callback = progress_callback

    def _update_progress(self, message):
        """Update progress via callback if available."""
        if self.progress_callback:
            self.progress_callback(message)
        
    def load_tsv_data(self, filepath, mix_type):
        """
        Load TSV data from Astral output.
        
        Parameters:
        -----------
        filepath : str or Path
            Path to the TSV file
        mix_type : str
            Either 'E100Y75' or 'E25Y150'
        
        Returns:
        --------
        pd.DataFrame
            Loaded data
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        data = pd.read_csv(filepath, sep='\t', low_memory=False)
        msg = f"Loaded {len(data)} entries from {filepath.name}"
        print(msg)
        self._update_progress(msg)

        if mix_type == 'E100Y75':
            self.mix1_data = data
        elif mix_type == 'E25Y150':
            self.mix2_data = data
        else:
            raise ValueError("mix_type must be 'E100Y75' or 'E25Y150'")

        return data
    
    def identify_organism(self, protein_id):
        """
        Identify organism from protein ID.
        
        Parameters:
        -----------
        protein_id : str
            Protein identifier
        
        Returns:
        --------
        str
            'HeLa', 'E.coli', 'Yeast', or 'Unknown'
        """
        if pd.isna(protein_id):
            return 'Unknown'
        
        protein_id = str(protein_id).upper()
        
        # Common patterns for each organism
        if any(x in protein_id for x in ['HUMAN', 'HOMO_SAPIENS', 'HS_']):
            return 'HeLa'
        elif any(x in protein_id for x in ['ECOLI', 'E.COLI', 'EC_', 'ESCHERICHIA']):
            return 'E.coli'
        elif any(x in protein_id for x in ['YEAST', 'SACCHAROMYCES', 'SC_', 'CEREVISIAE']):
            return 'Yeast'
        else:
            return 'Unknown'
    
    def annotate_organism(self, data, protein_id_column='Protein.IDs'):
        """
        Add organism annotation to dataframe.
        
        Parameters:
        -----------
        data : pd.DataFrame
            Input data
        protein_id_column : str
            Name of column containing protein IDs
        
        Returns:
        --------
        pd.DataFrame
            Data with added 'Organism' column
        """
        if protein_id_column not in data.columns:
            # Try to find a suitable column
            possible_columns = [col for col in data.columns if 'protein' in col.lower()]
            if possible_columns:
                protein_id_column = possible_columns[0]
                msg = f"Using column: {protein_id_column}"
                print(msg)
                self._update_progress(msg)
            else:
                raise ValueError(f"Could not find protein ID column. Available columns: {list(data.columns)}")

        data['Organism'] = data[protein_id_column].apply(self.identify_organism)
        return data
    
    def calculate_organism_stats(self, data):
        """
        Calculate statistics for each organism.
        
        Parameters:
        -----------
        data : pd.DataFrame
            Annotated data with 'Organism' column
        
        Returns:
        --------
        pd.DataFrame
            Summary statistics by organism
        """
        stats = data.groupby('Organism').agg({
            'Organism': 'count'
        }).rename(columns={'Organism': 'Protein_Count'})
        
        return stats
    
    def compare_protein_ids(self, protein_id_column='Protein.IDs'):
        """
        Compare protein IDs between the two mixes.
        
        Parameters:
        -----------
        protein_id_column : str
            Name of column containing protein IDs
        
        Returns:
        --------
        dict
            Comparison results with shared and unique proteins
        """
        if self.mix1_data is None or self.mix2_data is None:
            raise ValueError("Both mix datasets must be loaded first")
        
        # Annotate organisms
        mix1 = self.annotate_organism(self.mix1_data.copy(), protein_id_column)
        mix2 = self.annotate_organism(self.mix2_data.copy(), protein_id_column)
        
        # Get protein ID sets
        mix1_ids = set(mix1[protein_id_column].dropna())
        mix2_ids = set(mix2[protein_id_column].dropna())
        
        # Compare
        shared_ids = mix1_ids & mix2_ids
        mix1_only = mix1_ids - mix2_ids
        mix2_only = mix2_ids - mix1_ids
        
        # Organism-specific stats
        results = {
            'total_shared': len(shared_ids),
            'total_mix1_only': len(mix1_only),
            'total_mix2_only': len(mix2_only),
            'shared_proteins': shared_ids,
            'mix1_only_proteins': mix1_only,
            'mix2_only_proteins': mix2_only,
            'mix1_stats': self.calculate_organism_stats(mix1),
            'mix2_stats': self.calculate_organism_stats(mix2)
        }
        
        # Organism-specific comparison
        for organism in ['HeLa', 'E.coli', 'Yeast']:
            mix1_org = set(mix1[mix1['Organism'] == organism][protein_id_column].dropna())
            mix2_org = set(mix2[mix2['Organism'] == organism][protein_id_column].dropna())
            
            results[f'{organism}_shared'] = len(mix1_org & mix2_org)
            results[f'{organism}_mix1_only'] = len(mix1_org - mix2_org)
            results[f'{organism}_mix2_only'] = len(mix2_org - mix1_org)
        
        self.comparison_results = results
        return results
    
    def calculate_fold_changes(self, intensity_column='Intensity', protein_id_column='Protein.IDs'):
        """
        Calculate fold changes for proteins detected in both mixes.
        
        Parameters:
        -----------
        intensity_column : str
            Name of column containing intensity values
        protein_id_column : str
            Name of column containing protein IDs
        
        Returns:
        --------
        pd.DataFrame
            Fold change data with organism annotations
        """
        if self.mix1_data is None or self.mix2_data is None:
            raise ValueError("Both mix datasets must be loaded first")
        
        # Annotate organisms
        mix1 = self.annotate_organism(self.mix1_data.copy(), protein_id_column)
        mix2 = self.annotate_organism(self.mix2_data.copy(), protein_id_column)
        
        # Find suitable intensity column if not present
        if intensity_column not in mix1.columns:
            possible_cols = [col for col in mix1.columns if 'intensity' in col.lower() or 'abundance' in col.lower()]
            if possible_cols:
                intensity_column = possible_cols[0]
                msg = f"Using intensity column: {intensity_column}"
                print(msg)
                self._update_progress(msg)

        # Merge on protein IDs
        merged = mix1[[protein_id_column, intensity_column, 'Organism']].merge(
            mix2[[protein_id_column, intensity_column, 'Organism']],
            on=[protein_id_column, 'Organism'],
            suffixes=('_E100Y75', '_E25Y150')
        )
        
        # Calculate fold change (log2)
        merged['Log2FC'] = np.log2(
            (merged[f'{intensity_column}_E25Y150'] + 1) / 
            (merged[f'{intensity_column}_E100Y75'] + 1)
        )
        
        # Calculate absolute fold change
        merged['FoldChange'] = merged[f'{intensity_column}_E25Y150'] / merged[f'{intensity_column}_E100Y75']
        
        return merged
    
    def plot_protein_overlap(self, save_path=None):
        """
        Visualize protein ID overlap between mixes.
        
        Parameters:
        -----------
        save_path : str or Path, optional
            Path to save the plot
        """
        if self.comparison_results is None:
            raise ValueError("Run compare_protein_ids() first")
        
        results = self.comparison_results
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Overall comparison
        ax = axes[0]
        categories = ['Shared', 'E100Y75\nOnly', 'E25Y150\nOnly']
        values = [results['total_shared'], results['total_mix1_only'], results['total_mix2_only']]
        colors = ['#2ecc71', '#3498db', '#e74c3c']
        
        ax.bar(categories, values, color=colors, edgecolor='black', linewidth=1.5)
        ax.set_ylabel('Number of Proteins', fontsize=12)
        ax.set_title('Overall Protein ID Comparison', fontsize=14, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        for i, v in enumerate(values):
            ax.text(i, v + max(values)*0.02, str(v), ha='center', va='bottom', fontweight='bold')
        
        # Organism-specific comparison
        ax = axes[1]
        organisms = ['HeLa', 'E.coli', 'Yeast']
        shared = [results.get(f'{org}_shared', 0) for org in organisms]
        mix1_only = [results.get(f'{org}_mix1_only', 0) for org in organisms]
        mix2_only = [results.get(f'{org}_mix2_only', 0) for org in organisms]
        
        x = np.arange(len(organisms))
        width = 0.25
        
        ax.bar(x - width, shared, width, label='Shared', color='#2ecc71', edgecolor='black')
        ax.bar(x, mix1_only, width, label='E100Y75 Only', color='#3498db', edgecolor='black')
        ax.bar(x + width, mix2_only, width, label='E25Y150 Only', color='#e74c3c', edgecolor='black')
        
        ax.set_ylabel('Number of Proteins', fontsize=12)
        ax.set_xlabel('Organism', fontsize=12)
        ax.set_title('Organism-Specific Protein Overlap', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(organisms)
        ax.legend()
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")
        
        plt.show()
    
    def plot_fold_change_distribution(self, fold_change_df, save_path=None):
        """
        Plot fold change distributions by organism.
        
        Parameters:
        -----------
        fold_change_df : pd.DataFrame
            Fold change data from calculate_fold_changes()
        save_path : str or Path, optional
            Path to save the plot
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        organisms = ['HeLa', 'E.coli', 'Yeast']
        colors = {'HeLa': '#9b59b6', 'E.coli': '#e67e22', 'Yeast': '#16a085'}
        
        # Histogram of log2 fold changes by organism
        ax = axes[0, 0]
        for organism in organisms:
            org_data = fold_change_df[fold_change_df['Organism'] == organism]
            if len(org_data) > 0:
                ax.hist(org_data['Log2FC'].dropna(), bins=30, alpha=0.6, 
                       label=organism, color=colors.get(organism), edgecolor='black')
        
        ax.set_xlabel('Log2 Fold Change (E25Y150 / E100Y75)', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title('Log2 Fold Change Distribution by Organism', fontsize=12, fontweight='bold')
        ax.legend()
        ax.axvline(x=0, color='red', linestyle='--', linewidth=2, alpha=0.7)
        ax.grid(alpha=0.3)
        
        # Box plot by organism
        ax = axes[0, 1]
        data_to_plot = [fold_change_df[fold_change_df['Organism'] == org]['Log2FC'].dropna() 
                        for org in organisms]
        bp = ax.boxplot(data_to_plot, labels=organisms, patch_artist=True, showfliers=False)
        
        for patch, organism in zip(bp['boxes'], organisms):
            patch.set_facecolor(colors.get(organism))
            patch.set_alpha(0.6)
        
        ax.set_ylabel('Log2 Fold Change (E25Y150 / E100Y75)', fontsize=11)
        ax.set_title('Log2 Fold Change by Organism', fontsize=12, fontweight='bold')
        ax.axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.7)
        ax.grid(axis='y', alpha=0.3)
        
        # Scatter plot: Mix1 vs Mix2 intensities
        ax = axes[1, 0]
        for organism in organisms:
            org_data = fold_change_df[fold_change_df['Organism'] == organism]
            if len(org_data) > 0:
                int_col = [col for col in org_data.columns if 'ntensity' in col or 'bundance' in col][0]
                mix1_col = f'{int_col}_E100Y75'
                mix2_col = f'{int_col}_E25Y150'
                
                ax.scatter(np.log10(org_data[mix1_col] + 1), 
                          np.log10(org_data[mix2_col] + 1),
                          alpha=0.5, s=20, label=organism, color=colors.get(organism))
        
        # Add diagonal line
        max_val = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=2)
        
        ax.set_xlabel('Log10(Intensity) E100Y75', fontsize=11)
        ax.set_ylabel('Log10(Intensity) E25Y150', fontsize=11)
        ax.set_title('Intensity Comparison', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3)
        
        # Expected vs Observed fold changes for E.coli and Yeast
        ax = axes[1, 1]
        
        # Expected log2 fold changes
        # E.coli: 25/100 = 0.25, log2(0.25) = -2
        # Yeast: 150/75 = 2, log2(2) = 1
        expected = {'E.coli': -2, 'Yeast': 1, 'HeLa': 0}
        
        observed_medians = []
        organism_labels = []
        expected_values = []
        
        for organism in organisms:
            org_data = fold_change_df[fold_change_df['Organism'] == organism]['Log2FC'].dropna()
            if len(org_data) > 0:
                median_fc = org_data.median()
                observed_medians.append(median_fc)
                organism_labels.append(organism)
                expected_values.append(expected[organism])
        
        x = np.arange(len(organism_labels))
        width = 0.35
        
        ax.bar(x - width/2, expected_values, width, label='Expected', 
               color='gray', alpha=0.7, edgecolor='black')
        ax.bar(x + width/2, observed_medians, width, label='Observed (Median)', 
               alpha=0.7, edgecolor='black',
               color=[colors.get(org) for org in organism_labels])
        
        ax.set_ylabel('Log2 Fold Change', fontsize=11)
        ax.set_title('Expected vs Observed Fold Changes', fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(organism_labels)
        ax.legend()
        ax.axhline(y=0, color='red', linestyle='--', linewidth=1, alpha=0.5)
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")
        
        plt.show()
    
    def generate_qc_report(self, output_dir=None):
        """
        Generate comprehensive QC report.
        
        Parameters:
        -----------
        output_dir : str or Path, optional
            Directory to save report files
        """
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(exist_ok=True, parents=True)
        else:
            output_dir = Path('.')
        
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        print("\n" + "="*70)
        print("HEY MIX QC REPORT")
        print("="*70)
        print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("="*70)
        
        # Protein ID comparison
        if self.comparison_results:
            results = self.comparison_results
            
            print("\n--- PROTEIN IDENTIFICATION SUMMARY ---")
            print(f"Total shared proteins: {results['total_shared']}")
            print(f"E100Y75 unique proteins: {results['total_mix1_only']}")
            print(f"E25Y150 unique proteins: {results['total_mix2_only']}")
            
            print("\n--- ORGANISM-SPECIFIC IDENTIFICATION ---")
            print(f"{'Organism':<12} {'Shared':>10} {'E100Y75 Only':>15} {'E25Y150 Only':>15}")
            print("-" * 55)
            for org in ['HeLa', 'E.coli', 'Yeast']:
                shared = results.get(f'{org}_shared', 0)
                mix1 = results.get(f'{org}_mix1_only', 0)
                mix2 = results.get(f'{org}_mix2_only', 0)
                print(f"{org:<12} {shared:>10} {mix1:>15} {mix2:>15}")
            
            # Save detailed comparison
            comparison_file = output_dir / f'protein_comparison_{timestamp}.txt'
            with open(comparison_file, 'w') as f:
                f.write("HEY Mix Protein Comparison Report\n")
                f.write(f"Generated: {datetime.now()}\n\n")
                f.write(f"Total Shared: {results['total_shared']}\n")
                f.write(f"E100Y75 Only: {results['total_mix1_only']}\n")
                f.write(f"E25Y150 Only: {results['total_mix2_only']}\n")
            
            print(f"\nDetailed comparison saved to: {comparison_file}")
        
        # Fold change analysis
        try:
            fc_df = self.calculate_fold_changes()
            
            print("\n--- FOLD CHANGE ANALYSIS ---")
            for org in ['HeLa', 'E.coli', 'Yeast']:
                org_data = fc_df[fc_df['Organism'] == org]['Log2FC'].dropna()
                if len(org_data) > 0:
                    print(f"\n{org}:")
                    print(f"  Median Log2FC: {org_data.median():.3f}")
                    print(f"  Mean Log2FC: {org_data.mean():.3f}")
                    print(f"  Std Dev: {org_data.std():.3f}")
                    print(f"  Proteins analyzed: {len(org_data)}")
            
            # Save fold change data
            fc_file = output_dir / f'fold_changes_{timestamp}.csv'
            fc_df.to_csv(fc_file, index=False)
            print(f"\nFold change data saved to: {fc_file}")
            
            # Generate plots
            overlap_plot = output_dir / f'protein_overlap_{timestamp}.png'
            self.plot_protein_overlap(save_path=overlap_plot)
            
            fc_plot = output_dir / f'fold_change_distribution_{timestamp}.png'
            self.plot_fold_change_distribution(fc_df, save_path=fc_plot)
            
        except Exception as e:
            print(f"\nWarning: Could not complete fold change analysis: {e}")
        
        print("\n" + "="*70)
        print("QC REPORT COMPLETE")
        print("="*70 + "\n")


class HEYAnalyzerGUI:
    """Modern GUI for HEY Analyzer with dark mode support."""

    # Color schemes
    DARK_THEME = {
        'bg': '#1e1e1e',
        'fg': '#e0e0e0',
        'frame_bg': '#2d2d2d',
        'button_bg': '#0e639c',
        'button_fg': '#ffffff',
        'button_hover': '#1177bb',
        'entry_bg': '#3c3c3c',
        'entry_fg': '#e0e0e0',
        'text_bg': '#252526',
        'text_fg': '#d4d4d4',
        'highlight': '#007acc',
        'success': '#4ec9b0',
        'warning': '#ce9178',
        'error': '#f48771',
    }

    LIGHT_THEME = {
        'bg': '#f3f3f3',
        'fg': '#1e1e1e',
        'frame_bg': '#ffffff',
        'button_bg': '#0078d4',
        'button_fg': '#ffffff',
        'button_hover': '#106ebe',
        'entry_bg': '#ffffff',
        'entry_fg': '#1e1e1e',
        'text_bg': '#ffffff',
        'text_fg': '#1e1e1e',
        'highlight': '#0078d4',
        'success': '#107c10',
        'warning': '#ca5010',
        'error': '#e81123',
    }

    def __init__(self, root):
        self.root = root
        self.root.title("Astral HEY Analysis - QC Tool")
        self.root.geometry("1000x700")

        self.analyzer = None
        self.dark_mode = True
        self.theme = self.DARK_THEME

        # File paths
        self.mix1_path = tk.StringVar()
        self.mix2_path = tk.StringVar()
        self.output_dir = tk.StringVar(value="qc_results")

        # Column names
        self.protein_id_col = tk.StringVar(value="Protein.IDs")
        self.intensity_col = tk.StringVar(value="Intensity")

        self.setup_ui()
        self.apply_theme()

    def setup_ui(self):
        """Setup the user interface."""
        # Menu bar
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        view_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="View", menu=view_menu)
        view_menu.add_command(label="Toggle Dark/Light Mode", command=self.toggle_theme)

        # Main container
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(4, weight=1)

        # Title
        title_frame = ttk.Frame(main_frame)
        title_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=(0, 10))

        self.title_label = tk.Label(
            title_frame,
            text="🔬 Astral HEY QC Analysis",
            font=('Segoe UI', 18, 'bold')
        )
        self.title_label.pack()

        subtitle = tk.Label(
            title_frame,
            text="Orbitrap Astral + Vanquish Neo Quality Control",
            font=('Segoe UI', 10)
        )
        subtitle.pack()

        # File selection frame
        file_frame = ttk.LabelFrame(main_frame, text="Data Files", padding="10")
        file_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)
        file_frame.columnconfigure(1, weight=1)

        # Mix 1 (E100Y75)
        ttk.Label(file_frame, text="E100Y75 Mix:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.mix1_entry = ttk.Entry(file_frame, textvariable=self.mix1_path)
        self.mix1_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=5, pady=5)
        self.mix1_btn = ttk.Button(file_frame, text="Browse", command=lambda: self.browse_file(self.mix1_path))
        self.mix1_btn.grid(row=0, column=2, pady=5)

        # Mix 2 (E25Y150)
        ttk.Label(file_frame, text="E25Y150 Mix:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.mix2_entry = ttk.Entry(file_frame, textvariable=self.mix2_path)
        self.mix2_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=5, pady=5)
        self.mix2_btn = ttk.Button(file_frame, text="Browse", command=lambda: self.browse_file(self.mix2_path))
        self.mix2_btn.grid(row=1, column=2, pady=5)

        # Settings frame
        settings_frame = ttk.LabelFrame(main_frame, text="Analysis Settings", padding="10")
        settings_frame.grid(row=2, column=0, sticky=(tk.W, tk.E), pady=5)
        settings_frame.columnconfigure(1, weight=1)
        settings_frame.columnconfigure(3, weight=1)

        # Column settings
        ttk.Label(settings_frame, text="Protein ID Column:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))
        ttk.Entry(settings_frame, textvariable=self.protein_id_col, width=20).grid(row=0, column=1, sticky=tk.W, padx=5)

        ttk.Label(settings_frame, text="Intensity Column:").grid(row=0, column=2, sticky=tk.W, padx=(20, 5))
        ttk.Entry(settings_frame, textvariable=self.intensity_col, width=20).grid(row=0, column=3, sticky=tk.W, padx=5)

        # Output directory
        ttk.Label(settings_frame, text="Output Directory:").grid(row=1, column=0, sticky=tk.W, pady=(10, 0))
        ttk.Entry(settings_frame, textvariable=self.output_dir, width=30).grid(row=1, column=1, columnspan=2, sticky=(tk.W, tk.E), padx=5, pady=(10, 0))
        ttk.Button(settings_frame, text="Browse", command=self.browse_output_dir).grid(row=1, column=3, sticky=tk.W, pady=(10, 0))

        # Control buttons
        btn_frame = ttk.Frame(main_frame)
        btn_frame.grid(row=3, column=0, pady=10)

        self.run_btn = tk.Button(
            btn_frame,
            text="▶ Run Analysis",
            font=('Segoe UI', 11, 'bold'),
            command=self.run_analysis,
            padx=20,
            pady=8
        )
        self.run_btn.pack(side=tk.LEFT, padx=5)

        self.clear_btn = tk.Button(
            btn_frame,
            text="Clear Log",
            font=('Segoe UI', 10),
            command=self.clear_log,
            padx=15,
            pady=6
        )
        self.clear_btn.pack(side=tk.LEFT, padx=5)

        # Progress/Log frame
        log_frame = ttk.LabelFrame(main_frame, text="Analysis Log", padding="5")
        log_frame.grid(row=4, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)

        self.log_text = scrolledtext.ScrolledText(
            log_frame,
            wrap=tk.WORD,
            font=('Consolas', 9),
            height=15
        )
        self.log_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Status bar
        self.status_bar = ttk.Label(main_frame, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.grid(row=5, column=0, sticky=(tk.W, tk.E), pady=(5, 0))

        # Store references for theme application
        self.themed_widgets = {
            'title': self.title_label,
            'subtitle': subtitle,
            'buttons': [self.run_btn, self.clear_btn],
            'log': self.log_text
        }

    def apply_theme(self):
        """Apply current theme colors."""
        theme = self.theme

        # Root window
        self.root.configure(bg=theme['bg'])

        # Title and subtitle
        self.themed_widgets['title'].configure(bg=theme['bg'], fg=theme['fg'])
        self.themed_widgets['subtitle'].configure(bg=theme['bg'], fg=theme['fg'])

        # Buttons
        for btn in self.themed_widgets['buttons']:
            btn.configure(
                bg=theme['button_bg'],
                fg=theme['button_fg'],
                activebackground=theme['button_hover'],
                activeforeground=theme['button_fg'],
                relief=tk.FLAT,
                bd=0
            )

        # Log text
        self.log_text.configure(
            bg=theme['text_bg'],
            fg=theme['text_fg'],
            insertbackground=theme['fg'],
            selectbackground=theme['highlight']
        )

        # Configure ttk styles
        style = ttk.Style()
        style.theme_use('clam')

        style.configure('TFrame', background=theme['bg'])
        style.configure('TLabelframe', background=theme['bg'], foreground=theme['fg'])
        style.configure('TLabelframe.Label', background=theme['bg'], foreground=theme['fg'], font=('Segoe UI', 10, 'bold'))
        style.configure('TLabel', background=theme['bg'], foreground=theme['fg'])
        style.configure('TButton', background=theme['button_bg'], foreground=theme['fg'])
        style.map('TButton', background=[('active', theme['button_hover'])])
        style.configure('TEntry', fieldbackground=theme['entry_bg'], foreground=theme['entry_fg'])

    def toggle_theme(self):
        """Toggle between dark and light mode."""
        self.dark_mode = not self.dark_mode
        self.theme = self.DARK_THEME if self.dark_mode else self.LIGHT_THEME
        self.apply_theme()
        mode = "Dark" if self.dark_mode else "Light"
        self.log(f"Switched to {mode} mode")

    def browse_file(self, var):
        """Browse for TSV file."""
        filename = filedialog.askopenfilename(
            title="Select TSV File",
            filetypes=[("TSV files", "*.tsv"), ("Text files", "*.txt"), ("All files", "*.*")]
        )
        if filename:
            var.set(filename)
            self.log(f"Selected file: {Path(filename).name}")

    def browse_output_dir(self):
        """Browse for output directory."""
        dirname = filedialog.askdirectory(title="Select Output Directory")
        if dirname:
            self.output_dir.set(dirname)
            self.log(f"Output directory: {dirname}")

    def log(self, message, level='INFO'):
        """Add message to log."""
        timestamp = datetime.now().strftime('%H:%M:%S')
        color_tag = 'normal'

        if level == 'SUCCESS':
            color_tag = 'success'
        elif level == 'WARNING':
            color_tag = 'warning'
        elif level == 'ERROR':
            color_tag = 'error'

        self.log_text.insert(tk.END, f"[{timestamp}] {message}\n", color_tag)
        self.log_text.see(tk.END)

        # Configure tags with theme colors
        self.log_text.tag_config('success', foreground=self.theme['success'])
        self.log_text.tag_config('warning', foreground=self.theme['warning'])
        self.log_text.tag_config('error', foreground=self.theme['error'])
        self.log_text.tag_config('normal', foreground=self.theme['text_fg'])

        self.root.update_idletasks()

    def clear_log(self):
        """Clear the log text."""
        self.log_text.delete(1.0, tk.END)
        self.log("Log cleared")

    def update_status(self, message):
        """Update status bar."""
        self.status_bar.config(text=message)
        self.root.update_idletasks()

    def run_analysis(self):
        """Run the HEY analysis in a separate thread."""
        # Validate inputs
        if not self.mix1_path.get() or not self.mix2_path.get():
            messagebox.showerror("Error", "Please select both TSV files")
            return

        if not Path(self.mix1_path.get()).exists():
            messagebox.showerror("Error", f"File not found: {self.mix1_path.get()}")
            return

        if not Path(self.mix2_path.get()).exists():
            messagebox.showerror("Error", f"File not found: {self.mix2_path.get()}")
            return

        # Disable run button
        self.run_btn.config(state=tk.DISABLED)
        self.update_status("Running analysis...")

        # Run in thread to keep GUI responsive
        thread = threading.Thread(target=self._run_analysis_thread, daemon=True)
        thread.start()

    def _run_analysis_thread(self):
        """Thread worker for running analysis."""
        try:
            self.log("="*70)
            self.log("Starting HEY QC Analysis", 'SUCCESS')
            self.log("="*70)

            # Create analyzer with progress callback
            self.analyzer = HEYAnalyzer(progress_callback=self.log)

            # Load data
            self.log("\n>>> Loading Mix 1 (E100Y75)...")
            self.analyzer.load_tsv_data(self.mix1_path.get(), 'E100Y75')

            self.log("\n>>> Loading Mix 2 (E25Y150)...")
            self.analyzer.load_tsv_data(self.mix2_path.get(), 'E25Y150')

            # Compare protein IDs
            self.log("\n>>> Comparing protein identifications...")
            results = self.analyzer.compare_protein_ids(
                protein_id_column=self.protein_id_col.get()
            )

            self.log(f"Total shared proteins: {results['total_shared']}", 'SUCCESS')
            self.log(f"E100Y75 unique: {results['total_mix1_only']}")
            self.log(f"E25Y150 unique: {results['total_mix2_only']}")

            # Calculate fold changes
            self.log("\n>>> Calculating fold changes...")
            fc_df = self.analyzer.calculate_fold_changes(
                intensity_column=self.intensity_col.get(),
                protein_id_column=self.protein_id_col.get()
            )

            # Display organism-specific results
            self.log("\n>>> Organism-specific results:")
            for org in ['HeLa', 'E.coli', 'Yeast']:
                org_data = fc_df[fc_df['Organism'] == org]['Log2FC'].dropna()
                if len(org_data) > 0:
                    self.log(f"  {org}: Median Log2FC = {org_data.median():.3f} ({len(org_data)} proteins)")

            # Generate report
            self.log("\n>>> Generating QC report and plots...")
            output_path = Path(self.output_dir.get())
            output_path.mkdir(exist_ok=True, parents=True)

            self.analyzer.generate_qc_report(output_dir=output_path)

            self.log("\n" + "="*70, 'SUCCESS')
            self.log("Analysis complete! Results saved to:", 'SUCCESS')
            self.log(f"  {output_path.absolute()}", 'SUCCESS')
            self.log("="*70, 'SUCCESS')

            self.update_status("Analysis complete!")
            messagebox.showinfo("Success", f"Analysis complete!\n\nResults saved to:\n{output_path.absolute()}")

        except Exception as e:
            self.log(f"\n!!! ERROR: {str(e)}", 'ERROR')
            self.update_status("Analysis failed")
            messagebox.showerror("Analysis Error", f"An error occurred:\n\n{str(e)}")

        finally:
            # Re-enable run button
            self.run_btn.config(state=tk.NORMAL)


def main():
    """Launch the GUI application."""
    root = tk.Tk()
    HEYAnalyzerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
