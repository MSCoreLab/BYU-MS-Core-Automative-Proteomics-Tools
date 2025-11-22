#!/usr/bin/env python3
"""
Simple Protein Data Plotter
Quick visualization of protein IDs and relative abundances from TSV files
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


class SimpleProteinPlotter:
    """Simple GUI for plotting protein data."""
    
    # Class-level constants
    ORGANISMS = ['HeLa', 'E.coli', 'Yeast']
    COLORS = {'HeLa': '#9b59b6', 'E.coli': '#e67e22', 'Yeast': '#16a085', 'Unknown': '#95a5a6'}

    def __init__(self, root):
        self.root = root
        self.root.title("Simple Protein Plotter")
        self.root.geometry("450x600")
        
        self.files = []
        self.setup_ui()

    def setup_ui(self):
        """Setup the user interface."""
        main_frame = ttk.Frame(self.root, padding="20")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        
        # Title
        ttk.Label(main_frame, text="Protein Data Plotter", 
                 font=('Arial', 16, 'bold')).grid(row=0, column=0, columnspan=2, pady=(0, 20))
        
        # File selection
        ttk.Label(main_frame, text="Select TSV files:").grid(row=1, column=0, sticky=tk.W, pady=5)
        
        self.file_listbox = tk.Listbox(main_frame, height=6, width=50)
        self.file_listbox.grid(row=2, column=0, columnspan=2, pady=5, sticky=(tk.W, tk.E))
        
        # Buttons
        btn_frame = ttk.Frame(main_frame)
        btn_frame.grid(row=3, column=0, columnspan=2, pady=10)
        
        ttk.Button(btn_frame, text="Add Files", command=self.add_files).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Clear", command=self.clear_files).pack(side=tk.LEFT, padx=5)
        
        # Plot buttons
        plot_frame = ttk.LabelFrame(main_frame, text="Generate Plots", padding="10")
        plot_frame.grid(row=4, column=0, columnspan=2, pady=20, sticky=(tk.W, tk.E))
        
        ttk.Button(plot_frame, text="📊 Protein ID Bar Chart", 
                  command=self.plot_protein_ids, width=35).pack(pady=5, fill=tk.X, padx=10)
        
        # File selector for box plot
        selector_frame = ttk.Frame(plot_frame)
        selector_frame.pack(pady=10, fill=tk.X, padx=10)
        
        ttk.Label(selector_frame, text="Select file/sample:").pack(side=tk.LEFT, padx=(0, 5))
        self.selected_file = tk.StringVar()
        self.file_dropdown = ttk.Combobox(selector_frame, textvariable=self.selected_file, 
                                          state='readonly', width=35)
        self.file_dropdown.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        ttk.Button(plot_frame, text="📦 Abundance Box Plot (by Organism)", 
                  command=self.plot_abundances, width=35).pack(pady=5, fill=tk.X, padx=10)
        
        ttk.Button(plot_frame, text="📈 Both Plots", 
                  command=self.plot_both, width=35).pack(pady=5, fill=tk.X, padx=10)

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
        
        # Update dropdown
        self.update_file_dropdown()

    def clear_files(self):
        """Clear file list."""
        self.files.clear()
        self.file_listbox.delete(0, tk.END)
        self.file_dropdown['values'] = []
        self.selected_file.set('')
    
    def update_file_dropdown(self):
        """Update the file dropdown with current files."""
        file_names = [Path(f).stem for f in self.files]
        self.file_dropdown['values'] = file_names
        if file_names and not self.selected_file.get():
            self.selected_file.set(file_names[0])

    def identify_organism(self, protein_name):
        """Identify organism from protein name."""
        if pd.isna(protein_name):
            return 'Unknown'
        
        protein_name = str(protein_name).upper()
        
        # Human/HeLa proteins
        if any(x in protein_name for x in ['_HUMAN', 'HOMO_SAPIENS']):
            return 'HeLa'
        # E. coli proteins
        elif any(x in protein_name for x in ['_ECOLI', '_ECOL', '_ECO2', '_ECO5', '_ECO7', 
                                               '_SHIF', '_SHIB', '_SHIS', 'ESCHERICHIA']):
            return 'E.coli'
        # Yeast proteins
        elif any(x in protein_name for x in ['_YEAST', 'SACCHAROMYCES', 'CEREVISIAE']):
            return 'Yeast'
        else:
            return 'Unknown'
    
    def load_data(self):
        """Load data from selected files."""
        if not self.files:
            messagebox.showerror("Error", "Please add at least one TSV file")
            return None
        
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
                                   if col in df.columns), None)
                if not protein_col:
                    protein_cols = [col for col in df.columns if 'protein' in col.lower()]
                    protein_col = protein_cols[0] if protein_cols else None
                
                df['Organism'] = df[protein_col].apply(self.identify_organism) if protein_col else 'Unknown'
                all_data.append(df)
                
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load {Path(filepath).name}:\n{e}")
                return None
        
        return pd.concat(all_data, ignore_index=True)

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
    
    def _create_box_plot(self, data, selected, ax):
        """Helper: Create box plot for selected file."""
        file_data = data[data['Source_File'] == selected]
        if file_data.empty:
            messagebox.showerror("Error", f"No data found for {selected}")
            return False
        
        if selected not in self.file_to_raw_column:
            messagebox.showerror("Error", f"Could not find intensity column for {selected}")
            return False
        
        intensity_col = self.file_to_raw_column[selected]
        if intensity_col not in file_data.columns:
            messagebox.showerror("Error", f"Column {intensity_col} not found in filtered data")
            return False
        
        plot_data, labels, plot_colors = [], [], []
        for organism in self.ORGANISMS:
            org_data = pd.to_numeric(file_data[file_data['Organism'] == organism][intensity_col], errors='coerce')
            org_data = org_data[org_data > 0]
            if len(org_data) > 0:
                plot_data.append(np.log10(org_data))
                labels.append(organism)
                plot_colors.append(self.COLORS[organism])
        
        if not plot_data:
            messagebox.showerror("Error", "No organism data with intensity values found")
            return False
        
        bp = ax.boxplot(plot_data, tick_labels=labels, patch_artist=True, showfliers=False, widths=0.6)
        for patch, color in zip(bp['boxes'], plot_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
            patch.set_edgecolor('black')
            patch.set_linewidth(1.5)
        
        for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(bp[element], color='black', linewidth=1.5)
        
        ax.set_ylabel('Log10 Intensity', fontsize=12, fontweight='bold')
        ax.set_xlabel('Organism', fontsize=12, fontweight='bold')
        ax.set_title(f'Protein Abundance - {selected}', fontsize=14, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        for i, (label, data_arr) in enumerate(zip(labels, plot_data), 1):
            ax.text(i, ax.get_ylim()[1] * 0.95, f'n={len(data_arr)}', 
                   ha='center', fontsize=9, style='italic')
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

    def plot_abundances(self):
        """Plot relative abundances as box plots by organism for selected file."""
        if not self.selected_file.get():
            messagebox.showerror("Error", "Please select a file from the dropdown")
            return
        
        data = self.load_data()
        if data is None:
            return
        
        fig, ax = plt.subplots(figsize=(10, 7))
        if self._create_box_plot(data, self.selected_file.get(), ax):
            plt.tight_layout()
            plt.show()

    def plot_both(self):
        """Generate both plots simultaneously."""
        if not self.selected_file.get():
            messagebox.showerror("Error", "Please select a file from the dropdown for the box plot")
            return
        
        data = self.load_data()
        if data is None:
            return
        
        fig1, ax1 = plt.subplots(figsize=(12, 7))
        fig2, ax2 = plt.subplots(figsize=(10, 7))
        
        self._create_bar_chart(data, ax1)
        fig1.tight_layout()
        
        if self._create_box_plot(data, self.selected_file.get(), ax2):
            fig2.tight_layout()
            plt.show()


def main():
    """Launch the application."""
    root = tk.Tk()
    SimpleProteinPlotter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
