#!/usr/bin/env python3
"""
MSPP Data Plotter - Flask Backend API
Provides REST endpoints for MS proteomics data visualization
"""

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
from pathlib import Path
import tempfile
import base64
import io
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

app = Flask(__name__, static_folder='../frontend/dist', static_url_path='')
CORS(app)  # Enable CORS for React dev server

# Set dark mode for matplotlib
plt.style.use('dark_background')


class MSPPDataProcessor:
    """Core data processing logic extracted from GUI version."""
    
    ORGANISMS = ['HeLa', 'E.coli', 'Yeast']
    COLORS = {'HeLa': '#9b59b6', 'E.coli': '#e67e22', 'Yeast': '#16a085', 'Unknown': '#95a5a6'}
    
    ORGANISM_PATTERNS = {
        'HeLa': ['_HUMAN', 'HOMO_SAPIENS'],
        'E.coli': ['_ECOLI', '_ECOL', '_ECO2', '_ECO5', '_ECO7', '_SHIF', '_SHIB', '_SHIS', 'ESCHERICHIA'],
        'Yeast': ['_YEAST', 'SACCHAROMYCES', 'CEREVISIAE']
    }
    
    def __init__(self):
        self.file_to_raw_column = {}
        self._loaded_data = None
        self._loaded_files = None
    
    def _get_cache_key(self, file_paths):
        """Generate cache key from sorted file paths."""
        return tuple(sorted(file_paths))
    
    def identify_organism_vectorized(self, series):
        """Vectorized organism identification."""
        upper = series.fillna('').astype(str).str.upper()
        result = pd.Series('Unknown', index=series.index)
        
        for organism, patterns in self.ORGANISM_PATTERNS.items():
            mask = upper.str.contains('|'.join(patterns), regex=True)
            result = result.where(~mask, organism)
        
        return pd.Categorical(result, categories=self.ORGANISMS + ['Unknown'])
    
    def load_tsv_files(self, file_paths):
        """Load and process TSV files with caching."""
        # Check if we can use cached data
        cache_key = self._get_cache_key(file_paths)
        if self._loaded_files == cache_key and self._loaded_data is not None:
            return self._loaded_data
        
        # Reset file_to_raw_column mapping for new dataset
        self.file_to_raw_column = {}
        all_data = []
        
        for filepath in file_paths:
            # Use low_memory=False and specify dtypes for common columns to speed up parsing
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
            
            if protein_col:
                df['Organism'] = self.identify_organism_vectorized(df[protein_col])
            else:
                df['Organism'] = 'Unknown'
            all_data.append(df)
        
        # Cache the result
        self._loaded_data = pd.concat(all_data, ignore_index=True)
        self._loaded_files = cache_key
        
        return self._loaded_data
    
    def _get_organism_data(self, file_data, intensity_col, organism):
        """Extract positive numeric intensity data for an organism."""
        data = pd.to_numeric(
            file_data[file_data['Organism'] == organism][intensity_col],
            errors='coerce'
        )
        return data[data > 0].dropna()
    
    def _get_hela_median(self, file_data, intensity_col):
        """Get HeLa median for normalization."""
        hela_data = self._get_organism_data(file_data, intensity_col, 'HeLa')
        return hela_data.median() if len(hela_data) > 0 else 1.0
    
    def _calculate_protein_fold_changes(self, file_data, intensity_col):
        """Calculate per-protein log2 fold changes (E.coli/Yeast)."""
        ecoli_data = self._get_organism_data(file_data, intensity_col, 'E.coli')
        yeast_data = self._get_organism_data(file_data, intensity_col, 'Yeast')
        
        if len(ecoli_data) == 0 or len(yeast_data) == 0:
            return None
        
        hela_median = self._get_hela_median(file_data, intensity_col)
        ecoli_norm = ecoli_data / hela_median
        yeast_median = (yeast_data / hela_median).median()
        
        protein_fcs = np.log2(ecoli_norm / yeast_median)
        return protein_fcs.replace([np.inf, -np.inf], np.nan).dropna().values
    
    def _calculate_organism_vs_hela(self, file_data, intensity_col, organism):
        """Calculate per-protein log2 fold changes (Organism/HeLa median)."""
        org_data = self._get_organism_data(file_data, intensity_col, organism)
        
        if len(org_data) == 0:
            return None
        
        hela_median = self._get_hela_median(file_data, intensity_col)
        if hela_median == 1.0:
            return None
        
        protein_fcs = np.log2(org_data / hela_median)
        return protein_fcs.replace([np.inf, -np.inf], np.nan).dropna().values
    
    def generate_bar_chart(self, data):
        """Generate protein ID counts bar chart."""
        fig, ax = plt.subplots(figsize=(12, 7))
        
        counts = data.groupby(['Source_File', 'Organism'], observed=True).size().unstack(fill_value=0)
        org_order = self.ORGANISMS + ['Unknown']
        counts = counts.reindex(columns=[col for col in org_order if col in counts.columns], fill_value=0)
        plot_colors = [self.COLORS[col] for col in counts.columns]
        
        counts.plot(kind='bar', stacked=True, ax=ax, color=plot_colors, 
                   edgecolor='black', linewidth=0.5, alpha=0.8)
        ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Protein IDs', fontsize=12, fontweight='bold')
        ax.set_title('Protein ID Counts by Organism', fontsize=14, fontweight='bold')
        ax.legend(title='Organism', fontsize=10, loc='upper right')
        ax.grid(axis='y', alpha=0.3)
        ax.tick_params(axis='x', rotation=45, labelbottom=True)
        plt.setp(ax.xaxis.get_majorticklabels(), ha='right')
        plt.tight_layout()
        
        return self._fig_to_base64(fig)
    
    def generate_fold_change_plot(self, data):
        """Generate E.coli vs Yeast fold change plot."""
        fig, ax = plt.subplots(figsize=(12, 7))
        
        # Calculate fold changes per sample
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
                results.append((source_file, fcs))
        
        if not results:
            plt.close(fig)
            return None
        
        # Sort by median
        sorted_results = sorted(results, key=lambda x: np.median(x[1]), reverse=True)
        sample_names = [x[0] for x in sorted_results]
        sample_fcs = [x[1] for x in sorted_results]
        
        # Create box plots
        positions = np.arange(1, len(sample_fcs) + 1)
        bp = ax.boxplot(sample_fcs, positions=positions, widths=0.6,
                        patch_artist=True, showfliers=True, showmeans=True,
                        flierprops=dict(marker='o', markerfacecolor='#3498db', 
                                       markersize=3, alpha=0.4, markeredgecolor='none'),
                        meanprops=dict(marker='s', markerfacecolor='white', 
                                      markeredgecolor='white', markersize=4))
        
        self._style_boxplot(bp)
        
        ax.axhline(y=0, color='#f39c12', linestyle='--', linewidth=2, alpha=0.9)
        
        # Labels
        ax.set_ylabel('Log2 Abundance Ratio (E.coli / Yeast median)', fontsize=12, fontweight='bold')
        ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
        ax.set_title('Protein Fold Change (HeLa-Normalized)', fontsize=14, fontweight='bold')
        ax.set_xticks(positions)
        ax.set_xticklabels([s.replace('report.pg_matrix_', '') for s in sample_names], 
                          rotation=45, ha='right', fontsize=9)
        ax.grid(axis='y', alpha=0.3)
        
        # Legend
        legend_elements = [
            Line2D([0], [0], color='#f39c12', linestyle='--', linewidth=2, label='Expected 1:1'),
            Line2D([0], [0], color='#2c3e50', linewidth=2, label='Median')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
        
        plt.tight_layout()
        return self._fig_to_base64(fig)
    
    def generate_organisms_vs_hela(self, data):
        """Generate organisms vs HeLa comparison plot."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
        
        for ax, organism in [(ax1, 'E.coli'), (ax2, 'Yeast')]:
            sample_data = []
            sample_names = []
            
            for source_file in data['Source_File'].unique():
                if source_file not in self.file_to_raw_column:
                    continue
                
                intensity_col = self.file_to_raw_column[source_file]
                file_data = data[data['Source_File'] == source_file]
                
                fcs = self._calculate_organism_vs_hela(file_data, intensity_col, organism)
                if fcs is not None and len(fcs) > 0:
                    sample_data.append(fcs)
                    sample_names.append(source_file)
            
            if not sample_data:
                ax.text(0.5, 0.5, f'No {organism} data', ha='center', va='center', 
                       transform=ax.transAxes, fontsize=14)
                continue
            
            # Sort by median
            sorted_pairs = sorted(zip(sample_data, sample_names), 
                                 key=lambda x: np.median(x[0]), reverse=True)
            sample_data, sample_names = zip(*sorted_pairs)
            
            positions = np.arange(1, len(sample_data) + 1)
            bp = ax.boxplot(sample_data, positions=positions, widths=0.6,
                           patch_artist=True, showfliers=True, showmeans=True,
                           flierprops=dict(marker='o', markerfacecolor='#e74c3c', 
                                          markersize=3, alpha=0.4, markeredgecolor='none'),
                           meanprops=dict(marker='s', markerfacecolor='white', 
                                         markeredgecolor='white', markersize=4))
            
            self._style_boxplot(bp)
            
            ax.axhline(y=0, color='#f39c12', linestyle='--', linewidth=2, alpha=0.9)
            ax.set_ylabel(f'Log2 Ratio ({organism} / HeLa median)', fontsize=11, fontweight='bold')
            ax.set_xlabel('Sample', fontsize=11, fontweight='bold')
            ax.set_title(f'{organism} vs HeLa', fontsize=13, fontweight='bold')
            ax.set_xticks(positions)
            ax.set_xticklabels([s.replace('report.pg_matrix_', '') for s in sample_names], 
                              rotation=45, ha='right', fontsize=8)
            ax.grid(axis='y', alpha=0.3)
        
        plt.suptitle('Spike-in Validation', fontsize=15, fontweight='bold', y=0.98)
        plt.tight_layout()
        return self._fig_to_base64(fig)
    
    def generate_grouped_fold_change(self, data, pattern):
        """Generate grouped fold change plot by regex pattern."""
        import re
        
        try:
            regex = re.compile(pattern, re.IGNORECASE)
        except re.error as e:
            return {'error': f'Invalid regex pattern: {str(e)}'}
        
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
            available_files = ', '.join(list(data['Source_File'].unique())[:3])
            return {'error': f'No files matched pattern \'{pattern}\'. Available files: {available_files}...'}
        
        # Sort groups by median fold change
        sorted_groups = sorted(groups.keys(), key=lambda g: np.median(groups[g]), reverse=True)
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(12, 7))
        
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
        y_max = max([max(d) for d in box_data]) if box_data else 1
        for i, group in enumerate(sorted_groups):
            med = np.median(groups[group])
            ax.text(i + 1.3, med, f'{med:.2f}', fontsize=10, va='center', color='#f39c12', fontweight='bold')
            ax.text(i + 1, y_max + 0.3, f'{len(groups[group])} proteins', 
                   fontsize=9, ha='center', va='bottom', color='white')
        
        # Labels
        ax.set_ylabel('Log2 Abundance Ratio (E.coli / Yeast)', fontsize=12, fontweight='bold')
        ax.set_xlabel('Group', fontsize=12, fontweight='bold')
        ax.set_title('Per-Protein Fold Change by Group (HeLa-Normalized)', fontsize=14, fontweight='bold')
        ax.set_xticks(positions)
        ax.set_xticklabels(sorted_groups, rotation=45, ha='right', fontsize=10)
        ax.grid(axis='y', alpha=0.3)
        
        # Legend
        legend_elements = [
            Line2D([0], [0], color='#f39c12', linestyle='--', linewidth=2, label='Expected 1:1'),
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
        
        img_base64 = self._fig_to_base64(fig)
        
        return {
            'image': img_base64,
            'unmatched_count': len(unmatched),
            'group_count': len(sorted_groups)
        }
    
    def _style_boxplot(self, bp):
        """Apply consistent styling to boxplot elements."""
        for patch in bp['boxes']:
            patch.set_facecolor('#5dade2')
            patch.set_alpha(0.7)
            patch.set_edgecolor('white')
        for element in ['whiskers', 'caps']:
            plt.setp(bp[element], color='white', linewidth=1)
        plt.setp(bp['medians'], color='#2c3e50', linewidth=2)
    
    def _fig_to_base64(self, fig):
        """Convert matplotlib figure to base64 encoded PNG."""
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        plt.close(fig)
        return img_base64


# Global processor instance
processor = MSPPDataProcessor()
uploaded_files = {}  # Store uploaded files temporarily


@app.route('/')
def serve_react_app():
    """Serve React frontend."""
    return send_from_directory(app.static_folder, 'index.html')


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
        try:
            Path(filepath).unlink(missing_ok=True)
        except Exception:
            pass
    uploaded_files.clear()
    return jsonify({'message': 'All files cleared'})


@app.route('/api/plot/bar-chart', methods=['POST'])
def generate_bar_chart():
    """Generate protein ID bar chart."""
    if not uploaded_files:
        return jsonify({'error': 'No files uploaded'}), 400
    
    try:
        data = processor.load_tsv_files(list(uploaded_files.values()))
        img_base64 = processor.generate_bar_chart(data)
        return jsonify({'image': img_base64})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/plot/fold-change', methods=['POST'])
def generate_fold_change():
    """Generate E.coli vs Yeast fold change plot."""
    if not uploaded_files:
        return jsonify({'error': 'No files uploaded'}), 400
    
    try:
        data = processor.load_tsv_files(list(uploaded_files.values()))
        img_base64 = processor.generate_fold_change_plot(data)
        if img_base64 is None:
            return jsonify({'error': 'No valid data for fold change analysis'}), 400
        return jsonify({'image': img_base64})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/plot/organisms-vs-hela', methods=['POST'])
def generate_organisms_vs_hela():
    """Generate organisms vs HeLa comparison plot."""
    if not uploaded_files:
        return jsonify({'error': 'No files uploaded'}), 400
    
    try:
        data = processor.load_tsv_files(list(uploaded_files.values()))
        img_base64 = processor.generate_organisms_vs_hela(data)
        return jsonify({'image': img_base64})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/plot/grouped-fold-change', methods=['POST'])
def generate_grouped_fold_change():
    """Generate grouped fold change plot."""
    if not uploaded_files:
        return jsonify({'error': 'No files uploaded'}), 400
    
    try:
        request_data = request.get_json()
        pattern = request_data.get('pattern', r'(E\d+)')
        
        data = processor.load_tsv_files(list(uploaded_files.values()))
        result = processor.generate_grouped_fold_change(data, pattern)
        
        if 'error' in result:
            return jsonify(result), 400
        
        return jsonify(result)
    except Exception as e:
        return jsonify({'error': str(e)}), 500


if __name__ == '__main__':
    app.run(debug=True, port=5000)
