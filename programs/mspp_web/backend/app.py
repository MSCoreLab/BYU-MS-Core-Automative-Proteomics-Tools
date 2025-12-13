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
    
    def _extract_mix_identifier(self, filename):
        """Extract mix identifier from filename, excluding E25/E100 prefix.
        
        For 'report.pg_matrix_E25_30_4_440960_600.tsv', returns '30_4_440960_600'.
        This ensures E25 and E100 files from the same mix are grouped together.
        """
        import re
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
    
    def generate_sample_comparison(self, data):
        """Generate E25 vs E100 intensity comparison plot."""
        import re
        
        # Get all sample files
        sample_files = sorted(data['Source_File'].unique())
        
        if len(sample_files) == 0:
            return None
        
        # Group samples by mix identifier
        mix_groups = {}
        for source_file in sample_files:
            mix_id = self._extract_mix_identifier(source_file)
            if mix_id not in mix_groups:
                mix_groups[mix_id] = []
            mix_groups[mix_id].append(source_file)
        
        sorted_mixes = sorted(mix_groups.keys())
        
        # Calculate intensities for E25 and E100 separately
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
                # E.coli
                e25_ecoli = self._calculate_sample_intensities(
                    data[data['Source_File'] == e25_file], e25_file, 'E.coli'
                )
                e100_ecoli = self._calculate_sample_intensities(
                    data[data['Source_File'] == e100_file], e100_file, 'E.coli'
                )
                
                # Yeast
                e25_yeast = self._calculate_sample_intensities(
                    data[data['Source_File'] == e25_file], e25_file, 'Yeast'
                )
                e100_yeast = self._calculate_sample_intensities(
                    data[data['Source_File'] == e100_file], e100_file, 'Yeast'
                )
                
                if e25_ecoli is not None:
                    ecoli_results.append(('E25', e25_ecoli, mix_id))
                    current_position += 1
                if e100_ecoli is not None:
                    ecoli_results.append(('E100', e100_ecoli, mix_id))
                    current_position += 1
                
                if e25_yeast is not None:
                    yeast_results.append(('E25', e25_yeast, mix_id))
                if e100_yeast is not None:
                    yeast_results.append(('E100', e100_yeast, mix_id))
            
            if current_position > 0:
                mix_boundaries.append(current_position)
        
        if not ecoli_results and not yeast_results:
            return None
        
        # Create the plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))
        
        # E.coli plot
        if ecoli_results:
            ecoli_labels = [r[0] for r in ecoli_results]
            ecoli_data = [r[1] for r in ecoli_results]
            positions = np.arange(1, len(ecoli_data) + 1)
            
            bp1 = ax1.boxplot(
                ecoli_data, positions=positions, widths=0.6,
                patch_artist=True, showfliers=True, showmeans=True,
                flierprops=dict(marker='o', markerfacecolor='#e67e22', markersize=3, alpha=0.4, markeredgecolor='none'),
                meanprops=dict(marker='s', markerfacecolor='white', markeredgecolor='white', markersize=5)
            )
            
            for i, (patch, label) in enumerate(zip(bp1['boxes'], ecoli_labels)):
                if label == 'E25':
                    patch.set_facecolor('#e67e22')
                else:
                    patch.set_facecolor('#3498db')
                patch.set_alpha(0.7)
                patch.set_edgecolor('white')
                patch.set_linewidth(1.5)
            
            plt.setp(bp1['whiskers'], color='white', linewidth=1.5)
            plt.setp(bp1['caps'], color='white', linewidth=1.5)
            plt.setp(bp1['medians'], color='#2c3e50', linewidth=2.5)
            
            # Add median annotations
            for i, intensities in enumerate(ecoli_data):
                median_val = np.median(intensities)
                ax1.text(i + 1.35, median_val, f'{median_val:.2f}',
                        fontsize=9, va='center', color='white', fontweight='bold')
            
            # Add visual separators
            ylim = ax1.get_ylim()
            prev_boundary = 0
            for idx, boundary in enumerate(mix_boundaries[:-1]):
                ax1.axvline(x=boundary + 0.5, color='#555555', linestyle='-', linewidth=2, alpha=0.8)
                mid_point = (prev_boundary + boundary) / 2 + 0.5
                if idx < len(sorted_mixes):
                    ax1.text(mid_point, ylim[1] * 0.95, f'Mix: {sorted_mixes[idx]}',
                            ha='center', va='top', fontsize=9, color='#aaaaaa',
                            bbox=dict(boxstyle='round,pad=0.3', facecolor='#2c2c2c', edgecolor='#555555', alpha=0.8))
                prev_boundary = boundary
            
            if sorted_mixes:
                mid_point = (prev_boundary + len(ecoli_data)) / 2 + 0.5
                ax1.text(mid_point, ylim[1] * 0.95, f'Mix: {sorted_mixes[-1]}',
                        ha='center', va='top', fontsize=9, color='#aaaaaa',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='#2c2c2c', edgecolor='#555555', alpha=0.8))
            
            ax1.axhline(y=0, color='#f39c12', linestyle='--', linewidth=2, alpha=0.9, label='Reference (1:1)')
            ax1.set_ylabel('Log2 Intensity (HeLa-Normalized)', fontsize=12, fontweight='bold')
            ax1.set_xlabel('Sample', fontsize=12, fontweight='bold')
            ax1.set_title('E.coli Intensity Comparison', fontsize=14, fontweight='bold')
            ax1.set_xticks(positions)
            ax1.set_xticklabels(ecoli_labels, rotation=0, ha='center', fontsize=10)
            ax1.grid(axis='y', alpha=0.3)
            ax1.legend(fontsize=9)
        else:
            ax1.text(0.5, 0.5, 'No E.coli data', ha='center', va='center', 
                    transform=ax1.transAxes, fontsize=14)
        
        # Yeast plot
        if yeast_results:
            yeast_labels = [r[0] for r in yeast_results]
            yeast_data = [r[1] for r in yeast_results]
            positions = np.arange(1, len(yeast_data) + 1)
            
            bp2 = ax2.boxplot(
                yeast_data, positions=positions, widths=0.6,
                patch_artist=True, showfliers=True, showmeans=True,
                flierprops=dict(marker='o', markerfacecolor='#16a085', markersize=3, alpha=0.4, markeredgecolor='none'),
                meanprops=dict(marker='s', markerfacecolor='white', markeredgecolor='white', markersize=5)
            )
            
            for i, (patch, label) in enumerate(zip(bp2['boxes'], yeast_labels)):
                if label == 'E25':
                    patch.set_facecolor('#16a085')
                else:
                    patch.set_facecolor('#9b59b6')
                patch.set_alpha(0.7)
                patch.set_edgecolor('white')
                patch.set_linewidth(1.5)
            
            plt.setp(bp2['whiskers'], color='white', linewidth=1.5)
            plt.setp(bp2['caps'], color='white', linewidth=1.5)
            plt.setp(bp2['medians'], color='#2c3e50', linewidth=2.5)
            
            # Add median annotations
            for i, intensities in enumerate(yeast_data):
                median_val = np.median(intensities)
                ax2.text(i + 1.35, median_val, f'{median_val:.2f}',
                        fontsize=9, va='center', color='white', fontweight='bold')
            
            # Add visual separators
            ylim = ax2.get_ylim()
            prev_boundary = 0
            for idx, boundary in enumerate(mix_boundaries[:-1]):
                ax2.axvline(x=boundary + 0.5, color='#555555', linestyle='-', linewidth=2, alpha=0.8)
                mid_point = (prev_boundary + boundary) / 2 + 0.5
                if idx < len(sorted_mixes):
                    ax2.text(mid_point, ylim[1] * 0.95, f'Mix: {sorted_mixes[idx]}',
                            ha='center', va='top', fontsize=9, color='#aaaaaa',
                            bbox=dict(boxstyle='round,pad=0.3', facecolor='#2c2c2c', edgecolor='#555555', alpha=0.8))
                prev_boundary = boundary
            
            if sorted_mixes:
                mid_point = (prev_boundary + len(yeast_data)) / 2 + 0.5
                ax2.text(mid_point, ylim[1] * 0.95, f'Mix: {sorted_mixes[-1]}',
                        ha='center', va='top', fontsize=9, color='#aaaaaa',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='#2c2c2c', edgecolor='#555555', alpha=0.8))
            
            ax2.axhline(y=0, color='#f39c12', linestyle='--', linewidth=2, alpha=0.9, label='Reference (1:1)')
            ax2.set_ylabel('Log2 Intensity (HeLa-Normalized)', fontsize=12, fontweight='bold')
            ax2.set_xlabel('Sample', fontsize=12, fontweight='bold')
            ax2.set_title('Yeast Intensity Comparison', fontsize=14, fontweight='bold')
            ax2.set_xticks(positions)
            ax2.set_xticklabels(yeast_labels, rotation=0, ha='center', fontsize=10)
            ax2.grid(axis='y', alpha=0.3)
            ax2.legend(fontsize=9)
        else:
            ax2.text(0.5, 0.5, 'No Yeast data', ha='center', va='center',
                    transform=ax2.transAxes, fontsize=14)
        
        plt.suptitle('E25 vs E100 Intensity Comparison by Mix (HeLa-Normalized)',
                    fontsize=15, fontweight='bold', y=0.98)
        plt.tight_layout()
        return self._fig_to_base64(fig)
    
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


@app.route('/api/plot/sample-comparison', methods=['POST'])
def generate_sample_comparison():
    """Generate E25 vs E100 intensity comparison plot."""
    if not uploaded_files:
        return jsonify({'error': 'No files uploaded'}), 400
    
    try:
        data = processor.load_tsv_files(list(uploaded_files.values()))
        img_base64 = processor.generate_sample_comparison(data)
        if img_base64 is None:
            return jsonify({'error': 'No valid E25/E100 pairs found'}), 400
        return jsonify({'image': img_base64})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


if __name__ == '__main__':
    app.run(debug=True, port=5000)
