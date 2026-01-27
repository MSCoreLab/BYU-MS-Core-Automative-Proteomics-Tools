#!/usr/bin/env python3
"""
MSPP Data Plotter - Flask Backend API
Modular refactor: Logic is in logic.py, Routes are here.
"""

import contextlib
import io
import logging
import mimetypes
import os
import tempfile
from pathlib import Path

from flask import Flask, jsonify, request, send_file, send_from_directory
from flask_cors import CORS
from werkzeug.utils import secure_filename

# Import our custom logic
from .logic import DataProcessor, PlotGenerator, fig_to_base64

# Force correct MIME types
mimetypes.add_type('application/javascript', '.js')
mimetypes.add_type('text/css', '.css')

app = Flask(__name__, static_folder='../frontend/dist', static_url_path='')
CORS(app)

# Global instances
processor = DataProcessor()
plotter = PlotGenerator(processor)
uploaded_files = {}

@app.after_request
def add_security_headers(response):
    """Add security and performance headers for local and work environments."""
    # Ensure no old data is served
    response.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '0'
    # Strict but functional CSP
    response.headers['Content-Security-Policy'] = "default-src 'self' 'unsafe-inline' 'unsafe-eval' data: blob:;"
    return response

@app.route('/')
def serve_react_app():
    return send_from_directory(app.static_folder, 'index.html')

@app.route('/api/health')
def health_check():
    return jsonify({'status': 'ok'})

@app.route('/api/upload', methods=['POST'])
def upload_files():
    if 'files' not in request.files:
        return jsonify({'error': 'No files provided'}), 400

    files = request.files.getlist('files')
    temp_paths = []

    for file in files:
        if file and file.filename and file.filename.lower().endswith(('.tsv', '.txt')):
            # SECURITY: Use secure_filename to prevent path traversal
            safe_name = secure_filename(file.filename)
            temp_path = Path(tempfile.gettempdir()) / safe_name
            file.save(temp_path)
            uploaded_files[safe_name] = str(temp_path)
            temp_paths.append(safe_name)

    return jsonify({
        'message': f'{len(temp_paths)} files uploaded successfully',
        'files': temp_paths
    })

@app.route('/api/files', methods=['GET', 'DELETE'])
def manage_files():
    if request.method == 'DELETE':
        for path in uploaded_files.values():
            with contextlib.suppress(Exception):
                Path(path).unlink(missing_ok=True)
        uploaded_files.clear()
        processor.clear_cache()
        return jsonify({'message': 'Cleared'})
    return jsonify({'files': list(uploaded_files.keys())})

@app.route('/api/plot/<chart_type>', methods=['POST'])
def generate_plot(chart_type):
    if not uploaded_files:
        return jsonify({'error': 'No files uploaded'}), 400

    try:
        data = processor.load_data(list(uploaded_files.values()))
        if chart_type == 'bar-chart':
            fig = plotter.create_bar_chart_figure(data)
        elif chart_type == 'sample-comparison':
            fig = plotter.create_comparison_figure(data)
        else:
            return jsonify({'error': 'Invalid plot type'}), 400

        return jsonify({'image': fig_to_base64(fig)})
    except Exception as e:
        logging.exception(f"Plot generation failed: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/export/<chart_type>', methods=['POST'])
def export_plot(chart_type):
    if not uploaded_files:
        return jsonify({'error': 'No files uploaded'}), 400

    try:
        data = processor.load_data(list(uploaded_files.values()))
        if chart_type == 'bar-chart':
            fig = plotter.create_bar_chart_figure(data, figsize=(10, 6))
            name = 'protein_id_bar_chart.png'
        elif chart_type == 'sample-comparison':
            fig = plotter.create_comparison_figure(data, figsize=(18, 16))
            name = 'intensity_ratio_comparison.png'
        else:
            return jsonify({'error': 'Invalid plot type'}), 400

        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
        buf.seek(0)
        return send_file(buf, mimetype='image/png', as_attachment=True, download_name=name)
    except Exception as e:
        logging.exception(f"Export failed: {e}")
        return jsonify({'error': str(e)}), 500

if __name__ == "__main__":
    # For standalone run (dev); enable debug only if FLASK_DEBUG=1
    debug_mode = bool(int(os.getenv("FLASK_DEBUG", "0")))
    app.run(port=8050, debug=debug_mode)
