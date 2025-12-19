# MS Protein & Peptide Data Plotter (Web App)

Modern web application for analyzing and visualizing mass spectrometry proteomics data from Orbitrap Astral MS (DIA-NN output).

## Features

### 📊 Visualizations
- **Protein Count Bar Chart**: Stacked bar chart showing organism composition across samples
- **Intensity Ratio Comparison**: 3-panel box plots (HeLa, E.coli, Yeast) showing log2 intensity ratios

### 📥 File Management
- Drag-and-drop or click to upload TSV files
- Multi-file support (E25 and E100 samples)
- Real-time file list display
- Clear all files with one click

### 💾 Export Options
- Export individual plots (300 DPI PNG)
- Export all plots as ZIP file
- High-quality images suitable for publication

### 🔬 Analysis Features
- Automatic organism identification (HeLa, E.coli, Yeast)
- Consensus protein calculation (proteins in both E25 and E100)
- Log2 intensity ratio calculations
- Quality control validation

## Quick Start

### Prerequisites
- Python 3.14+ (or compatible)
- Node.js 18+ and npm
- Modern web browser (Chrome, Firefox, Edge, Safari)

### Installation

1. **Set up backend**
   ```bash
   cd backend
   pip install -r requirements.txt
   ```

2. **Set up frontend**
   ```bash
   cd frontend
   npm install
   ```

### Running the Application

#### Development Mode

**Terminal 1 - Backend:**
```bash
cd backend
python app.py
```
Backend runs on `http://localhost:5000`

**Terminal 2 - Frontend:**
```bash
cd frontend
npm run dev
```
Frontend runs on `http://localhost:5173`

#### Production Mode

1. **Build frontend:**
   ```bash
   cd frontend
   npm run build
   ```

2. **Run backend (serves built frontend):**
   ```bash
   cd backend
   python app.py
   ```

3. **Access application:**
   Open `http://localhost:5000` in your browser

## Usage

### 1. Upload Files
- Click the upload zone or drag-and-drop TSV files
- Upload E25 and E100 sample files (e.g., `report.pg_matrix_E25_*.tsv` and `report.pg_matrix_E100_*.tsv`)
- Files appear in the list below the upload zone

### 2. Generate Plots
- Click **"Generate Protein Count Plot"** for the bar chart
- Click **"Generate Intensity Ratio Plot"** for the comparison plots
- Plots appear in the main viewing area

### 3. Export Results
- **Export Current Plot**: Click the export button below the displayed plot (300 DPI PNG)
- **Export All Plots (ZIP)**: Click the green "Export All Plots (ZIP)" button to download both plots in one archive

### 4. Clear Files
- Click **"Clear All Files"** to remove uploaded files and start over

## Expected File Format

Input TSV files should follow DIA-NN protein group matrix format:

```
Protein.Ids	Protein.Names	Genes	First.Protein.Description	Proteotypic	Protein.Q.Value	Global.Q.Value	[Sample columns...]
```

**File Naming Convention:**
- E25 samples: `report.pg_matrix_E25_*`
- E100 samples: `report.pg_matrix_E100_*`

## Algorithm Details

### Organism Identification
- **HeLa**: Proteins matching `sp\|` (UniProt SwissProt)
- **E.coli**: Proteins containing `ECOLI`
- **Yeast**: Proteins containing `YEAST`
- **Unknown**: Everything else

### Intensity Ratio Calculation
1. Identify consensus proteins (present in both E25 and E100)
2. Calculate log2(E25/E100) for each protein
3. Group by organism
4. Display as box plots with expected ratio lines:
   - HeLa: 0 (constant concentration)
   - E.coli: -2 (log2(25/100))
   - Yeast: +1 (log2(150/75))

## API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/api/upload` | Upload TSV file |
| GET | `/api/files` | List uploaded files |
| DELETE | `/api/files/<filename>` | Delete specific file |
| POST | `/api/generate-plot` | Generate protein count plot |
| POST | `/api/generate-comparison` | Generate intensity ratio plot |
| POST | `/api/export/bar` | Export bar chart (PNG) |
| POST | `/api/export/comparison` | Export comparison plot (PNG) |
| POST | `/api/export/all` | Export both plots (ZIP) |

## Tech Stack

### Backend
- **Flask**: Web framework and REST API
- **Pandas**: Data processing and analysis
- **NumPy**: Numerical operations
- **Matplotlib**: Visualization generation

### Frontend
- **React 18**: UI framework
- **TypeScript**: Type-safe JavaScript
- **Vite**: Build tool and dev server
- **Lucide React**: Icon library

## Troubleshooting

### Port Already in Use
```bash
# Find process using port 5000
netstat -ano | findstr :5000

# Kill the process (replace PID)
taskkill /PID <PID> /F
```

### CORS Errors
Ensure Flask backend has CORS enabled:
```python
from flask_cors import CORS
CORS(app)
```

### File Upload Fails
- Check file format (must be TSV)
- Verify file size limits in Flask config
- Ensure files contain required columns

### Plot Not Generating
- Verify at least one E25 and one E100 file uploaded
- Check browser console for error messages
- Ensure files contain protein intensity data

## Documentation

- **[ARCHITECTURE.md](./ARCHITECTURE.md)**: Technical architecture and design patterns
- **[Backend README](./backend/README.md)**: Flask API documentation (if exists)
- **[Frontend README](./frontend/README.md)**: React component documentation

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md) in the root repository.

## License

See [LICENSE](../../LICENSE) in the root repository.

## Features

- ✨ Modern, dark-themed UI
- 📤 Drag-and-drop file upload
- 📊 Two visualization modes:
  - Protein ID bar chart
  - Sample Intensity Comparison (E25 vs E100)
- 🖥️ Native desktop application
- 🚀 Can be deployed as web server (Flask standalone)
- 🔍 Automatic mix detection and grouping
- 🎨 Color-coded E25/E100 comparison plots

## API Endpoints

- `GET /api/health` - Health check
- `POST /api/upload` - Upload TSV files
- `GET /api/files` - List uploaded files
- `DELETE /api/files` - Clear all files
- `POST /api/plot/bar-chart` - Generate protein ID bar chart
- `POST /api/plot/sample-comparison` - Generate E25 vs E100 intensity comparison

See [API_CHANGES.md](./API_CHANGES.md) for detailed migration information.

## Tech Stack

**Frontend:**
- React 18
- TypeScript 5
- Vite (build tool)
- Lucide React (icons)
- Modern CSS

**Backend:**
- Flask (REST API)
- Pandas (data processing)
- Matplotlib (plotting)
- NumPy (numerical operations)

**Desktop:**
- PyWebView (native window)
- Threading (Flask background server)

## Converting to Web App

To deploy as a standalone web server instead of desktop app:

1. Remove PyWebView dependency
2. Run Flask in production mode (e.g., with Gunicorn)
3. Serve React build from Flask static folder
4. Deploy to server (Heroku, DigitalOcean, AWS, etc.)

## File Structure

```
mspp_web/
├── backend/
│   ├── app.py              # Flask API
│   └── requirements.txt
├── frontend/
│   ├── src/
│   │   ├── App.tsx         # Main React component
│   │   ├── App.css         # Styles
│   │   ├── main.tsx        # Entry point
│   │   └── index.css       # Global styles
│   ├── index.html
│   ├── package.json
│   ├── tsconfig.json
│   └── vite.config.ts
├── desktop_app.py          # Desktop launcher
└── README.md
```

## Expected File Format

Files should follow the naming convention:
- `report.pg_matrix_E25_<mix_params>.tsv`
- `report.pg_matrix_E100_<mix_params>.tsv`

Example:
- `report.pg_matrix_E25_30_4_440960_600.tsv`
- `report.pg_matrix_E100_30_4_440960_600.tsv`

The application automatically:
- Detects E25/E100 patterns (flexible: E25, E-25, E_25, etc.)
- Groups samples by mix identifier
- Separates different mixes with visual boundaries

## Next Steps

- Add plot export (PNG/SVG download)
- Add data table view
- Implement plot zoom/pan
- Add loading progress indicators
- Create installer (PyInstaller + NSIS)
- Add statistics summary panel
