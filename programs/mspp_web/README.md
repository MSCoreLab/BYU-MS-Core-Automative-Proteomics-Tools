# MSPP Data Plotter - React + TypeScript + Python

Modern desktop application for MS proteomics data visualization with a React/TypeScript frontend and Python backend.

## Architecture

```
┌─────────────────────────────────────────┐
│   Desktop App (PyWebView)              │
│  ┌───────────────────────────────────┐ │
│  │  React + TypeScript Frontend      │ │
│  │  (Vite, Modern UI)                │ │
│  └───────────────┬───────────────────┘ │
│                  │ HTTP/REST API        │
│  ┌───────────────▼───────────────────┐ │
│  │  Flask Backend                    │ │
│  │  (Pandas, Matplotlib, NumPy)      │ │
│  └───────────────────────────────────┘ │
└─────────────────────────────────────────┘
```

## Setup

### 1. Backend (Python)

```powershell
cd backend
pip install -r requirements.txt
```

### 2. Frontend (React + TypeScript)

```powershell
cd frontend
npm install
```

## Development

### Run in Development Mode

**Terminal 1 - Backend:**
```powershell
cd backend
python app.py
```

**Terminal 2 - Frontend:**
```powershell
cd frontend
npm run dev
```

**Terminal 3 - Desktop App:**
```powershell
python desktop_app.py
```

The desktop app will open a native window showing the React UI at http://localhost:3000, which proxies API calls to the Flask backend at http://localhost:5000.

## Production Build

### 1. Build Frontend
```powershell
cd frontend
npm run build
```

### 2. Update desktop_app.py
Change the webview URL from `http://localhost:3000` to:
```python
url = str(Path(__file__).parent / 'frontend' / 'dist' / 'index.html')
```

### 3. Run Desktop App
```powershell
python desktop_app.py
```

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
