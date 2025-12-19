# MSPP Web Application Architecture

## Overview
The MS Protein & Peptide (MSPP) Data Plotter is a full-stack web application for analyzing and visualizing mass spectrometry proteomics data from Orbitrap Astral MS (DIA-NN output).

## Technology Stack

### Backend
- **Framework**: Flask (Python 3.14+)
- **Data Processing**: pandas, numpy
- **Visualization**: matplotlib (Agg backend)
- **CORS**: flask-cors for cross-origin requests

### Frontend
- **Framework**: React 18 + TypeScript
- **Build Tool**: Vite
- **Icons**: lucide-react
- **Styling**: CSS3 with CSS variables

## Architecture Pattern

### Backend (Flask REST API)

#### Class Structure
```
DataProcessor
├── Handles all data loading and processing
├── Caching mechanism for performance
├── Organism identification (vectorized)
└── Intensity ratio calculations

PlotGenerator
├── Creates matplotlib visualizations
├── _create_bar_chart_figure() - Reusable bar chart generator
├── _create_comparison_figure() - Reusable comparison plot generator
├── create_bar_chart() - Returns base64 for web display
├── create_sample_comparison_plot() - Returns base64 for web display
└── _fig_to_base64() - Converts matplotlib figures to base64
```

#### Key Design Decisions
1. **Separation of Concerns**: Data processing logic separated from plotting logic
2. **DRY Principle**: Reusable `_create_*_figure()` methods eliminate code duplication
3. **Caching**: File data cached to avoid redundant reads
4. **Vectorization**: Pandas vectorized operations for performance
5. **Consensus Proteins**: Only proteins with valid intensities in both E25 and E100
6. **Flexible Sizing**: Figure size parameterized for different use cases (display vs. export)

#### Refactoring (December 2024)
- **Eliminated 200+ lines of duplicated code** by extracting reusable plot generation methods
- **Reduced backend from 857 to 714 lines** (~17% reduction)
- **Created `_create_bar_chart_figure()`**: Shared by display, individual export, and ZIP export
- **Created `_create_comparison_figure()`**: Shared by display, individual export, and ZIP export
- **Benefits**: Single source of truth, easier maintenance, consistent visualizations

#### API Endpoints

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/api/upload` | POST | Upload TSV files |
| `/api/files` | GET | List uploaded files |
| `/api/files` | DELETE | Clear all files |
| `/api/plot/bar-chart` | POST | Generate protein count bar chart (base64) |
| `/api/plot/sample-comparison` | POST | Generate intensity ratio plots (base64) |
| `/api/export/bar-chart` | POST | Export bar chart as PNG (300 DPI) |
| `/api/export/sample-comparison` | POST | Export intensity plots as PNG (300 DPI) |
| `/api/export/all` | POST | Export all plots as ZIP |

### Frontend (React + TypeScript)

#### Component Structure
```
App.tsx (Main Component)
├── File Upload Section
│   ├── File browser
│   ├── Upload button
│   └── File list display
├── Plot Controls Section
│   ├── Generate Bar Chart button
│   ├── Generate Comparison button
│   └── Export All Plots button
└── Plot Display Section
    ├── Loading state
    ├── Plot image display
    ├── Export Plot button (per plot)
    └── Empty state
```

#### State Management
- `files`: Selected files from file input
- `uploadedFiles`: Files successfully uploaded to backend
- `loading`: Loading state for API calls
- `plotImage`: Base64 image data from backend
- `currentPlotType`: Tracks which plot is displayed
- `error`: Error messages

#### Key Features
1. **File Management**: Multi-file upload with TSV/TXT filtering
2. **Real-time Feedback**: Loading indicators and error messages
3. **Export Options**: Individual plot export or ZIP of all plots
4. **Responsive Design**: Grid layout with mobile support

## Data Flow

### Upload Flow
```
User selects files → Frontend validates → POST /api/upload
→ Backend saves to temp → Returns file list → Frontend updates UI
```

### Plot Generation Flow
```
User clicks plot button → POST /api/plot/{type}
→ Backend loads data → Processes → Generates plot
→ Returns base64 PNG → Frontend displays image
```

### Export Flow
```
User clicks export → POST /api/export/{type}
→ Backend generates high-res plot → Returns binary PNG/ZIP
→ Frontend triggers download
```

## Algorithm: Intensity Ratio Calculation

### Purpose
Calculate log2(E25/E100) intensity ratios for quality control assessment.

### Steps
1. **Identify Samples**: Parse filenames to detect E25 and E100 pairs
2. **Find Consensus Proteins**: Proteins with valid intensities in BOTH samples
3. **Calculate Ratios**: `log2(E25_intensity / E100_intensity)` per protein
4. **Filter**: Remove inf/NaN values
5. **Return**: Array of ratios (one per protein)

### Expected Ratios
- **HeLa**: ~0 (constant concentration)
- **E.coli**: ~-2 (log2(25/100) = -2, 4-fold dilution)
- **Yeast**: ~+1 (log2(150/75) = 1, 2-fold concentration)

### Validation
- HeLa median within ±0.5 of 0
- E.coli median between -2.5 and -1.5
- Yeast median between 0.5 and 1.5

## Jupyter Notebook vs Web App

### Feature Parity
Both tools provide:
- ✓ File upload/loading
- ✓ Organism identification
- ✓ Protein count bar charts
- ✓ Intensity ratio box plots (3-panel)
- ✓ Statistical summaries
- ✓ Validation checks
- ✓ High-quality plot export

### Differences

| Feature | Notebook | Web App |
|---------|----------|---------|
| Interface | Jupyter cells | React UI |
| File Loading | tkinter dialog | Browser upload |
| Plot Display | Inline matplotlib | Base64 images |
| Export | Manual save | Download buttons |
| Portability | Requires Python | Browser-based |
| Use Case | Research/Development | Production/Routine |

## Performance Optimizations

### Backend
1. **Caching**: Loaded data cached to avoid re-reading files
2. **Vectorization**: Pandas operations instead of row-by-row processing
3. **Lazy Loading**: Data only loaded when needed
4. **Memory Management**: Figures closed after conversion to avoid leaks

### Frontend
1. **Build Optimization**: Vite tree-shaking and minification
2. **Code Splitting**: React lazy loading (future enhancement)
3. **Asset Optimization**: CSS minification

## Security Considerations

1. **File Upload**: Only TSV/TXT files accepted
2. **Temp Files**: Stored in system temp directory
3. **CORS**: Enabled for development, should be restricted in production
4. **Input Validation**: File existence and format checks

## Future Enhancements

### Planned Features
- [ ] Multiple file comparison (>2 replicates)
- [ ] Custom organism patterns
- [ ] CSV export of statistics
- [ ] Batch processing mode
- [ ] User authentication
- [ ] Database storage for results

### Technical Debt
- [ ] Add comprehensive unit tests
- [ ] Implement logging configuration
- [ ] Add API rate limiting
- [ ] Improve error handling granularity
- [ ] Add TypeScript strict mode
- [ ] Implement proper state management (Redux/Zustand)

## Development Workflow

### Backend Development
```bash
cd programs/mspp_web/backend
python app.py  # Runs on http://localhost:5000
```

### Frontend Development
```bash
cd programs/mspp_web/frontend
npm run dev    # Runs on http://localhost:5173
```

### Production Build
```bash
cd programs/mspp_web/frontend
npm run build  # Outputs to dist/
# Backend serves from dist/ folder
```

## Testing

### Manual Testing Checklist
- [ ] Upload multiple TSV files
- [ ] Generate bar chart
- [ ] Generate intensity comparison
- [ ] Export individual plots
- [ ] Export all plots as ZIP
- [ ] Clear files and re-upload
- [ ] Test with invalid files
- [ ] Test with missing E25 or E100

### Example Test Files
Located in: `tests/Sample Files/Hey Astral/`
- `report.pg_matrix_E25_8_4_20_0.tsv`
- `report.pg_matrix_E100_8_4_20_0.tsv`

## Maintenance

### Common Issues

**Issue**: Plots not displaying
**Solution**: Rebuild frontend with `npm run build`

**Issue**: FileNotFoundError in notebook
**Solution**: Use file upload dialog instead of hardcoded paths

**Issue**: No data in plots
**Solution**: Verify E25/E100 files contain matching proteins

**Issue**: Memory errors
**Solution**: Clear cache with `/api/files DELETE` endpoint

## Contributing

See `CONTRIBUTING.md` for guidelines on:
- Code style (PEP 8 for Python, ESLint for TypeScript)
- Commit message format
- Pull request process
- Testing requirements
