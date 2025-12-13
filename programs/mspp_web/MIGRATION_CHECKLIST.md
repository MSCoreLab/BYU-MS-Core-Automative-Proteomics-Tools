# MSPP Web App Migration Checklist

## ✅ Completed Changes

### Backend (app.py)
- [x] Removed `Line2D` import
- [x] Removed `_calculate_protein_fold_changes()` method
- [x] Removed `_calculate_organism_vs_hela()` method  
- [x] Removed `generate_fold_change_plot()` method
- [x] Removed `generate_organisms_vs_hela()` method
- [x] Removed `generate_grouped_fold_change()` method
- [x] Added `_extract_mix_identifier()` method
- [x] Added `_calculate_sample_intensities()` method
- [x] Added `generate_sample_comparison()` method
- [x] Updated endpoint: `/api/plot/fold-change` → `/api/plot/sample-comparison`
- [x] Removed endpoint: `/api/plot/organisms-vs-hela`
- [x] Removed endpoint: `/api/plot/grouped-fold-change`

### Frontend (App.tsx)
- [x] Removed unused imports: `Activity`, `GitBranch`
- [x] Removed `PlotResult` interface properties: `unmatched_count`, `group_count`
- [x] Removed `groupPattern` state variable
- [x] Updated button: "E.coli vs Yeast Fold Change" → "Sample Intensity Comparison (E25 vs E100)"
- [x] Updated API call: `fold-change` → `sample-comparison`
- [x] Removed "Organisms vs HeLa" button
- [x] Removed "Grouped Fold Change" section
- [x] Removed pattern input UI
- [x] Simplified plot result handling (removed group count logic)

### Styling (App.css)
- [x] Removed `.group-pattern-section` styles
- [x] Removed `.pattern-input` styles
- [x] Removed `.pattern-field` styles
- [x] Removed `.pattern-hint` styles
- [x] Simplified `.plot-buttons` styles

### Documentation
- [x] Created `API_CHANGES.md` with migration guide
- [x] Updated `README.md` features section
- [x] Updated `README.md` API endpoints section
- [x] Added expected file format documentation
- [x] Updated Next Steps roadmap

## 🧪 Testing Checklist

### Backend Testing
- [ ] Test health check: `GET /api/health`
- [ ] Test file upload with E25/E100 files
- [ ] Test bar chart generation
- [ ] Test sample comparison with valid E25/E100 pairs
- [ ] Test error handling for missing E25 or E100 files
- [ ] Test mix identifier extraction with various formats
- [ ] Verify removed endpoints return 404

### Frontend Testing
- [ ] Upload multiple TSV files
- [ ] Generate Protein ID bar chart
- [ ] Generate Sample Intensity Comparison plot
- [ ] Verify error messages display correctly
- [ ] Test file clearing functionality
- [ ] Verify layout with no uploaded files
- [ ] Test responsive design at different screen sizes

### Integration Testing
- [ ] Test complete workflow: upload → bar chart → sample comparison
- [ ] Test with actual E25/E100 data files
- [ ] Verify plot displays correctly with base64 image
- [ ] Test with multiple mixes (e.g., Mix A: _600, Mix B: _800)
- [ ] Verify mix boundaries and labels appear correctly
- [ ] Test color coding (E25: orange/teal, E100: blue/purple)

## 🚀 Deployment Steps

1. **Install Dependencies**
   ```bash
   cd backend
   pip install -r requirements.txt
   
   cd ../frontend
   npm install
   ```

2. **Build Frontend**
   ```bash
   cd frontend
   npm run build
   ```

3. **Test Backend**
   ```bash
   cd backend
   python app.py
   # Verify at http://localhost:5000/api/health
   ```

4. **Test Frontend Dev Mode**
   ```bash
   cd frontend
   npm run dev
   # Verify at http://localhost:3000
   ```

5. **Test Desktop App**
   ```bash
   python desktop_app.py
   ```

## 📝 Known Changes

### Breaking Changes
- Old `/api/plot/fold-change` endpoint no longer exists
- Old `/api/plot/organisms-vs-hela` endpoint no longer exists  
- Old `/api/plot/grouped-fold-change` endpoint no longer exists
- Pattern matching UI removed from frontend

### New Features
- Automatic E25/E100 detection from filenames
- Mix-based grouping and visual bracketing
- Side-by-side E.coli and Yeast comparison plots
- Color-coded box plots for easy comparison
- Flexible pattern matching (E25, E-25, E_25, etc.)

### Feature Parity
✅ Both GUI and Web versions now have identical functionality:
- Protein ID bar chart
- E25 vs E100 intensity comparison with mix grouping
