# MSPP Web API Changes

## Summary
The backend API has been updated to match the new GUI functionality. The focus has shifted from E.coli vs Yeast fold change analysis to E25 vs E100 intensity comparison within experimental mixes.

## Breaking Changes

### Removed Endpoints
1. **`POST /api/plot/fold-change`** - Removed E.coli vs Yeast fold change analysis
2. **`POST /api/plot/organisms-vs-hela`** - Removed spike-in validation plots
3. **`POST /api/plot/grouped-fold-change`** - Removed grouped fold change by regex pattern

### New Endpoints
1. **`POST /api/plot/sample-comparison`** - Generate E25 vs E100 intensity comparison
   - **Request**: No body required
   - **Response**: 
     ```json
     {
       "image": "base64_encoded_png_string"
     }
     ```
   - **Error Response** (400):
     ```json
     {
       "error": "No valid E25/E100 pairs found"
     }
     ```

## Unchanged Endpoints
- `GET /api/health` - Health check
- `POST /api/upload` - Upload TSV files
- `GET /api/files` - List uploaded files
- `DELETE /api/files` - Clear all files
- `POST /api/plot/bar-chart` - Generate protein ID bar chart

## New Features
The sample comparison plot now:
- Automatically groups samples by mix identifier (extracted from filename)
- Identifies E25 and E100 samples using flexible pattern matching (E25, E-25, E_25, etc.)
- Shows side-by-side box plots for E.coli and Yeast
- Color codes E25 (orange/teal) vs E100 (blue/purple)
- Adds visual separators and labels between different mixes
- Displays HeLa-normalized log2 intensities

## Frontend Updates Needed

### Update Component Names
- Rename "Fold Change" button/component to "Sample Intensity Comparison"
- Remove "Organisms vs HeLa" button/component
- Remove "Grouped Fold Change" feature and pattern input field

### Update API Calls
Replace:
```javascript
// Old
fetch('/api/plot/fold-change', { method: 'POST' })

// New
fetch('/api/plot/sample-comparison', { method: 'POST' })
```

Remove:
```javascript
fetch('/api/plot/organisms-vs-hela', { method: 'POST' })
fetch('/api/plot/grouped-fold-change', { 
  method: 'POST',
  body: JSON.stringify({ pattern: '...' })
})
```

### Expected File Format
The backend expects files with E25/E100 in the filename:
- `report.pg_matrix_E25_30_4_440960_600.tsv`
- `report.pg_matrix_E100_30_4_440960_600.tsv`
- Mix ID extracted: `30_4_440960_600`

Pattern matching is flexible and case-insensitive:
- ✅ E25, E-25, E_25, e25
- ✅ E100, E-100, E_100, e100

## Testing
Test with sample files that include:
1. Both E25 and E100 for each mix
2. Multiple mixes (different experimental parameters)
3. Protein data with HeLa, E.coli, and Yeast organisms

## Migration Guide
1. Update frontend routing to remove old endpoints
2. Update button labels and component names
3. Test with sample data files
4. Update user documentation
