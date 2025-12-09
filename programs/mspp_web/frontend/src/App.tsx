import { useState } from 'react'
import { Upload, BarChart3, TrendingUp, Activity, X, Loader2, GitBranch } from 'lucide-react'
import './App.css'

interface PlotResult {
  image: string
  error?: string
  unmatched_count?: number
  group_count?: number
}

function App() {
  const [files, setFiles] = useState<File[]>([])
  const [uploadedFiles, setUploadedFiles] = useState<string[]>([])
  const [loading, setLoading] = useState(false)
  const [plotImage, setPlotImage] = useState<string | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [groupPattern, setGroupPattern] = useState<string>(String.raw`(E\d+)`)

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files) {
      setFiles(Array.from(e.target.files))
    }
  }

  const uploadFiles = async () => {
    if (files.length === 0) {
      setError('Please select TSV files first')
      return
    }

    setLoading(true)
    setError(null)

    const formData = new FormData()
    files.forEach(file => formData.append('files', file))

    try {
      const response = await fetch('/api/upload', {
        method: 'POST',
        body: formData,
      })

      if (!response.ok) throw new Error('Upload failed')

      const data = await response.json()
      setUploadedFiles(data.files)
      setFiles([])
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Upload failed')
    } finally {
      setLoading(false)
    }
  }

  const clearFiles = async () => {
    try {
      await fetch('/api/files', { method: 'DELETE' })
      setUploadedFiles([])
      setPlotImage(null)
    } catch (err) {
      setError('Failed to clear files')
    }
  }

  const generatePlot = async (endpoint: string, requestBody?: object) => {
    if (uploadedFiles.length === 0) {
      setError('Please upload files first')
      return
    }

    setLoading(true)
    setError(null)
    setPlotImage(null)

    try {
      const response = await fetch(`/api/plot/${endpoint}`, {
        method: 'POST',
        headers: requestBody ? { 'Content-Type': 'application/json' } : undefined,
        body: requestBody ? JSON.stringify(requestBody) : undefined,
      })

      const data: PlotResult = await response.json()

      if (!response.ok) throw new Error(data.error || 'Plot generation failed')

      setPlotImage(data.image)
      
      // Show info about grouped results
      if (data.group_count !== undefined) {
        const msg = `Generated plot with ${data.group_count} group(s)`
        if (data.unmatched_count && data.unmatched_count > 0) {
          setError(`${msg} (${data.unmatched_count} files did not match pattern)`)
        }
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Plot generation failed')
    } finally {
      setLoading(false)
    }
  }

  return (
    <div className="app">
      <header className="header">
        <h1>MS Protein & Peptide Data Plotter</h1>
        <p>Modern interface for mass spectrometry data visualization</p>
      </header>

      <main className="main">
        <div className="controls-section">
          <div className="file-upload-zone">
            <div className="upload-area">
              <Upload size={32} />
              <input
                type="file"
                multiple
                accept=".tsv,.txt"
                onChange={handleFileChange}
                className="file-input"
                id="file-input"
              />
              <label htmlFor="file-input" className="upload-label">
                {files.length > 0
                  ? `${files.length} file(s) selected`
                  : 'Click to select TSV files'}
              </label>
            </div>

            <div className="upload-actions">
              <button
                onClick={uploadFiles}
                disabled={files.length === 0 || loading}
                className="btn btn-primary"
              >
                {loading ? (
                  <>
                    <Loader2 className="spin" size={16} />
                    Uploading...
                  </>
                ) : (
                  <>
                    <Upload size={16} />
                    Upload Files
                  </>
                )}
              </button>

              {uploadedFiles.length > 0 && (
                <button onClick={clearFiles} className="btn btn-secondary">
                  <X size={16} />
                  Clear All
                </button>
              )}
            </div>

            {uploadedFiles.length > 0 && (
              <div className="uploaded-files">
                <h3>Uploaded Files ({uploadedFiles.length})</h3>
                <ul>
                  {uploadedFiles.map((file, idx) => (
                    <li key={idx}>{file}</li>
                  ))}
                </ul>
              </div>
            )}
          </div>

          <div className="plot-controls">
            <h3>Generate Plots</h3>
            <div className="plot-buttons">
              <button
                onClick={() => generatePlot('bar-chart')}
                disabled={uploadedFiles.length === 0 || loading}
                className="btn btn-plot"
              >
                <BarChart3 size={20} />
                Protein ID Bar Chart
              </button>

              <button
                onClick={() => generatePlot('fold-change')}
                disabled={uploadedFiles.length === 0 || loading}
                className="btn btn-plot"
              >
                <TrendingUp size={20} />
                E.coli vs Yeast Fold Change
              </button>

              <button
                onClick={() => generatePlot('organisms-vs-hela')}
                disabled={uploadedFiles.length === 0 || loading}
                className="btn btn-plot"
              >
                <Activity size={20} />
                Organisms vs HeLa
              </button>
            </div>

            <div className="group-pattern-section">
              <h4>Group Files by Pattern</h4>
              <div className="pattern-input">
                <label htmlFor="pattern">Regex Pattern:</label>
                <input
                  id="pattern"
                  type="text"
                  value={groupPattern}
                  onChange={(e) => setGroupPattern(e.target.value)}
                  placeholder="(E\d+)"
                  className="pattern-field"
                />
              </div>
              <p className="pattern-hint">
                Default: (E\d+) groups by E25, E100, etc.
              </p>
              <button
                onClick={() => generatePlot('grouped-fold-change', { pattern: groupPattern })}
                disabled={uploadedFiles.length === 0 || loading || !groupPattern}
                className="btn btn-plot"
              >
                <GitBranch size={20} />
                Grouped Fold Change
              </button>
            </div>
          </div>

          {error && (
            <div className="error-message">
              <strong>Error:</strong> {error}
            </div>
          )}
        </div>

        <div className="plot-display">
          {loading && !plotImage && (
            <div className="loading">
              <Loader2 className="spin" size={48} />
              <p>Generating plot...</p>
            </div>
          )}

          {plotImage && (
            <div className="plot-container">
              <img
                src={`data:image/png;base64,${plotImage}`}
                alt="Generated plot"
                className="plot-image"
              />
            </div>
          )}

          {!plotImage && !loading && (
            <div className="empty-state">
              <BarChart3 size={64} opacity={0.3} />
              <p>Upload TSV files and select a plot type to visualize your data</p>
            </div>
          )}
        </div>
      </main>
    </div>
  )
}

export default App
