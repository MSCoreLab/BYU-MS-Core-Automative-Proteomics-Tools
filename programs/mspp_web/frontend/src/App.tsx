import { useState } from 'react'
import { Upload, BarChart3, TrendingUp, X, Loader2, Download, Activity } from 'lucide-react'
import './App.css'

interface PlotResult {
  image: string
  error?: string
}

function App() {
  // TSV file state
  const [files, setFiles] = useState<File[]>([])
  const [uploadedFiles, setUploadedFiles] = useState<string[]>([])
  
  // mzML file state
  const [mzmlFiles, setMzmlFiles] = useState<File[]>([])
  const [uploadedMzmlFiles, setUploadedMzmlFiles] = useState<string[]>([])
  
  // UI state
  const [loading, setLoading] = useState(false)
  const [plotImage, setPlotImage] = useState<string | null>(null)
  const [currentPlotType, setCurrentPlotType] = useState<string | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [activeTab, setActiveTab] = useState<'protein' | 'tic'>('protein')

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files) {
      setFiles(Array.from(e.target.files))
    }
  }

  const handleMzmlFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files) {
      setMzmlFiles(Array.from(e.target.files))
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

  const uploadMzmlFiles = async () => {
    if (mzmlFiles.length === 0) {
      setError('Please select mzML files first')
      return
    }

    setLoading(true)
    setError(null)

    const formData = new FormData()
    mzmlFiles.forEach(file => formData.append('files', file))

    try {
      const response = await fetch('/api/upload/mzml', {
        method: 'POST',
        body: formData,
      })

      if (!response.ok) throw new Error('mzML upload failed')

      const data = await response.json()
      setUploadedMzmlFiles(data.files)
      setMzmlFiles([])
    } catch (err) {
      setError(err instanceof Error ? err.message : 'mzML upload failed')
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

  const clearMzmlFiles = async () => {
    try {
      await fetch('/api/mzml/files', { method: 'DELETE' })
      setUploadedMzmlFiles([])
      setPlotImage(null)
    } catch (err) {
      setError('Failed to clear mzML files')
    }
  }

  const generatePlot = async (endpoint: string, requestBody?: object) => {
    // Check appropriate files based on endpoint
    if (endpoint === 'tic') {
      if (uploadedMzmlFiles.length === 0) {
        setError('Please upload mzML files first')
        return
      }
    } else {
      if (uploadedFiles.length === 0) {
        setError('Please upload TSV files first')
        return
      }
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
      setCurrentPlotType(endpoint)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Plot generation failed')
    } finally {
      setLoading(false)
    }
  }

  const exportPlot = async () => {
    if (!currentPlotType) return

    try {
      const response = await fetch(`/api/export/${currentPlotType}`, {
        method: 'POST',
      })

      if (!response.ok) throw new Error('Export failed')

      const blob = await response.blob()
      const url = window.URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      
      // Determine filename based on plot type
      let filename = 'plot.png'
      if (currentPlotType === 'bar-chart') {
        filename = 'protein_id_bar_chart.png'
      } else if (currentPlotType === 'sample-comparison') {
        filename = 'intensity_ratio_comparison.png'
      } else if (currentPlotType === 'tic') {
        filename = 'tic_comparison.png'
      }
      
      a.download = filename
      document.body.appendChild(a)
      a.click()
      window.URL.revokeObjectURL(url)
      document.body.removeChild(a)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Export failed')
    }
  }

  const exportAllPlots = async () => {
    try {
      const response = await fetch('/api/export/all', {
        method: 'POST',
      })

      if (!response.ok) throw new Error('Export failed')

      const blob = await response.blob()
      const url = window.URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = 'mspp_plots.zip'
      document.body.appendChild(a)
      a.click()
      window.URL.revokeObjectURL(url)
      document.body.removeChild(a)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Export all plots failed')
    }
  }

  return (
    <div className="app">
      <header className="header">
        <h1>MS Protein & Peptide Data Plotter</h1>
        <p>Modern interface for mass spectrometry data visualization</p>
      </header>

      <main className="main">
        {/* Tab Navigation */}
        <div className="tab-navigation">
          <button
            className={`tab-button ${activeTab === 'protein' ? 'active' : ''}`}
            onClick={() => setActiveTab('protein')}
          >
            <BarChart3 size={20} />
            Protein Analysis
          </button>
          <button
            className={`tab-button ${activeTab === 'tic' ? 'active' : ''}`}
            onClick={() => setActiveTab('tic')}
          >
            <Activity size={20} />
            TIC Analysis
          </button>
        </div>

        {/* Protein Analysis Tab */}
        {activeTab === 'protein' && (
          <>
            <div className="main-content">
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
                onClick={() => generatePlot('sample-comparison')}
                disabled={uploadedFiles.length === 0 || loading}
                className="btn btn-plot"
              >
                <TrendingUp size={20} />
                Sample Intensity Comparison (E25 vs E100)
              </button>

              <button
                onClick={exportAllPlots}
                disabled={uploadedFiles.length === 0 || loading}
                className="btn btn-export-all"
                title="Download all plots as ZIP"
              >
                <Download size={20} />
                Export All Plots (ZIP)
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
                    <div className="plot-header">
                      <button
                        onClick={exportPlot}
                        className="btn btn-export"
                        title="Download plot as PNG"
                      >
                        <Download size={16} />
                        Export Plot
                      </button>
                    </div>
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
            </div>
          </>
        )}

        {/* TIC Analysis Tab */}
        {activeTab === 'tic' && (
          <>
            <div className="main-content">
              <div className="controls-section">
            <div className="file-upload-zone">
              <div className="upload-area">
                <Upload size={32} />
                <input
                  type="file"
                  multiple
                  accept=".mzML"
                  onChange={handleMzmlFileChange}
                  className="file-input"
                  id="mzml-file-input"
                />
                <label htmlFor="mzml-file-input" className="upload-label">
                  {mzmlFiles.length > 0
                    ? `${mzmlFiles.length} mzML file(s) selected`
                    : 'Click to select mzML files'}
                </label>
              </div>

              <div className="upload-actions">
                <button
                  onClick={uploadMzmlFiles}
                  disabled={mzmlFiles.length === 0 || loading}
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
                      Upload mzML Files
                    </>
                  )}
                </button>

                {uploadedMzmlFiles.length > 0 && (
                  <button onClick={clearMzmlFiles} className="btn btn-secondary">
                    <X size={16} />
                    Clear All
                  </button>
                )}
              </div>

              {uploadedMzmlFiles.length > 0 && (
                <div className="uploaded-files">
                  <h3>Uploaded mzML Files ({uploadedMzmlFiles.length})</h3>
                  <ul>
                    {uploadedMzmlFiles.map((file, idx) => (
                      <li key={idx}>{file}</li>
                    ))}
                  </ul>
                </div>
              )}
            </div>

            <div className="plot-controls">
              <h3>Generate TIC Plot</h3>
              <div className="plot-buttons">
                <button
                  onClick={() => generatePlot('tic')}
                  disabled={uploadedMzmlFiles.length === 0 || loading}
                  className="btn btn-plot"
                >
                  <Activity size={20} />
                  Total Ion Current Comparison
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
                    <p>Processing mzML files and generating TIC plot...</p>
                  </div>
                )}

                {plotImage && (
                  <div className="plot-container">
                    <div className="plot-header">
                      <button
                        onClick={exportPlot}
                        className="btn btn-export"
                        title="Download plot as PNG"
                      >
                        <Download size={16} />
                        Export TIC Plot
                      </button>
                    </div>
                    <img
                      src={`data:image/png;base64,${plotImage}`}
                      alt="TIC Comparison Plot"
                      className="plot-image"
                    />
                  </div>
                )}

                {!plotImage && !loading && (
                  <div className="empty-state">
                    <Activity size={64} opacity={0.3} />
                    <p>Upload mzML files to generate Total Ion Current comparison plots</p>
                  </div>
                )}
              </div>
            </div>
          </>
        )}
      </main>
    </div>
  )
}

export default App
