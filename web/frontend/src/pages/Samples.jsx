import { useState, useEffect } from 'react'
import { API_BASE } from '../config'
import { useToast } from '../context/ToastContext'
import './Samples.css'

function Samples() {
    const toast = useToast()
    const [datasets, setDatasets] = useState([])
    const [selectedDataset, setSelectedDataset] = useState(null)
    const [loading, setLoading] = useState(true)
    const [uploading, setUploading] = useState(false)
    const [validating, setValidating] = useState(false)
    const [validationResult, setValidationResult] = useState(null)
    const [formData, setFormData] = useState({
        name: '',
        samples: []
    })

    useEffect(() => {
        fetchDatasets()
    }, [])

    const fetchDatasets = async () => {
        try {
            const res = await fetch(`${API_BASE}/samples/datasets`)
            if (res.ok) {
                setDatasets(await res.json())
            }
        } catch (error) {
            console.error('Failed to fetch datasets:', error)
        } finally {
            setLoading(false)
        }
    }

    const fetchDataset = async (id) => {
        try {
            const res = await fetch(`${API_BASE}/samples/datasets/${id}`)
            if (res.ok) {
                setSelectedDataset(await res.json())
                setValidationResult(null) // Clear validation when switching datasets
            }
        } catch (error) {
            console.error('Failed to fetch dataset:', error)
        }
    }

    const validateDataset = async () => {
        if (!selectedDataset) return

        setValidating(true)
        try {
            const res = await fetch(`${API_BASE}/samples/datasets/${selectedDataset.id}/validate`)
            if (res.ok) {
                const result = await res.json()
                setValidationResult(result)
                if (result.pass_rate >= 0.9) {
                    toast.success(`Data quality check passed: ${(result.pass_rate * 100).toFixed(0)}% valid`)
                } else if (result.pass_rate >= 0.7) {
                    toast.warning(`Data quality issues found: ${(result.pass_rate * 100).toFixed(0)}% valid`)
                } else {
                    toast.error(`Significant data quality issues: ${(result.pass_rate * 100).toFixed(0)}% valid`)
                }
            } else {
                toast.error('Failed to validate dataset')
            }
        } catch (error) {
            console.error('Validation error:', error)
            toast.error('Failed to validate dataset')
        } finally {
            setValidating(false)
        }
    }

    const handleFileUpload = async (e) => {
        const file = e.target.files[0]
        if (!file) return

        // Validate file type
        const validExtensions = ['.csv', '.json']
        const fileExt = file.name.toLowerCase().slice(file.name.lastIndexOf('.'))
        if (!validExtensions.includes(fileExt)) {
            toast.warning('Please upload a CSV or JSON file')
            return
        }

        setUploading(true)
        const formData = new FormData()
        formData.append('file', file)

        try {
            const res = await fetch(`${API_BASE}/samples/upload-file`, {
                method: 'POST',
                body: formData
            })

            if (res.ok) {
                const result = await res.json()
                await fetchDatasets()
                toast.success(`${result.file_type} file uploaded successfully! ${result.sample_count} samples loaded.`)
            } else {
                const error = await res.json()
                toast.error(`Upload failed: ${error.detail}`)
            }
        } catch (error) {
            console.error('Upload error:', error)
            toast.error('Upload failed. Please try again.')
        } finally {
            setUploading(false)
            // Reset file input
            e.target.value = ''
        }
    }

    const deleteDataset = async (id) => {
        if (!confirm('Are you sure you want to delete this dataset?')) return

        try {
            const res = await fetch(`${API_BASE}/samples/datasets/${id}`, {
                method: 'DELETE'
            })
            if (res.ok) {
                fetchDatasets()
                if (selectedDataset?.id === id) {
                    setSelectedDataset(null)
                }
            }
        } catch (error) {
            console.error('Delete error:', error)
        }
    }

    return (
        <div className="samples-page">
            <header className="page-header">
                <div>
                    <h1 className="page-title">Water Samples</h1>
                    <p className="page-subtitle">
                        Manage and analyze groundwater chemistry data
                    </p>
                </div>
            </header>

            <div className="samples-layout">
                {/* Sidebar - Dataset List */}
                <aside className="datasets-panel card">
                    <div className="panel-header">
                        <h3>Datasets</h3>
                        <label className="btn btn-primary upload-btn">
                            <input
                                type="file"
                                accept=".csv,.json"
                                onChange={handleFileUpload}
                                disabled={uploading}
                                hidden
                            />
                            {uploading ? 'Uploading...' : 'Upload File'}
                        </label>
                    </div>

                    <div className="datasets-list">
                        {loading ? (
                            <div className="loading-state">
                                <div className="spinner"></div>
                                <p>Loading datasets...</p>
                            </div>
                        ) : datasets.length === 0 ? (
                            <div className="empty-state">
                                <span className="empty-icon">📂</span>
                                <p>No datasets yet</p>
                                <span className="empty-hint">Upload a CSV or JSON file to get started</span>
                            </div>
                        ) : (
                            datasets.map((ds) => (
                                <div
                                    key={ds.id}
                                    className={`dataset-item ${selectedDataset?.id === ds.id ? 'active' : ''}`}
                                    onClick={() => fetchDataset(ds.id)}
                                >
                                    <div className="dataset-info">
                                        <span className="dataset-name">{ds.name}</span>
                                        <span className="dataset-count">{ds.sample_count} samples</span>
                                    </div>
                                    <button
                                        className="delete-btn"
                                        onClick={(e) => {
                                            e.stopPropagation()
                                            deleteDataset(ds.id)
                                        }}
                                    >
                                        🗑️
                                    </button>
                                </div>
                            ))
                        )}
                    </div>
                </aside>

                {/* Main Content - Sample Details */}
                <main className="samples-content">
                    {selectedDataset ? (
                        <>
                            <div className="card">
                                <div className="card-header">
                                    <h3 className="card-title">{selectedDataset.name}</h3>
                                    <span className="badge badge-primary">{selectedDataset.sample_count} Samples</span>
                                </div>
                                {selectedDataset.description && (
                                    <p className="dataset-description">{selectedDataset.description}</p>
                                )}
                            </div>

                            <div className="card mt-lg">
                                <div className="card-header">
                                    <h3 className="card-title">Sample Data</h3>
                                </div>
                                <div className="table-container">
                                    <table className="samples-table">
                                        <thead>
                                            <tr>
                                                <th>ID</th>
                                                <th>Location</th>
                                                <th>Date</th>
                                                <th>Ca²⁺</th>
                                                <th>Mg²⁺</th>
                                                <th>Na⁺</th>
                                                <th>HCO₃⁻</th>
                                                <th>SO₄²⁻</th>
                                                <th>Cl⁻</th>
                                                <th>pH</th>
                                                <th>EC</th>
                                                <th>δ¹⁸O</th>
                                                <th>δ²H</th>
                                            </tr>
                                        </thead>
                                        <tbody>
                                            {selectedDataset.samples.map((sample, idx) => (
                                                <tr key={idx}>
                                                    <td>{sample.sample_id || '-'}</td>
                                                    <td>{sample.location_id || '-'}</td>
                                                    <td>{sample.date || '-'}</td>
                                                    <td>{sample.ca?.toFixed(1) || '-'}</td>
                                                    <td>{sample.mg?.toFixed(1) || '-'}</td>
                                                    <td>{sample.na?.toFixed(1) || '-'}</td>
                                                    <td>{sample.hco3?.toFixed(1) || '-'}</td>
                                                    <td>{sample.so4?.toFixed(1) || '-'}</td>
                                                    <td>{sample.cl?.toFixed(1) || '-'}</td>
                                                    <td>{sample.ph?.toFixed(2) || '-'}</td>
                                                    <td>{sample.ec || '-'}</td>
                                                    <td>{sample.d18o?.toFixed(2) || '-'}</td>
                                                    <td>{sample.d2h?.toFixed(2) || '-'}</td>
                                                </tr>
                                            ))}
                                        </tbody>
                                    </table>
                                </div>
                            </div>

                            {/* Data Quality Check */}
                            <div className="card mt-lg">
                                <div className="card-header">
                                    <h3 className="card-title">Data Quality Check</h3>
                                    <button
                                        className="btn btn-secondary"
                                        onClick={validateDataset}
                                        disabled={validating}
                                    >
                                        {validating ? 'Validating...' : 'Run QC Check'}
                                    </button>
                                </div>

                                {validationResult ? (
                                    <div className="validation-results">
                                        <div className="validation-summary">
                                            <div className={`pass-rate ${
                                                validationResult.summary.pass_rate >= 0.9 ? 'pass-rate-good' :
                                                validationResult.summary.pass_rate >= 0.7 ? 'pass-rate-warning' : 'pass-rate-error'
                                            }`}>
                                                <span className="rate-value">
                                                    {(validationResult.summary.pass_rate * 100).toFixed(0)}%
                                                </span>
                                                <span className="rate-label">Pass Rate</span>
                                            </div>
                                            <div className="validation-stats">
                                                <div className="stat-item">
                                                    <span className="stat-value">{validationResult.valid_samples}</span>
                                                    <span className="stat-label">Valid</span>
                                                </div>
                                                <div className="stat-item">
                                                    <span className="stat-value">{validationResult.invalid_samples}</span>
                                                    <span className="stat-label">Invalid</span>
                                                </div>
                                                <div className="stat-item">
                                                    <span className="stat-value">{validationResult.total_samples}</span>
                                                    <span className="stat-label">Total</span>
                                                </div>
                                            </div>
                                        </div>

                                        {/* Common Issues */}
                                        {Object.keys(validationResult.summary.common_issues || {}).length > 0 && (
                                            <div className="common-issues">
                                                <h4>Common Issues</h4>
                                                <div className="issues-list">
                                                    {Object.entries(validationResult.summary.common_issues).map(([issue, count]) => (
                                                        <div key={issue} className="issue-item">
                                                            <span className="issue-name">{issue.replace(/_/g, ' ')}</span>
                                                            <span className="issue-count">{count} samples</span>
                                                        </div>
                                                    ))}
                                                </div>
                                            </div>
                                        )}

                                        {/* Invalid Samples */}
                                        {validationResult.results.filter(r => !r.is_valid).length > 0 && (
                                            <div className="invalid-samples">
                                                <h4>Invalid Samples</h4>
                                                <div className="invalid-list">
                                                    {validationResult.results
                                                        .filter(r => !r.is_valid)
                                                        .slice(0, 10)
                                                        .map((r, idx) => (
                                                            <div key={idx} className="invalid-item">
                                                                <span className="sample-id">{r.sample_id}</span>
                                                                <span className="sample-flags">
                                                                    {r.flags.join(', ')}
                                                                </span>
                                                            </div>
                                                        ))}
                                                    {validationResult.results.filter(r => !r.is_valid).length > 10 && (
                                                        <div className="more-samples">
                                                            ...and {validationResult.results.filter(r => !r.is_valid).length - 10} more
                                                        </div>
                                                    )}
                                                </div>
                                            </div>
                                        )}
                                    </div>
                                ) : (
                                    <div className="validation-empty">
                                        <span className="empty-icon">🔍</span>
                                        <p>Run a quality check to validate your data</p>
                                        <span className="empty-hint">Checks charge balance, negative values, and required fields</span>
                                    </div>
                                )}
                            </div>
                        </>
                    ) : (
                        <div className="card empty-content">
                            <div className="empty-state large">
                                <span className="empty-icon">🧪</span>
                                <h3>Select a Dataset</h3>
                                <p>Choose a dataset from the sidebar to view sample details</p>
                            </div>
                        </div>
                    )}
                </main>
            </div>
        </div>
    )
}

export default Samples
