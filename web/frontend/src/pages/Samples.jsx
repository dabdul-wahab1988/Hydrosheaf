import { useState, useEffect } from 'react'
import './Samples.css'

const API_BASE = 'http://localhost:8000/api'

function Samples() {
    const [datasets, setDatasets] = useState([])
    const [selectedDataset, setSelectedDataset] = useState(null)
    const [loading, setLoading] = useState(true)
    const [uploading, setUploading] = useState(false)
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
            }
        } catch (error) {
            console.error('Failed to fetch dataset:', error)
        }
    }

    const handleFileUpload = async (e) => {
        const file = e.target.files[0]
        if (!file) return

        // Validate file type
        const validExtensions = ['.csv', '.json']
        const fileExt = file.name.toLowerCase().slice(file.name.lastIndexOf('.'))
        if (!validExtensions.includes(fileExt)) {
            alert('Please upload a CSV or JSON file')
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
                alert(`${result.file_type} file uploaded successfully! ${result.sample_count} samples loaded.`)
            } else {
                const error = await res.json()
                alert(`Upload failed: ${error.detail}`)
            }
        } catch (error) {
            console.error('Upload error:', error)
            alert('Upload failed. Please try again.')
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
