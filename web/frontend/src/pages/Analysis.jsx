import { useState, useEffect } from 'react'
import { Link } from 'react-router-dom'
import { useProject } from '../context/ProjectContext'
import './Analysis.css'

const API_BASE = 'http://localhost:8000/api'

function Analysis() {
    const { currentProject, saveResultToProject } = useProject()
    const [datasets, setDatasets] = useState([])
    const [formData, setFormData] = useState({
        name: '',
        analysis_type: 'full_pipeline',
        dataset_id: '',
        config: {
            lasso_penalty: 0.1,
            enable_phreeqc: false,
            enable_isotopes: true,
            enable_uncertainty: true,
            bootstrap_iterations: 1000
        }
    })
    const [running, setRunning] = useState(false)
    const [currentJob, setCurrentJob] = useState(null)
    const [progress, setProgress] = useState(0)

    useEffect(() => {
        fetchDatasets()
    }, [])

    useEffect(() => {
        let interval
        if (currentJob && currentJob.status !== 'completed' && currentJob.status !== 'failed') {
            interval = setInterval(checkJobStatus, 1000)
        }
        return () => clearInterval(interval)
    }, [currentJob])

    const fetchDatasets = async () => {
        try {
            const res = await fetch(`${API_BASE}/samples/datasets`)
            if (res.ok) {
                setDatasets(await res.json())
            }
        } catch (error) {
            console.error('Failed to fetch datasets:', error)
        }
    }

    const checkJobStatus = async () => {
        if (!currentJob) return

        try {
            const res = await fetch(`${API_BASE}/analysis/status/${currentJob.job_id}`)
            if (res.ok) {
                const job = await res.json()
                setCurrentJob(job)

                if (job.status === 'running') {
                    setProgress(prev => Math.min(prev + 10, 90))
                } else if (job.status === 'completed') {
                    setProgress(100)
                    setRunning(false)

                    // Save results to project if a project is active
                    if (currentProject && job.results) {
                        const resultToSave = {
                            name: formData.name,
                            analysis_type: formData.analysis_type,
                            job_id: currentJob.job_id,
                            ...job.results
                        }
                        await saveResultToProject(resultToSave)
                    }
                } else if (job.status === 'failed') {
                    setRunning(false)
                }
            }
        } catch (error) {
            console.error('Status check error:', error)
        }
    }

    const startAnalysis = async () => {
        if (!formData.name || !formData.dataset_id) {
            alert('Please provide an analysis name and select a dataset')
            return
        }

        setRunning(true)
        setProgress(0)

        try {
            const res = await fetch(`${API_BASE}/samples/datasets/${formData.dataset_id}`)
            const datasetData = await res.json()

            const analysisRes = await fetch(`${API_BASE}/analysis/run`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    name: formData.name,
                    analysis_type: formData.analysis_type,
                    samples: datasetData.samples,
                    edges: [],
                    config: formData.config
                })
            })

            if (analysisRes.ok) {
                const job = await analysisRes.json()
                setCurrentJob(job)
                setProgress(10)
            }
        } catch (error) {
            console.error('Analysis error:', error)
            setRunning(false)
            alert('Failed to start analysis')
        }
    }

    return (
        <div className="analysis-page">
            <header className="page-header">
                <div>
                    <h1 className="page-title">Run Analysis</h1>
                    <p className="page-subtitle">
                        Configure and execute hydrogeochemical inverse modeling
                    </p>
                </div>
            </header>

            {/* Project Warning */}
            {!currentProject && (
                <div className="no-project-alert">
                    <div className="alert-icon">⚠️</div>
                    <div className="alert-content">
                        <h4>No Project Selected</h4>
                        <p>Your results will <strong>not be saved</strong>. Create or select a project to save your analysis results.</p>
                    </div>
                    <Link to="/projects" className="btn btn-primary">
                        Create Project
                    </Link>
                </div>
            )}

            {/* Project Active Indicator */}
            {currentProject && (
                <div className="project-active-alert">
                    <div className="alert-icon">✅</div>
                    <div className="alert-content">
                        <h4>Saving to: {currentProject.name}</h4>
                        <p>Analysis results will be automatically saved to your project.</p>
                    </div>
                </div>
            )}

            <div className="analysis-layout">
                {/* Configuration Form */}
                <div className="card config-card">
                    <div className="card-header">
                        <h3 className="card-title">Analysis Configuration</h3>
                    </div>

                    <div className="form-group">
                        <label className="form-label">Analysis Name</label>
                        <input
                            type="text"
                            className="form-input"
                            value={formData.name}
                            onChange={(e) => setFormData({ ...formData, name: e.target.value })}
                            placeholder="e.g., Site-A Transport Analysis"
                            disabled={running}
                        />
                    </div>

                    <div className="form-group">
                        <label className="form-label">Analysis Type</label>
                        <select
                            className="form-select"
                            value={formData.analysis_type}
                            onChange={(e) => setFormData({ ...formData, analysis_type: e.target.value })}
                            disabled={running}
                        >
                            <option value="full_pipeline">Full Pipeline (Recommended)</option>
                            <option value="transport">Transport Modeling Only</option>
                            <option value="reaction_path">Reaction Path Only</option>
                            <option value="network_inference">Network Inference Only</option>
                            <option value="nitrate_source">Nitrate Source Discrimination</option>
                        </select>
                    </div>

                    <div className="form-group">
                        <label className="form-label">Select Dataset</label>
                        <select
                            className="form-select"
                            value={formData.dataset_id}
                            onChange={(e) => setFormData({ ...formData, dataset_id: e.target.value })}
                            disabled={running}
                        >
                            <option value="">-- Select a dataset --</option>
                            {datasets.map(ds => (
                                <option key={ds.id} value={ds.id}>
                                    {ds.name} ({ds.sample_count} samples)
                                </option>
                            ))}
                        </select>
                    </div>

                    <div className="divider"></div>

                    <h4 className="section-title">Advanced Settings</h4>

                    <div className="form-group">
                        <label className="form-label">
                            LASSO Penalty (λ): {formData.config.lasso_penalty}
                        </label>
                        <input
                            type="range"
                            className="form-range"
                            min="0"
                            max="1"
                            step="0.01"
                            value={formData.config.lasso_penalty}
                            onChange={(e) => setFormData({
                                ...formData,
                                config: { ...formData.config, lasso_penalty: parseFloat(e.target.value) }
                            })}
                            disabled={running}
                        />
                        <span className="range-hint">Higher = more sparse reaction set</span>
                    </div>

                    <div className="toggle-group">
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_isotopes}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_isotopes: e.target.checked }
                                })}
                                disabled={running}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">Enable Isotope Analysis</span>
                        </label>
                    </div>

                    <div className="toggle-group">
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_uncertainty}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_uncertainty: e.target.checked }
                                })}
                                disabled={running}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">Enable Uncertainty Quantification</span>
                        </label>
                    </div>

                    <div className="toggle-group">
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_phreeqc}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_phreeqc: e.target.checked }
                                })}
                                disabled={running}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">Enable PHREEQC Constraints</span>
                        </label>
                    </div>

                    {formData.config.enable_uncertainty && (
                        <div className="form-group">
                            <label className="form-label">Bootstrap Iterations</label>
                            <input
                                type="number"
                                className="form-input"
                                value={formData.config.bootstrap_iterations}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, bootstrap_iterations: parseInt(e.target.value) }
                                })}
                                min="100"
                                max="10000"
                                disabled={running}
                            />
                        </div>
                    )}

                    <button
                        className="btn btn-primary w-full mt-lg"
                        onClick={startAnalysis}
                        disabled={running}
                    >
                        {running ? (
                            <>
                                <span className="spinner-small"></span>
                                Running Analysis...
                            </>
                        ) : (
                            <>
                                <span>⚗️</span> Start Analysis
                            </>
                        )}
                    </button>

                    {!currentProject && (
                        <p className="save-warning">⚠️ Results will not be saved without a project</p>
                    )}
                </div>

                {/* Progress Panel */}
                <div className="progress-panel">
                    <div className="card">
                        <div className="card-header">
                            <h3 className="card-title">Analysis Progress</h3>
                            {currentJob && (
                                <span className={`badge badge-${currentJob.status === 'completed' ? 'success' :
                                        currentJob.status === 'failed' ? 'error' : 'warning'
                                    }`}>
                                    {currentJob.status}
                                </span>
                            )}
                        </div>

                        {currentJob ? (
                            <div className="progress-content">
                                <div className="progress-header">
                                    <span className="job-name">{currentJob.name}</span>
                                    <span className="job-type">{currentJob.analysis_type}</span>
                                </div>

                                <div className="progress-bar-container">
                                    <div className="progress-bar">
                                        <div
                                            className="progress-fill"
                                            style={{ width: `${progress}%` }}
                                        ></div>
                                    </div>
                                    <span className="progress-text">{progress}%</span>
                                </div>

                                <div className="progress-steps">
                                    <div className={`step ${progress >= 10 ? 'completed' : progress > 0 ? 'active' : ''}`}>
                                        <span className="step-icon">📥</span>
                                        <span className="step-label">Loading Data</span>
                                    </div>
                                    <div className={`step ${progress >= 30 ? 'completed' : progress > 20 ? 'active' : ''}`}>
                                        <span className="step-icon">🔄</span>
                                        <span className="step-label">Transport Modeling</span>
                                    </div>
                                    <div className={`step ${progress >= 50 ? 'completed' : progress > 40 ? 'active' : ''}`}>
                                        <span className="step-icon">⚛️</span>
                                        <span className="step-label">Reaction Fitting</span>
                                    </div>
                                    <div className={`step ${progress >= 70 ? 'completed' : progress > 60 ? 'active' : ''}`}>
                                        <span className="step-icon">📊</span>
                                        <span className="step-label">Uncertainty Analysis</span>
                                    </div>
                                    <div className={`step ${progress >= 100 ? 'completed' : progress > 90 ? 'active' : ''}`}>
                                        <span className="step-icon">✅</span>
                                        <span className="step-label">Complete</span>
                                    </div>
                                </div>

                                {currentJob.status === 'completed' && (
                                    <div className="completion-actions">
                                        {currentProject ? (
                                            <div className="saved-notice">
                                                <span>✅ Results saved to "{currentProject.name}"</span>
                                            </div>
                                        ) : (
                                            <div className="not-saved-notice">
                                                <span>⚠️ Results not saved (no project)</span>
                                            </div>
                                        )}
                                        <Link to="/results" className="btn btn-primary w-full">
                                            View Results
                                        </Link>
                                    </div>
                                )}
                            </div>
                        ) : (
                            <div className="empty-state">
                                <span className="empty-icon">⚗️</span>
                                <p>No active analysis</p>
                                <span className="empty-hint">Configure and start an analysis to see progress</span>
                            </div>
                        )}
                    </div>

                    {/* Analysis Info */}
                    <div className="card mt-lg">
                        <div className="card-header">
                            <h3 className="card-title">Analysis Modules</h3>
                        </div>
                        <div className="modules-info">
                            <div className="module-item">
                                <span className="module-icon">🔀</span>
                                <div className="module-text">
                                    <h4>Transport Modeling</h4>
                                    <p>Weighted least squares for evaporation vs mixing</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">⚛️</span>
                                <div className="module-text">
                                    <h4>Reaction Path</h4>
                                    <p>Sparse LASSO regression for mineral reactions</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">🔗</span>
                                <div className="module-text">
                                    <h4>Network Inference</h4>
                                    <p>Probabilistic flow direction from head gradients</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">📊</span>
                                <div className="module-text">
                                    <h4>Uncertainty</h4>
                                    <p>Bayesian MCMC and bootstrap confidence intervals</p>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    )
}

export default Analysis
