import { useState, useEffect } from 'react'
import { Link } from 'react-router-dom'
import { useProject } from '../context/ProjectContext'
import './Analysis.css'

const API_BASE = 'http://localhost:8000/api'

function Analysis() {
    const { currentProject, saveResultToProject } = useProject()
    const [datasets, setDatasets] = useState([])
    const [capabilities, setCapabilities] = useState(null)
    const [capabilitiesLoading, setCapabilitiesLoading] = useState(false)
    const [formData, setFormData] = useState({
        name: '',
        analysis_type: 'full_pipeline',
        dataset_id: '',
        config: {
            // Core settings
            lasso_penalty: 0.1,

            // Core modules
            enable_phreeqc: false,
            enable_isotopes: true,
            enable_uncertainty: false,

            // Uncertainty settings
            uncertainty_method: 'bootstrap',
            bootstrap_iterations: 1000,
            bayesian_samples: 5000,
            bayesian_chains: 4,

            // Nitrate source
            enable_nitrate_source: false,
            nitrate_isotope_mixing: true,

            // Temporal analysis
            enable_temporal: false,
            temporal_window_days: 365,
            residence_time_method: 'cross_correlation',

            // 3D network
            enable_3d_network: false,
            vertical_flow_enabled: true,
            vertical_anisotropy: 0.1,
            enable_layer_system: false,

            // Reactive transport
            enable_reactive_transport: false,
            rt_simulator: 'phreeqc_kinetic',
            rt_time_steps: 100,

            // Vadose zone
            enable_vadose_zone: false,

            // CoDA
            enable_coda: false,

            // Edge inference
            edge_radius_km: 10.0,
            edge_max_neighbors: 3,
            edge_p_min: 0.75,
            edge_head_inference: 'heuristic',

            // Gibbs & Exchange
            enable_gibbs: true,
            gibbs_weight: 0.5,
            enable_exchange: true
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

    const fetchCapabilities = async (datasetId) => {
        if (!datasetId) {
            setCapabilities(null)
            return
        }
        setCapabilitiesLoading(true)
        try {
            const res = await fetch(`${API_BASE}/samples/datasets/${datasetId}/capabilities`)
            if (res.ok) {
                const caps = await res.json()
                setCapabilities(caps)
                // Auto-disable modules that aren't available
                setFormData(prev => ({
                    ...prev,
                    config: {
                        ...prev.config,
                        enable_isotopes: caps.available_analyses.isotope_analysis && prev.config.enable_isotopes,
                        enable_phreeqc: caps.available_analyses.phreeqc && prev.config.enable_phreeqc,
                        enable_nitrate_source: caps.available_analyses.nitrate_source && prev.config.enable_nitrate_source,
                        enable_temporal: caps.available_analyses.temporal && prev.config.enable_temporal,
                        enable_3d_network: caps.available_analyses.network_3d && prev.config.enable_3d_network,
                    }
                }))
            }
        } catch (error) {
            console.error('Failed to fetch capabilities:', error)
        } finally {
            setCapabilitiesLoading(false)
        }
    }

    // Fetch capabilities when dataset changes
    useEffect(() => {
        fetchCapabilities(formData.dataset_id)
    }, [formData.dataset_id])

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
                    config: formData.config,
                    project_id: currentProject?.id || null,
                    dataset_id: formData.dataset_id,
                    dataset_name: datasetData.name || '',
                    dataset_description: datasetData.description || ''
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

                    {/* Dataset Capabilities Info */}
                    {capabilitiesLoading && (
                        <div className="capabilities-loading">
                            <span className="spinner-small"></span> Analyzing dataset...
                        </div>
                    )}

                    {capabilities && !capabilitiesLoading && (
                        <div className="capabilities-panel">
                            <div className="capabilities-header">
                                <h5>Dataset Analysis</h5>
                                <span className="sample-count">{capabilities.sample_count} samples</span>
                            </div>

                            <div className="data-fields">
                                {capabilities.available_fields.ions.length > 0 && (
                                    <div className="field-group">
                                        <span className="field-label">Ions:</span>
                                        <span className="field-values">{capabilities.available_fields.ions.join(', ')}</span>
                                    </div>
                                )}
                                {capabilities.available_fields.isotopes.length > 0 && (
                                    <div className="field-group">
                                        <span className="field-label">Isotopes:</span>
                                        <span className="field-values">{capabilities.available_fields.isotopes.join(', ')}</span>
                                    </div>
                                )}
                                {capabilities.available_fields.spatial.length > 0 && (
                                    <div className="field-group">
                                        <span className="field-label">Spatial:</span>
                                        <span className="field-values">{capabilities.available_fields.spatial.join(', ')}</span>
                                    </div>
                                )}
                                {capabilities.available_fields.temporal.length > 0 && (
                                    <div className="field-group">
                                        <span className="field-label">Temporal:</span>
                                        <span className="field-values">{capabilities.available_fields.temporal.join(', ')}</span>
                                    </div>
                                )}
                            </div>

                            {capabilities.warnings.length > 0 && (
                                <div className="capabilities-warnings">
                                    {capabilities.warnings.map((warning, idx) => (
                                        <div key={idx} className="warning-item">
                                            <span className="warning-icon">⚠️</span>
                                            <span>{warning}</span>
                                        </div>
                                    ))}
                                </div>
                            )}
                        </div>
                    )}

                    <div className="divider"></div>

                    <h4 className="section-title">Core Settings</h4>

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

                    <div className={`toggle-group ${capabilities && !capabilities.available_analyses.isotope_analysis ? 'toggle-disabled' : ''}`}>
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_isotopes}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_isotopes: e.target.checked }
                                })}
                                disabled={running || (capabilities && !capabilities.available_analyses.isotope_analysis)}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">
                                Enable Isotope Analysis (δ18O/δ2H)
                                {capabilities && !capabilities.available_analyses.isotope_analysis && (
                                    <span className="toggle-warning" title="Requires d18o and d2h fields">⚠️ No isotope data</span>
                                )}
                            </span>
                        </label>
                    </div>

                    <div className={`toggle-group ${capabilities && !capabilities.available_analyses.phreeqc ? 'toggle-disabled' : ''}`}>
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_phreeqc}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_phreeqc: e.target.checked }
                                })}
                                disabled={running || (capabilities && !capabilities.available_analyses.phreeqc)}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">
                                Enable PHREEQC Constraints
                                {capabilities && !capabilities.available_analyses.phreeqc && (
                                    <span className="toggle-warning" title="Requires pH and temperature fields">⚠️ No pH/temp data</span>
                                )}
                            </span>
                        </label>
                    </div>

                    <div className="toggle-group">
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_gibbs}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_gibbs: e.target.checked }
                                })}
                                disabled={running}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">Enable Gibbs Diagram Analysis</span>
                        </label>
                    </div>

                    <div className="toggle-group">
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_exchange}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_exchange: e.target.checked }
                                })}
                                disabled={running}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">Enable Ion Exchange Detection</span>
                        </label>
                    </div>

                    <div className="divider"></div>
                    <h4 className="section-title">Uncertainty Quantification</h4>

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
                            <span className="toggle-label">Enable Uncertainty Analysis</span>
                        </label>
                    </div>

                    {formData.config.enable_uncertainty && (
                        <>
                            <div className="form-group">
                                <label className="form-label">Uncertainty Method</label>
                                <select
                                    className="form-select"
                                    value={formData.config.uncertainty_method}
                                    onChange={(e) => setFormData({
                                        ...formData,
                                        config: { ...formData.config, uncertainty_method: e.target.value }
                                    })}
                                    disabled={running}
                                >
                                    <option value="bootstrap">Bootstrap (Recommended)</option>
                                    <option value="bayesian">Bayesian MCMC</option>
                                    <option value="monte_carlo">Monte Carlo</option>
                                </select>
                            </div>

                            {formData.config.uncertainty_method === 'bootstrap' && (
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

                            {formData.config.uncertainty_method === 'bayesian' && (
                                <>
                                    <div className="form-group">
                                        <label className="form-label">Bayesian Samples</label>
                                        <input
                                            type="number"
                                            className="form-input"
                                            value={formData.config.bayesian_samples}
                                            onChange={(e) => setFormData({
                                                ...formData,
                                                config: { ...formData.config, bayesian_samples: parseInt(e.target.value) }
                                            })}
                                            min="1000"
                                            max="50000"
                                            disabled={running}
                                        />
                                    </div>
                                    <div className="form-group">
                                        <label className="form-label">MCMC Chains</label>
                                        <input
                                            type="number"
                                            className="form-input"
                                            value={formData.config.bayesian_chains}
                                            onChange={(e) => setFormData({
                                                ...formData,
                                                config: { ...formData.config, bayesian_chains: parseInt(e.target.value) }
                                            })}
                                            min="1"
                                            max="8"
                                            disabled={running}
                                        />
                                    </div>
                                </>
                            )}
                        </>
                    )}

                    <div className="divider"></div>
                    <h4 className="section-title">Nitrate Source Discrimination</h4>

                    <div className={`toggle-group ${capabilities && !capabilities.available_analyses.nitrate_source ? 'toggle-disabled' : ''}`}>
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_nitrate_source}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_nitrate_source: e.target.checked }
                                })}
                                disabled={running || (capabilities && !capabilities.available_analyses.nitrate_source)}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">
                                Enable Nitrate Source Analysis
                                {capabilities && !capabilities.available_analyses.nitrate_source && (
                                    <span className="toggle-warning" title="Requires NO3, d15n, and d18o_no3 fields">⚠️ Missing NO3 isotopes</span>
                                )}
                            </span>
                        </label>
                    </div>

                    {formData.config.enable_nitrate_source && (
                        <div className="toggle-group">
                            <label className="toggle">
                                <input
                                    type="checkbox"
                                    checked={formData.config.nitrate_isotope_mixing}
                                    onChange={(e) => setFormData({
                                        ...formData,
                                        config: { ...formData.config, nitrate_isotope_mixing: e.target.checked }
                                    })}
                                    disabled={running}
                                />
                                <span className="toggle-slider"></span>
                                <span className="toggle-label">Use Dual-Isotope Mixing (δ15N/δ18O-NO3)</span>
                            </label>
                        </div>
                    )}

                    <div className="divider"></div>
                    <h4 className="section-title">Temporal Analysis</h4>

                    <div className={`toggle-group ${capabilities && !capabilities.available_analyses.temporal ? 'toggle-disabled' : ''}`}>
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_temporal}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_temporal: e.target.checked }
                                })}
                                disabled={running || (capabilities && !capabilities.available_analyses.temporal)}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">
                                Enable Temporal/Residence Time Analysis
                                {capabilities && !capabilities.available_analyses.temporal && (
                                    <span className="toggle-warning" title="Requires date field with time series data">⚠️ No temporal data</span>
                                )}
                            </span>
                        </label>
                    </div>

                    {formData.config.enable_temporal && (
                        <>
                            <div className="form-group">
                                <label className="form-label">Residence Time Method</label>
                                <select
                                    className="form-select"
                                    value={formData.config.residence_time_method}
                                    onChange={(e) => setFormData({
                                        ...formData,
                                        config: { ...formData.config, residence_time_method: e.target.value }
                                    })}
                                    disabled={running}
                                >
                                    <option value="cross_correlation">Cross-Correlation (Recommended)</option>
                                    <option value="gradient">Hydraulic Gradient</option>
                                    <option value="bayesian_lag">Bayesian Lag</option>
                                    <option value="ttd">Transit Time Distribution</option>
                                    <option value="tracer_decay">Tracer Decay</option>
                                </select>
                            </div>
                            <div className="form-group">
                                <label className="form-label">Temporal Window (days)</label>
                                <input
                                    type="number"
                                    className="form-input"
                                    value={formData.config.temporal_window_days}
                                    onChange={(e) => setFormData({
                                        ...formData,
                                        config: { ...formData.config, temporal_window_days: parseInt(e.target.value) }
                                    })}
                                    min="30"
                                    max="3650"
                                    disabled={running}
                                />
                            </div>
                        </>
                    )}

                    <div className="divider"></div>
                    <h4 className="section-title">3D Network / Layered Aquifer</h4>

                    <div className={`toggle-group ${capabilities && !capabilities.available_analyses.network_3d ? 'toggle-disabled' : ''}`}>
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_3d_network}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_3d_network: e.target.checked }
                                })}
                                disabled={running || (capabilities && !capabilities.available_analyses.network_3d)}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">
                                Enable 3D Flow Network
                                {capabilities && !capabilities.available_analyses.network_3d && (
                                    <span className="toggle-warning" title="Requires x, y, and z/elevation coordinates">⚠️ No 3D coordinates</span>
                                )}
                            </span>
                        </label>
                    </div>

                    {formData.config.enable_3d_network && (
                        <>
                            <div className="toggle-group">
                                <label className="toggle">
                                    <input
                                        type="checkbox"
                                        checked={formData.config.vertical_flow_enabled}
                                        onChange={(e) => setFormData({
                                            ...formData,
                                            config: { ...formData.config, vertical_flow_enabled: e.target.checked }
                                        })}
                                        disabled={running}
                                    />
                                    <span className="toggle-slider"></span>
                                    <span className="toggle-label">Allow Vertical Flow</span>
                                </label>
                            </div>
                            <div className="form-group">
                                <label className="form-label">
                                    Vertical Anisotropy (Kh/Kv): {formData.config.vertical_anisotropy}
                                </label>
                                <input
                                    type="range"
                                    className="form-range"
                                    min="0.001"
                                    max="1"
                                    step="0.01"
                                    value={formData.config.vertical_anisotropy}
                                    onChange={(e) => setFormData({
                                        ...formData,
                                        config: { ...formData.config, vertical_anisotropy: parseFloat(e.target.value) }
                                    })}
                                    disabled={running}
                                />
                            </div>
                            <div className="toggle-group">
                                <label className="toggle">
                                    <input
                                        type="checkbox"
                                        checked={formData.config.enable_layer_system}
                                        onChange={(e) => setFormData({
                                            ...formData,
                                            config: { ...formData.config, enable_layer_system: e.target.checked }
                                        })}
                                        disabled={running}
                                    />
                                    <span className="toggle-slider"></span>
                                    <span className="toggle-label">Enable Layer System</span>
                                </label>
                            </div>
                        </>
                    )}

                    <div className="divider"></div>
                    <h4 className="section-title">Reactive Transport Validation</h4>

                    <div className="toggle-group">
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_reactive_transport}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_reactive_transport: e.target.checked }
                                })}
                                disabled={running}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">Enable Reactive Transport Validation</span>
                        </label>
                    </div>

                    {formData.config.enable_reactive_transport && (
                        <>
                            <div className="form-group">
                                <label className="form-label">RT Simulator</label>
                                <select
                                    className="form-select"
                                    value={formData.config.rt_simulator}
                                    onChange={(e) => setFormData({
                                        ...formData,
                                        config: { ...formData.config, rt_simulator: e.target.value }
                                    })}
                                    disabled={running}
                                >
                                    <option value="phreeqc_kinetic">PHREEQC Kinetic</option>
                                    <option value="mt3dms">MT3DMS</option>
                                </select>
                            </div>
                            <div className="form-group">
                                <label className="form-label">Time Steps</label>
                                <input
                                    type="number"
                                    className="form-input"
                                    value={formData.config.rt_time_steps}
                                    onChange={(e) => setFormData({
                                        ...formData,
                                        config: { ...formData.config, rt_time_steps: parseInt(e.target.value) }
                                    })}
                                    min="10"
                                    max="1000"
                                    disabled={running}
                                />
                            </div>
                        </>
                    )}

                    <div className="divider"></div>
                    <h4 className="section-title">Additional Modules</h4>

                    <div className="toggle-group">
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_vadose_zone}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_vadose_zone: e.target.checked }
                                })}
                                disabled={running}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">Enable Vadose Zone Modeling</span>
                        </label>
                    </div>

                    <div className="toggle-group">
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_coda}
                                onChange={(e) => setFormData({
                                    ...formData,
                                    config: { ...formData.config, enable_coda: e.target.checked }
                                })}
                                disabled={running}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">Enable CoDA (Compositional Data Analysis)</span>
                        </label>
                    </div>

                    <div className="divider"></div>
                    <h4 className="section-title">Edge Inference Settings</h4>

                    <div className="form-group">
                        <label className="form-label">Head Inference Method</label>
                        <select
                            className="form-select"
                            value={formData.config.edge_head_inference}
                            onChange={(e) => setFormData({
                                ...formData,
                                config: { ...formData.config, edge_head_inference: e.target.value }
                            })}
                            disabled={running}
                        >
                            <option value="heuristic">Heuristic (Fast)</option>
                            <option value="bayesian">Bayesian Linear</option>
                            <option value="bayesian_mcmc">Bayesian MCMC</option>
                        </select>
                    </div>

                    <div className="form-group">
                        <label className="form-label">
                            Search Radius: {formData.config.edge_radius_km} km
                        </label>
                        <input
                            type="range"
                            className="form-range"
                            min="0.5"
                            max="50"
                            step="0.5"
                            value={formData.config.edge_radius_km}
                            onChange={(e) => setFormData({
                                ...formData,
                                config: { ...formData.config, edge_radius_km: parseFloat(e.target.value) }
                            })}
                            disabled={running}
                        />
                    </div>

                    <div className="form-group">
                        <label className="form-label">Max Neighbors: {formData.config.edge_max_neighbors}</label>
                        <input
                            type="range"
                            className="form-range"
                            min="1"
                            max="10"
                            step="1"
                            value={formData.config.edge_max_neighbors}
                            onChange={(e) => setFormData({
                                ...formData,
                                config: { ...formData.config, edge_max_neighbors: parseInt(e.target.value) }
                            })}
                            disabled={running}
                        />
                    </div>

                    <div className="form-group">
                        <label className="form-label">
                            Min Edge Probability: {formData.config.edge_p_min}
                        </label>
                        <input
                            type="range"
                            className="form-range"
                            min="0"
                            max="1"
                            step="0.05"
                            value={formData.config.edge_p_min}
                            onChange={(e) => setFormData({
                                ...formData,
                                config: { ...formData.config, edge_p_min: parseFloat(e.target.value) }
                            })}
                            disabled={running}
                        />
                    </div>

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
                            <h3 className="card-title">Available Modules</h3>
                        </div>
                        <div className="modules-info">
                            <div className="module-item">
                                <span className="module-icon">🔀</span>
                                <div className="module-text">
                                    <h4>Transport Modeling</h4>
                                    <p>Evaporation vs mixing selection</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">⚛️</span>
                                <div className="module-text">
                                    <h4>Reaction Path</h4>
                                    <p>LASSO sparse mineral reactions</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">🔗</span>
                                <div className="module-text">
                                    <h4>Network Inference</h4>
                                    <p>Probabilistic flow direction</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">📊</span>
                                <div className="module-text">
                                    <h4>Uncertainty</h4>
                                    <p>Bootstrap/Bayesian MCMC</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">🌡️</span>
                                <div className="module-text">
                                    <h4>Isotopes</h4>
                                    <p>δ18O/δ2H fractionation</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">🧪</span>
                                <div className="module-text">
                                    <h4>PHREEQC</h4>
                                    <p>Saturation index constraints</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">🌾</span>
                                <div className="module-text">
                                    <h4>Nitrate Source</h4>
                                    <p>Dual-isotope discrimination</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">⏱️</span>
                                <div className="module-text">
                                    <h4>Temporal</h4>
                                    <p>Residence time estimation</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">📐</span>
                                <div className="module-text">
                                    <h4>3D Network</h4>
                                    <p>Layered aquifer systems</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">🔄</span>
                                <div className="module-text">
                                    <h4>Reactive Transport</h4>
                                    <p>Kinetic validation</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">💧</span>
                                <div className="module-text">
                                    <h4>Vadose Zone</h4>
                                    <p>Unsaturated flow modeling</p>
                                </div>
                            </div>
                            <div className="module-item">
                                <span className="module-icon">📈</span>
                                <div className="module-text">
                                    <h4>CoDA</h4>
                                    <p>Compositional data analysis</p>
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
