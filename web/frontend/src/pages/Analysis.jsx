import { useReducer, useEffect, useCallback } from 'react'
import { Link } from 'react-router-dom'
import { useProject } from '../context/ProjectContext'
import { useToast } from '../context/ToastContext'
import { useWebSocket } from '../hooks/useWebSocket'
import { API_BASE } from '../config'
import './Analysis.css'

// Initial configuration defaults
const initialConfig = {
    // Core settings
    lasso_penalty: 0.1,

    // Core modules
    enable_phreeqc: true,  // Always enabled by default - critical for preventing modeling errors
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

// Initial state for the reducer
const initialState = {
    // Data
    datasets: [],
    capabilities: null,
    capabilitiesLoading: false,

    // Form
    formData: {
        name: '',
        analysis_type: 'full_pipeline',
        dataset_id: '',
        config: initialConfig
    },

    // Analysis execution
    running: false,
    currentJob: null,
    progress: 0,
    currentStep: '',

    // Validation
    validationErrors: {}
}

// Action types
const ActionTypes = {
    SET_DATASETS: 'SET_DATASETS',
    SET_CAPABILITIES: 'SET_CAPABILITIES',
    SET_CAPABILITIES_LOADING: 'SET_CAPABILITIES_LOADING',
    UPDATE_FORM_FIELD: 'UPDATE_FORM_FIELD',
    UPDATE_CONFIG: 'UPDATE_CONFIG',
    UPDATE_CONFIG_FROM_CAPABILITIES: 'UPDATE_CONFIG_FROM_CAPABILITIES',
    START_ANALYSIS: 'START_ANALYSIS',
    SET_CURRENT_JOB: 'SET_CURRENT_JOB',
    UPDATE_PROGRESS: 'UPDATE_PROGRESS',
    ANALYSIS_COMPLETE: 'ANALYSIS_COMPLETE',
    ANALYSIS_FAILED: 'ANALYSIS_FAILED',
    SET_VALIDATION_ERRORS: 'SET_VALIDATION_ERRORS',
    RESET_ANALYSIS: 'RESET_ANALYSIS'
}

// Reducer function
function analysisReducer(state, action) {
    switch (action.type) {
        case ActionTypes.SET_DATASETS:
            return { ...state, datasets: action.payload }

        case ActionTypes.SET_CAPABILITIES:
            return { ...state, capabilities: action.payload }

        case ActionTypes.SET_CAPABILITIES_LOADING:
            return { ...state, capabilitiesLoading: action.payload }

        case ActionTypes.UPDATE_FORM_FIELD:
            return {
                ...state,
                formData: { ...state.formData, [action.field]: action.value }
            }

        case ActionTypes.UPDATE_CONFIG:
            return {
                ...state,
                formData: {
                    ...state.formData,
                    config: { ...state.formData.config, [action.field]: action.value }
                }
            }

        case ActionTypes.UPDATE_CONFIG_FROM_CAPABILITIES: {
            const caps = action.payload
            return {
                ...state,
                formData: {
                    ...state.formData,
                    config: {
                        ...state.formData.config,
                        enable_isotopes: caps.available_analyses.isotope_analysis && state.formData.config.enable_isotopes,
                        // PHREEQC is always kept enabled - it's critical for model accuracy
                        // Users can manually disable if needed, but we don't auto-disable
                        enable_nitrate_source: caps.available_analyses.nitrate_source && state.formData.config.enable_nitrate_source,
                        enable_temporal: caps.available_analyses.temporal && state.formData.config.enable_temporal,
                        enable_3d_network: caps.available_analyses.network_3d && state.formData.config.enable_3d_network,
                    }
                }
            }
        }

        case ActionTypes.START_ANALYSIS:
            return { ...state, running: true, progress: 0, currentStep: '' }

        case ActionTypes.SET_CURRENT_JOB:
            return { ...state, currentJob: action.payload }

        case ActionTypes.UPDATE_PROGRESS:
            return {
                ...state,
                progress: action.progress,
                currentStep: action.step !== undefined ? action.step : state.currentStep
            }

        case ActionTypes.ANALYSIS_COMPLETE:
            return {
                ...state,
                running: false,
                progress: 100,
                currentJob: { ...state.currentJob, status: 'completed', results: action.results }
            }

        case ActionTypes.ANALYSIS_FAILED:
            return {
                ...state,
                running: false,
                currentJob: { ...state.currentJob, status: 'failed', error: action.error }
            }

        case ActionTypes.SET_VALIDATION_ERRORS:
            return { ...state, validationErrors: action.payload }

        case ActionTypes.RESET_ANALYSIS:
            return {
                ...state,
                running: false,
                currentJob: null,
                progress: 0,
                currentStep: ''
            }

        default:
            return state
    }
}

function Analysis() {
    const { currentProject, saveResultToProject } = useProject()
    const toast = useToast()
    const [state, dispatch] = useReducer(analysisReducer, initialState)

    // Destructure state for easier access
    const {
        datasets,
        capabilities,
        capabilitiesLoading,
        formData,
        running,
        currentJob,
        progress,
        currentStep,
        validationErrors
    } = state

    // Validation function for config values
    const validateConfig = useCallback(() => {
        const errors = {}
        const config = formData.config

        // Bootstrap iterations: 10-10,000
        if (config.enable_uncertainty && config.uncertainty_method === 'bootstrap') {
            if (config.bootstrap_iterations < 10 || config.bootstrap_iterations > 10000) {
                errors.bootstrap_iterations = 'Must be between 10 and 10,000'
            }
        }

        // Bayesian samples: 100-50,000
        if (config.enable_uncertainty && config.uncertainty_method === 'bayesian') {
            if (config.bayesian_samples < 100 || config.bayesian_samples > 50000) {
                errors.bayesian_samples = 'Must be between 100 and 50,000'
            }
            if (config.bayesian_chains < 1 || config.bayesian_chains > 8) {
                errors.bayesian_chains = 'Must be between 1 and 8'
            }
        }

        // Edge radius: 0.1-100 km
        if (config.edge_radius_km < 0.1 || config.edge_radius_km > 100) {
            errors.edge_radius_km = 'Must be between 0.1 and 100 km'
        }

        // LASSO penalty: 0-10
        if (config.lasso_penalty < 0 || config.lasso_penalty > 10) {
            errors.lasso_penalty = 'Must be between 0 and 10'
        }

        // p_min: 0-1
        if (config.edge_p_min < 0 || config.edge_p_min > 1) {
            errors.edge_p_min = 'Must be between 0 and 1'
        }

        // Temporal window: 30-3650 days
        if (config.enable_temporal) {
            if (config.temporal_window_days < 30 || config.temporal_window_days > 3650) {
                errors.temporal_window_days = 'Must be between 30 and 3,650 days'
            }
        }

        // RT time steps: 10-1000
        if (config.enable_reactive_transport) {
            if (config.rt_time_steps < 10 || config.rt_time_steps > 1000) {
                errors.rt_time_steps = 'Must be between 10 and 1,000'
            }
        }

        dispatch({ type: ActionTypes.SET_VALIDATION_ERRORS, payload: errors })
        return Object.keys(errors).length === 0
    }, [formData.config])

    // Handle WebSocket progress updates
    const handleProgress = useCallback((data) => {
        dispatch({ type: ActionTypes.UPDATE_PROGRESS, progress: data.progress, step: data.step })
    }, [])

    // Handle WebSocket completion
    const handleComplete = useCallback(async (data) => {
        dispatch({ type: ActionTypes.UPDATE_PROGRESS, progress: data.progress })
        if (data.status === 'completed') {
            // Fetch full results
            try {
                const res = await fetch(`${API_BASE}/analysis/results/${currentJob?.job_id}`)
                if (res.ok) {
                    const results = await res.json()
                    dispatch({ type: ActionTypes.ANALYSIS_COMPLETE, results })

                    // Save results to project if a project is active
                    if (currentProject && results) {
                        const resultToSave = {
                            name: formData.name,
                            analysis_type: formData.analysis_type,
                            job_id: currentJob?.job_id,
                            ...results
                        }
                        await saveResultToProject(resultToSave)
                    }
                }
            } catch (error) {
                console.error('Failed to fetch results:', error)
            }
        } else if (data.status === 'failed') {
            dispatch({ type: ActionTypes.ANALYSIS_FAILED, error: data.error })
            toast.error(`Analysis failed: ${data.error || 'Unknown error'}`)
        }
    }, [currentJob?.job_id, currentProject, formData.name, formData.analysis_type, saveResultToProject, toast])

    // WebSocket connection for real-time progress
    const { isConnected: wsConnected } = useWebSocket(
        currentJob?.job_id,
        {
            enabled: running && !!currentJob?.job_id,
            onProgress: handleProgress,
            onComplete: handleComplete
        }
    )

    // Fetch datasets on mount
    const fetchDatasets = useCallback(async () => {
        try {
            const res = await fetch(`${API_BASE}/samples/datasets`)
            if (res.ok) {
                dispatch({ type: ActionTypes.SET_DATASETS, payload: await res.json() })
            }
        } catch (error) {
            console.error('Failed to fetch datasets:', error)
        }
    }, [])

    useEffect(() => {
        fetchDatasets()
    }, [fetchDatasets])

    // Fetch capabilities for selected dataset
    const fetchCapabilities = useCallback(async (datasetId) => {
        if (!datasetId) {
            dispatch({ type: ActionTypes.SET_CAPABILITIES, payload: null })
            return
        }
        dispatch({ type: ActionTypes.SET_CAPABILITIES_LOADING, payload: true })
        try {
            const res = await fetch(`${API_BASE}/samples/datasets/${datasetId}/capabilities`)
            if (res.ok) {
                const caps = await res.json()
                dispatch({ type: ActionTypes.SET_CAPABILITIES, payload: caps })
                // Auto-disable modules that aren't available
                dispatch({ type: ActionTypes.UPDATE_CONFIG_FROM_CAPABILITIES, payload: caps })
            }
        } catch (error) {
            console.error('Failed to fetch capabilities:', error)
        } finally {
            dispatch({ type: ActionTypes.SET_CAPABILITIES_LOADING, payload: false })
        }
    }, [])

    // Fetch capabilities when dataset changes
    useEffect(() => {
        fetchCapabilities(formData.dataset_id)
    }, [formData.dataset_id, fetchCapabilities])

    // Helper functions for form updates
    const updateFormField = useCallback((field, value) => {
        dispatch({ type: ActionTypes.UPDATE_FORM_FIELD, field, value })
    }, [])

    const updateConfig = useCallback((field, value) => {
        dispatch({ type: ActionTypes.UPDATE_CONFIG, field, value })
    }, [])

    // Start analysis
    const startAnalysis = useCallback(async () => {
        if (!formData.name || !formData.dataset_id) {
            toast.warning('Please provide an analysis name and select a dataset')
            return
        }

        // Validate config values
        if (!validateConfig()) {
            toast.error('Please fix validation errors before running analysis')
            return
        }

        dispatch({ type: ActionTypes.START_ANALYSIS })

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
                dispatch({ type: ActionTypes.SET_CURRENT_JOB, payload: job })
                dispatch({ type: ActionTypes.UPDATE_PROGRESS, progress: 10 })
            }
        } catch (error) {
            console.error('Analysis error:', error)
            dispatch({ type: ActionTypes.RESET_ANALYSIS })
            toast.error('Failed to start analysis. Please try again.')
        }
    }, [formData, currentProject, validateConfig, toast])

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
                            onChange={(e) => updateFormField('name', e.target.value)}
                            placeholder="e.g., Site-A Transport Analysis"
                            disabled={running}
                        />
                    </div>

                    <div className="form-group">
                        <label className="form-label">Analysis Type</label>
                        <select
                            className="form-select"
                            value={formData.analysis_type}
                            onChange={(e) => updateFormField('analysis_type', e.target.value)}
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
                            onChange={(e) => updateFormField('dataset_id', e.target.value)}
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
                            onChange={(e) => updateConfig('lasso_penalty', parseFloat(e.target.value))}
                            disabled={running}
                        />
                        <span className="range-hint">Higher = more sparse reaction set</span>
                    </div>

                    <div className={`toggle-group ${capabilities && !capabilities.available_analyses.isotope_analysis ? 'toggle-disabled' : ''}`}>
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_isotopes}
                                onChange={(e) => updateConfig('enable_isotopes', e.target.checked)}
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

                    <div className="toggle-group">
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_phreeqc}
                                onChange={(e) => updateConfig('enable_phreeqc', e.target.checked)}
                                disabled={running}
                            />
                            <span className="toggle-slider"></span>
                            <span className="toggle-label">
                                Enable PHREEQC Constraints (Recommended)
                                {capabilities && !capabilities.available_analyses.phreeqc && (
                                    <span className="toggle-warning" title="Dataset lacks pH/temperature but PHREEQC can still provide valuable constraints">⚠️ No pH/temp data</span>
                                )}
                            </span>
                        </label>
                    </div>

                    <div className="toggle-group">
                        <label className="toggle">
                            <input
                                type="checkbox"
                                checked={formData.config.enable_gibbs}
                                onChange={(e) => updateConfig('enable_gibbs', e.target.checked)}
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
                                onChange={(e) => updateConfig('enable_exchange', e.target.checked)}
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
                                onChange={(e) => updateConfig('enable_uncertainty', e.target.checked)}
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
                                    onChange={(e) => updateConfig('uncertainty_method', e.target.value)}
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
                                        className={`form-input ${validationErrors.bootstrap_iterations ? 'input-error' : ''}`}
                                        value={formData.config.bootstrap_iterations}
                                        onChange={(e) => updateConfig('bootstrap_iterations', parseInt(e.target.value) || 0)}
                                        min="10"
                                        max="10000"
                                        disabled={running}
                                    />
                                    {validationErrors.bootstrap_iterations && (
                                        <span className="validation-error">{validationErrors.bootstrap_iterations}</span>
                                    )}
                                </div>
                            )}

                            {formData.config.uncertainty_method === 'bayesian' && (
                                <>
                                    <div className="form-group">
                                        <label className="form-label">Bayesian Samples</label>
                                        <input
                                            type="number"
                                            className={`form-input ${validationErrors.bayesian_samples ? 'input-error' : ''}`}
                                            value={formData.config.bayesian_samples}
                                            onChange={(e) => updateConfig('bayesian_samples', parseInt(e.target.value) || 0)}
                                            min="100"
                                            max="50000"
                                            disabled={running}
                                        />
                                        {validationErrors.bayesian_samples && (
                                            <span className="validation-error">{validationErrors.bayesian_samples}</span>
                                        )}
                                    </div>
                                    <div className="form-group">
                                        <label className="form-label">MCMC Chains</label>
                                        <input
                                            type="number"
                                            className={`form-input ${validationErrors.bayesian_chains ? 'input-error' : ''}`}
                                            value={formData.config.bayesian_chains}
                                            onChange={(e) => updateConfig('bayesian_chains', parseInt(e.target.value) || 0)}
                                            min="1"
                                            max="8"
                                            disabled={running}
                                        />
                                        {validationErrors.bayesian_chains && (
                                            <span className="validation-error">{validationErrors.bayesian_chains}</span>
                                        )}
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
                                onChange={(e) => updateConfig('enable_nitrate_source', e.target.checked)}
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
                                    onChange={(e) => updateConfig('nitrate_isotope_mixing', e.target.checked)}
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
                                onChange={(e) => updateConfig('enable_temporal', e.target.checked)}
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
                                    onChange={(e) => updateConfig('residence_time_method', e.target.value)}
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
                                    onChange={(e) => updateConfig('temporal_window_days', parseInt(e.target.value))}
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
                                onChange={(e) => updateConfig('enable_3d_network', e.target.checked)}
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
                                        onChange={(e) => updateConfig('vertical_flow_enabled', e.target.checked)}
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
                                    onChange={(e) => updateConfig('vertical_anisotropy', parseFloat(e.target.value))}
                                    disabled={running}
                                />
                            </div>
                            <div className="toggle-group">
                                <label className="toggle">
                                    <input
                                        type="checkbox"
                                        checked={formData.config.enable_layer_system}
                                        onChange={(e) => updateConfig('enable_layer_system', e.target.checked)}
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
                                onChange={(e) => updateConfig('enable_reactive_transport', e.target.checked)}
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
                                    onChange={(e) => updateConfig('rt_simulator', e.target.value)}
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
                                    onChange={(e) => updateConfig('rt_time_steps', parseInt(e.target.value))}
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
                                onChange={(e) => updateConfig('enable_vadose_zone', e.target.checked)}
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
                                onChange={(e) => updateConfig('enable_coda', e.target.checked)}
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
                            onChange={(e) => updateConfig('edge_head_inference', e.target.value)}
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
                            onChange={(e) => updateConfig('edge_radius_km', parseFloat(e.target.value))}
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
                            onChange={(e) => updateConfig('edge_max_neighbors', parseInt(e.target.value))}
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
                            onChange={(e) => updateConfig('edge_p_min', parseFloat(e.target.value))}
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

                                {currentStep && (
                                    <div className="current-step-display">
                                        <span className="step-indicator">Current Step:</span>
                                        <span className="step-name">{currentStep}</span>
                                        {wsConnected && <span className="ws-indicator" title="Real-time updates">Live</span>}
                                    </div>
                                )}

                                <div className="progress-steps">
                                    <div className={`step ${progress >= 20 ? 'completed' : progress >= 10 ? 'active' : ''}`}>
                                        <span className="step-icon">📥</span>
                                        <span className="step-label">Initializing</span>
                                    </div>
                                    <div className={`step ${progress >= 40 ? 'completed' : progress >= 30 ? 'active' : ''}`}>
                                        <span className="step-icon">🔄</span>
                                        <span className="step-label">Converting Samples</span>
                                    </div>
                                    <div className={`step ${progress >= 60 ? 'completed' : progress >= 50 ? 'active' : ''}`}>
                                        <span className="step-icon">🔗</span>
                                        <span className="step-label">Building Network</span>
                                    </div>
                                    <div className={`step ${progress >= 80 ? 'completed' : progress >= 70 ? 'active' : ''}`}>
                                        <span className="step-icon">⚛️</span>
                                        <span className="step-label">Running Pipeline</span>
                                    </div>
                                    <div className={`step ${progress >= 100 ? 'completed' : progress >= 90 ? 'active' : ''}`}>
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
