import { Link } from 'react-router-dom'

const progressSteps = [
    { threshold: 10, icon: '📥', label: 'Loading Data' },
    { threshold: 30, icon: '🔄', label: 'Transport Modeling' },
    { threshold: 50, icon: '⚛️', label: 'Reaction Fitting' },
    { threshold: 70, icon: '📊', label: 'Uncertainty Analysis' },
    { threshold: 100, icon: '✅', label: 'Complete' },
]

function AnalysisProgress({ currentJob, progress, currentProject }) {
    const getStepStatus = (threshold, index) => {
        const prevThreshold = index > 0 ? progressSteps[index - 1].threshold : 0
        if (progress >= threshold) return 'completed'
        if (progress > prevThreshold) return 'active'
        return ''
    }

    return (
        <div className="card">
            <div className="card-header">
                <h3 className="card-title">Analysis Progress</h3>
                {currentJob && (
                    <span className={`badge badge-${
                        currentJob.status === 'completed' ? 'success' :
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
                            />
                        </div>
                        <span className="progress-text">{progress}%</span>
                    </div>

                    <div className="progress-steps">
                        {progressSteps.map((step, idx) => (
                            <div key={idx} className={`step ${getStepStatus(step.threshold, idx)}`}>
                                <span className="step-icon">{step.icon}</span>
                                <span className="step-label">{step.label}</span>
                            </div>
                        ))}
                    </div>

                    {currentJob.status === 'completed' && (
                        <div className="completion-actions">
                            {currentProject ? (
                                <div className="saved-notice">
                                    <span>Results saved to "{currentProject.name}"</span>
                                </div>
                            ) : (
                                <div className="not-saved-notice">
                                    <span>Results not saved (no project)</span>
                                </div>
                            )}
                            <Link to="/results" className="btn btn-primary w-full">
                                View Results
                            </Link>
                        </div>
                    )}

                    {currentJob.status === 'failed' && currentJob.error && (
                        <div className="error-notice">
                            <strong>Error:</strong> {currentJob.error}
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
    )
}

export default AnalysisProgress
