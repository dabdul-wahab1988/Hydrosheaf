import { useState, useEffect } from 'react'
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, RadarChart, PolarGrid, PolarAngleAxis, PolarRadiusAxis, Radar, Legend } from 'recharts'
import { API_BASE } from '../config'
import ResultsGraph from '../components/ResultsGraph'
import './Results.css'

function Results() {
    const [jobs, setJobs] = useState([])
    const [selectedJob, setSelectedJob] = useState(null)
    const [results, setResults] = useState(null)
    const [loading, setLoading] = useState(true)

    useEffect(() => {
        fetchJobs()
    }, [])

    const fetchJobs = async () => {
        try {
            const res = await fetch(`${API_BASE}/analysis/jobs`)
            if (res.ok) {
                const data = await res.json()
                setJobs(data)
                if (data.length > 0) {
                    const completed = data.find(j => j.status === 'completed')
                    if (completed) {
                        selectJob(completed.job_id)
                    }
                }
            }
        } catch (error) {
            console.error('Failed to fetch jobs:', error)
        } finally {
            setLoading(false)
        }
    }

    const selectJob = async (jobId) => {
        setSelectedJob(jobId)
        try {
            const res = await fetch(`${API_BASE}/analysis/results/${jobId}`)
            if (res.ok) {
                setResults(await res.json())
            }
        } catch (error) {
            console.error('Failed to fetch results:', error)
        }
    }

    const deleteJob = async (jobId) => {
        if (!confirm('Delete this analysis job?')) return
        try {
            await fetch(`${API_BASE}/analysis/jobs/${jobId}`, { method: 'DELETE' })
            fetchJobs()
            if (selectedJob === jobId) {
                setSelectedJob(null)
                setResults(null)
            }
        } catch (error) {
            console.error('Delete error:', error)
        }
    }

    // Chart data from results
    const reactionData = results?.reactions?.map(r => ({
        name: r.mineral,
        rate: Math.abs(r.rate_mmol_L),
        direction: r.direction
    })) || []

    const mixingData = results?.transport_model?.mixing_fractions ?
        Object.entries(results.transport_model.mixing_fractions).map(([key, value]) => ({
            subject: key.replace('_', ' '),
            A: value * 100,
            fullMark: 100
        })) : []

    const handleDownload = async (format) => {
        if (!selectedJob) return
        try {
            const res = await fetch(`${API_BASE}/analysis/export/${selectedJob}?format=${format}`)
            if (res.ok) {
                const data = await res.json()
                const blob = new Blob([data.content], { type: format === 'csv' ? 'text/csv' : 'application/json' })
                const url = window.URL.createObjectURL(blob)
                const a = document.createElement('a')
                a.href = url
                a.download = data.filename
                document.body.appendChild(a)
                a.click()
                window.URL.revokeObjectURL(url)
                document.body.removeChild(a)
            }
        } catch (error) {
            console.error('Download failed:', error)
        }
    }

    return (
        <div className="results-page">
            <header className="page-header">
                <div>
                    <h1 className="page-title">Analysis Results</h1>
                    <p className="page-subtitle">
                        View and interpret hydrogeochemical modeling outputs
                    </p>
                </div>
                {selectedJob && (
                    <div className="header-actions">
                        <button className="btn btn-secondary" onClick={() => handleDownload('json')}>
                            📥 JSON
                        </button>
                        <button className="btn btn-secondary" onClick={() => handleDownload('csv')}>
                            📥 CSV
                        </button>
                    </div>
                )}
            </header>

            <div className="results-layout">
                {/* Jobs List */}
                <aside className="jobs-panel card">
                    <div className="panel-header">
                        <h3>Completed Analyses</h3>
                    </div>

                    {loading ? (
                        <div className="loading-state">
                            <div className="spinner"></div>
                        </div>
                    ) : jobs.filter(j => j.status === 'completed').length === 0 ? (
                        <div className="empty-state">
                            <span className="empty-icon">📊</span>
                            <p>No completed analyses</p>
                        </div>
                    ) : (
                        <div className="jobs-list">
                            {jobs.filter(j => j.status === 'completed').map(job => (
                                <div
                                    key={job.job_id}
                                    className={`job-item ${selectedJob === job.job_id ? 'active' : ''}`}
                                    onClick={() => selectJob(job.job_id)}
                                >
                                    <div className="job-info">
                                        <span className="job-name">{job.name}</span>
                                        <span className="job-type">{job.analysis_type}</span>
                                    </div>
                                    <button
                                        className="delete-btn"
                                        onClick={(e) => {
                                            e.stopPropagation()
                                            deleteJob(job.job_id)
                                        }}
                                    >
                                        🗑️
                                    </button>
                                </div>
                            ))}
                        </div>
                    )}
                </aside>

                {/* Results Content */}
                <main className="results-content">
                    {results ? (
                        <>
                            {/* Summary Cards */}
                            <div className="summary-grid">
                                <div className="summary-card">
                                    <span className="summary-icon">🧪</span>
                                    <div className="summary-info">
                                        <span className="summary-value">{results.summary?.total_samples || 0}</span>
                                        <span className="summary-label">Samples Analyzed</span>
                                    </div>
                                </div>
                                <div className="summary-card">
                                    <span className="summary-icon">🔀</span>
                                    <div className="summary-info">
                                        <span className="summary-value highlight">
                                            {results.transport_model?.dominant_process || 'N/A'}
                                        </span>
                                        <span className="summary-label">Dominant Process</span>
                                    </div>
                                </div>
                                <div className="summary-card">
                                    <span className="summary-icon">⚛️</span>
                                    <div className="summary-info">
                                        <span className="summary-value">{results.reactions?.length || 0}</span>
                                        <span className="summary-label">Active Reactions</span>
                                    </div>
                                </div>
                                <div className="summary-card">
                                    <span className="summary-icon">⏳</span>
                                    <div className="summary-info">
                                        <span className="summary-value">
                                            {results.nuclear_aging?.average_age_years ? 
                                                `${results.nuclear_aging.average_age_years.toFixed(1)}y` : 'N/A'}
                                        </span>
                                        <span className="summary-label">Average Age</span>
                                    </div>
                                </div>
                                <div className="summary-card">
                                    <span className="summary-icon">📈</span>
                                    <div className="summary-info">
                                        <span className="summary-value">{results.uncertainty?.confidence_level ? `${(results.uncertainty.confidence_level * 100).toFixed(0)}%` : 'N/A'}</span>
                                        <span className="summary-label">Confidence Level</span>
                                    </div>
                                </div>
                            </div>

                            {/* Scientific Plots Gallery */}
                            {results.plots && Object.keys(results.plots).length > 0 && (
                                <div className="card mt-lg plots-gallery">
                                    <div className="card-header">
                                        <h3 className="card-title">Scientific Visualization</h3>
                                        <span className="badge badge-primary">Generated</span>
                                    </div>
                                    <div className="plots-grid">
                                        {Object.entries(results.plots).map(([key, url]) => (
                                            <div key={key} className="plot-item">
                                                <h4>{key.replace('_', ' ').toUpperCase()}</h4>
                                                <a href={`${API_BASE}${url}`} target="_blank" rel="noopener noreferrer">
                                                    <img 
                                                        src={`${API_BASE}${url}`} 
                                                        alt={`${key} plot`} 
                                                        className="scientific-plot"
                                                        onError={(e) => {e.target.style.display='none'}}
                                                    />
                                                </a>
                                            </div>
                                        ))}
                                    </div>
                                </div>
                            )}

                             {/* Transport Model */}
                             <div className="card mt-lg">
                                 {/* ... existing transport model code ... */}
                             </div>

                             {/* Nuclear Dating Results */}
                             {results.nuclear_aging && (
                                 <div className="card mt-lg">
                                     <div className="card-header">
                                         <h3 className="card-title">Network-Enhanced Nuclear Dating</h3>
                                         <span className="badge badge-success">Bayesian MCMC</span>
                                     </div>
                                     <div className="nuclear-results-grid">
                                         <div className="table-container">
                                             <table className="nuclear-table">
                                                 <thead>
                                                     <tr>
                                                         <th>Node</th>
                                                         <th>Mean Age (y)</th>
                                                         <th>95% CI (y)</th>
                                                         <th>Modern Prob.</th>
                                                         <th>Status</th>
                                                     </tr>
                                                 </thead>
                                                 <tbody>
                                                     {Object.entries(results.nuclear_aging.nodes || {}).map(([id, data]) => (
                                                         <tr key={id}>
                                                             <td className="node-id">{id}</td>
                                                             <td>{data.mean_age_years.toFixed(1)}</td>
                                                             <td className="text-muted text-small">
                                                                 [{data.age_95_low.toFixed(1)}, {data.age_95_high.toFixed(1)}]
                                                             </td>
                                                             <td>
                                                                 <div className="prob-pill-container">
                                                                     <div 
                                                                         className="prob-pill-fill" 
                                                                         style={{ 
                                                                             width: `${data.p_modern * 100}%`,
                                                                             backgroundColor: data.p_modern > 0.7 ? '#22c55e' : (data.p_modern > 0.3 ? '#eab308' : '#ef4444')
                                                                         }}
                                                                     ></div>
                                                                     <span className="prob-text">{(data.p_modern * 100).toFixed(0)}%</span>
                                                                 </div>
                                                             </td>
                                                             <td>
                                                                 <span className={`badge badge-${data.p_modern > 0.5 ? 'success' : 'warning'}`}>
                                                                     {data.p_modern > 0.5 ? 'Modern' : 'Sub-modern'}
                                                                 </span>
                                                             </td>
                                                         </tr>
                                                     ))}
                                                 </tbody>
                                             </table>
                                         </div>
                                         {results.nuclear_aging.diagnostics && (
                                             <div className="nuclear-diagnostics mt-md">
                                                 <h4>Inferred Global Parameters</h4>
                                                 <div className="diag-metrics">
                                                     <div className="diag-item">
                                                         <span className="label">Input Scaling</span>
                                                         <span className="value">{results.nuclear_aging.diagnostics.inferred_input_scale.toFixed(2)}</span>
                                                     </div>
                                                     <div className="diag-item">
                                                         <span className="label">Mean Velocity</span>
                                                         <span className="value">{results.nuclear_aging.diagnostics.inferred_velocity_m_y.toFixed(1)} m/y</span>
                                                     </div>
                                                     {results.nuclear_aging.diagnostics.inferred_c14_a0_global > 0 && (
                                                         <div className="diag-item">
                                                             <span className="label">C14 A0 (Global)</span>
                                                             <span className="value">{results.nuclear_aging.diagnostics.inferred_c14_a0_global.toFixed(1)} pmc</span>
                                                         </div>
                                                     )}
                                                 </div>
                                             </div>
                                         )}
                                     </div>
                                 </div>
                             )}


                            {/* Edge Network Visualization */}
                            {results.edge_results && results.edge_results.length > 0 && (
                                <div className="card mt-lg">
                                    <div className="card-header">
                                        <h3 className="card-title">Flow Network</h3>
                                        <span className="badge badge-primary">
                                            {results.edge_results.length} Edges
                                        </span>
                                    </div>
                                    <ResultsGraph
                                        edgeResults={results.edge_results}
                                        width={700}
                                        height={400}
                                    />
                                </div>
                            )}

                            {/* Reactions Chart */}
                            <div className="card mt-lg">
                                <div className="card-header">
                                    <h3 className="card-title">Reaction Rates</h3>
                                    <span className="badge badge-success">LASSO Fit</span>
                                </div>
                                <div className="chart-container">
                                    <ResponsiveContainer width="100%" height={300}>
                                        <BarChart data={reactionData} layout="vertical">
                                            <CartesianGrid strokeDasharray="3 3" stroke="rgba(148, 163, 184, 0.1)" />
                                            <XAxis type="number" stroke="#94a3b8" fontSize={12} />
                                            <YAxis type="category" dataKey="name" stroke="#94a3b8" fontSize={12} width={80} />
                                            <Tooltip
                                                contentStyle={{
                                                    background: 'rgba(17, 24, 39, 0.95)',
                                                    border: '1px solid rgba(148, 163, 184, 0.2)',
                                                    borderRadius: '8px',
                                                    color: '#f8fafc'
                                                }}
                                                formatter={(value, name, props) => [
                                                    `${value.toFixed(3)} mmol/L (${props.payload.direction})`,
                                                    'Rate'
                                                ]}
                                            />
                                            <Bar
                                                dataKey="rate"
                                                fill="#0ea5e9"
                                                radius={[0, 4, 4, 0]}
                                            />
                                        </BarChart>
                                    </ResponsiveContainer>
                                </div>
                            </div>

                            {/* Reactions Table */}
                            <div className="card mt-lg">
                                <div className="card-header">
                                    <h3 className="card-title">Mineral Reactions</h3>
                                </div>
                                <div className="table-container">
                                    <table>
                                        <thead>
                                            <tr>
                                                <th>Mineral</th>
                                                <th>Rate (mmol/L)</th>
                                                <th>Direction</th>
                                                <th>Uncertainty</th>
                                            </tr>
                                        </thead>
                                        <tbody>
                                            {results.reactions?.map((reaction, idx) => (
                                                <tr key={idx}>
                                                    <td>{reaction.mineral}</td>
                                                    <td className={reaction.rate_mmol_L >= 0 ? 'text-success' : 'text-error'}>
                                                        {reaction.rate_mmol_L >= 0 ? '+' : ''}{reaction.rate_mmol_L.toFixed(3)}
                                                    </td>
                                                    <td>
                                                        <span className={`badge badge-${reaction.direction === 'dissolution' ? 'success' : 'warning'}`}>
                                                            {reaction.direction}
                                                        </span>
                                                    </td>
                                                    <td>
                                                        {results.uncertainty?.reaction_uncertainties?.[reaction.mineral] ? (
                                                            <span className="uncertainty-range">
                                                                [{results.uncertainty.reaction_uncertainties[reaction.mineral].lower.toFixed(3)},
                                                                {results.uncertainty.reaction_uncertainties[reaction.mineral].upper.toFixed(3)}]
                                                            </span>
                                                        ) : '-'}
                                                    </td>
                                                </tr>
                                            ))}
                                        </tbody>
                                    </table>
                                </div>
                            </div>

                            {/* Network Flow */}
                            {results.network_inference && (
                                <div className="card mt-lg">
                                    <div className="card-header">
                                        <h3 className="card-title">Flow Network Probabilities</h3>
                                    </div>
                                    <div className="flow-list">
                                        {results.network_inference.flow_probabilities?.map((flow, idx) => (
                                            <div key={idx} className="flow-item">
                                                <div className="flow-path">
                                                    <span className="flow-node">{flow.from}</span>
                                                    <span className="flow-arrow">→</span>
                                                    <span className="flow-node">{flow.to}</span>
                                                </div>
                                                <div className="flow-probability">
                                                    <div className="prob-bar">
                                                        <div
                                                            className="prob-fill"
                                                            style={{ width: `${flow.probability * 100}%` }}
                                                        ></div>
                                                    </div>
                                                    <span className="prob-value">{(flow.probability * 100).toFixed(1)}%</span>
                                                </div>
                                            </div>
                                        ))}
                                    </div>
                                </div>
                            )}
                        </>
                    ) : (
                        <div className="card empty-content">
                            <div className="empty-state large">
                                <span className="empty-icon">📊</span>
                                <h3>Select an Analysis</h3>
                                <p>Choose a completed analysis from the sidebar to view results</p>
                            </div>
                        </div>
                    )}
                </main>
            </div>
        </div>
    )
}

export default Results
