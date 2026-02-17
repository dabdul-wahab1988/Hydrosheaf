import { useState, useEffect } from 'react'
import { Link } from 'react-router-dom'
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, LineChart, Line, PieChart, Pie, Cell } from 'recharts'
import { API_BASE } from '../config'
import { useToast } from '../context/ToastContext'
import { StatCardSkeleton, ChartSkeleton } from '../components/Skeleton'
import './Dashboard.css'

function Dashboard() {
    const [datasets, setDatasets] = useState([])
    const [jobs, setJobs] = useState([])
    const [loading, setLoading] = useState(true)
    const toast = useToast()

    useEffect(() => {
        fetchData()
    }, [])

    const fetchData = async () => {
        try {
            const [datasetsRes, jobsRes] = await Promise.all([
                fetch(`${API_BASE}/samples/datasets`),
                fetch(`${API_BASE}/analysis/jobs`)
            ])

            if (datasetsRes.ok) {
                setDatasets(await datasetsRes.json())
            }
            if (jobsRes.ok) {
                setJobs(await jobsRes.json())
            }
        } catch (error) {
            toast.error('Failed to load dashboard data')
        } finally {
            setLoading(false)
        }
    }

    const totalSamples = datasets.reduce((sum, ds) => sum + ds.sample_count, 0)
    const completedJobs = jobs.filter(j => j.status === 'completed').length

    // Demo chart data
    const ionData = [
        { name: 'Ca²⁺', value: 85 },
        { name: 'Mg²⁺', value: 32 },
        { name: 'Na⁺', value: 45 },
        { name: 'K⁺', value: 5 },
        { name: 'HCO₃⁻', value: 245 },
        { name: 'SO₄²⁻', value: 78 },
        { name: 'Cl⁻', value: 52 },
    ]

    const isotopeData = [
        { d18O: -6.5, d2H: -45 },
        { d18O: -5.8, d2H: -40 },
        { d18O: -5.2, d2H: -35 },
        { d18O: -4.8, d2H: -32 },
        { d18O: -5.5, d2H: -37 },
        { d18O: -4.5, d2H: -30 },
    ]

    const processData = [
        { name: 'Dissolution', value: 45 },
        { name: 'Mixing', value: 30 },
        { name: 'Evaporation', value: 15 },
        { name: 'Ion Exchange', value: 10 },
    ]

    const COLORS = ['#0ea5e9', '#06b6d4', '#a855f7', '#22c55e']

    return (
        <div className="dashboard">
            {/* Header */}
            <header className="page-header">
                <div>
                    <h1 className="page-title">Dashboard</h1>
                    <p className="page-subtitle">
                        Groundwater Hydrogeochemistry Analysis Platform
                    </p>
                </div>
                <Link to="/analysis" className="btn btn-primary">
                    <span>⚗️</span> New Analysis
                </Link>
            </header>

            {/* Stats Grid */}
            <div className="stats-grid">
                {loading ? (
                    <>
                        <StatCardSkeleton />
                        <StatCardSkeleton />
                        <StatCardSkeleton />
                        <StatCardSkeleton />
                    </>
                ) : (
                    <>
                        <div className="stat-card animate-fadeIn">
                            <div className="stat-icon">🧪</div>
                            <div className="stat-value">{totalSamples}</div>
                            <div className="stat-label">Water Samples</div>
                        </div>
                        <div className="stat-card animate-fadeIn" style={{ animationDelay: '0.1s' }}>
                            <div className="stat-icon">📂</div>
                            <div className="stat-value">{datasets.length}</div>
                            <div className="stat-label">Datasets</div>
                        </div>
                        <div className="stat-card animate-fadeIn" style={{ animationDelay: '0.2s' }}>
                            <div className="stat-icon">⚗️</div>
                            <div className="stat-value">{jobs.length}</div>
                            <div className="stat-label">Analysis Jobs</div>
                        </div>
                        <div className="stat-card animate-fadeIn" style={{ animationDelay: '0.3s' }}>
                            <div className="stat-icon">✅</div>
                            <div className="stat-value">{completedJobs}</div>
                            <div className="stat-label">Completed</div>
                        </div>
                    </>
                )}
            </div>

            {/* Charts Grid */}
            <div className="charts-grid">
                {/* Ion Concentrations */}
                <div className="card chart-card">
                    <div className="card-header">
                        <h3 className="card-title">Major Ion Distribution</h3>
                        <span className="badge badge-primary">mg/L</span>
                    </div>
                    <div className="chart-container">
                        <ResponsiveContainer width="100%" height={300}>
                            <BarChart data={ionData}>
                                <CartesianGrid strokeDasharray="3 3" stroke="rgba(148, 163, 184, 0.1)" />
                                <XAxis dataKey="name" stroke="#94a3b8" fontSize={12} />
                                <YAxis stroke="#94a3b8" fontSize={12} />
                                <Tooltip
                                    contentStyle={{
                                        background: 'rgba(17, 24, 39, 0.95)',
                                        border: '1px solid rgba(148, 163, 184, 0.2)',
                                        borderRadius: '8px',
                                        color: '#f8fafc'
                                    }}
                                />
                                <Bar
                                    dataKey="value"
                                    fill="url(#colorGradient)"
                                    radius={[4, 4, 0, 0]}
                                />
                                <defs>
                                    <linearGradient id="colorGradient" x1="0" y1="0" x2="0" y2="1">
                                        <stop offset="0%" stopColor="#0ea5e9" />
                                        <stop offset="100%" stopColor="#06b6d4" />
                                    </linearGradient>
                                </defs>
                            </BarChart>
                        </ResponsiveContainer>
                    </div>
                </div>

                {/* Isotope Plot */}
                <div className="card chart-card">
                    <div className="card-header">
                        <h3 className="card-title">Stable Isotopes (δ¹⁸O vs δ²H)</h3>
                        <span className="badge badge-success">LMWL</span>
                    </div>
                    <div className="chart-container">
                        <ResponsiveContainer width="100%" height={300}>
                            <LineChart data={isotopeData}>
                                <CartesianGrid strokeDasharray="3 3" stroke="rgba(148, 163, 184, 0.1)" />
                                <XAxis
                                    dataKey="d18O"
                                    stroke="#94a3b8"
                                    fontSize={12}
                                    label={{ value: 'δ¹⁸O (‰)', position: 'bottom', fill: '#94a3b8' }}
                                />
                                <YAxis
                                    dataKey="d2H"
                                    stroke="#94a3b8"
                                    fontSize={12}
                                    label={{ value: 'δ²H (‰)', angle: -90, position: 'left', fill: '#94a3b8' }}
                                />
                                <Tooltip
                                    contentStyle={{
                                        background: 'rgba(17, 24, 39, 0.95)',
                                        border: '1px solid rgba(148, 163, 184, 0.2)',
                                        borderRadius: '8px',
                                        color: '#f8fafc'
                                    }}
                                />
                                <Line
                                    type="monotone"
                                    dataKey="d2H"
                                    stroke="#a855f7"
                                    strokeWidth={2}
                                    dot={{ fill: '#a855f7', strokeWidth: 2 }}
                                />
                            </LineChart>
                        </ResponsiveContainer>
                    </div>
                </div>

                {/* Processes Pie Chart */}
                <div className="card chart-card">
                    <div className="card-header">
                        <h3 className="card-title">Dominant Processes</h3>
                    </div>
                    <div className="chart-container">
                        <ResponsiveContainer width="100%" height={300}>
                            <PieChart>
                                <Pie
                                    data={processData}
                                    cx="50%"
                                    cy="50%"
                                    innerRadius={60}
                                    outerRadius={100}
                                    paddingAngle={5}
                                    dataKey="value"
                                    label={({ name, percent }) => `${name} ${(percent * 100).toFixed(0)}%`}
                                    labelLine={{ stroke: '#94a3b8' }}
                                >
                                    {processData.map((entry, index) => (
                                        <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                                    ))}
                                </Pie>
                                <Tooltip
                                    contentStyle={{
                                        background: 'rgba(17, 24, 39, 0.95)',
                                        border: '1px solid rgba(148, 163, 184, 0.2)',
                                        borderRadius: '8px',
                                        color: '#f8fafc'
                                    }}
                                />
                            </PieChart>
                        </ResponsiveContainer>
                    </div>
                </div>
            </div>

            {/* Quick Actions */}
            <div className="quick-actions">
                <h2>Quick Actions</h2>
                <div className="actions-grid">
                    <Link to="/samples" className="action-card">
                        <span className="action-icon">📤</span>
                        <span className="action-title">Upload Samples</span>
                        <span className="action-desc">Import water chemistry data</span>
                    </Link>
                    <Link to="/network" className="action-card">
                        <span className="action-icon">🔗</span>
                        <span className="action-title">Create Network</span>
                        <span className="action-desc">Define flow connections</span>
                    </Link>
                    <Link to="/analysis" className="action-card">
                        <span className="action-icon">⚗️</span>
                        <span className="action-title">Run Analysis</span>
                        <span className="action-desc">Start inverse modeling</span>
                    </Link>
                    <Link to="/results" className="action-card">
                        <span className="action-icon">📈</span>
                        <span className="action-title">View Results</span>
                        <span className="action-desc">Explore analysis outputs</span>
                    </Link>
                </div>
            </div>

            {/* Recent Jobs */}
            {jobs.length > 0 && (
                <div className="card mt-lg">
                    <div className="card-header">
                        <h3 className="card-title">Recent Analysis Jobs</h3>
                        <Link to="/results" className="btn btn-secondary">View All</Link>
                    </div>
                    <div className="table-container">
                        <table>
                            <thead>
                                <tr>
                                    <th>Name</th>
                                    <th>Type</th>
                                    <th>Status</th>
                                    <th>Created</th>
                                </tr>
                            </thead>
                            <tbody>
                                {jobs.slice(0, 5).map((job) => (
                                    <tr key={job.job_id}>
                                        <td>{job.name}</td>
                                        <td>{job.analysis_type}</td>
                                        <td>
                                            <span className={`badge badge-${job.status === 'completed' ? 'success' : job.status === 'failed' ? 'error' : 'warning'}`}>
                                                {job.status}
                                            </span>
                                        </td>
                                        <td>{new Date(job.created_at).toLocaleString()}</td>
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>
                </div>
            )}
        </div>
    )
}

export default Dashboard
