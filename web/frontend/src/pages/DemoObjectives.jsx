import { useState, useEffect } from 'react'
import {
    LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer,
    BarChart, Bar, ScatterChart, Scatter, ReferenceLine
} from 'recharts'
import './DemoObjectives.css'

const API_BASE = 'http://localhost:8000/api/demo'

function DemoObjectives() {
    const [activeTab, setActiveTab] = useState(1)
    const [loading, setLoading] = useState(false)
    const [data, setData] = useState(null)
    const [error, setError] = useState(null)

    useEffect(() => {
        fetchData(activeTab)
    }, [activeTab])

    const fetchData = async (objectiveId) => {
        setLoading(true)
        setError(null)
        setData(null)
        
        let endpoint = ''
        switch (objectiveId) {
            case 1: endpoint = '/objectives/1/vadose'; break;
            case 2: endpoint = '/objectives/2/sources'; break;
            case 3: endpoint = '/objectives/3/recharge'; break;
            case 4: endpoint = '/objectives/4/transport'; break;
            default: endpoint = '';
        }

        try {
            const response = await fetch(`${API_BASE}${endpoint}`)
            if (!response.ok) {
                throw new Error(`Error fetching data: ${response.statusText}`)
            }
            const result = await response.json()
            setData(result)
        } catch (err) {
            setError(err.message)
        } finally {
            setLoading(false)
        }
    }

    const renderObjective1 = () => {
        if (!data) return null
        return (
            <div className="objective-content">
                <div className="objective-header">
                    <h2>Vadose Zone Nitrate Transport</h2>
                    <p>Quantifying nitrate leaching dynamics and seasonal patterns.</p>
                </div>

                <div className="objective-grid">
                    <div className="card">
                        <h3 className="card-title">Seasonal Statistics</h3>
                        <div className="stat-grid">
                            {data.seasonal_stats.map((stat, idx) => (
                                <div key={idx} className="stat-item">
                                    <span className="stat-value">{stat.mean.toFixed(1)}</span>
                                    <span className="stat-label">{stat.season} Season (mg/L)</span>
                                </div>
                            ))}
                        </div>
                        <div style={{ marginTop: '1rem', height: 300 }}>
                            <ResponsiveContainer width="100%" height="100%">
                                <BarChart data={data.seasonal_stats}>
                                    <CartesianGrid strokeDasharray="3 3" />
                                    <XAxis dataKey="season" />
                                    <YAxis label={{ value: 'NO3 (mg/L)', angle: -90, position: 'insideLeft' }} />
                                    <Tooltip />
                                    <Bar dataKey="mean" fill="#8884d8" name="Mean NO3" />
                                </BarChart>
                            </ResponsiveContainer>
                        </div>
                    </div>

                    <div className="card">
                        <h3 className="card-title">Time Series Analysis</h3>
                        <div style={{ height: 300 }}>
                            <ResponsiveContainer width="100%" height="100%">
                                <LineChart data={data.time_series}>
                                    <CartesianGrid strokeDasharray="3 3" />
                                    <XAxis dataKey="date" tickFormatter={(str) => str.substring(5)} />
                                    <YAxis />
                                    <Tooltip />
                                    <Legend />
                                    <Line type="monotone" dataKey="no3" stroke="#82ca9d" name="NO3 Concentration" dot={{ r: 4 }} />
                                </LineChart>
                            </ResponsiveContainer>
                        </div>
                    </div>
                </div>
            </div>
        )
    }

    const renderObjective2 = () => {
        if (!data) return null
        return (
            <div className="objective-content">
                <div className="objective-header">
                    <h2>Nitrate Source Discrimination</h2>
                    <p>Identifying contamination sources using dual isotopes (δ¹⁵N, δ¹⁸O).</p>
                </div>

                <div className="objective-grid">
                    <div className="card">
                        <h3 className="card-title">Endmember Mixing Plot</h3>
                        <div style={{ height: 400 }}>
                            <ResponsiveContainer width="100%" height="100%">
                                <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
                                    <CartesianGrid strokeDasharray="3 3" />
                                    <XAxis type="number" dataKey="d15N" name="δ15N" unit="‰" domain={[-5, 25]} label={{ value: 'δ15N-NO3 (‰)', position: 'bottom', offset: 0 }} />
                                    <YAxis type="number" dataKey="d18O" name="δ18O" unit="‰" domain={[-5, 70]} label={{ value: 'δ18O-NO3 (‰)', angle: -90, position: 'insideLeft' }} />
                                    <Tooltip cursor={{ strokeDasharray: '3 3' }} />
                                    <Legend />
                                    
                                    {/* Sources */}
                                    <Scatter name="Sources" data={data.sources} fill="#ff7300" shape="star" />
                                    
                                    {/* Samples */}
                                    <Scatter name="Samples" data={data.samples} fill="#387908" />
                                </ScatterChart>
                            </ResponsiveContainer>
                        </div>
                    </div>

                    <div className="card">
                        <h3 className="card-title">Source Contributions</h3>
                        <table className="data-table">
                            <thead>
                                <tr>
                                    <th>Station</th>
                                    <th>Dominant Source</th>
                                    <th>Probability</th>
                                </tr>
                            </thead>
                            <tbody>
                                {data.samples.map((sample, idx) => {
                                    // Find max probability source
                                    const bestSource = Object.entries(sample.probabilities).reduce((a, b) => a[1] > b[1] ? a : b);
                                    return (
                                        <tr key={idx}>
                                            <td>{sample.station}</td>
                                            <td>{bestSource[0]}</td>
                                            <td>{(bestSource[1] * 100).toFixed(1)}%</td>
                                        </tr>
                                    )
                                })}
                            </tbody>
                        </table>
                        <div style={{ marginTop: '1rem', fontSize: '0.9rem', color: '#666' }}>
                            MCMC Analysis Available: {data.mcmc_available ? '✅ Yes' : '❌ No'}
                        </div>
                    </div>
                </div>
            </div>
        )
    }

    const renderObjective3 = () => {
        if (!data) return null
        return (
            <div className="objective-content">
                <div className="objective-header">
                    <h2>Groundwater Recharge Tracing</h2>
                    <p>Characterizing recharge mechanisms via stable water isotopes.</p>
                </div>

                <div className="objective-grid">
                    <div className="card" style={{ gridColumn: 'span 2' }}>
                        <h3 className="card-title">Dual Isotope Plot (δ²H vs δ¹⁸O)</h3>
                        <div style={{ height: 500 }}>
                            <ResponsiveContainer width="100%" height="100%">
                                <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
                                    <CartesianGrid strokeDasharray="3 3" />
                                    <XAxis type="number" dataKey="d18O" name="δ18O" unit="‰" domain={[-10, 5]} label={{ value: 'δ18O (‰)', position: 'bottom', offset: 0 }} />
                                    <YAxis type="number" dataKey="d2H" name="δ2H" unit="‰" domain={[-60, 20]} label={{ value: 'δ2H (‰)', angle: -90, position: 'insideLeft' }} />
                                    <Tooltip cursor={{ strokeDasharray: '3 3' }} />
                                    <Legend />
                                    
                                    {/* Reference Lines - drawn as scatters for simplicity or using ReferenceLine if fixed */}
                                    <Scatter name="Samples" data={data.samples} fill="#8884d8" />
                                    <Line dataKey="y" data={data.gmwl} stroke="#ff7300" name="GMWL" dot={false} isAnimationActive={false} />
                                    <Line dataKey="y" data={data.lmwl} stroke="#82ca9d" name="LMWL" dot={false} isAnimationActive={false} />
                                </ScatterChart>
                            </ResponsiveContainer>
                        </div>
                    </div>

                    <div className="card">
                        <h3 className="card-title">Recharge Metrics</h3>
                        <div className="stat-grid">
                            <div className="stat-item">
                                <span className="stat-value">{data.samples.filter(s => s.is_evaporated).length}</span>
                                <span className="stat-label">Evaporated Samples</span>
                            </div>
                            <div className="stat-item">
                                <span className="stat-value">{data.samples.length}</span>
                                <span className="stat-label">Total Samples</span>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        )
    }

    const renderObjective4 = () => {
        if (!data) return null
        
        if (!data.available) {
            return (
                <div className="error-message">
                    FloPy is not installed or configured on the server. Please install it to view transport modeling results.
                </div>
            )
        }

        const chartData = data.output_times.map((t, i) => ({
            time: t,
            concentration: data.output_conc[i]
        }));

        return (
            <div className="objective-content">
                <div className="objective-header">
                    <h2>Saturated Zone Transport Modeling</h2>
                    <p>Coupled Vadose-Saturated transport simulation using FloPy/MT3DMS.</p>
                </div>

                <div className="objective-grid">
                    <div className="card">
                        <h3 className="card-title">Breakthrough Curve</h3>
                        <div style={{ height: 300 }}>
                            <ResponsiveContainer width="100%" height="100%">
                                <LineChart data={chartData}>
                                    <CartesianGrid strokeDasharray="3 3" />
                                    <XAxis dataKey="time" label={{ value: 'Time (days)', position: 'bottom', offset: 0 }} />
                                    <YAxis label={{ value: 'Concentration (mg/L)', angle: -90, position: 'insideLeft' }} />
                                    <Tooltip />
                                    <Line type="monotone" dataKey="concentration" stroke="#8884d8" dot={false} strokeWidth={2} />
                                </LineChart>
                            </ResponsiveContainer>
                        </div>
                    </div>

                    <div className="card">
                        <h3 className="card-title">Model Metrics</h3>
                        <div className="stat-grid">
                            <div className="stat-item">
                                <span className="stat-value">{data.metrics.peak_conc.toFixed(1)}</span>
                                <span className="stat-label">Peak Conc. (mg/L)</span>
                            </div>
                            <div className="stat-item">
                                <span className="stat-value">{data.metrics.travel_time.toFixed(0)}</span>
                                <span className="stat-label">Travel Time (days)</span>
                            </div>
                            <div className="stat-item">
                                <span className="stat-value">{(data.metrics.attenuation * 100).toFixed(1)}%</span>
                                <span className="stat-label">Attenuation</span>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        )
    }

    return (
        <div className="demo-container">
            <h1>Research Objectives Demo</h1>
            
            <div className="demo-tabs">
                <button className={`demo-tab ${activeTab === 1 ? 'active' : ''}`} onClick={() => setActiveTab(1)}>
                    1. Vadose Transport
                </button>
                <button className={`demo-tab ${activeTab === 2 ? 'active' : ''}`} onClick={() => setActiveTab(2)}>
                    2. Source ID
                </button>
                <button className={`demo-tab ${activeTab === 3 ? 'active' : ''}`} onClick={() => setActiveTab(3)}>
                    3. Recharge Tracing
                </button>
                <button className={`demo-tab ${activeTab === 4 ? 'active' : ''}`} onClick={() => setActiveTab(4)}>
                    4. Transport Modeling
                </button>
            </div>

            {error && <div className="error-message">{error}</div>}
            
            {loading ? (
                <div className="loading-spinner">Loading analysis data...</div>
            ) : (
                <>
                    {activeTab === 1 && renderObjective1()}
                    {activeTab === 2 && renderObjective2()}
                    {activeTab === 3 && renderObjective3()}
                    {activeTab === 4 && renderObjective4()}
                </>
            )}
        </div>
    )
}

export default DemoObjectives
