import { useState } from 'react'
import './Network.css'

const API_BASE = 'http://localhost:8000/api'

function Network() {
    const [nodes, setNodes] = useState([
        { id: 'well_a', name: 'Well A', x: 100, y: 150, hydraulic_head: 250 },
        { id: 'well_b', name: 'Well B', x: 300, y: 100, hydraulic_head: 245 },
        { id: 'well_c', name: 'Well C', x: 250, y: 300, hydraulic_head: 238 },
    ])
    const [edges, setEdges] = useState([
        { source: 'well_a', target: 'well_b' },
        { source: 'well_b', target: 'well_c' },
    ])
    const [newNode, setNewNode] = useState({ id: '', name: '', hydraulic_head: '' })
    const [selectedNode, setSelectedNode] = useState(null)
    const [inferenceResult, setInferenceResult] = useState(null)

    const addNode = () => {
        if (!newNode.id || !newNode.name) return

        setNodes([...nodes, {
            ...newNode,
            x: 200 + Math.random() * 200,
            y: 150 + Math.random() * 150,
            hydraulic_head: parseFloat(newNode.hydraulic_head) || null
        }])
        setNewNode({ id: '', name: '', hydraulic_head: '' })
    }

    const removeNode = (id) => {
        setNodes(nodes.filter(n => n.id !== id))
        setEdges(edges.filter(e => e.source !== id && e.target !== id))
        if (selectedNode === id) setSelectedNode(null)
    }

    const addEdge = (source, target) => {
        if (!edges.find(e => e.source === source && e.target === target)) {
            setEdges([...edges, { source, target }])
        }
    }

    const runInference = async () => {
        try {
            // Create network first
            await fetch(`${API_BASE}/network/create`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ name: 'flow_network', nodes, edges })
            })

            // Run inference
            const res = await fetch(`${API_BASE}/network/flow_network/infer-flow`, {
                method: 'POST'
            })

            if (res.ok) {
                setInferenceResult(await res.json())
            }
        } catch (error) {
            console.error('Inference error:', error)
        }
    }

    return (
        <div className="network-page">
            <header className="page-header">
                <div>
                    <h1 className="page-title">Flow Network</h1>
                    <p className="page-subtitle">
                        Define aquifer connections and infer flow directions
                    </p>
                </div>
                <button className="btn btn-primary" onClick={runInference}>
                    <span>🔮</span> Infer Flow Directions
                </button>
            </header>

            <div className="network-layout">
                {/* Network Canvas */}
                <div className="card network-canvas-card">
                    <div className="card-header">
                        <h3 className="card-title">Network Visualization</h3>
                        <span className="badge badge-primary">{nodes.length} Nodes</span>
                    </div>

                    <div className="network-canvas">
                        <svg width="100%" height="400" viewBox="0 0 600 400">
                            {/* Grid background */}
                            <defs>
                                <pattern id="grid" width="50" height="50" patternUnits="userSpaceOnUse">
                                    <path d="M 50 0 L 0 0 0 50" fill="none" stroke="rgba(148,163,184,0.1)" strokeWidth="1" />
                                </pattern>
                            </defs>
                            <rect width="100%" height="100%" fill="url(#grid)" />

                            {/* Edges */}
                            {edges.map((edge, idx) => {
                                const sourceNode = nodes.find(n => n.id === edge.source)
                                const targetNode = nodes.find(n => n.id === edge.target)
                                if (!sourceNode || !targetNode) return null

                                const inference = inferenceResult?.inferred_edges?.find(
                                    e => e.source === edge.source && e.target === edge.target
                                )

                                return (
                                    <g key={idx}>
                                        <line
                                            x1={sourceNode.x}
                                            y1={sourceNode.y}
                                            x2={targetNode.x}
                                            y2={targetNode.y}
                                            stroke={inference ?
                                                inference.flow_direction === 'forward' ? '#22c55e' :
                                                    inference.flow_direction === 'reverse' ? '#ef4444' : '#94a3b8'
                                                : '#0ea5e9'}
                                            strokeWidth="3"
                                            strokeOpacity="0.7"
                                        />
                                        {inference && (
                                            <text
                                                x={(sourceNode.x + targetNode.x) / 2}
                                                y={(sourceNode.y + targetNode.y) / 2 - 10}
                                                fill="#94a3b8"
                                                fontSize="12"
                                                textAnchor="middle"
                                            >
                                                P={inference.flow_probability?.toFixed(2) || '?'}
                                            </text>
                                        )}
                                    </g>
                                )
                            })}

                            {/* Nodes */}
                            {nodes.map((node) => (
                                <g
                                    key={node.id}
                                    className={`network-node ${selectedNode === node.id ? 'selected' : ''}`}
                                    onClick={() => setSelectedNode(selectedNode === node.id ? null : node.id)}
                                >
                                    <circle
                                        cx={node.x}
                                        cy={node.y}
                                        r="25"
                                        fill="url(#nodeGradient)"
                                        stroke={selectedNode === node.id ? '#a855f7' : '#0ea5e9'}
                                        strokeWidth="2"
                                        style={{ cursor: 'pointer' }}
                                    />
                                    <text
                                        x={node.x}
                                        y={node.y + 5}
                                        textAnchor="middle"
                                        fill="#f8fafc"
                                        fontSize="12"
                                        fontWeight="600"
                                    >
                                        {node.name.slice(0, 5)}
                                    </text>
                                    {node.hydraulic_head && (
                                        <text
                                            x={node.x}
                                            y={node.y + 45}
                                            textAnchor="middle"
                                            fill="#94a3b8"
                                            fontSize="10"
                                        >
                                            H={node.hydraulic_head}m
                                        </text>
                                    )}
                                </g>
                            ))}

                            <defs>
                                <radialGradient id="nodeGradient" cx="30%" cy="30%">
                                    <stop offset="0%" stopColor="#1e3a5f" />
                                    <stop offset="100%" stopColor="#0f172a" />
                                </radialGradient>
                            </defs>
                        </svg>
                    </div>
                </div>

                {/* Controls Panel */}
                <div className="controls-panel">
                    {/* Add Node */}
                    <div className="card">
                        <div className="card-header">
                            <h3 className="card-title">Add Node</h3>
                        </div>
                        <div className="form-group">
                            <label className="form-label">Node ID</label>
                            <input
                                type="text"
                                className="form-input"
                                value={newNode.id}
                                onChange={(e) => setNewNode({ ...newNode, id: e.target.value })}
                                placeholder="e.g., well_d"
                            />
                        </div>
                        <div className="form-group">
                            <label className="form-label">Name</label>
                            <input
                                type="text"
                                className="form-input"
                                value={newNode.name}
                                onChange={(e) => setNewNode({ ...newNode, name: e.target.value })}
                                placeholder="e.g., Well D"
                            />
                        </div>
                        <div className="form-group">
                            <label className="form-label">Hydraulic Head (m)</label>
                            <input
                                type="number"
                                className="form-input"
                                value={newNode.hydraulic_head}
                                onChange={(e) => setNewNode({ ...newNode, hydraulic_head: e.target.value })}
                                placeholder="e.g., 235"
                            />
                        </div>
                        <button className="btn btn-primary w-full" onClick={addNode}>
                            Add Node
                        </button>
                    </div>

                    {/* Node List */}
                    <div className="card mt-lg">
                        <div className="card-header">
                            <h3 className="card-title">Nodes</h3>
                        </div>
                        <div className="nodes-list">
                            {nodes.map((node) => (
                                <div key={node.id} className={`node-item ${selectedNode === node.id ? 'selected' : ''}`}>
                                    <div className="node-info">
                                        <span className="node-name">{node.name}</span>
                                        <span className="node-head">H: {node.hydraulic_head || '-'}m</span>
                                    </div>
                                    <button
                                        className="delete-btn"
                                        onClick={() => removeNode(node.id)}
                                    >
                                        ✕
                                    </button>
                                </div>
                            ))}
                        </div>
                    </div>

                    {/* Inference Results */}
                    {inferenceResult && (
                        <div className="card mt-lg inference-results">
                            <div className="card-header">
                                <h3 className="card-title">Inference Results</h3>
                                <span className="badge badge-success">Computed</span>
                            </div>
                            <div className="results-list">
                                {inferenceResult.inferred_edges.map((edge, idx) => (
                                    <div key={idx} className="result-item">
                                        <span className="edge-label">{edge.source} → {edge.target}</span>
                                        <span className={`flow-badge ${edge.flow_direction}`}>
                                            {edge.flow_direction}
                                        </span>
                                        <span className="prob-value">P: {edge.flow_probability?.toFixed(2)}</span>
                                    </div>
                                ))}
                            </div>
                        </div>
                    )}
                </div>
            </div>
        </div>
    )
}

export default Network
