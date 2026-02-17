import { useState, useCallback } from 'react'
import { API_BASE } from '../config'
import { useToast } from '../context/ToastContext'
import NetworkGraph from '../components/NetworkGraph'
import './Network.css'

function Network() {
    const toast = useToast()
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
    const [newEdge, setNewEdge] = useState({ source: '', target: '' })
    const [selectedNode, setSelectedNode] = useState(null)
    const [inferenceResult, setInferenceResult] = useState(null)
    const [inferring, setInferring] = useState(false)

    const addNode = () => {
        if (!newNode.id || !newNode.name) {
            toast.warning('Please provide node ID and name')
            return
        }

        if (nodes.find(n => n.id === newNode.id)) {
            toast.warning('Node ID already exists')
            return
        }

        setNodes([...nodes, {
            ...newNode,
            x: 200 + Math.random() * 200,
            y: 150 + Math.random() * 150,
            hydraulic_head: parseFloat(newNode.hydraulic_head) || null
        }])
        setNewNode({ id: '', name: '', hydraulic_head: '' })
        toast.success(`Added node "${newNode.name}"`)
    }

    const removeNode = (id) => {
        const node = nodes.find(n => n.id === id)
        setNodes(nodes.filter(n => n.id !== id))
        setEdges(edges.filter(e => e.source !== id && e.target !== id))
        if (selectedNode === id) setSelectedNode(null)
        toast.info(`Removed node "${node?.name || id}"`)
    }

    const addEdge = () => {
        if (!newEdge.source || !newEdge.target) {
            toast.warning('Please select source and target nodes')
            return
        }

        if (newEdge.source === newEdge.target) {
            toast.warning('Source and target must be different')
            return
        }

        if (edges.find(e => e.source === newEdge.source && e.target === newEdge.target)) {
            toast.warning('Edge already exists')
            return
        }

        setEdges([...edges, { source: newEdge.source, target: newEdge.target }])
        setNewEdge({ source: '', target: '' })
        toast.success('Added edge')
    }

    const removeEdge = (source, target) => {
        setEdges(edges.filter(e => !(e.source === source && e.target === target)))
        toast.info('Removed edge')
    }

    const handleNodeMove = useCallback((nodeId, x, y) => {
        setNodes(prev => prev.map(n =>
            n.id === nodeId ? { ...n, x, y } : n
        ))
    }, [])

    const runInference = async () => {
        if (nodes.length < 2) {
            toast.warning('Add at least 2 nodes before running inference')
            return
        }

        if (edges.length < 1) {
            toast.warning('Add at least 1 edge before running inference')
            return
        }

        setInferring(true)
        try {
            // Create network first
            const createRes = await fetch(`${API_BASE}/network/create`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ name: 'flow_network', nodes, edges })
            })

            if (!createRes.ok) {
                throw new Error('Failed to create network')
            }

            // Run inference
            const res = await fetch(`${API_BASE}/network/flow_network/infer-flow`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    method: 'bayesian',
                    radius_km: 10.0,
                    max_neighbors: 3
                })
            })

            if (res.ok) {
                const result = await res.json()
                setInferenceResult(result)
                toast.success(result.hydrosheaf_used
                    ? 'Flow inference completed using Hydrosheaf'
                    : 'Flow inference completed using gradient method'
                )
            } else {
                throw new Error('Inference request failed')
            }
        } catch (error) {
            toast.error(`Inference error: ${error.message}`)
        } finally {
            setInferring(false)
        }
    }

    const clearNetwork = () => {
        setNodes([])
        setEdges([])
        setSelectedNode(null)
        setInferenceResult(null)
        toast.info('Network cleared')
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
                <div className="header-actions">
                    <button
                        className="btn btn-secondary"
                        onClick={clearNetwork}
                        disabled={nodes.length === 0}
                    >
                        <span>🗑️</span> Clear
                    </button>
                    <button
                        className="btn btn-primary"
                        onClick={runInference}
                        disabled={inferring || nodes.length < 2}
                    >
                        {inferring ? (
                            <>
                                <span className="spinner"></span> Inferring...
                            </>
                        ) : (
                            <>
                                <span>🔮</span> Infer Flow Directions
                            </>
                        )}
                    </button>
                </div>
            </header>

            <div className="network-layout">
                {/* Network Canvas */}
                <div className="card network-canvas-card">
                    <div className="card-header">
                        <h3 className="card-title">Network Visualization</h3>
                        <div className="header-badges">
                            <span className="badge badge-primary">{nodes.length} Nodes</span>
                            <span className="badge badge-secondary">{edges.length} Edges</span>
                        </div>
                    </div>

                    <NetworkGraph
                        nodes={nodes}
                        edges={edges}
                        inferenceResult={inferenceResult}
                        selectedNode={selectedNode}
                        onNodeSelect={setSelectedNode}
                        onNodeMove={handleNodeMove}
                        width={600}
                        height={400}
                    />

                    {/* Flow direction legend */}
                    {inferenceResult && (
                        <div className="flow-legend">
                            <div className="legend-item">
                                <span className="legend-color forward"></span>
                                <span>Forward flow</span>
                            </div>
                            <div className="legend-item">
                                <span className="legend-color reverse"></span>
                                <span>Reverse flow</span>
                            </div>
                            <div className="legend-item">
                                <span className="legend-color uncertain"></span>
                                <span>Uncertain</span>
                            </div>
                        </div>
                    )}
                </div>

                {/* Controls Panel */}
                <div className="controls-panel">
                    {/* Add Node */}
                    <div className="card">
                        <div className="card-header">
                            <h3 className="card-title">Add Node</h3>
                        </div>
                        <div className="form-group">
                            <label className="form-label" htmlFor="node-id">Node ID</label>
                            <input
                                id="node-id"
                                type="text"
                                className="form-input"
                                value={newNode.id}
                                onChange={(e) => setNewNode({ ...newNode, id: e.target.value })}
                                placeholder="e.g., well_d"
                            />
                        </div>
                        <div className="form-group">
                            <label className="form-label" htmlFor="node-name">Name</label>
                            <input
                                id="node-name"
                                type="text"
                                className="form-input"
                                value={newNode.name}
                                onChange={(e) => setNewNode({ ...newNode, name: e.target.value })}
                                placeholder="e.g., Well D"
                            />
                        </div>
                        <div className="form-group">
                            <label className="form-label" htmlFor="node-head">Hydraulic Head (m)</label>
                            <input
                                id="node-head"
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

                    {/* Add Edge */}
                    <div className="card mt-lg">
                        <div className="card-header">
                            <h3 className="card-title">Add Edge</h3>
                        </div>
                        <div className="form-group">
                            <label className="form-label" htmlFor="edge-source">Source Node</label>
                            <select
                                id="edge-source"
                                className="form-select"
                                value={newEdge.source}
                                onChange={(e) => setNewEdge({ ...newEdge, source: e.target.value })}
                            >
                                <option value="">Select source...</option>
                                {nodes.map(n => (
                                    <option key={n.id} value={n.id}>{n.name}</option>
                                ))}
                            </select>
                        </div>
                        <div className="form-group">
                            <label className="form-label" htmlFor="edge-target">Target Node</label>
                            <select
                                id="edge-target"
                                className="form-select"
                                value={newEdge.target}
                                onChange={(e) => setNewEdge({ ...newEdge, target: e.target.value })}
                            >
                                <option value="">Select target...</option>
                                {nodes.filter(n => n.id !== newEdge.source).map(n => (
                                    <option key={n.id} value={n.id}>{n.name}</option>
                                ))}
                            </select>
                        </div>
                        <button
                            className="btn btn-secondary w-full"
                            onClick={addEdge}
                            disabled={!newEdge.source || !newEdge.target}
                        >
                            Add Edge
                        </button>
                    </div>

                    {/* Node List */}
                    <div className="card mt-lg">
                        <div className="card-header">
                            <h3 className="card-title">Nodes</h3>
                        </div>
                        <div className="nodes-list">
                            {nodes.length === 0 ? (
                                <div className="empty-list">No nodes added</div>
                            ) : (
                                nodes.map((node) => (
                                    <div
                                        key={node.id}
                                        className={`node-item ${selectedNode === node.id ? 'selected' : ''}`}
                                        onClick={() => setSelectedNode(selectedNode === node.id ? null : node.id)}
                                    >
                                        <div className="node-info">
                                            <span className="node-name">{node.name}</span>
                                            <span className="node-head">H: {node.hydraulic_head || '-'}m</span>
                                        </div>
                                        <button
                                            className="delete-btn"
                                            onClick={(e) => {
                                                e.stopPropagation()
                                                removeNode(node.id)
                                            }}
                                            aria-label={`Remove ${node.name}`}
                                        >
                                            ✕
                                        </button>
                                    </div>
                                ))
                            )}
                        </div>
                    </div>

                    {/* Edge List */}
                    <div className="card mt-lg">
                        <div className="card-header">
                            <h3 className="card-title">Edges</h3>
                        </div>
                        <div className="edges-list">
                            {edges.length === 0 ? (
                                <div className="empty-list">No edges added</div>
                            ) : (
                                edges.map((edge, idx) => {
                                    const sourceNode = nodes.find(n => n.id === edge.source)
                                    const targetNode = nodes.find(n => n.id === edge.target)
                                    return (
                                        <div key={idx} className="edge-item">
                                            <span className="edge-label">
                                                {sourceNode?.name || edge.source} → {targetNode?.name || edge.target}
                                            </span>
                                            <button
                                                className="delete-btn"
                                                onClick={() => removeEdge(edge.source, edge.target)}
                                                aria-label={`Remove edge ${edge.source} to ${edge.target}`}
                                            >
                                                ✕
                                            </button>
                                        </div>
                                    )
                                })
                            )}
                        </div>
                    </div>

                    {/* Inference Results */}
                    {inferenceResult && (
                        <div className="card mt-lg inference-results">
                            <div className="card-header">
                                <h3 className="card-title">Inference Results</h3>
                                <span className="badge badge-success">
                                    {inferenceResult.hydrosheaf_used ? 'Hydrosheaf' : 'Gradient'}
                                </span>
                            </div>
                            <div className="results-info">
                                <span className="method-label">
                                    Method: {inferenceResult.method}
                                </span>
                            </div>
                            <div className="results-list">
                                {inferenceResult.inferred_edges?.map((edge, idx) => (
                                    <div key={idx} className="result-item">
                                        <span className="edge-label">{edge.source} → {edge.target}</span>
                                        <span className={`flow-badge ${edge.flow_direction}`}>
                                            {edge.flow_direction}
                                        </span>
                                        <span className="prob-value">
                                            P: {edge.flow_probability?.toFixed(2)}
                                        </span>
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
