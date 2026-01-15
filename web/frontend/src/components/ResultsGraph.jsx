import { useEffect, useRef, useState } from 'react'
import * as d3 from 'd3'
import './ResultsGraph.css'

/**
 * D3.js force-directed graph visualization for analysis results.
 * Shows nodes (sampling locations) and edges (flow connections) with
 * colors indicating transport model type.
 */
function ResultsGraph({ edgeResults = [], width = 600, height = 400 }) {
    const svgRef = useRef(null)
    const [selectedEdge, setSelectedEdge] = useState(null)

    useEffect(() => {
        if (!svgRef.current || edgeResults.length === 0) return

        // Clear previous content
        d3.select(svgRef.current).selectAll('*').remove()

        // Extract unique nodes from edges
        const nodesSet = new Set()
        edgeResults.forEach(edge => {
            const parts = (edge.edge_id || '').split('->')
            if (parts.length >= 2) {
                nodesSet.add(parts[0].trim())
                nodesSet.add(parts[1].trim())
            }
        })

        const nodes = Array.from(nodesSet).map(id => ({ id, name: id }))

        // Create links from edge results
        const links = edgeResults.map(edge => {
            const parts = (edge.edge_id || '').split('->')
            return {
                source: parts[0]?.trim() || '',
                target: parts[1]?.trim() || '',
                gamma: edge.gamma || 0,
                f: edge.f || 0,
                model: edge.transport_model || 'unknown',
                edge_id: edge.edge_id
            }
        }).filter(link => link.source && link.target)

        if (nodes.length === 0 || links.length === 0) return

        const svg = d3.select(svgRef.current)
            .attr('width', width)
            .attr('height', height)

        // Add zoom behavior
        const g = svg.append('g')

        const zoom = d3.zoom()
            .scaleExtent([0.5, 3])
            .on('zoom', (event) => {
                g.attr('transform', event.transform)
            })

        svg.call(zoom)

        // Define arrow markers for different edge types
        const defs = svg.append('defs')

        // Evaporation arrow (red)
        defs.append('marker')
            .attr('id', 'arrow-evap')
            .attr('viewBox', '0 -5 10 10')
            .attr('refX', 20)
            .attr('refY', 0)
            .attr('markerWidth', 6)
            .attr('markerHeight', 6)
            .attr('orient', 'auto')
            .append('path')
            .attr('fill', '#ef4444')
            .attr('d', 'M0,-5L10,0L0,5')

        // Mixing arrow (green)
        defs.append('marker')
            .attr('id', 'arrow-mix')
            .attr('viewBox', '0 -5 10 10')
            .attr('refX', 20)
            .attr('refY', 0)
            .attr('markerWidth', 6)
            .attr('markerHeight', 6)
            .attr('orient', 'auto')
            .append('path')
            .attr('fill', '#22c55e')
            .attr('d', 'M0,-5L10,0L0,5')

        // Unknown arrow (gray)
        defs.append('marker')
            .attr('id', 'arrow-unknown')
            .attr('viewBox', '0 -5 10 10')
            .attr('refX', 20)
            .attr('refY', 0)
            .attr('markerWidth', 6)
            .attr('markerHeight', 6)
            .attr('orient', 'auto')
            .append('path')
            .attr('fill', '#94a3b8')
            .attr('d', 'M0,-5L10,0L0,5')

        // Create force simulation
        const simulation = d3.forceSimulation(nodes)
            .force('link', d3.forceLink(links).id(d => d.id).distance(100))
            .force('charge', d3.forceManyBody().strength(-300))
            .force('center', d3.forceCenter(width / 2, height / 2))
            .force('collision', d3.forceCollide().radius(30))

        // Get edge color based on transport model
        const getEdgeColor = (model) => {
            if (model === 'evaporation' || model === 'evap') return '#ef4444'
            if (model === 'mixing' || model === 'mix') return '#22c55e'
            return '#94a3b8'
        }

        const getArrowMarker = (model) => {
            if (model === 'evaporation' || model === 'evap') return 'url(#arrow-evap)'
            if (model === 'mixing' || model === 'mix') return 'url(#arrow-mix)'
            return 'url(#arrow-unknown)'
        }

        // Draw links (edges)
        const link = g.append('g')
            .attr('class', 'links')
            .selectAll('line')
            .data(links)
            .join('line')
            .attr('stroke', d => getEdgeColor(d.model))
            .attr('stroke-width', d => Math.max(2, Math.min(6, Math.abs(d.gamma) * 10 + 2)))
            .attr('stroke-opacity', 0.8)
            .attr('marker-end', d => getArrowMarker(d.model))
            .style('cursor', 'pointer')
            .on('click', (event, d) => {
                setSelectedEdge(d)
            })
            .on('mouseover', function(event, d) {
                d3.select(this)
                    .attr('stroke-width', d => Math.max(4, Math.min(8, Math.abs(d.gamma) * 10 + 4)))
                    .attr('stroke-opacity', 1)
            })
            .on('mouseout', function(event, d) {
                d3.select(this)
                    .attr('stroke-width', d => Math.max(2, Math.min(6, Math.abs(d.gamma) * 10 + 2)))
                    .attr('stroke-opacity', 0.8)
            })

        // Draw nodes
        const node = g.append('g')
            .attr('class', 'nodes')
            .selectAll('g')
            .data(nodes)
            .join('g')
            .call(d3.drag()
                .on('start', dragstarted)
                .on('drag', dragged)
                .on('end', dragended))

        // Node circles
        node.append('circle')
            .attr('r', 15)
            .attr('fill', 'url(#nodeGradient)')
            .attr('stroke', '#0ea5e9')
            .attr('stroke-width', 2)

        // Add gradient for nodes
        const nodeGradient = defs.append('radialGradient')
            .attr('id', 'nodeGradient')

        nodeGradient.append('stop')
            .attr('offset', '0%')
            .attr('stop-color', '#38bdf8')

        nodeGradient.append('stop')
            .attr('offset', '100%')
            .attr('stop-color', '#0284c7')

        // Node labels
        node.append('text')
            .attr('dy', 4)
            .attr('text-anchor', 'middle')
            .attr('fill', 'white')
            .attr('font-size', '10px')
            .attr('font-weight', 'bold')
            .text(d => d.id.substring(0, 3))

        // Node tooltips
        node.append('title')
            .text(d => d.id)

        // Update positions on tick
        simulation.on('tick', () => {
            link
                .attr('x1', d => d.source.x)
                .attr('y1', d => d.source.y)
                .attr('x2', d => d.target.x)
                .attr('y2', d => d.target.y)

            node
                .attr('transform', d => `translate(${d.x},${d.y})`)
        })

        // Drag functions
        function dragstarted(event, d) {
            if (!event.active) simulation.alphaTarget(0.3).restart()
            d.fx = d.x
            d.fy = d.y
        }

        function dragged(event, d) {
            d.fx = event.x
            d.fy = event.y
        }

        function dragended(event, d) {
            if (!event.active) simulation.alphaTarget(0)
            d.fx = null
            d.fy = null
        }

        // Cleanup
        return () => {
            simulation.stop()
        }
    }, [edgeResults, width, height])

    return (
        <div className="results-graph">
            <svg ref={svgRef}></svg>

            {/* Legend */}
            <div className="graph-legend">
                <div className="legend-item">
                    <span className="legend-line evap"></span>
                    <span>Evaporation</span>
                </div>
                <div className="legend-item">
                    <span className="legend-line mix"></span>
                    <span>Mixing</span>
                </div>
                <div className="legend-item">
                    <span className="legend-line unknown"></span>
                    <span>Unknown</span>
                </div>
            </div>

            {/* Selected Edge Info */}
            {selectedEdge && (
                <div className="edge-info-panel">
                    <div className="edge-info-header">
                        <h4>{selectedEdge.edge_id}</h4>
                        <button onClick={() => setSelectedEdge(null)} className="close-btn">x</button>
                    </div>
                    <div className="edge-info-content">
                        <div className="info-row">
                            <span className="info-label">Transport Model:</span>
                            <span className={`info-value model-${selectedEdge.model}`}>
                                {selectedEdge.model}
                            </span>
                        </div>
                        <div className="info-row">
                            <span className="info-label">Gamma (evap):</span>
                            <span className="info-value">{selectedEdge.gamma?.toFixed(4) || 'N/A'}</span>
                        </div>
                        <div className="info-row">
                            <span className="info-label">f (mixing):</span>
                            <span className="info-value">{selectedEdge.f?.toFixed(4) || 'N/A'}</span>
                        </div>
                    </div>
                </div>
            )}
        </div>
    )
}

export default ResultsGraph
