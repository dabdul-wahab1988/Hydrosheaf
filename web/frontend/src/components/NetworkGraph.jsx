import { useEffect, useRef, useCallback } from 'react'
import * as d3 from 'd3'
import './NetworkGraph.css'

/**
 * Interactive D3.js force-directed network graph
 * Features: draggable nodes, zoom/pan, color-coded edges by flow direction
 */
function NetworkGraph({
    nodes,
    edges,
    inferenceResult,
    selectedNode,
    onNodeSelect,
    onNodeMove,
    width = 600,
    height = 400
}) {
    const svgRef = useRef(null)
    const simulationRef = useRef(null)
    const containerRef = useRef(null)

    // Get edge color based on flow direction
    const getEdgeColor = useCallback((edge) => {
        if (!inferenceResult?.inferred_edges) return '#0ea5e9'

        const inference = inferenceResult.inferred_edges.find(
            e => e.source === edge.source && e.target === edge.target
        )

        if (!inference) return '#0ea5e9'

        switch (inference.flow_direction) {
            case 'forward': return '#22c55e'
            case 'reverse': return '#ef4444'
            default: return '#94a3b8'
        }
    }, [inferenceResult])

    // Get flow probability for edge label
    const getFlowProbability = useCallback((edge) => {
        if (!inferenceResult?.inferred_edges) return null

        const inference = inferenceResult.inferred_edges.find(
            e => e.source === edge.source && e.target === edge.target
        )

        return inference?.flow_probability
    }, [inferenceResult])

    useEffect(() => {
        if (!svgRef.current || !nodes.length) return

        const svg = d3.select(svgRef.current)
        const container = svg.select('.graph-container')

        // Clear previous content
        container.selectAll('*').remove()

        // Create zoom behavior
        const zoom = d3.zoom()
            .scaleExtent([0.3, 3])
            .on('zoom', (event) => {
                container.attr('transform', event.transform)
            })

        svg.call(zoom)

        // Prepare data for D3
        const nodeData = nodes.map(n => ({
            ...n,
            x: n.x || width / 2,
            y: n.y || height / 2
        }))

        const linkData = edges.map(e => ({
            ...e,
            sourceId: e.source,
            targetId: e.target
        }))

        // Create force simulation
        const simulation = d3.forceSimulation(nodeData)
            .force('link', d3.forceLink(linkData)
                .id(d => d.id)
                .distance(120)
                .strength(0.5))
            .force('charge', d3.forceManyBody()
                .strength(-300)
                .distanceMax(250))
            .force('center', d3.forceCenter(width / 2, height / 2))
            .force('collision', d3.forceCollide().radius(40))

        simulationRef.current = simulation

        // Create gradient for nodes
        const defs = svg.select('defs')
        if (defs.select('#nodeGradient').empty()) {
            defs.append('radialGradient')
                .attr('id', 'nodeGradient')
                .attr('cx', '30%')
                .attr('cy', '30%')
                .selectAll('stop')
                .data([
                    { offset: '0%', color: '#1e3a5f' },
                    { offset: '100%', color: '#0f172a' }
                ])
                .join('stop')
                .attr('offset', d => d.offset)
                .attr('stop-color', d => d.color)
        }

        // Create arrow markers for directed edges
        defs.selectAll('marker.arrow')
            .data(['forward', 'reverse', 'uncertain', 'default'])
            .join('marker')
            .attr('class', 'arrow')
            .attr('id', d => `arrow-${d}`)
            .attr('viewBox', '0 -5 10 10')
            .attr('refX', 35)
            .attr('refY', 0)
            .attr('markerWidth', 6)
            .attr('markerHeight', 6)
            .attr('orient', 'auto')
            .append('path')
            .attr('fill', d => {
                switch(d) {
                    case 'forward': return '#22c55e'
                    case 'reverse': return '#ef4444'
                    case 'uncertain': return '#94a3b8'
                    default: return '#0ea5e9'
                }
            })
            .attr('d', 'M0,-5L10,0L0,5')

        // Create links (edges)
        const link = container.append('g')
            .attr('class', 'links')
            .selectAll('g')
            .data(linkData)
            .join('g')
            .attr('class', 'link-group')

        // Edge lines
        link.append('line')
            .attr('class', 'link')
            .attr('stroke', d => getEdgeColor(d))
            .attr('stroke-width', 3)
            .attr('stroke-opacity', 0.7)
            .attr('marker-end', d => {
                const inference = inferenceResult?.inferred_edges?.find(
                    e => e.source === d.sourceId && e.target === d.targetId
                )
                const direction = inference?.flow_direction || 'default'
                return `url(#arrow-${direction})`
            })

        // Edge probability labels
        link.append('text')
            .attr('class', 'link-label')
            .attr('fill', '#94a3b8')
            .attr('font-size', '11px')
            .attr('text-anchor', 'middle')
            .attr('dy', -8)
            .text(d => {
                const prob = getFlowProbability(d)
                return prob != null ? `P=${prob.toFixed(2)}` : ''
            })

        // Create nodes
        const node = container.append('g')
            .attr('class', 'nodes')
            .selectAll('g')
            .data(nodeData)
            .join('g')
            .attr('class', d => `node ${selectedNode === d.id ? 'selected' : ''}`)
            .attr('cursor', 'grab')

        // Drag behavior
        const drag = d3.drag()
            .on('start', (event, d) => {
                if (!event.active) simulation.alphaTarget(0.3).restart()
                d.fx = d.x
                d.fy = d.y
                d3.select(event.sourceEvent.target.parentNode).attr('cursor', 'grabbing')
            })
            .on('drag', (event, d) => {
                d.fx = event.x
                d.fy = event.y
            })
            .on('end', (event, d) => {
                if (!event.active) simulation.alphaTarget(0)
                // Keep node fixed at new position
                d3.select(event.sourceEvent.target.parentNode).attr('cursor', 'grab')
                if (onNodeMove) {
                    onNodeMove(d.id, d.x, d.y)
                }
            })

        node.call(drag)

        // Node circles
        node.append('circle')
            .attr('r', 25)
            .attr('fill', 'url(#nodeGradient)')
            .attr('stroke', d => selectedNode === d.id ? '#a855f7' : '#0ea5e9')
            .attr('stroke-width', 2)
            .on('click', (event, d) => {
                event.stopPropagation()
                if (onNodeSelect) {
                    onNodeSelect(selectedNode === d.id ? null : d.id)
                }
            })

        // Node labels
        node.append('text')
            .attr('class', 'node-label')
            .attr('text-anchor', 'middle')
            .attr('dy', 5)
            .attr('fill', '#f8fafc')
            .attr('font-size', '12px')
            .attr('font-weight', '600')
            .attr('pointer-events', 'none')
            .text(d => d.name.slice(0, 5))

        // Hydraulic head labels
        node.append('text')
            .attr('class', 'head-label')
            .attr('text-anchor', 'middle')
            .attr('dy', 45)
            .attr('fill', '#94a3b8')
            .attr('font-size', '10px')
            .attr('pointer-events', 'none')
            .text(d => d.hydraulic_head ? `H=${d.hydraulic_head}m` : '')

        // Glow effect on hover
        node.on('mouseenter', function() {
            d3.select(this).select('circle')
                .transition()
                .duration(200)
                .attr('filter', 'url(#glow)')
        }).on('mouseleave', function() {
            d3.select(this).select('circle')
                .transition()
                .duration(200)
                .attr('filter', null)
        })

        // Create glow filter
        if (defs.select('#glow').empty()) {
            const filter = defs.append('filter')
                .attr('id', 'glow')
                .attr('x', '-50%')
                .attr('y', '-50%')
                .attr('width', '200%')
                .attr('height', '200%')

            filter.append('feGaussianBlur')
                .attr('stdDeviation', '3')
                .attr('result', 'coloredBlur')

            const feMerge = filter.append('feMerge')
            feMerge.append('feMergeNode').attr('in', 'coloredBlur')
            feMerge.append('feMergeNode').attr('in', 'SourceGraphic')
        }

        // Update positions on simulation tick
        simulation.on('tick', () => {
            link.selectAll('line')
                .attr('x1', d => d.source.x)
                .attr('y1', d => d.source.y)
                .attr('x2', d => d.target.x)
                .attr('y2', d => d.target.y)

            link.selectAll('text')
                .attr('x', d => (d.source.x + d.target.x) / 2)
                .attr('y', d => (d.source.y + d.target.y) / 2)

            node.attr('transform', d => `translate(${d.x},${d.y})`)
        })

        // Double-click to reset zoom
        svg.on('dblclick.zoom', () => {
            svg.transition()
                .duration(500)
                .call(zoom.transform, d3.zoomIdentity)
        })

        // Cleanup
        return () => {
            simulation.stop()
        }
    }, [nodes, edges, inferenceResult, selectedNode, width, height, getEdgeColor, getFlowProbability, onNodeSelect, onNodeMove])

    // Update selected node appearance
    useEffect(() => {
        if (!svgRef.current) return

        const svg = d3.select(svgRef.current)

        svg.selectAll('.node')
            .classed('selected', d => d.id === selectedNode)
            .select('circle')
            .attr('stroke', d => d.id === selectedNode ? '#a855f7' : '#0ea5e9')
    }, [selectedNode])

    return (
        <div className="network-graph-container" ref={containerRef}>
            <svg
                ref={svgRef}
                width="100%"
                height={height}
                viewBox={`0 0 ${width} ${height}`}
                className="network-graph-svg"
            >
                <defs>
                    {/* Grid pattern */}
                    <pattern id="grid" width="50" height="50" patternUnits="userSpaceOnUse">
                        <path
                            d="M 50 0 L 0 0 0 50"
                            fill="none"
                            stroke="rgba(148,163,184,0.1)"
                            strokeWidth="1"
                        />
                    </pattern>
                </defs>
                {/* Background with grid */}
                <rect width="100%" height="100%" fill="url(#grid)" />
                {/* Main graph container (for zoom/pan) */}
                <g className="graph-container"></g>
            </svg>
            <div className="graph-controls">
                <span className="control-hint">Drag nodes to reposition</span>
                <span className="control-hint">Scroll to zoom</span>
                <span className="control-hint">Double-click to reset</span>
            </div>
        </div>
    )
}

export default NetworkGraph
