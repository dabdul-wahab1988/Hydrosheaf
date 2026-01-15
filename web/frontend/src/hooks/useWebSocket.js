import { useState, useEffect, useRef, useCallback } from 'react'

/**
 * Custom hook for WebSocket connection to analysis progress updates.
 *
 * @param {string} jobId - The analysis job ID to track
 * @param {Object} options - Configuration options
 * @param {boolean} options.enabled - Whether to enable the WebSocket connection
 * @param {function} options.onProgress - Callback for progress updates
 * @param {function} options.onComplete - Callback for completion
 * @param {function} options.onError - Callback for errors
 * @returns {Object} WebSocket state and controls
 */
export function useWebSocket(jobId, options = {}) {
    const {
        enabled = true,
        onProgress,
        onComplete,
        onError
    } = options

    const [isConnected, setIsConnected] = useState(false)
    const [progress, setProgress] = useState(0)
    const [step, setStep] = useState('')
    const [status, setStatus] = useState('pending')
    const [error, setError] = useState(null)

    const wsRef = useRef(null)
    const reconnectTimeoutRef = useRef(null)
    const pingIntervalRef = useRef(null)

    // Determine WebSocket URL based on current location
    const getWsUrl = useCallback(() => {
        const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:'
        const host = import.meta.env.VITE_API_HOST || window.location.host

        // In development, use the backend port
        const wsHost = host.includes(':5173') || host.includes(':5174')
            ? host.replace(/:\d+/, ':8000')
            : host

        return `${protocol}//${wsHost}/ws/analysis/${jobId}`
    }, [jobId])

    const connect = useCallback(() => {
        if (!jobId || !enabled || wsRef.current?.readyState === WebSocket.OPEN) {
            return
        }

        const url = getWsUrl()

        try {
            wsRef.current = new WebSocket(url)

            wsRef.current.onopen = () => {
                setIsConnected(true)
                setError(null)

                // Start ping interval to keep connection alive
                pingIntervalRef.current = setInterval(() => {
                    if (wsRef.current?.readyState === WebSocket.OPEN) {
                        wsRef.current.send('ping')
                    }
                }, 30000)
            }

            wsRef.current.onmessage = (event) => {
                try {
                    const data = JSON.parse(event.data)

                    if (data.type === 'pong') {
                        // Pong response, connection is alive
                        return
                    }

                    if (data.type === 'progress') {
                        setProgress(data.progress)
                        setStep(data.step)
                        setStatus(data.status)
                        onProgress?.(data)
                    } else if (data.type === 'completion') {
                        setProgress(data.progress)
                        setStatus(data.status)
                        if (data.error) {
                            setError(data.error)
                            onError?.(data.error)
                        }
                        onComplete?.(data)
                    }
                } catch (err) {
                    console.error('Failed to parse WebSocket message:', err)
                }
            }

            wsRef.current.onerror = (event) => {
                console.error('WebSocket error:', event)
                setError('WebSocket connection error')
            }

            wsRef.current.onclose = (event) => {
                setIsConnected(false)
                clearInterval(pingIntervalRef.current)

                // Reconnect after delay if not completed
                if (status !== 'completed' && status !== 'failed' && enabled) {
                    reconnectTimeoutRef.current = setTimeout(() => {
                        connect()
                    }, 3000)
                }
            }
        } catch (err) {
            console.error('Failed to create WebSocket:', err)
            setError(err.message)
        }
    }, [jobId, enabled, getWsUrl, onProgress, onComplete, onError, status])

    const disconnect = useCallback(() => {
        clearTimeout(reconnectTimeoutRef.current)
        clearInterval(pingIntervalRef.current)

        if (wsRef.current) {
            wsRef.current.close()
            wsRef.current = null
        }
        setIsConnected(false)
    }, [])

    // Connect when jobId changes
    useEffect(() => {
        if (jobId && enabled) {
            connect()
        }

        return () => {
            disconnect()
        }
    }, [jobId, enabled, connect, disconnect])

    // Cleanup on unmount
    useEffect(() => {
        return () => {
            disconnect()
        }
    }, [disconnect])

    return {
        isConnected,
        progress,
        step,
        status,
        error,
        connect,
        disconnect
    }
}

export default useWebSocket
