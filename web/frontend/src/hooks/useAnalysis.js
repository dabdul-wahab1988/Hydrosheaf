import { useState, useEffect, useCallback } from 'react'
import { API_BASE } from '../config'
import { useToast } from '../context/ToastContext'

/**
 * Custom hook for analysis-related data fetching and state management
 */
export function useAnalysis(currentProject, saveResultToProject) {
    const [datasets, setDatasets] = useState([])
    const [capabilities, setCapabilities] = useState(null)
    const [capabilitiesLoading, setCapabilitiesLoading] = useState(false)
    const [currentJob, setCurrentJob] = useState(null)
    const [running, setRunning] = useState(false)
    const [progress, setProgress] = useState(0)
    const toast = useToast()

    // Fetch available datasets
    const fetchDatasets = useCallback(async () => {
        try {
            const res = await fetch(`${API_BASE}/samples/datasets`)
            if (res.ok) {
                setDatasets(await res.json())
            }
        } catch (error) {
            toast.error('Failed to fetch datasets')
        }
    }, [toast])

    // Fetch dataset capabilities
    const fetchCapabilities = useCallback(async (datasetId) => {
        if (!datasetId) {
            setCapabilities(null)
            return null
        }
        setCapabilitiesLoading(true)
        try {
            const res = await fetch(`${API_BASE}/samples/datasets/${datasetId}/capabilities`)
            if (res.ok) {
                const caps = await res.json()
                setCapabilities(caps)
                return caps
            }
        } catch (error) {
            toast.error('Failed to analyze dataset capabilities')
        } finally {
            setCapabilitiesLoading(false)
        }
        return null
    }, [toast])

    // Check job status
    const checkJobStatus = useCallback(async (formData) => {
        if (!currentJob) return

        try {
            const res = await fetch(`${API_BASE}/analysis/status/${currentJob.job_id}`)
            if (res.ok) {
                const job = await res.json()
                setCurrentJob(job)

                if (job.status === 'running') {
                    setProgress(prev => Math.min(prev + 10, 90))
                } else if (job.status === 'completed') {
                    setProgress(100)
                    setRunning(false)
                    toast.success('Analysis completed successfully!')

                    // Save results to project if active
                    if (currentProject && job.results) {
                        const resultToSave = {
                            name: formData.name,
                            analysis_type: formData.analysis_type,
                            job_id: currentJob.job_id,
                            ...job.results
                        }
                        await saveResultToProject(resultToSave)
                    }
                } else if (job.status === 'failed') {
                    setRunning(false)
                    toast.error(`Analysis failed: ${job.error || 'Unknown error'}`)
                }
            }
        } catch (error) {
            toast.error('Failed to check analysis status')
        }
    }, [currentJob, currentProject, saveResultToProject, toast])

    // Start analysis
    const startAnalysis = useCallback(async (formData) => {
        if (!formData.name || !formData.dataset_id) {
            toast.warning('Please provide an analysis name and select a dataset')
            return false
        }

        setRunning(true)
        setProgress(0)

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
                setCurrentJob(job)
                setProgress(10)
                toast.info('Analysis started...')
                return true
            } else {
                throw new Error('Failed to start analysis')
            }
        } catch (error) {
            toast.error('Failed to start analysis')
            setRunning(false)
            return false
        }
    }, [currentProject, toast])

    // Initial fetch
    useEffect(() => {
        fetchDatasets()
    }, [fetchDatasets])

    // Poll job status
    useEffect(() => {
        let interval
        if (currentJob && currentJob.status !== 'completed' && currentJob.status !== 'failed') {
            interval = setInterval(() => checkJobStatus({}), 1000)
        }
        return () => clearInterval(interval)
    }, [currentJob, checkJobStatus])

    return {
        datasets,
        capabilities,
        capabilitiesLoading,
        currentJob,
        running,
        progress,
        fetchCapabilities,
        startAnalysis,
        checkJobStatus,
    }
}
