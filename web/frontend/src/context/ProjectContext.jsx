import { createContext, useContext, useState, useEffect } from 'react'
import { API_BASE } from '../config'
import { useToast } from './ToastContext'

const ProjectContext = createContext(null)

export function ProjectProvider({ children }) {
    const [currentProject, setCurrentProject] = useState(null)
    const [projects, setProjects] = useState([])
    const [loading, setLoading] = useState(true)
    const toast = useToast()

    useEffect(() => {
        fetchProjects()
    }, [])

    const fetchProjects = async () => {
        try {
            const res = await fetch(`${API_BASE}/projects/list`)
            if (res.ok) {
                setProjects(await res.json())
            }
        } catch (error) {
            toast.error('Failed to fetch projects. Check your connection.')
        } finally {
            setLoading(false)
        }
    }

    const createProject = async (name, description = '') => {
        try {
            const res = await fetch(`${API_BASE}/projects/create`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ name, description })
            })

            if (res.ok) {
                const data = await res.json()
                await fetchProjects()

                // Set the newly created project as current
                const newProject = {
                    id: data.project_id,
                    name: name,
                    description: description,
                }
                setCurrentProject(newProject)
                toast.success(`Project "${name}" created successfully`)
                return newProject
            }
            return null
        } catch (error) {
            toast.error('Failed to create project. Please try again.')
            return null
        }
    }

    const selectProject = async (projectId) => {
        if (!projectId) {
            setCurrentProject(null)
            return
        }

        try {
            const res = await fetch(`${API_BASE}/projects/${projectId}`)
            if (res.ok) {
                const project = await res.json()
                setCurrentProject(project)
                return project
            }
        } catch (error) {
            toast.error('Failed to select project')
        }
        return null
    }

    const saveResultToProject = async (result) => {
        if (!currentProject) {
            toast.warning('No project selected - results will not be saved')
            return false
        }

        try {
            const res = await fetch(`${API_BASE}/projects/${currentProject.id}/add-result`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(result)
            })

            if (res.ok) {
                // Refresh project data
                await selectProject(currentProject.id)
                toast.success('Result saved to project')
                return true
            }
        } catch (error) {
            toast.error('Failed to save result to project')
        }
        return false
    }

    const downloadProjectReport = async () => {
        if (!currentProject) return

        try {
            const res = await fetch(`${API_BASE}/projects/${currentProject.id}/download`)
            if (res.ok) {
                const blob = await res.blob()
                const url = window.URL.createObjectURL(blob)
                const a = document.createElement('a')
                a.href = url
                a.download = `${currentProject.name.replace(/\s+/g, '_')}_report.txt`
                document.body.appendChild(a)
                a.click()
                window.URL.revokeObjectURL(url)
                a.remove()
                toast.success('Report downloaded')
            }
        } catch (error) {
            toast.error('Failed to download report')
        }
    }

    const downloadCompleteProject = async () => {
        if (!currentProject) return

        try {
            const res = await fetch(`${API_BASE}/projects/${currentProject.id}/download-complete`)
            if (res.ok) {
                const blob = await res.blob()
                const url = window.URL.createObjectURL(blob)
                const a = document.createElement('a')
                a.href = url
                a.download = `${currentProject.name.replace(/\s+/g, '_')}.zip`
                document.body.appendChild(a)
                a.click()
                window.URL.revokeObjectURL(url)
                a.remove()
                toast.success('Project exported successfully')
            }
        } catch (error) {
            toast.error('Failed to download project export')
        }
    }

    const clearProject = () => {
        setCurrentProject(null)
    }

    const value = {
        currentProject,
        projects,
        loading,
        createProject,
        selectProject,
        clearProject,
        saveResultToProject,
        downloadProjectReport,
        downloadCompleteProject,
        fetchProjects,
    }

    return (
        <ProjectContext.Provider value={value}>
            {children}
        </ProjectContext.Provider>
    )
}

export function useProject() {
    const context = useContext(ProjectContext)
    if (!context) {
        throw new Error('useProject must be used within a ProjectProvider')
    }
    return context
}

export default ProjectContext
