import { createContext, useContext, useState, useEffect } from 'react'

const API_BASE = 'http://localhost:8000/api'

const ProjectContext = createContext(null)

export function ProjectProvider({ children }) {
    const [currentProject, setCurrentProject] = useState(null)
    const [projects, setProjects] = useState([])
    const [loading, setLoading] = useState(true)

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
            console.error('Failed to fetch projects:', error)
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
                return newProject
            }
            return null
        } catch (error) {
            console.error('Failed to create project:', error)
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
            console.error('Failed to select project:', error)
        }
        return null
    }

    const saveResultToProject = async (result) => {
        if (!currentProject) {
            console.warn('No project selected - results will not be saved')
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
                return true
            }
        } catch (error) {
            console.error('Failed to save result:', error)
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
            }
        } catch (error) {
            console.error('Failed to download report:', error)
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
