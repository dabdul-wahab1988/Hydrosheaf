import { useState, useEffect } from 'react'
import { useNavigate } from 'react-router-dom'
import { useProject } from '../context/ProjectContext'
import './Projects.css'

function Projects() {
    const navigate = useNavigate()
    const {
        projects,
        currentProject,
        createProject,
        selectProject,
        clearProject,
        downloadProjectReport,
        downloadCompleteProject,
        fetchProjects
    } = useProject()

    const [showCreateModal, setShowCreateModal] = useState(false)
    const [newProjectName, setNewProjectName] = useState('')
    const [newProjectDesc, setNewProjectDesc] = useState('')
    const [creating, setCreating] = useState(false)

    const handleCreateProject = async () => {
        if (!newProjectName.trim()) return

        setCreating(true)
        const project = await createProject(newProjectName, newProjectDesc)
        setCreating(false)

        if (project) {
            setShowCreateModal(false)
            setNewProjectName('')
            setNewProjectDesc('')
        }
    }

    const handleSelectProject = async (projectId) => {
        await selectProject(projectId)
    }

    const handleStartAnalysis = () => {
        if (currentProject) {
            navigate('/analysis')
        }
    }

    return (
        <div className="projects-page">
            <header className="page-header">
                <div>
                    <h1 className="page-title">Projects</h1>
                    <p className="page-subtitle">
                        Create a project to save your analysis results
                    </p>
                </div>
                <button
                    className="btn btn-primary"
                    onClick={() => setShowCreateModal(true)}
                >
                    <span>➕</span> New Project
                </button>
            </header>

            {/* Current Project Banner */}
            {currentProject ? (
                <div className="current-project-card">
                    <div className="project-active-indicator">
                        <span className="pulse-dot"></span>
                        <span>Active Project</span>
                    </div>
                    <div className="project-details">
                        <h2>{currentProject.name}</h2>
                        {currentProject.description && (
                            <p className="project-description">{currentProject.description}</p>
                        )}
                        <div className="project-stats">
                            <span className="stat">
                                <strong>{currentProject.analysis_jobs?.length || 0}</strong> Analyses
                            </span>
                            <span className="stat">
                                <strong>{currentProject.results?.length || 0}</strong> Results Saved
                            </span>
                        </div>
                    </div>
                    <div className="project-actions">
                        <button className="btn btn-primary" onClick={handleStartAnalysis}>
                            Run Analysis
                        </button>
                        <button className="btn btn-secondary" onClick={downloadCompleteProject}>
                            📦 Export Project ZIP
                        </button>
                        <button className="btn btn-secondary" onClick={downloadProjectReport}>
                            📥 Download Report
                        </button>
                        <button className="btn btn-secondary" onClick={clearProject}>
                            Close Project
                        </button>
                    </div>
                </div>
            ) : (
                <div className="no-project-warning">
                    <div className="warning-icon">⚠️</div>
                    <div className="warning-content">
                        <h3>No Project Selected</h3>
                        <p>Create or select a project to save your analysis results. Without a project, your results will not be saved.</p>
                    </div>
                </div>
            )}

            {/* Projects Grid */}
            <div className="projects-section">
                <h2>Your Projects</h2>
                {projects.length === 0 ? (
                    <div className="empty-projects">
                        <span className="empty-icon">📁</span>
                        <h3>No Projects Yet</h3>
                        <p>Create your first project to start saving analysis results</p>
                        <button
                            className="btn btn-primary"
                            onClick={() => setShowCreateModal(true)}
                        >
                            Create Your First Project
                        </button>
                    </div>
                ) : (
                    <div className="projects-grid">
                        {projects.map((project) => (
                            <div
                                key={project.id}
                                className={`project-card ${currentProject?.id === project.id ? 'active' : ''}`}
                                onClick={() => handleSelectProject(project.id)}
                            >
                                <div className="project-card-header">
                                    <h3>{project.name}</h3>
                                    {currentProject?.id === project.id && (
                                        <span className="active-badge">Active</span>
                                    )}
                                </div>
                                {project.description && (
                                    <p className="project-card-desc">{project.description}</p>
                                )}
                                <div className="project-card-meta">
                                    <span>{project.analysis_count || 0} analyses</span>
                                    <span>Created {new Date(project.created_at).toLocaleDateString()}</span>
                                </div>
                                {project.has_results && (
                                    <span className="has-results-badge">Has Results</span>
                                )}
                            </div>
                        ))}
                    </div>
                )}
            </div>

            {/* Create Project Modal */}
            {showCreateModal && (
                <div className="modal-overlay" onClick={() => setShowCreateModal(false)}>
                    <div className="modal" onClick={e => e.stopPropagation()}>
                        <div className="modal-header">
                            <h2>Create New Project</h2>
                            <button
                                className="modal-close"
                                onClick={() => setShowCreateModal(false)}
                            >
                                ✕
                            </button>
                        </div>
                        <div className="modal-body">
                            <div className="form-group">
                                <label className="form-label">Project Name *</label>
                                <input
                                    type="text"
                                    className="form-input"
                                    value={newProjectName}
                                    onChange={(e) => setNewProjectName(e.target.value)}
                                    placeholder="e.g., Accra Basin Study 2024"
                                    autoFocus
                                />
                            </div>
                            <div className="form-group">
                                <label className="form-label">Description (optional)</label>
                                <textarea
                                    className="form-textarea"
                                    value={newProjectDesc}
                                    onChange={(e) => setNewProjectDesc(e.target.value)}
                                    placeholder="Brief description of this project..."
                                    rows={3}
                                />
                            </div>
                            <div className="info-box">
                                <span className="info-icon">ℹ️</span>
                                <p>All analysis results will be saved to this project and can be downloaded as a complete report.</p>
                            </div>
                        </div>
                        <div className="modal-footer">
                            <button
                                className="btn btn-secondary"
                                onClick={() => setShowCreateModal(false)}
                            >
                                Cancel
                            </button>
                            <button
                                className="btn btn-primary"
                                onClick={handleCreateProject}
                                disabled={!newProjectName.trim() || creating}
                            >
                                {creating ? 'Creating...' : 'Create Project'}
                            </button>
                        </div>
                    </div>
                </div>
            )}
        </div>
    )
}

export default Projects
