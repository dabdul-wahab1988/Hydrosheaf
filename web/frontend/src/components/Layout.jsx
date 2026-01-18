import { useState, useEffect } from 'react'
import { NavLink, useLocation } from 'react-router-dom'
import { useProject } from '../context/ProjectContext'
import './Layout.css'

const navItems = [
    { path: '/', label: 'Dashboard', icon: '📊' },
    { path: '/projects', label: 'Projects', icon: '📁' },
    { path: '/samples', label: 'Samples', icon: '🧪' },
    { path: '/network', label: 'Network', icon: '🔗' },
    { path: '/analysis', label: 'Analysis', icon: '⚗️' },
    { path: '/results', label: 'Results', icon: '📈' },
    { path: '/demo', label: 'Research Demo', icon: '🎓' },
]

function Layout({ children }) {
    const { currentProject } = useProject()
    const [sidebarOpen, setSidebarOpen] = useState(false)
    const location = useLocation()

    // Close sidebar on route change (mobile)
    useEffect(() => {
        setSidebarOpen(false)
    }, [location.pathname])

    // Close sidebar on Escape key
    useEffect(() => {
        const handleEscape = (e) => {
            if (e.key === 'Escape') setSidebarOpen(false)
        }
        window.addEventListener('keydown', handleEscape)
        return () => window.removeEventListener('keydown', handleEscape)
    }, [])

    return (
        <div className="layout">
            {/* Mobile hamburger button */}
            <button
                className="sidebar-toggle"
                onClick={() => setSidebarOpen(!sidebarOpen)}
                aria-label={sidebarOpen ? 'Close menu' : 'Open menu'}
                aria-expanded={sidebarOpen}
            >
                <span className={`hamburger-icon ${sidebarOpen ? 'open' : ''}`} />
            </button>

            {/* Mobile overlay */}
            {sidebarOpen && (
                <div
                    className="sidebar-overlay"
                    onClick={() => setSidebarOpen(false)}
                    aria-hidden="true"
                />
            )}

            {/* Sidebar */}
            <aside className={`sidebar ${sidebarOpen ? 'open' : ''}`}>
                <div className="sidebar-header">
                    <div className="logo">
                        <span className="logo-icon">💧</span>
                        <div className="logo-text">
                            <h1>Hydrosheaf</h1>
                            <span className="logo-subtitle">Hydrogeochemistry Framework</span>
                        </div>
                    </div>
                </div>

                {/* Current Project Status */}
                {currentProject ? (
                    <div className="current-project-badge">
                        <span className="project-status-dot active"></span>
                        <div className="project-status-info">
                            <span className="project-status-label">Active Project</span>
                            <span className="project-status-name">{currentProject.name}</span>
                        </div>
                    </div>
                ) : (
                    <NavLink to="/projects" className="no-project-badge">
                        <span className="project-status-dot inactive"></span>
                        <div className="project-status-info">
                            <span className="project-status-label">No Project</span>
                            <span className="project-status-name">Results not saved</span>
                        </div>
                    </NavLink>
                )}

                <nav className="sidebar-nav">
                    {navItems.map((item) => (
                        <NavLink
                            key={item.path}
                            to={item.path}
                            className={({ isActive }) =>
                                `nav-item ${isActive ? 'active' : ''}`
                            }
                        >
                            <span className="nav-icon">{item.icon}</span>
                            <span className="nav-label">{item.label}</span>
                        </NavLink>
                    ))}
                </nav>

                <div className="sidebar-footer">
                    <div className="creators-section">
                        <span className="creators-label">Created by</span>
                        <div className="creators-list">
                            <span className="creator-name">Dickson Abdul-Wahab</span>
                            <span className="creator-name">Ebenezer Aquisman Asare</span>
                            <span className="creator-name">Abdul Rashid Dickson</span>
                        </div>
                    </div>
                    <div className="version-info">
                        <span className="version-badge">v1.0.0</span>
                        <span className="version-text">Sheaf-Theoretic Method</span>
                    </div>
                </div>
            </aside>

            {/* Main Content */}
            <main className="main-content">
                <div className="content-wrapper">
                    {children}
                </div>
            </main>
        </div>
    )
}

export default Layout
