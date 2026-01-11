import { NavLink } from 'react-router-dom'
import { useProject } from '../context/ProjectContext'
import './Layout.css'

const navItems = [
    { path: '/', label: 'Dashboard', icon: '📊' },
    { path: '/projects', label: 'Projects', icon: '📁' },
    { path: '/samples', label: 'Samples', icon: '🧪' },
    { path: '/network', label: 'Network', icon: '🔗' },
    { path: '/analysis', label: 'Analysis', icon: '⚗️' },
    { path: '/results', label: 'Results', icon: '📈' },
]

function Layout({ children }) {
    const { currentProject } = useProject()

    return (
        <div className="layout">
            {/* Sidebar */}
            <aside className="sidebar">
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
