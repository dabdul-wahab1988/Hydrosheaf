import './Toast.css'

const icons = {
    success: '✓',
    error: '✕',
    warning: '⚠',
    info: 'ℹ',
}

function Toast({ toasts, onRemove }) {
    if (toasts.length === 0) return null

    return (
        <div className="toast-container" role="region" aria-label="Notifications">
            {toasts.map((toast) => (
                <div
                    key={toast.id}
                    className={`toast toast-${toast.type}`}
                    role="alert"
                    aria-live="polite"
                >
                    <span className="toast-icon" aria-hidden="true">
                        {icons[toast.type]}
                    </span>
                    <span className="toast-message">{toast.message}</span>
                    <button
                        className="toast-close"
                        onClick={() => onRemove(toast.id)}
                        aria-label="Dismiss notification"
                    >
                        ✕
                    </button>
                </div>
            ))}
        </div>
    )
}

export default Toast
