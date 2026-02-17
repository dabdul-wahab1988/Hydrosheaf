import { createContext, useContext, useState, useCallback } from 'react'
import Toast from '../components/Toast'

const ToastContext = createContext(null)

export function ToastProvider({ children }) {
    const [toasts, setToasts] = useState([])

    const addToast = useCallback((message, type = 'info', duration = 5000) => {
        const id = Date.now() + Math.random()
        setToasts(prev => [...prev, { id, message, type }])

        // Auto-remove after duration
        setTimeout(() => {
            setToasts(prev => prev.filter(t => t.id !== id))
        }, duration)
    }, [])

    const removeToast = useCallback((id) => {
        setToasts(prev => prev.filter(t => t.id !== id))
    }, [])

    const success = useCallback((msg, duration) => addToast(msg, 'success', duration), [addToast])
    const error = useCallback((msg, duration) => addToast(msg, 'error', duration || 7000), [addToast])
    const warning = useCallback((msg, duration) => addToast(msg, 'warning', duration), [addToast])
    const info = useCallback((msg, duration) => addToast(msg, 'info', duration), [addToast])

    return (
        <ToastContext.Provider value={{ success, error, warning, info }}>
            {children}
            <Toast toasts={toasts} onRemove={removeToast} />
        </ToastContext.Provider>
    )
}

export function useToast() {
    const context = useContext(ToastContext)
    if (!context) {
        throw new Error('useToast must be used within a ToastProvider')
    }
    return context
}
