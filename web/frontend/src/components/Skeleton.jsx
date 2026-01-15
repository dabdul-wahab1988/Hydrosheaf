import './Skeleton.css'

/**
 * Base skeleton component with shimmer animation
 */
export function Skeleton({ width, height, variant = 'rect', className = '' }) {
    return (
        <div
            className={`skeleton skeleton-${variant} ${className}`}
            style={{ width, height }}
            aria-hidden="true"
        />
    )
}

/**
 * Skeleton for stat cards on dashboard
 */
export function StatCardSkeleton() {
    return (
        <div className="skeleton-stat-card">
            <Skeleton width="40px" height="40px" variant="circle" />
            <div className="skeleton-stat-content">
                <Skeleton width="60%" height="14px" />
                <Skeleton width="40%" height="24px" />
            </div>
        </div>
    )
}

/**
 * Skeleton for table rows
 */
export function TableRowSkeleton({ cols = 5 }) {
    return (
        <div className="skeleton-table-row">
            {Array(cols).fill(0).map((_, i) => (
                <Skeleton key={i} width={`${Math.random() * 30 + 50}%`} height="16px" />
            ))}
        </div>
    )
}

/**
 * Skeleton for a full table
 */
export function TableSkeleton({ rows = 5, cols = 5 }) {
    return (
        <div className="skeleton-table">
            <div className="skeleton-table-header">
                {Array(cols).fill(0).map((_, i) => (
                    <Skeleton key={i} width="80%" height="14px" />
                ))}
            </div>
            {Array(rows).fill(0).map((_, i) => (
                <TableRowSkeleton key={i} cols={cols} />
            ))}
        </div>
    )
}

/**
 * Skeleton for content cards
 */
export function CardSkeleton() {
    return (
        <div className="skeleton-card">
            <div className="skeleton-card-header">
                <Skeleton width="60%" height="20px" />
                <Skeleton width="80px" height="28px" />
            </div>
            <div className="skeleton-card-body">
                <Skeleton width="100%" height="12px" />
                <Skeleton width="90%" height="12px" />
                <Skeleton width="75%" height="12px" />
            </div>
        </div>
    )
}

/**
 * Skeleton for chart placeholders
 */
export function ChartSkeleton({ height = '200px' }) {
    return (
        <div className="skeleton-chart" style={{ height }}>
            <div className="skeleton-chart-bars">
                {Array(7).fill(0).map((_, i) => (
                    <Skeleton
                        key={i}
                        width="12%"
                        height={`${Math.random() * 60 + 20}%`}
                        className="skeleton-bar"
                    />
                ))}
            </div>
        </div>
    )
}

/**
 * Full page loading skeleton
 */
export function PageSkeleton() {
    return (
        <div className="skeleton-page" role="status" aria-label="Loading content">
            <Skeleton width="200px" height="32px" className="skeleton-title" />
            <Skeleton width="400px" height="16px" className="skeleton-subtitle" />
            <div className="skeleton-grid">
                <StatCardSkeleton />
                <StatCardSkeleton />
                <StatCardSkeleton />
                <StatCardSkeleton />
            </div>
            <CardSkeleton />
        </div>
    )
}

export default Skeleton
