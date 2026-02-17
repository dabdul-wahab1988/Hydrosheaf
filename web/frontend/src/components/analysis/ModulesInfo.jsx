import { memo } from 'react'

const modules = [
    { icon: '🔀', title: 'Transport Modeling', desc: 'Evaporation vs mixing selection' },
    { icon: '⚛️', title: 'Reaction Path', desc: 'LASSO sparse mineral reactions' },
    { icon: '🔗', title: 'Network Inference', desc: 'Probabilistic flow direction' },
    { icon: '📊', title: 'Uncertainty', desc: 'Bootstrap/Bayesian MCMC' },
    { icon: '🌡️', title: 'Isotopes', desc: 'δ18O/δ2H fractionation' },
    { icon: '🧪', title: 'PHREEQC', desc: 'Saturation index constraints' },
    { icon: '🌾', title: 'Nitrate Source', desc: 'Dual-isotope discrimination' },
    { icon: '⏱️', title: 'Temporal', desc: 'Residence time estimation' },
    { icon: '📐', title: '3D Network', desc: 'Layered aquifer systems' },
    { icon: '🔄', title: 'Reactive Transport', desc: 'Kinetic validation' },
    { icon: '💧', title: 'Vadose Zone', desc: 'Unsaturated flow modeling' },
    { icon: '📈', title: 'CoDA', desc: 'Compositional data analysis' },
]

function ModulesInfo() {
    return (
        <div className="card mt-lg">
            <div className="card-header">
                <h3 className="card-title">Available Modules</h3>
            </div>
            <div className="modules-info">
                {modules.map((module, idx) => (
                    <div key={idx} className="module-item">
                        <span className="module-icon">{module.icon}</span>
                        <div className="module-text">
                            <h4>{module.title}</h4>
                            <p>{module.desc}</p>
                        </div>
                    </div>
                ))}
            </div>
        </div>
    )
}

export default memo(ModulesInfo)
