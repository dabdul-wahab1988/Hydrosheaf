import { Routes, Route } from 'react-router-dom'
import Layout from './components/Layout'
import Dashboard from './pages/Dashboard'
import Samples from './pages/Samples'
import Network from './pages/Network'
import Analysis from './pages/Analysis'
import Results from './pages/Results'
import Projects from './pages/Projects'

function App() {
  return (
    <Layout>
      <Routes>
        <Route path="/" element={<Dashboard />} />
        <Route path="/projects" element={<Projects />} />
        <Route path="/samples" element={<Samples />} />
        <Route path="/network" element={<Network />} />
        <Route path="/analysis" element={<Analysis />} />
        <Route path="/results" element={<Results />} />
      </Routes>
    </Layout>
  )
}

export default App
