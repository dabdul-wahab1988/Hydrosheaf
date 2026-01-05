from hydrosheaf import Config, fit_network

samples = [
    {
        "sample_id": "s1",
        "site_id": "A",
        "Ca": 1.0,
        "Mg": 1.0,
        "Na": 0.5,
        "HCO3": 2.0,
        "Cl": 0.5,
        "SO4": 1.0,
        "NO3": 0.2,
        "F": 0.1,
        "EC": 100.0,
        "TDS": 200.0,
        "pH": 7.0,
    },
    {
        "sample_id": "s2",
        "site_id": "B",
        "Ca": 1.2,
        "Mg": 1.1,
        "Na": 0.6,
        "HCO3": 2.4,
        "Cl": 0.6,
        "SO4": 1.2,
        "NO3": 0.2,
        "F": 0.1,
        "EC": 105.0,
        "TDS": 210.0,
        "pH": 7.1,
    },
]

edges = [("A", "B")]

config = Config()
results = fit_network(samples, edges, config)

for result in results:
    print(result)
