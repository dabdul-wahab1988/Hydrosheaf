"""
Starter YAML configuration templates for all supported Hydrosheaf calibration types.
"""

KINETIC_TEMPLATE = """# Hydrosheaf Calibration Config: Kinetic
calibration:
  type: kinetic
  settings:
    n_workers: 2
    max_iterations: 30
    output_dir: calibration_results
    engine: internal # "internal" (PESTGLM) or "pestpp-glm" or "pestpp-ies"
  model:
    minerals:
      - calcite
      - dolomite
    observations_file: observations.csv
    temp_default_c: 25.0
    residence_time_days: 10.0
    initial_solution:
      pH: 7.0
      Ca: 0.1
      Mg: 0.05
    fit_parameters:
      - "calcite:k:global"
      - "dolomite:k:global"
"""

TRANSPORT_TEMPLATE = """# Hydrosheaf Calibration Config: Transport
calibration:
  type: transport
  settings:
    n_workers: 1
    max_iterations: 50
    output_dir: calibration_results
    engine: internal
  parameters:
    # Fixed parameter example:
    # - name: decay
    #   initial: 0.0
    #   bounds: [0.0, 1.0]
    #   fixed: true
    # Tied parameter example:
    # - name: velocity
    #   initial: 0.1
    #   bounds: [0.01, 1.0]
    #   tied_to: dispersivity
    - name: dispersivity
      initial: 1.0
      bounds: [0.1, 10.0]
      log: true
  model:
    observations_file: transport_obs.csv
    distance_m: 10.0
    velocity: 0.1
    dispersivity: 1.0
    decay: 0.0
    fit_parameters:
      - dispersivity
      - decay
"""

VADOSE_TEMPLATE = """# Hydrosheaf Calibration Config: Vadose
calibration:
  type: vadose
  settings:
    n_workers: 1
    max_iterations: 50
    output_dir: calibration_results
    engine: internal
  model:
    config_file: profile.json
    forcing_file: forcing.csv
    observations_file: vadose_obs.csv
    dz: 0.1
    kc: 1.0
    layers_to_fit:
      - 0
    fit_parameters:
      - ks_L0
      - alpha_L0
      - n_L0
      - kc
      - preferential_flow_fraction
      - ttd_cv
"""

NITRATE_TEMPLATE = """# Hydrosheaf Calibration Config: Nitrate
calibration:
  type: nitrate
  settings:
    n_workers: 1
    max_iterations: 50
    output_dir: calibration_results
    engine: internal
  model:
    samples_file: nodes.csv
    target_file: targets.csv
    fit_parameters:
      - prior_logit
      - w1_no3_cl
      - w2_no3_k
      - w3_po4
      - w4_fe
      - w5_denitrif
      - w6_alk_coupling
      - nitrate_source_min_mg_L
      - denitrification_min_extent
      - min_top_probability
      - min_top_margin
"""

AGE_TEMPORAL_TEMPLATE = """# Hydrosheaf Calibration Config: Age/Temporal
calibration:
  type: age_temporal
  settings:
    n_workers: 1
    max_iterations: 50
    output_dir: calibration_results
    engine: internal
  model:
    observations_file: age_obs.csv
    porosity: 0.2
    velocity: 0.1
    decay: 0.001
    ttd_cv: 0.7
    fit_parameters:
      - tau
      - decay
      - velocity
      - porosity
      - ttd_cv
"""

TOPOLOGY_TEMPLATE = """# Hydrosheaf Calibration Config: Flow Topology
calibration:
  type: topology
  settings:
    n_workers: 1
    max_iterations: 50
    output_dir: calibration_results
    engine: internal
  model:
    edges_file: candidate_edges.csv
    observations_file: topology_obs.csv
    prior_sigma: 2.0
    normalize_by_upstream: false
"""

COMPOSITE_TEMPLATE = """# Hydrosheaf Calibration Config: Composite Joint Calibration
calibration:
  type: composite
  settings:
    n_workers: 2
    max_iterations: 50
    output_dir: calibration_results
    engine: internal
  parameters:
    - name: dispersivity
      initial: 1.0
      bounds: [0.1, 10.0]
      log: true
    - name: decay
      initial: 0.0
      bounds: [0.0, 1.0]
      fixed: true               # do not adjust during calibration
    - name: ks_multiplier
      initial: 1.0
      bounds: [0.5, 2.0]
      tied_to: dispersivity     # tied to another parameter
  sub_models:
    - type: transport
      id: sub_transport
      observations_file: transport_obs.csv
      distance_m: 10.0
      fit_parameters:
        - dispersivity
    - type: vadose
      id: sub_vadose
      profile_file: profile.json
      forcing_file: forcing.csv
      observations_file: vadose_obs.csv
      layers_to_fit:
        - 0
      fit_parameters:
        - ks_multiplier
"""

TEMPLATES = {
    "kinetic": KINETIC_TEMPLATE,
    "transport": TRANSPORT_TEMPLATE,
    "vadose": VADOSE_TEMPLATE,
    "nitrate": NITRATE_TEMPLATE,
    "age": AGE_TEMPORAL_TEMPLATE,
    "temporal": AGE_TEMPORAL_TEMPLATE,
    "age_temporal": AGE_TEMPORAL_TEMPLATE,
    "topology": TOPOLOGY_TEMPLATE,
    "composite": COMPOSITE_TEMPLATE,
}
