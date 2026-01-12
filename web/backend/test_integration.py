"""
Test script to verify Hydrosheaf integration with web backend
Run this after starting the backend server to verify everything works
"""

import requests
import pytest
import time
import json
import sys

API_BASE = "http://localhost:8000/api"


def _backend_available() -> bool:
    try:
        res = requests.get(f"{API_BASE}/health", timeout=2)
        return res.status_code < 500
    except requests.RequestException:
        return False


def print_section(title):
    """Print a section header"""
    print("\n" + "=" * 60)
    print(f"  {title}")
    print("=" * 60)


def run_health_checks():
    """Run health check endpoints"""
    print_section("Testing Health Checks")

    # Main health check
    try:
        res = requests.get(f"{API_BASE}/health")
        print(f"   [OK] Main API health: {res.status_code}")
        print(f"  {json.dumps(res.json(), indent=2)}")
    except Exception as e:
        print(f"   [X] Main API health failed: {e}")
        return False

    # Analysis health check
    try:
        res = requests.get(f"{API_BASE}/analysis/health")
        health = res.json()
        print(f"\n   [OK] Analysis router health: {res.status_code}")
        print(f"  Hydrosheaf available: {health.get('hydrosheaf_available')}")
        print(f"  Status: {health.get('status')}")
        print(f"  Message: {health.get('message')}")

        if not health.get('hydrosheaf_available'):
            print("\n   [WARN] WARNING: Hydrosheaf not installed!")
            print("   Install with: cd web/backend && pip install ../../")
            return False
    except Exception as e:
        print(f"   [X] Analysis health failed: {e}")
        return False

    # Network health check
    try:
        res = requests.get(f"{API_BASE}/network/health")
        health = res.json()
        print(f"\n   [OK] Network router health: {res.status_code}")
        print(f"  Hydrosheaf available: {health.get('hydrosheaf_available')}")
    except Exception as e:
        print(f"   [X] Network health failed: {e}")

    return True


def run_full_analysis_workflow():
    """Run complete analysis workflow with real Hydrosheaf"""
    print_section("Testing Full Analysis Workflow")

    # Prepare test samples (mg/L units, will be converted to mmol/L by adapter)
    test_samples = [
        {
            "sample_id": "W1",
            "location_id": "Well_1",
            "ca": 80.0,
            "mg": 30.0,
            "na": 40.0,
            "k": 5.0,
            "hco3": 250.0,
            "so4": 70.0,
            "cl": 50.0,
            "no3": 15.0,
            "f": 1.0,
            "fe": 0.1,
            "po4": 0.05,
            "ph": 7.2,
        },
        {
            "sample_id": "W2",
            "location_id": "Well_2",
            "ca": 95.0,
            "mg": 35.0,
            "na": 50.0,
            "k": 6.0,
            "hco3": 280.0,
            "so4": 85.0,
            "cl": 65.0,
            "no3": 22.0,
            "f": 1.2,
            "fe": 0.15,
            "po4": 0.08,
            "ph": 7.0,
        },
        {
            "sample_id": "W3",
            "location_id": "Well_3",
            "ca": 110.0,
            "mg": 40.0,
            "na": 60.0,
            "k": 7.0,
            "hco3": 310.0,
            "so4": 100.0,
            "cl": 80.0,
            "no3": 28.0,
            "f": 1.5,
            "fe": 0.2,
            "po4": 0.1,
            "ph": 6.8,
        },
    ]

    print(f"\n1. Uploading test dataset with {len(test_samples)} samples...")
    try:
        upload_res = requests.post(
            f"{API_BASE}/samples/upload",
            json={"name": "Integration Test Dataset", "samples": test_samples}
        )
        upload_res.raise_for_status()
        dataset = upload_res.json()
        print(f"   [OK] Dataset uploaded: {dataset['name']} ({dataset['sample_count']} samples)")
    except Exception as e:
        print(f"   [X] Upload failed: {e}")
        return False

    print("\n2. Starting Hydrosheaf analysis...")
    analysis_req = {
        "name": "Integration Test Analysis",
        "analysis_type": "full_pipeline",
        "samples": test_samples,
        "config": {
            "lasso_penalty": 0.1,
            "enable_phreeqc": False,  # Disable PHREEQC for simplicity
            "enable_isotopes": False,  # No isotope data
            "enable_uncertainty": False,  # Disable for faster testing
        }
    }

    try:
        analysis_res = requests.post(
            f"{API_BASE}/analysis/run",
            json=analysis_req
        )
        analysis_res.raise_for_status()
        job = analysis_res.json()
        job_id = job['job_id']
        print(f"   [OK] Analysis job created: {job_id}")
        print(f"   Hydrosheaf available: {job.get('hydrosheaf_available')}")

        if not job.get('hydrosheaf_available'):
            print("   [X] Cannot proceed without Hydrosheaf!")
            return False

    except Exception as e:
        print(f"   [X] Analysis start failed: {e}")
        return False

    print("\n3. Polling for completion (max 60 seconds)...")
    for i in range(60):
        time.sleep(1)
        try:
            status_res = requests.get(f"{API_BASE}/analysis/status/{job_id}")
            status = status_res.json()
            status_str = status['status']

            if i % 5 == 0:  # Print every 5 seconds
                print(f"   [{i+1}s] Status: {status_str}")

            if status_str == 'completed':
                print(f"\n   [OK] Analysis completed in {i+1} seconds!")
                break
            elif status_str == 'failed':
                print(f"\n   [X] Analysis failed: {status.get('error')}")
                if 'error_preview' in status:
                    print(f"\n   Error preview:\n{status['error_preview']}")
                return False
        except Exception as e:
            print(f"   [X] Status check failed: {e}")
            return False
    else:
        print("\n   [X] Analysis timed out after 60 seconds")
        return False

    print("\n4. Fetching results...")
    try:
        results_res = requests.get(f"{API_BASE}/analysis/results/{job_id}")
        results_res.raise_for_status()
        results = results_res.json()

        print(f"   [OK] Results retrieved successfully!")
        print(f"\n   Summary:")
        print(f"     Total samples: {results['summary']['total_samples']}")
        print(f"     Total edges: {results['summary']['total_edges']}")

        print(f"\n   Transport Model:")
        tm = results['transport_model']
        print(f"     Dominant process: {tm['dominant_process']}")
        print(f"     Average gamma: {tm.get('average_gamma', 'N/A')}")

        print(f"\n   Reactions found: {len(results['reactions'])}")
        for i, rxn in enumerate(results['reactions'][:5], 1):  # Show top 5
            print(f"     {i}. {rxn['mineral']}: {rxn['rate_mmol_L']:.4f} mmol/L ({rxn['direction']})")

        # Check metadata to verify real Hydrosheaf was used
        if 'metadata' in results:
            metadata = results['metadata']
            print(f"\n   Metadata:")
            print(f"     Engine: {metadata.get('analysis_engine')}")
            print(f"     Mock data: {metadata.get('mock_data')}")

            if metadata.get('mock_data', True):
                print("\n   [WARN] WARNING: Results are still mock data!")
                return False

        print("\n   [OK] Real Hydrosheaf analysis confirmed!")

    except Exception as e:
        print(f"   [X] Results fetch failed: {e}")
        return False

    return True


def run_network_inference():
    """Run network flow inference"""
    print_section("Testing Network Flow Inference")

    # Create test network
    network_data = {
        "name": "Test Flow Network",
        "nodes": [
            {"id": "A", "name": "Well A", "x": 0, "y": 0, "hydraulic_head": 50.0},
            {"id": "B", "name": "Well B", "x": 100, "y": 50, "hydraulic_head": 48.0},
            {"id": "C", "name": "Well C", "x": 150, "y": 100, "hydraulic_head": 45.0},
        ],
        "edges": [
            {"source": "A", "target": "B"},
            {"source": "B", "target": "C"},
        ]
    }

    print("\n1. Creating network...")
    try:
        create_res = requests.post(
            f"{API_BASE}/network/create",
            json=network_data
        )
        create_res.raise_for_status()
        network = create_res.json()
        network_id = network['network_id']
        print(f"   [OK] Network created: {network_id}")
    except Exception as e:
        print(f"   [X] Network creation failed: {e}")
        return False

    print("\n2. Running flow inference...")
    try:
        infer_res = requests.post(
            f"{API_BASE}/network/{network_id}/infer-flow"
        )
        infer_res.raise_for_status()
        inference = infer_res.json()

        print(f"   [OK] Inference complete!")
        print(f"     Method: {inference['method']}")
        print(f"     Hydrosheaf used: {inference.get('hydrosheaf_used', False)}")

        print(f"\n     Inferred edges:")
        for edge in inference['inferred_edges']:
            print(f"       {edge['source']} → {edge['target']}: "
                  f"P={edge.get('flow_probability', 'N/A'):.2f} "
                  f"({edge.get('flow_direction', 'unknown')})")

    except Exception as e:
        print(f"   [X] Inference failed: {e}")
        return False

    return True


def test_health_checks():
    """Test health check endpoints"""
    if not _backend_available():
        pytest.skip("Backend not reachable; start the server to run integration tests.")
    assert run_health_checks()


def test_full_analysis_workflow():
    """Test complete analysis workflow with real Hydrosheaf"""
    if not _backend_available():
        pytest.skip("Backend not reachable; start the server to run integration tests.")
    assert run_full_analysis_workflow()


def test_network_inference():
    """Test network flow inference"""
    if not _backend_available():
        pytest.skip("Backend not reachable; start the server to run integration tests.")
    assert run_network_inference()


def main():
    """Run all tests"""
    print("\n" + "[SEARCH]" * 30)
    print("  HYDROSHEAF WEB BACKEND INTEGRATION TEST")
    print("[SEARCH]" * 30)

    print("\nMake sure the backend server is running:")
    print("  cd web/backend")
    print("  uvicorn app.main:app --reload")

    input("\nPress Enter to start tests...")

    # Run tests
    results = []

    results.append(("Health Checks", run_health_checks()))

    if results[-1][1]:  # Only continue if Hydrosheaf is available
        results.append(("Analysis Workflow", run_full_analysis_workflow()))
        results.append(("Network Inference", run_network_inference()))
    else:
        print("\n[WARN]  Skipping remaining tests (Hydrosheaf not available)")

    # Print summary
    print_section("Test Summary")
    for name, passed in results:
        status = "[OK] PASSED" if passed else "[X] FAILED"
        print(f"  {status}: {name}")

    all_passed = all(r[1] for r in results)

    if all_passed:
        print("\n[OK] ALL TESTS PASSED! [OK]")
        print("\nThe Hydrosheaf engine is successfully integrated!")
        return 0
    else:
        print("\n[X] SOME TESTS FAILED")
        print("\nCheck the errors above and ensure:")
        print("  1. Backend server is running")
        print("  2. Hydrosheaf is installed: pip install ../../")
        print("  3. All dependencies are installed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
