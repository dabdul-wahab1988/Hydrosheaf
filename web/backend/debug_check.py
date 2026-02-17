"""
Comprehensive debugging script for Hydrosheaf web backend integration
Run this to identify issues BEFORE starting the server
"""

import sys
import traceback
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))


def print_section(title, char="="):
    """Print a formatted section header"""
    print(f"\n{char * 70}", flush=True)
    print(f"  {title}", flush=True)
    print(f"{char * 70}\n", flush=True)


def test_imports():
    """Test all required imports"""
    print_section("Testing Imports")

    issues = []

    # Test Hydrosheaf core
    print("1. Testing Hydrosheaf core imports...", flush=True)
    try:
        from hydrosheaf import Config, fit_network_pipeline

        print("   [OK] Hydrosheaf core imported successfully")
    except ImportError as e:
        print(f"   [X] Hydrosheaf core import failed: {e}")
        issues.append(("Hydrosheaf core", str(e)))

    # Test Hydrosheaf graph types
    print("\n2. Testing Hydrosheaf graph types...")
    try:
        from hydrosheaf.graph.types import Edge

        print("   [OK] Hydrosheaf Edge imported successfully")
        print(f"   Edge class: {Edge}")

        # Test creating an Edge
        test_edge = Edge(edge_id="A->B", u="A", v="B", attrs={"weight": 1.0})
        print(f"   [OK] Test edge created: {test_edge.u} -> {test_edge.v}")
    except Exception as e:
        print(f"   [X] Hydrosheaf Edge import/creation failed: {e}")
        issues.append(("Hydrosheaf Edge", str(e)))

    # Test adapter module
    print("\n3. Testing adapter module...")
    try:
        from app.hydrosheaf_adapter import ConfigAdapter, SampleAdapter, ResultAdapter

        print("   [OK] Adapter module imported successfully")
        print(f"   ConfigAdapter: {ConfigAdapter}")
        print(f"   SampleAdapter: {SampleAdapter}")
        print(f"   ResultAdapter: {ResultAdapter}")
    except Exception as e:
        print(f"   [X] Adapter module import failed: {e}")
        print(f"\n   Traceback:")
        traceback.print_exc()
        issues.append(("Adapter module", str(e)))

    # Test routers
    print("\n4. Testing analysis router...")
    try:
        from app.routers.analysis import router as analysis_router

        print(f"   [OK] Analysis router imported successfully")
        print(f"   Routes: {len(analysis_router.routes)} endpoints")
    except Exception as e:
        print(f"   [X] Analysis router import failed: {e}")
        print(f"\n   Traceback:")
        traceback.print_exc()
        issues.append(("Analysis router", str(e)))

    print("\n5. Testing network router...")
    try:
        from app.routers.network import router as network_router

        print(f"   [OK] Network router imported successfully")
        print(f"   Routes: {len(network_router.routes)} endpoints")
    except Exception as e:
        print(f"   [X] Network router import failed: {e}")
        print(f"\n   Traceback:")
        traceback.print_exc()
        issues.append(("Network router", str(e)))

    return issues


def test_adapter_functions():
    """Test adapter functions with sample data"""
    print_section("Testing Adapter Functions")

    issues = []

    try:
        from app.hydrosheaf_adapter import ConfigAdapter, SampleAdapter, ResultAdapter
        from hydrosheaf import Config

        # Test ConfigAdapter
        print("1. Testing ConfigAdapter...")
        frontend_config = {
            "lasso_penalty": 0.1,
            "enable_phreeqc": False,
            "enable_isotopes": True,
        }
        try:
            hydrosheaf_config = ConfigAdapter.frontend_to_hydrosheaf(frontend_config)
            print(f"   [OK] Config conversion successful")
            print(f"     lambda_l1: {hydrosheaf_config.lambda_l1}")
            print(f"     phreeqc_enabled: {hydrosheaf_config.phreeqc_enabled}")
            print(f"     isotope_enabled: {hydrosheaf_config.isotope_enabled}")

            # Verify types
            assert isinstance(hydrosheaf_config, Config)
            assert hydrosheaf_config.lambda_l1 == 0.1
            print("   [OK] Config validation passed")
        except Exception as e:
            print(f"   [X] ConfigAdapter failed: {e}")
            traceback.print_exc()
            issues.append(("ConfigAdapter", str(e)))

        # Test SampleAdapter
        print("\n2. Testing SampleAdapter...")
        frontend_samples = [
            {
                "sample_id": "W1",
                "location_id": "Well_1",
                "ca": 80.0,  # mg/L
                "mg": 30.0,
                "na": 40.0,
                "hco3": 250.0,
                "cl": 50.0,
                "so4": 70.0,
                "no3": 15.0,
                "f": 1.0,
                "fe": 0.1,
                "po4": 0.05,
                "ph": 7.2,
            }
        ]
        try:
            hydrosheaf_samples = SampleAdapter.frontend_to_hydrosheaf(frontend_samples)
            print(f"   [OK] Sample conversion successful")
            print(f"     Samples converted: {len(hydrosheaf_samples)}")

            sample = hydrosheaf_samples[0]
            print(f"     site_id: {sample.get('site_id')}")
            print(
                f"     Ca (mmol/L): {sample.get('Ca', 'N/A'):.4f}"
            )  # Should be converted from mg/L
            print(f"     pH: {sample.get('ph')}")

            # Verify conversion
            assert "site_id" in sample
            assert "Ca" in sample
            assert sample["Ca"] < 10  # Should be in mmol/L, not mg/L
            print("   [OK] Sample validation passed")
        except Exception as e:
            print(f"   [X] SampleAdapter failed: {e}")
            traceback.print_exc()
            issues.append(("SampleAdapter", str(e)))

        # Test ResultAdapter with mock EdgeResult
        print("\n3. Testing ResultAdapter...")
        try:
            from hydrosheaf.inference.edge_fit import EdgeResult

            # Create mock edge results
            mock_results = [
                EdgeResult(
                    edge_id="W1->W2",
                    u="W1",
                    v="W2",
                    transport_model="evap",
                    gamma=1.15,
                    f=None,
                    z_labels=["Calcite", "Gypsum"],
                    z_extents=[0.084, -0.012],
                    transport_residual_norm=0.023,
                    anomaly_norm=0.015,
                )
            ]

            frontend_results = ResultAdapter.hydrosheaf_to_frontend(mock_results)
            print(f"   [OK] Result conversion successful")
            print(
                f"     Dominant process: {frontend_results['transport_model']['dominant_process']}"
            )
            print(f"     Reactions found: {len(frontend_results['reactions'])}")
            print(f"     Edges: {len(frontend_results['edges'])}")

            # Verify structure
            assert "transport_model" in frontend_results
            assert "reactions" in frontend_results
            assert "edges" in frontend_results
            assert frontend_results["edges"][0]["source"] == "W1"
            assert frontend_results["edges"][0]["target"] == "W2"
            print("   [OK] Result validation passed")
        except Exception as e:
            print(f"   [X] ResultAdapter failed: {e}")
            traceback.print_exc()
            issues.append(("ResultAdapter", str(e)))

    except ImportError as e:
        print(f"[X] Could not import adapters: {e}")
        issues.append(("Adapter import", str(e)))

    return issues


def test_hydrosheaf_integration():
    """Test actual Hydrosheaf function calls"""
    print_section("Testing Hydrosheaf Integration")

    issues = []

    try:
        from hydrosheaf import (
            Config,
            fit_network_pipeline,
            auto_disable_missing_modules,
        )
        from hydrosheaf.graph.types import Edge
        from app.hydrosheaf_adapter import SampleAdapter, ConfigAdapter

        print("1. Preparing test data...")
        # Create test samples
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
        ]

        # Convert to Hydrosheaf format
        hydrosheaf_samples = SampleAdapter.frontend_to_hydrosheaf(test_samples)
        print(f"   [OK] Converted {len(hydrosheaf_samples)} samples")

        # Create config
        config = Config(
            lambda_l1=0.1,
            phreeqc_enabled=False,
            isotope_enabled=False,
        )
        print(f"   [OK] Config created")

        # Auto-disable modules
        config = auto_disable_missing_modules(hydrosheaf_samples, config)
        print(f"   [OK] Auto-disable completed")
        print(f"     PHREEQC: {config.phreeqc_enabled}")
        print(f"     Isotopes: {config.isotope_enabled}")

        # Create edges
        edges = [
            Edge(
                edge_id="Well_1->Well_2", u="Well_1", v="Well_2", attrs={"weight": 1.0}
            )
        ]
        print(f"   [OK] Created {len(edges)} edges")

        # Run actual Hydrosheaf analysis!
        print("\n2. Running Hydrosheaf analysis...")
        try:
            edge_results, extras = fit_network_pipeline(
                hydrosheaf_samples,
                edges,
                config,
                auto_disable_missing=True,
            )
            print(f"   [OK] Analysis completed successfully!")
            print(f"     Edge results: {len(edge_results)}")

            if edge_results:
                result = edge_results[0]
                print(f"     Edge: {result.edge_id}")
                print(f"     Transport model: {result.transport_model}")
                print(f"     Gamma: {result.gamma}")
                print(f"     Reactions: {len(result.z_labels)}")

                # Verify result structure
                assert hasattr(result, "u")
                assert hasattr(result, "v")
                assert hasattr(result, "z_labels")
                assert hasattr(result, "z_extents")
                print("   [OK] Result structure validated")
        except Exception as e:
            print(f"   [X] Hydrosheaf analysis failed: {e}")
            traceback.print_exc()
            issues.append(("Hydrosheaf analysis", str(e)))

    except Exception as e:
        print(f"[X] Integration test failed: {e}")
        traceback.print_exc()
        issues.append(("Hydrosheaf integration", str(e)))

    return issues


def test_api_models():
    """Test API model compatibility"""
    print_section("Testing API Models")

    issues = []

    try:
        from app.routers.analysis import AnalysisRequest, AnalysisConfig
        from app.routers.network import NetworkData, Node, EdgeAPI

        print("1. Testing AnalysisRequest model...")
        try:
            request_data = {
                "name": "Test Analysis",
                "analysis_type": "full_pipeline",
                "samples": [{"sample_id": "W1", "location_id": "Well_1", "ca": 80.0}],
                "config": {
                    "lasso_penalty": 0.1,
                    "enable_phreeqc": False,
                },
            }
            request = AnalysisRequest(**request_data)
            print(f"   [OK] AnalysisRequest created successfully")
            print(f"     Name: {request.name}")
            print(f"     Type: {request.analysis_type}")
        except Exception as e:
            print(f"   [X] AnalysisRequest failed: {e}")
            issues.append(("AnalysisRequest", str(e)))

        print("\n2. Testing NetworkData model...")
        try:
            network_data = {
                "name": "Test Network",
                "nodes": [{"id": "A", "name": "Well A", "x": 0, "y": 0}],
                "edges": [{"source": "A", "target": "B"}],
            }
            network = NetworkData(**network_data)
            print(f"   [OK] NetworkData created successfully")
            print(f"     Name: {network.name}")
            print(f"     Nodes: {len(network.nodes)}")
        except Exception as e:
            print(f"   [X] NetworkData failed: {e}")
            issues.append(("NetworkData", str(e)))

    except Exception as e:
        print(f"[X] API model test failed: {e}")
        traceback.print_exc()
        issues.append(("API models", str(e)))

    return issues


def main():
    """Run all debugging checks"""
    print("\n" + "[SEARCH]" * 35, flush=True)
    print("  HYDROSHEAF WEB BACKEND - COMPREHENSIVE DEBUG CHECK")
    print("[SEARCH]" * 5)

    all_issues = []

    # Run all tests
    all_issues.extend(test_imports())
    all_issues.extend(test_adapter_functions())
    all_issues.extend(test_hydrosheaf_integration())
    all_issues.extend(test_api_models())

    # Print summary
    print_section("Debug Summary", "=")

    if not all_issues:
        print("All CHECKS PASSED!")
        print("\nNo issues found. The integration is working correctly!")
        print("\nNext steps:")
        print("  1. Start the backend: uvicorn app.main:app --reload")
        print("  2. Run integration tests: python test_integration.py")
        return 0
    else:
        print(f"[X] FOUND {len(all_issues)} ISSUE(S)")
        print("\nIssues found:")
        for i, (component, error) in enumerate(all_issues, 1):
            print(f"\n{i}. {component}:")
            print(f"   {error[:200]}...")  # Truncate long errors

        print("\nRecommendations:")
        if any(
            "Hydrosheaf" in issue[0] and "import" in issue[1].lower()
            for issue in all_issues
        ):
            print("  - Install Hydrosheaf: pip install ../../")
        if any("Adapter" in issue[0] for issue in all_issues):
            print("  - Check adapter module for syntax errors")
        if any("router" in issue[0].lower() for issue in all_issues):
            print("  - Check router modules for import conflicts")

        return 1


if __name__ == "__main__":
    sys.exit(main())
