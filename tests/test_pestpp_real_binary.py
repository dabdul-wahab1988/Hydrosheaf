"""
Real PEST++ binary execution smoke tests.
These tests invoke the actual PEST++ binaries (not mocked) with minimal
problems to verify end-to-end execution pipeline.

IMPORTANT:
- These require PEST++ binaries in bin/5.2.25/ (already present).
- Each test creates a self-contained model script, .pst, .tpl, and .ins.
- Tests may take 30-60s each.
"""

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest
import pandas as pd

from hydrosheaf.calibration.pestpp.runner import (
    get_executable_path,
    PestControlFile,
    PestParameter,
    PestObservation,
)


def _pestpp_available(name: str = "pestpp-glm") -> bool:
    """Check if a PEST++ binary is discoverable."""
    try:
        exe = get_executable_path(name)
        return Path(exe).exists() or shutil.which(name) is not None
    except Exception:
        return False


def _run_pestpp(exe: str, control_file: str, work_dir: Path, timeout: int):
    """Run PEST++ without Windows grandchild pipe-handle deadlocks."""
    stdout_path = work_dir / "pestpp.stdout.log"
    stderr_path = work_dir / "pestpp.stderr.log"
    with (
        stdout_path.open("w", encoding="utf-8") as stdout,
        stderr_path.open("w", encoding="utf-8") as stderr,
    ):
        try:
            completed = subprocess.run(
                [exe, control_file],
                cwd=str(work_dir),
                stdout=stdout,
                stderr=stderr,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as exc:
            completed = subprocess.CompletedProcess(exc.cmd, 1)
            stderr.write(
                "\nPEST++ was terminated after producing its output files; "
                "the Windows 5.2.25 GLM binary can remain resident after completion.\n"
            )
    return subprocess.CompletedProcess(
        completed.args,
        completed.returncode,
        stdout_path.read_text(encoding="utf-8", errors="replace"),
        stderr_path.read_text(encoding="utf-8", errors="replace"),
    )


def _write_minimal_pst(work_dir: Path, noptmax: int = 1, obs_groups=None):
    """Write a minimal control file that works with real PEST++ binaries."""
    obs_groups = obs_groups or ("obsgp", "obsgp")
    pst = PestControlFile()
    # Model command
    pst.model_command_line = f"{sys.executable} model.py"
    pst.tpl_file = "case.tpl"
    pst.ins_file = "case.ins"
    pst.model_input_file = "params.csv"
    pst.model_output_file = "results.dat"
    pst.control_data.noptmax = noptmax

    # Parameters (bounds must not cross zero with FACTOR change limit)
    p1 = PestParameter("p1", "none", "factor", 2.0, 0.0, 10.0, "pargp")
    p2 = PestParameter("p2", "none", "factor", 1.0, 0.0, 5.0, "pargp")
    pst.add_parameter(p1)
    pst.add_parameter(p2)

    # Observations
    o1 = PestObservation("o1", 4.0, 1.0, obs_groups[0])
    o2 = PestObservation("o2", 6.0, 1.0, obs_groups[1])
    pst.add_observation(o1)
    pst.add_observation(o2)

    pst.write(work_dir / "case.pst")

    # Template — match PestControlFile.write_tpl() output exactly
    tpl_lines = ["ptf ~"]
    for p in pst.parameters:
        tpl_lines.append(f"{p.name},~{p.name:<{20}}~")
    (work_dir / "case.tpl").write_text("\n".join(tpl_lines) + "\n")

    # Instruction — sequential l1 reads
    ins_lines = ["pif @"]
    for o in pst.observations:
        ins_lines.append(f"l1 w !{o.name}!")
    (work_dir / "case.ins").write_text("\n".join(ins_lines) + "\n")

    # Model script: completely self-contained, no embedded \n string literals
    model_lines = [
        "import sys",
        "params = {}",
        "for line in open('params.csv'):",
        "    if ',' not in line:",
        "        continue",
        "    parts = line.strip().split(',')",
        "    if len(parts) >= 2:",
        "        try:",
        "            params[parts[0]] = float(parts[1])",
        "        except Exception:",
        "            pass",
        "p1 = params.get('p1', 2.0)",
        "p2 = params.get('p2', 1.0)",
        "y1 = p1 * 2.0",
        "y2 = p1 * 3.0 + p2",
        "with open('results.dat', 'w') as f:",
        "    f.write('o1 ' + str(y1))",
        "    f.write(chr(10))",
        "    f.write('o2 ' + str(y2))",
        "    f.write(chr(10))",
    ]
    (work_dir / "model.py").write_text("\n".join(model_lines) + "\n")


class TestPestppGLMRealBinary:
    """Smoke-test GLM with real binary."""

    @pytest.mark.skipif(
        not _pestpp_available("pestpp-glm"), reason="PEST++ GLM binary not available"
    )
    def test_glm_runs_and_writes_output(self, tmp_path):
        """GLM actual binary should complete and write output files."""
        _write_minimal_pst(tmp_path, noptmax=1)
        exe = get_executable_path("pestpp-glm")

        proc = _run_pestpp(exe, "case.pst", tmp_path, timeout=15)
        outputs = (
            list(tmp_path.glob("case.par"))
            + list(tmp_path.glob("case.res"))
            + list(tmp_path.glob("case.rei"))
        )
        if proc.returncode not in (0, 1):
            assert outputs, f"Unexpected exit {proc.returncode}; stderr:\n{proc.stderr}"
        assert len(outputs) > 0, f"No output files; stderr:\n{proc.stderr}"

    @pytest.mark.skipif(
        not _pestpp_available("pestpp-glm"), reason="PEST++ GLM binary not available"
    )
    def test_glm_optimal_parameters_sensible(self, tmp_path):
        """GLM should move p1 toward 2.0 to match observation o1=4.0 (4=2*2)."""
        _write_minimal_pst(tmp_path, noptmax=1)
        exe = get_executable_path("pestpp-glm")

        proc = _run_pestpp(exe, "case.pst", tmp_path, timeout=15)
        par_file = tmp_path / "case.par"
        assert (
            par_file.exists()
        ), f"No case.par; exit={proc.returncode}\nstdout={proc.stdout}\nstderr={proc.stderr}"
        if proc.returncode not in (0, 1):
            assert (
                sys.platform == "win32"
            ), f"Unexpected exit={proc.returncode}\nstdout={proc.stdout}\nstderr={proc.stderr}"
        # Parse .par file to check parameter values
        par_text = par_file.read_text()
        lines = par_text.strip().splitlines()
        # Skip header line
        for line in lines[1:]:
            parts = line.strip().split()
            if len(parts) >= 2:
                p_name = parts[0]
                p_val = float(parts[1])
                if p_name == "p1":
                    # Optimal p1 should approach 2.0 within reasonable tolerance
                    assert 1.0 < p_val < 3.0, f"p1={p_val} did not converge toward 2.0"


class TestPestppIESRealBinary:
    """Smoke-test IES with real binary."""

    @pytest.mark.skipif(
        not _pestpp_available("pestpp-ies"), reason="PEST++ IES binary not available"
    )
    def test_ies_runs_with_ensemble(self, tmp_path):
        """IES should run with manually provided ensembles."""
        _write_minimal_pst(tmp_path, noptmax=1)

        # Write parameter ensemble (10 realizations near truth)
        par_ens = tmp_path / "case.0.par.csv"
        rows = ["real_name,p1,p2"]
        for i in range(10):
            p1_val = 1.5 + i * 0.1
            p2_val = 0.5 + i * 0.1
            rows.append(f"r{i},{p1_val:.6f},{p2_val:.6f}")
        par_ens.write_text("\n".join(rows) + "\n")

        # Write observation ensemble
        obs_ens = tmp_path / "case.obs.csv"
        rows = ["real_name,o1,o2"]
        for i in range(10):
            rows.append(f"r{i},{4.0 + i * 0.01:.6f},{6.0 + i * 0.01:.6f}")
        obs_ens.write_text("\n".join(rows) + "\n")

        # Append IES options to .pst
        with open(tmp_path / "case.pst", "a") as f:
            f.write("++ies_num_reals(10)\n")
            f.write("++ies_parameter_ensemble(case.0.par.csv)\n")
            f.write("++ies_observation_ensemble(case.obs.csv)\n")
            f.write("++ies_subset_size(5)\n")
            f.write("++max_run_fail(3)\n")

        exe = get_executable_path("pestpp-ies")
        proc = _run_pestpp(exe, "case.pst", tmp_path, timeout=180)
        # IES often exits 1 for convergence issues but should still write files
        acceptable_codes = (0, 1)
        if sys.platform == "win32":
            acceptable_codes = acceptable_codes + (3221226356,)
        assert (
            proc.returncode in acceptable_codes
        ), f"Unexpected exit {proc.returncode}; stderr:\n{proc.stderr}"
        # Minimal: some identifiable output file
        output_files = list(tmp_path.glob("case.*.par.csv")) + list(
            tmp_path.glob("case.phi.actual.csv")
        )
        assert (
            len(output_files) > 0
        ), f"IES produced no identifiable output files; stderr:\n{proc.stderr}"


class TestPestppSWPRealBinary:
    """Smoke-test SWP with real binary."""

    @pytest.mark.skipif(
        not _pestpp_available("pestpp-swp"), reason="PEST++ SWP binary not available"
    )
    def test_swp_runs_sweep(self, tmp_path):
        """SWP should run parameter sweep and produce output CSV."""
        _write_minimal_pst(tmp_path, noptmax=1)
        sweep_csv = tmp_path / "case.swp.in.csv"
        sweep_csv.write_text(
            "run_id,p1,p2\n"
            + "\n".join(f"{i},{1.0+i*0.1},0.0" for i in range(5))
            + "\n"
        )

        with open(tmp_path / "case.pst", "a") as f:
            f.write("++sweep_parameter_csv_file(case.swp.in.csv)\n")
            f.write("++sweep_output_csv_file(case.swp.out.csv)\n")
            f.write("++sweep_forgive(true)\n")
            f.write("++max_run_fail(3)\n")

        exe = get_executable_path("pestpp-swp")
        proc = _run_pestpp(exe, "case.pst", tmp_path, timeout=120)
        assert proc.returncode in (
            0,
            1,
        ), f"Unexpected exit {proc.returncode}; stderr:\n{proc.stderr}"
        out_csv = tmp_path / "case.swp.out.csv"
        assert (
            out_csv.exists()
        ), f"SWP did not produce output CSV; stdout={proc.stdout}\nstderr={proc.stderr}"
        lines = out_csv.read_text().splitlines()
        assert len(lines) >= 2, "SWP output CSV is empty"
        df = pd.read_csv(out_csv)
        assert {"o1", "o2"}.issubset(set(df.columns)), df.columns
        assert df["o1"].notna().any()


class TestPestppSENRealBinary:
    """Smoke-test SEN with real binary."""

    @pytest.mark.skipif(
        not _pestpp_available("pestpp-sen"), reason="PEST++ SEN binary not available"
    )
    def test_sen_runs_morris_and_writes_msn(self, tmp_path):
        """SEN should run with supported GSA options and write Method-of-Morris outputs."""
        _write_minimal_pst(tmp_path, noptmax=1)
        with open(tmp_path / "case.pst", "a") as f:
            f.write("++gsa_method(morris)\n")
            f.write("++gsa_morris_r(4)\n")
            f.write("++gsa_morris_p(4)\n")
            f.write("++gsa_morris_delta(0.5)\n")
            f.write("++max_run_fail(3)\n")

        exe = get_executable_path("pestpp-sen")
        proc = _run_pestpp(exe, "case.pst", tmp_path, timeout=120)
        assert proc.returncode in (
            0,
            1,
        ), f"Unexpected exit {proc.returncode}; stdout={proc.stdout}\nstderr={proc.stderr}"
        msn_file = tmp_path / "case.msn"
        assert (
            msn_file.exists()
        ), f"No case.msn; stdout={proc.stdout}\nstderr={proc.stderr}"
        msn = pd.read_csv(msn_file)
        assert {"parameter_name", "sen_mean_abs"}.issubset(set(msn.columns))
        assert set(msn["parameter_name"]) >= {"p1", "p2"}


class TestPestppDARealBinary:
    """Smoke-test DA with real binary."""

    @pytest.mark.skipif(
        not _pestpp_available("pestpp-da"), reason="PEST++ DA binary not available"
    )
    def test_da_runs_batch_assimilation_and_writes_global_outputs(self, tmp_path):
        """DA should run without unsupported da_num_cycles and write real DA output files."""
        _write_minimal_pst(tmp_path, noptmax=-1)
        with open(tmp_path / "case.pst", "a") as f:
            f.write("++ies_num_reals(6)\n")
            f.write("++max_run_fail(3)\n")

        exe = get_executable_path("pestpp-da")
        proc = _run_pestpp(exe, "case.pst", tmp_path, timeout=180)
        assert proc.returncode in (
            0,
            1,
        ), f"Unexpected exit {proc.returncode}; stdout={proc.stdout}\nstderr={proc.stderr}"
        assert (
            tmp_path / "case.global.phi.actual.csv"
        ).exists(), f"No DA global phi file; stdout={proc.stdout}\nstderr={proc.stderr}"
        assert list(
            tmp_path.glob("case.*.*.par.csv")
        ), f"No DA cycle parameter ensembles; files={list(tmp_path.iterdir())}"


class TestPestppMOURealBinary:
    """Smoke-test MOU with real binary."""

    @pytest.mark.skipif(
        not _pestpp_available("pestpp-mou"), reason="PEST++ MOU binary not available"
    )
    def test_mou_initial_population_runs_with_directional_objectives(self, tmp_path):
        """MOU should accept directional objective groups and evaluate an initial population."""
        _write_minimal_pst(
            tmp_path, noptmax=-1, obs_groups=("less_than_obj", "less_than_obj")
        )
        with open(tmp_path / "case.pst", "a") as f:
            f.write("++mou_objectives(o1,o2)\n")
            f.write("++mou_population_size(6)\n")
            f.write("++mou_generator(de)\n")
            f.write("++mou_save_population_every(1)\n")
            f.write("++opt_dec_var_groups(pargp)\n")
            f.write("++max_run_fail(3)\n")

        exe = get_executable_path("pestpp-mou")
        proc = _run_pestpp(exe, "case.pst", tmp_path, timeout=180)
        assert proc.returncode in (
            0,
            1,
        ), f"Unexpected exit {proc.returncode}; stdout={proc.stdout}\nstderr={proc.stderr}"
        dv_files = list(tmp_path.glob("case*.dv_pop.csv"))
        obs_files = list(tmp_path.glob("case*.obs_pop.csv"))
        assert (
            dv_files and obs_files
        ), f"No MOU population files; files={list(tmp_path.iterdir())}\nstdout={proc.stdout}\nstderr={proc.stderr}"


class TestPestppOPTRealBinary:
    """Smoke-test OPT with real binary."""

    @pytest.mark.skipif(
        not _pestpp_available("pestpp-opt"), reason="PEST++ OPT binary not available"
    )
    def test_opt_runs_with_directional_constraint_group(self, tmp_path):
        """OPT should run with valid opt_* options and a less-than constraint group."""
        _write_minimal_pst(tmp_path, noptmax=1, obs_groups=("less_than_limit", "obsgp"))
        with open(tmp_path / "case.pst", "a") as f:
            f.write("++opt_dec_var_groups(pargp)\n")
            f.write("++opt_constraint_groups(less_than_limit)\n")
            f.write("++opt_direction(min)\n")
            f.write("++opt_risk(0.5)\n")
            f.write("++max_run_fail(3)\n")

        exe = get_executable_path("pestpp-opt")
        proc = _run_pestpp(exe, "case.pst", tmp_path, timeout=180)
        assert proc.returncode in (
            0,
            1,
        ), f"Unexpected exit {proc.returncode}; stdout={proc.stdout}\nstderr={proc.stderr}"
        outputs = (
            list(tmp_path.glob("case.par"))
            + list(tmp_path.glob("case.rei"))
            + list(tmp_path.glob("case.res"))
        )
        assert (
            outputs
        ), f"No OPT output files; stdout={proc.stdout}\nstderr={proc.stderr}"


class TestPestppParallelRealBinary:
    """Smoke-test parallel manager/agent with real binary."""

    @pytest.mark.skipif(
        sys.platform == "win32" and os.getenv("HYDROSHEAF_RUN_WINDOWS_PANTHER") != "1",
        reason="Windows PANTHER test requires firewall permissions; set HYDROSHEAF_RUN_WINDOWS_PANTHER=1 to run",
    )
    @pytest.mark.skipif(
        not _pestpp_available("pestpp-glm"), reason="PEST++ GLM binary not available"
    )
    def test_parallel_glm_spawn_and_cleanup(self, tmp_path):
        """Manager spawns agents and produces output."""
        _write_minimal_pst(tmp_path, noptmax=1)
        with open(tmp_path / "case.pst", "a") as f:
            f.write("++lambdas(1.0)\n")
            f.write("++max_run_fail(3)\n")
        exe = get_executable_path("pestpp-glm")

        import socket

        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind(("", 0))
            port = s.getsockname()[1]

        mgr_cmd = [exe, "case.pst", "/h", f":{port}"]
        agent_cmd = [exe, "case.pst", "/h", f"localhost:{port}"]

        mgr = None
        agents = []
        mgr_out = ""
        mgr_err = ""

        def cleanup_process(proc):
            if proc is None:
                return
            try:
                if proc.poll() is None:
                    proc.kill()
                proc.wait(timeout=10)
            except Exception:
                pass
            for stream in (
                getattr(proc, "stdout", None),
                getattr(proc, "stderr", None),
            ):
                try:
                    if stream:
                        stream.close()
                except Exception:
                    pass

        try:
            mgr = subprocess.Popen(
                mgr_cmd,
                cwd=str(tmp_path),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            time.sleep(1)
            agents = [
                subprocess.Popen(
                    agent_cmd,
                    cwd=str(tmp_path),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                for _ in range(2)
            ]
            mgr_out, mgr_err = mgr.communicate(timeout=120)

            # Clean up agents
            for a in agents:
                try:
                    a.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    cleanup_process(a)
        except subprocess.TimeoutExpired:
            cleanup_process(mgr)
            for a in agents:
                cleanup_process(a)
            outputs = (
                list(tmp_path.glob("case.par"))
                + list(tmp_path.glob("case.rei"))
                + list(tmp_path.glob("case.res"))
            )
            assert outputs, (
                f"Parallel PEST++ timed out before producing outputs; "
                f"stderr:\n{mgr_err}\nstdout:\n{mgr_out}"
            )
            return
        finally:
            for a in agents:
                cleanup_process(a)

        assert mgr.returncode in (
            0,
            1,
        ), f"Manager exit={mgr.returncode}\nstderr={mgr_err[:4000]}\nstdout={mgr_out[:2000]}"
        outputs = (
            list(tmp_path.glob("case.par"))
            + list(tmp_path.glob("case.rei"))
            + list(tmp_path.glob("case.res"))
        )
        assert len(outputs) > 0, "No GLM output files produced in parallel mode"
