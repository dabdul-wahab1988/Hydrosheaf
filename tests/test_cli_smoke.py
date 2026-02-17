import subprocess
import sys
import shutil

def test_hydrosheaf_cli_help():
    """Test that the main CLI entry point runs and shows help."""
    # We run it as a subprocess to simulate real usage
    cmd = [sys.executable, "-m", "hydrosheaf.cli", "--help"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    assert result.returncode == 0
    # Check for common help output indicators, case-insensitive
    stdout_lower = result.stdout.lower()
    assert "usage:" in stdout_lower or "options:" in stdout_lower or "arguments:" in stdout_lower

def test_hydrosheaf_cal_cli_help():
    """Test that the calibration CLI entry point runs and shows help."""
    cmd = [sys.executable, "-m", "hydrosheaf.calibration.cli", "--help"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    assert result.returncode == 0
    stdout_lower = result.stdout.lower()
    assert "usage:" in stdout_lower or "options:" in stdout_lower or "arguments:" in stdout_lower
