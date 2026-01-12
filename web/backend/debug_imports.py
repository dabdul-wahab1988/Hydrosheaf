
import sys
import os

print(f"Python executable: {sys.executable}")
print(f"CWD: {os.getcwd()}")
print("Attempting to import hydrosheaf...")

try:
    import hydrosheaf
    print(f"SUCCESS: hydrosheaf imported from {hydrosheaf.__file__}")
except ImportError as e:
    print(f"FAILURE: Could not import hydrosheaf. Error: {e}")
    print("sys.path:")
    for p in sys.path:
        print(f"  {p}")

print("Attempting to import app.hydrosheaf_adapter...")
try:
    sys.path.append(os.getcwd())
    from app import hydrosheaf_adapter
    print(f"SUCCESS: hydrosheaf_adapter imported")
except ImportError as e:
    print(f"FAILURE: Could not import hydrosheaf_adapter. Error: {e}")
