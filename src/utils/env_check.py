import subprocess
import sys
import shutil

def check_python_dependencies():
    print("--- Checking Python Dependencies ---")
    
    # Try to fix torch DLL issue on Windows
    if sys.platform == "win32":
        import os
        torch_lib_path = os.path.join(os.path.dirname(sys.executable), "Lib", "site-packages", "torch", "lib")
        if os.path.exists(torch_lib_path):
            os.add_dll_directory(torch_lib_path)
            print(f"[INFO] Added DLL directory: {torch_lib_path}")

    dependencies = [
        "rdkit", "torch", "torch_geometric", "cclib", "h5py", "numpy", "pandas"
    ]
    # xtb-python is excluded here and will be handled via subprocess if needed
    all_found = True
    for dep in dependencies:
        try:
            __import__(dep.replace("-", "_"))
            print(f"[OK] {dep} is installed.")
        except ImportError as e:
            if dep == "torch" or dep == "torch_geometric":
                print(f"[WARNING] {dep} failed to load (DLL issue): {e}")
                # We can still proceed with Phase 1/2 without torch
            else:
                print(f"[ERROR] {dep} is NOT installed or failed to load: {e}")
                all_found = False
    return all_found

def check_external_binaries():
    print("\n--- Checking External Binaries ---")
    local_bin_root = r"c:\workspace\222_cc_project\orca_bin"
    xtb_path = os.path.join(local_bin_root, "xtb-6.7.1pre", "xtb.exe")
    orca_path = os.path.join(local_bin_root, "orca.exe")
    
    binaries = {
        "orca": orca_path,
        "xtb": xtb_path
    }
    
    all_found = True
    for bin_name, bin_path in binaries.items():
        if os.path.exists(bin_path):
            print(f"[OK] {bin_name} found at: {bin_path}")
            try:
                if bin_name == "xtb":
                    result = subprocess.run([bin_path, "--version"], capture_output=True, text=True)
                    print(f"     Version: {result.stdout.splitlines()[0]}")
                elif bin_name == "orca":
                    print(f"     {bin_name} is found.")
            except Exception as e:
                print(f"     [WARNING] Could not run {bin_name}: {e}")
        else:
            # Check if in PATH
            path = shutil.which(bin_name)
            if path:
                print(f"[OK] {bin_name} found in PATH at: {path}")
            else:
                print(f"[ERROR] {bin_name} NOT found at {bin_path} or in PATH.")
                all_found = False
    return all_found

def main():
    print("=== Delta-Ensemble Environment Checker ===\n")
    py_ok = check_python_dependencies()
    bin_ok = check_external_binaries()
    
    if py_ok and bin_ok:
        print("\n[SUCCESS] Environment is ready for Project ALCHEMIST.")
    else:
        print("\n[FAILED] Environment setup is incomplete.")
        sys.exit(1)

if __name__ == "__main__":
    main()
