import os
import subprocess
from src.utils.config import Config

class QuantumEngine:
    def __init__(self, orca_path=None, work_dir=None):
        self.orca_path = orca_path or Config.ORCA_PATH
        self.work_dir = work_dir or Config.get_work_dir("orca_work")
        os.makedirs(self.work_dir, exist_ok=True)

    def write_orca_input(self, xyz_coords, name, charge=0, multiplicity=1, solvent=None):
        """Generates an ORCA input file for single-point energy with optional solvation."""
        input_path = os.path.join(self.work_dir, f"{name}.inp")
        
        # Determine Solvation Method (CPCM for water/common solvents)
        solv_str = ""
        if solvent:
            solv_str = f" CPCM({solvent})"
            
        with open(input_path, "w") as f:
            f.write(f"! B3LYP 6-31G* TightSCF SP{solv_str}\n")
            # f.write(f"%pal nprocs 4 end\n") # Commented out due to mpiexec issues
            f.write(f"* xyz {charge} {multiplicity}\n")
            f.write(xyz_coords)
            f.write(f"*\n")
        return input_path

    def run_orca(self, input_path):
        """Runs ORCA and handles basic recovery."""
        output_path = input_path.replace(".inp", ".out")
        print(f"Running ORCA for {os.path.basename(input_path)}...")
        
        try:
            with open(output_path, "w") as out:
                subprocess.run([self.orca_path, input_path], stdout=out, stderr=subprocess.STDOUT, check=True)
            return output_path
        except subprocess.CalledProcessError:
            print(f"ORCA failed for {input_path}. Analyzing output for recovery...")
            return self.recover_orca(input_path, output_path)

    def recover_orca(self, input_path, output_path):
        """Simple recovery: Retry with SlowConv if SCF failed."""
        with open(output_path, "r") as f:
            content = f.read()
        
        if "SCF ITERATIONS HAVE NOT CONVERGED" in content:
            print("SCF non-convergence detected. Retrying with SlowConv...")
            new_input = input_path.replace(".inp", "_retry.inp")
            with open(input_path, "r") as f_in, open(new_input, "w") as f_out:
                lines = f_in.readlines()
                for line in lines:
                    if line.startswith("!"):
                        f_out.write(line.strip() + " SlowConv\n")
                    else:
                        f_out.write(line)
            return self.run_orca(new_input)
        return None

if __name__ == "__main__":
    # xyz for water as a test
    h2o_xyz = "O 0.000 0.000 0.000\nH 0.000 0.757 0.586\nH 0.000 -0.757 0.586"
    engine = QuantumEngine()
    inp = engine.write_orca_input(h2o_xyz, "test_h2o")
    # Not running ORCA here to avoid long wait in demo, but structure is ready.
    print(f"ORCA input generated at: {inp}")
