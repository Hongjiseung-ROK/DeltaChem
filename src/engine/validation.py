from src.engine.quantum import QuantumEngine
from src.utils.parser import ResultParser
import os

class ValidationEngine:
    def __init__(self, work_dir="data/validation"):
        self.quantum = QuantumEngine(work_dir=work_dir)
        self.parser = ResultParser()
        self.work_dir = work_dir

    def run_full_validation(self, xyz_coords, name, xtb_energy_ev, solvent=None):
        """Runs ORCA for a molecule and compares it with xTB energy."""
        print(f"--- Quantum Validation: {name} (Solvent: {solvent or 'Gas'}) ---")
        
        # 1. ORCA Calculation
        inp_path = self.quantum.write_orca_input(xyz_coords, name, solvent=solvent)
        out_path = self.quantum.run_orca(inp_path)
        
        if out_path and os.path.exists(out_path):
            # 2. Parse Results
            orca_results = self.parser.parse_orca(out_path)
            
            if orca_results:
                e_dft = orca_results["energy_ev"]
                delta_e = e_dft - xtb_energy_ev
                
                # Report error metrics in kcal/mol for better visibility
                error_kcal = abs(delta_e) * 23.0605
                
                print(f"[RESULTS] {name}")
                print(f"Energy (xTB): {xtb_energy_ev:.4f} eV")
                print(f"Energy (DFT): {e_dft:.4f} eV")
                print(f"Delta-Energy: {delta_e:.4f} eV ({error_kcal:.2f} kcal/mol)")
                
                return {
                    "name": name,
                    "e_xtb": xtb_energy_ev,
                    "e_dft": e_dft,
                    "delta_e": delta_e,
                    "error_kcal": error_kcal
                }
        return None

if __name__ == "__main__":
    # Test with H2O
    h2o_xyz = "O 0.0000 0.0000 0.1173\nH 0.0000 0.7572 -0.4692\nH 0.0000 -0.7572 -0.4692"
    val = ValidationEngine()
    # Mock xTB energy for H2O (actual is approx -76 Hartrees, around -2000 eV)
    # We'll just run ORCA to verify the connection works.
    res = val.run_full_validation(h2o_xyz, "H2O_Val", xtb_energy_ev=-2080.5)
