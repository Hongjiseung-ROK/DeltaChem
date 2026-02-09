import os
import pandas as pd
import numpy as np
from src.engine.validation import ValidationEngine
from src.engine.sampling import SamplingEngine
from rdkit import Chem
from src.utils.cleanup import CleanupManager

class AccuracyStudy:
    def __init__(self):
        self.validator = ValidationEngine(work_dir="data/accuracy_study/quantum")
        self.sampler = SamplingEngine(output_dir="data/accuracy_study/sampling")
        self.cleanup = CleanupManager()
        self.results = []
        
        # 1 kcal/mol = 0.04336 eV
        self.EV_TO_KCAL = 23.0605

    def run_study(self, molecules):
        print("=== DeltaChem Accuracy Verification Study ===")
        print(f"Targeting {len(molecules)} molecules for baseline MAE analysis.\n")
        
        for name, smiles in molecules:
            print(f"[*] Processing: {name}")
            try:
                # 1. Generate 1 optimized conformer
                mol = Chem.MolFromSmiles(smiles)
                mol = Chem.AddHs(mol)
                ids = Chem.AllChem.EmbedMultipleConfs(mol, numConfs=1)
                
                xtb_energy_ev = 0.0
                try:
                    opt_mol, _ = self.sampler.optimize_conformers(mol, ids, name)
                    if opt_mol.GetNumConformers() == 0:
                        print(f"    [!] xTB optimization returned 0 confs, falling back to RDKit geometry.")
                        opt_mol = mol
                    else:
                        # Success: Parse xTB energy
                        if opt_mol.HasProp("REMARK"):
                            import re
                            match = re.search(r"energy:\s+(-?\d+\.\d+)", opt_mol.GetProp("REMARK"))
                            if match:
                                xtb_energy_ev = float(match.group(1)) * 27.2114
                except Exception as e:
                    print(f"    [!] xTB failed: {e}. Using RDKit geometry.")
                    opt_mol = mol
                
                # 2. Run ORCA Validation
                conf = opt_mol.GetConformer()
                xyz = ""
                for atom in opt_mol.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    xyz += f"{atom.GetSymbol()} {pos.x:.3f} {pos.y:.3f} {pos.z:.3f}\n"
                
                val_result = self.validator.run_full_validation(xyz, name, xtb_energy_ev)
                
                if val_result:
                    diff_ev = abs(val_result["delta_e"])
                    diff_kcal = diff_ev * self.EV_TO_KCAL
                    
                    self.results.append({
                        "Molecule": name,
                        "E_xTB (eV)": xtb_energy_ev,
                        "E_DFT (eV)": val_result["e_dft"],
                        "Delta_E (eV)": val_result["delta_e"],
                        "MAE (kcal/mol)": diff_kcal
                    })
                    print(f"    -> Delta_E: {val_result['delta_e']:.4f} eV | Error: {diff_kcal:.2f} kcal/mol")
                
                # Cleanup huge ORCA files
                self.cleanup.cleanup_orca_temporary_files("data/accuracy_study/quantum")
                
            except Exception as e:
                print(f"    [!] Error processing {name}: {e}")

    def report(self):
        if not self.results:
            print("No results to report.")
            return
            
        df = pd.DataFrame(self.results)
        df.to_csv("data/accuracy_study/accuracy_report.csv", index=False)
        
        avg_mae = df["MAE (kcal/mol)"].mean()
        std_mae = df["MAE (kcal/mol)"].std()
        
        print("\n" + "="*40)
        print("FINAL ACCURACY VERIFICATION REPORT")
        print("="*40)
        print(df[["Molecule", "MAE (kcal/mol)"]].to_string(index=False))
        print("-"*40)
        print(f"BASELINE MAE (xTB vs DFT): {avg_mae:.2f} ± {std_mae:.2f} kcal/mol")
        print("="*40)
        print("Interpretation: A Delta-Learning model must correct this ~50-100 kcal/mol gap")
        print("to reach Chemical Accuracy (< 1.0 kcal/mol).")
        
        return avg_mae

if __name__ == "__main__":
    study = AccuracyStudy()
    test_set = [
        ("Methane", "C"),
        ("Water", "O"),
        ("Ethane", "CC"),
        ("Methanol", "CO"),
        ("Formaldehyde", "C=O")
    ]
    study.run_study(test_set)
    study.report()
