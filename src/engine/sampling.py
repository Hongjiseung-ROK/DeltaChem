import os
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from tqdm import tqdm

from src.utils.config import Config

class SamplingEngine:
    def __init__(self, output_dir=None, xtb_path=None):
        self.output_dir = output_dir or Config.get_work_dir("sampling")
        self.xtb_path = xtb_path or Config.XTB_PATH
        os.makedirs(self.output_dir, exist_ok=True)

    def generate_conformers(self, smiles, name, num_confs=50, rmsd_threshold=0.5):
        """Generates conformers using ETKDG and applies a diversity filter."""
        print(f"[{name}] Generating {num_confs} conformers...")
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 0xf00d
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=ps)
        
        # Initial MMFF optimization for speed
        AllChem.MMFFOptimizeMoleculeConfs(mol)
        
        # Diversity Filter (RMSD-based)
        unique_conf_ids = []
        if conf_ids:
            unique_conf_ids.append(conf_ids[0])
            for i in range(1, len(conf_ids)):
                is_unique = True
                for j in unique_conf_ids:
                    rmsd = rdMolAlign.GetBestRMS(mol, mol, i, j)
                    if rmsd < rmsd_threshold:
                        is_unique = False
                        break
                if is_unique:
                    unique_conf_ids.append(i)
        
        print(f"[{name}] Diversity filter: {len(conf_ids)} -> {len(unique_conf_ids)} unique conformers.")
        return mol, unique_conf_ids

    def optimize_conformers(self, mol, conf_ids, name, solvent=None):
        """Optimizes selected conformers using GFN2-xTB with optional solvation."""
        optimized_mol = Chem.Mol(mol)
        optimized_mol.RemoveAllConformers()
        
        work_dir = os.path.join(self.output_dir, f"{name}_xtb")
        os.makedirs(work_dir, exist_ok=True)
        
        # Build solvation command
        solv_cmd = []
        if solvent:
            solv_cmd = ["--gbsa", solvent]
            
        final_confs = []
        for i, conf_id in enumerate(tqdm(conf_ids, desc=f"Optimizing {name} with xTB")):
            conf_file = os.path.join(work_dir, f"conf_{i}.sdf")
            writer = Chem.SDWriter(conf_file)
            writer.write(mol, confId=conf_id)
            writer.close()
            
            # Run xTB optimization
            try:
                # --opt for geometry optimization
                cmd = [self.xtb_path, conf_file, "--opt", "tight"] + solv_cmd
                subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True, check=True)
                
                # Load optimized geometry (xtbopt.sdf)
                opt_file = os.path.join(work_dir, "xtbopt.sdf")
                if os.path.exists(opt_file):
                    supplier = Chem.SDMolSupplier(opt_file, removeHs=False)
                    opt_mol = next(supplier)
                    if opt_mol:
                        optimized_mol.AddConformer(opt_mol.GetConformer(), assignId=True)
                        # Optionally extract energy from xtb output or sdf properties
                        final_confs.append(i)
            except Exception as e:
                print(f"[{name}] Optimization failed for conformer {i}: {e}")
                
        output_path = os.path.join(self.output_dir, f"{name}_optimized.sdf")
        writer = Chem.SDWriter(output_path)
        for conf in optimized_mol.GetConformers():
            writer.write(optimized_mol, confId=conf.GetId())
        writer.close()
        
        print(f"[{name}] Saved {optimized_mol.GetNumConformers()} optimized conformers to {output_path}")
        return optimized_mol, output_path

if __name__ == "__main__":
    engine = SamplingEngine()
    # Test with Caffeine
    mol, ids = engine.generate_conformers("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Caffeine", num_confs=5)
    engine.optimize_conformers(mol, ids, "Caffeine")
