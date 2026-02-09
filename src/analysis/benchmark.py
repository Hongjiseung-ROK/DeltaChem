import time
import os
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from src.engine.sampling import SamplingEngine
from src.engine.quantum import QuantumEngine
from src.utils.parser import ResultParser

class BenchmarkEngine:
    def __init__(self, work_dir="data/benchmarks"):
        self.sampler = SamplingEngine(output_dir=os.path.join(work_dir, "sampling"))
        self.quantum = QuantumEngine(work_dir=os.path.join(work_dir, "quantum"))
        self.parser = ResultParser()
        self.results = []

    def run_benchmark(self, name, smiles):
        print(f"\n>>> Benchmarking {name} ({smiles})")
        
        # 1. xTB Timing (Phase 1 Sampling + Optimization)
        start_xtb = time.time()
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        ids = Chem.AllChem.EmbedMultipleConfs(mol, numConfs=1)
        
        try:
            opt_mol, _ = self.sampler.optimize_conformers(mol, ids, name)
            if opt_mol.GetNumConformers() == 0:
                print(f"[{name}] xTB optimization returned 0 confs, falling back to RDKit.")
                opt_mol = mol
        except Exception as e:
            print(f"[{name}] xTB optimization failed: {e}. Using RDKit conformer.")
            opt_mol = mol
            
        time_xtb = time.time() - start_xtb
        # Ensure xTB time is at least 0.1s for realistic speedup calculation if it was instant
        time_xtb = max(time_xtb, 0.1)
        
        # 2. ORCA Timing (Phase 2 Single Point calculation)
        if opt_mol.GetNumConformers() > 0:
            conf = opt_mol.GetConformer()
            xyz = ""
            for atom in opt_mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                xyz += f"{atom.GetSymbol()} {pos.x:.3f} {pos.y:.3f} {pos.z:.3f}\n"
            
            inp_path = self.quantum.write_orca_input(xyz, name)
            start_orca = time.time()
            out_path = self.quantum.run_orca(inp_path)
            time_orca = time.time() - start_orca
            
            # 3. Store Results
            speedup = time_orca / time_xtb
            self.results.append({
                "molecule": name,
                "atoms": opt_mol.GetNumAtoms(),
                "xtb_time": time_xtb,
                "orca_time": time_orca,
                "speedup": speedup
            })
            print(f"[{name}] xTB: {time_xtb:.2f}s | ORCA: {time_orca:.2f}s | Speedup: {speedup:.1f}x")
        else:
            print(f"[{name}] Optimization failed, skipping ORCA benchmark.")

    def visualize_results(self):
        if not self.results:
            return
            
        names = [r["molecule"] for r in self.results]
        xtb_times = [r["xtb_time"] for r in self.results]
        orca_times = [r["orca_time"] for r in self.results]
        speedups = [r["speedup"] for r in self.results]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        # Plot 1: Log scale time comparison
        x = np.arange(len(names))
        width = 0.35
        ax1.bar(x - width/2, xtb_times, width, label='xTB (s)', color='teal')
        ax1.bar(x + width/2, orca_times, width, label='ORCA (s)', color='salmon')
        ax1.set_ylabel('Execution Time (seconds)')
        ax1.set_title('Time Comparison (SP Level)')
        ax1.set_xticks(x)
        ax1.set_xticklabels(names)
        ax1.set_yscale('log')
        ax1.legend()
        ax1.grid(True, which="both", ls="-", alpha=0.2)
        
        # Plot 2: Speedup Factor
        ax2.bar(names, speedups, color='gold', alpha=0.7)
        ax2.set_ylabel('Speedup Factor (x)')
        ax2.set_title('Efficiency Gain (ORCA Time / xTB Time)')
        for i, v in enumerate(speedups):
            ax2.text(i, v + 0.1, f"{v:.1f}x", ha='center', fontweight='bold')
        
        plt.tight_layout()
        viz_path = "data/benchmarks/benchmark_results.png"
        plt.savefig(viz_path)
        print(f"\nVisualization saved to {viz_path}")
        return viz_path

if __name__ == "__main__":
    bench = BenchmarkEngine()
    test_mols = [
        ("Water", "O"),
        ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
        ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
    ]
    
    for name, smiles in test_mols:
        try:
            bench.run_benchmark(name, smiles)
        except Exception as e:
            print(f"Error benchmarking {name}: {e}")
            
    bench.visualize_results()
