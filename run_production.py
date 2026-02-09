import os
import time
from src.engine.sampling import SamplingEngine
from src.engine.quantum import QuantumEngine
from src.engine.data_loader import QM9SubsetLoader
from src.engine.inference import CleanupManager

def run_qm9_production():
    print("=== PROJECT ALCHEMIST: QM9 Production Launch ===")
    
    # 1. Load Data
    loader = QM9SubsetLoader()
    molecules = loader.get_smiles_subset()
    print(f"Loaded {len(molecules)} molecules from QM9 subset.\n")
    
    # 2. Initialize Engines
    sampler = SamplingEngine(output_dir="data/production/sampling")
    quantum = QuantumEngine(work_dir="data/production/orca")
    
    results = []
    
    # 3. Batch Processing Loop
    for name, smiles in molecules:
        print(f"\n>>> Processing {name} ({smiles})")
        start_time = time.time()
        
        try:
            # Phase 1: Sampling & Diversity Filtering
            # Lowering num_confs for production demo speed
            mol, unique_ids = sampler.generate_conformers(smiles, name, num_confs=10, rmsd_threshold=0.5)
            
            # Phase 2: Quantum Input Generation
            if unique_ids:
                conf = mol.GetConformer(unique_ids[0])
                xyz = ""
                for atom in mol.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    xyz += f"{atom.GetSymbol()} {pos.x:.3f} {pos.y:.3f} {pos.z:.3f}\n"
                
                inp_path = quantum.write_orca_input(xyz, f"{name}_top_conf")
                
                processing_time = time.time() - start_time
                results.append({
                    "name": name,
                    "status": "SUCCESS",
                    "unique_conformers": len(unique_ids),
                    "orca_input": inp_path,
                    "time": f"{processing_time:.2f}s"
                })
            else:
                results.append({"name": name, "status": "FAILED", "reason": "No conformers generated"})
                
        except Exception as e:
            print(f"Error processing {name}: {e}")
            results.append({"name": name, "status": "ERROR", "error": str(e)})

    # 4. Final Report
    print("\n" + "="*50)
    print("QM9 PRODUCTION SUMMARY REPORT")
    print("="*50)
    print(f"{'Molecule':<15} | {'Status':<8} | {'Confs':<6} | {'Time':<6}")
    print("-" * 50)
    for res in results:
        status = res.get("status", "N/A")
        confs = res.get("unique_conformers", 0)
        p_time = res.get("time", "N/A")
        print(f"{res['name']:<15} | {status:<8} | {confs:<6} | {p_time:<6}")
    print("="*50)
    
    # 5. Cleanup
    CleanupManager.cleanup_orca("data/production/orca")
    print("\n[FINISH] Production run completed.")

if __name__ == "__main__":
    run_qm9_production()
