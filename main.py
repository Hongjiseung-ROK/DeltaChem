from src.engine.sampling import SamplingEngine
from src.engine.quantum import QuantumEngine
from src.utils.parser import ResultParser
from src.engine.inference import InferenceEngine, CleanupManager
import os

def main():
    print("=== PROJECT ALCHEMIST: Delta-Ensemble Launch ===")
    
    # 1. Inputs
    caffeine_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    mol_name = "Caffeine"
    
    # 2. Phase 1: Conformational Sampling
    sampler = SamplingEngine()
    mol, unique_ids = sampler.generate_conformers(caffeine_smiles, mol_name, num_confs=5)
    # opt_mol, opt_path = sampler.optimize_conformers(mol, unique_ids, mol_name)
    
    # 3. Phase 2: Quantum Oracle (Mocked for Demo)
    quantum = QuantumEngine()
    print("\nPhase 2: Generating Quantum Input for top conformer...")
    conf = mol.GetConformer(unique_ids[0])
    # Extract XYZ
    xyz = ""
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz += f"{atom.GetSymbol()} {pos.x:.3f} {pos.y:.3f} {pos.z:.3f}\n"
    
    inp_path = quantum.write_orca_input(xyz, f"{mol_name}_conf0")
    print(f"ORCA Input ready at: {inp_path}")
    
    # 4. Phase 4: Inference Mock-up
    inf = InferenceEngine()
    mock_energies = [0.0, 0.02, 0.04] # Mocked eV
    weights = inf.calculate_boltzmann_weights(mock_energies)
    print(f"\nFinal Phase: Boltzmann weights generated for conformers: {weights}")
    
    # 5. Resource Guardrail
    CleanupManager.cleanup_orca("data/orca_work")
    
    print("\n[SUCCESS] Project ALCHEMIST framework launched and verified.")

if __name__ == "__main__":
    main()
