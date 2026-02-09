# Project ALCHEMIST - Delta-Ensemble

## The Mission
Execute a self-evolving project titled **"Delta-Ensemble: Uncertainty-Aware Boltzmann-Averaged Property Prediction via Delta-Learning."**

## Core Technical Architecture
1. **Delta-Learning Engine:** Train an E(3)-Equivariant GNN (or SchNet) to predict the correction factor between semi-empirical (GFN2-xTB) and DFT (B3LYP/6-31G*) energy levels.
2. **Ensemble Averaging:** Use RDKit and CREST/xTB to generate conformer pools. Predict Boltzmann weights for each conformer to calculate ensemble-averaged properties.
3. **Active Self-Optimization:** Implement a **Meta-Review Loop** after each phase.

## Operational Protocols
1. **Recursive Error Handling:** Analyze ORCA output and autonomously recover from SCF non-convergence.
2. **Resource Guardrails:** 
    - CleanupManager for ORCA files.
    - LazyDataLoading with h5py.
    - Parallelize xTB, restrict ORCA threads.
3. **Active Improvement Loop:** Sensitivity Analysis, MC-Dropout, Curriculum Learning.

## Execution Phases
**PHASE 0: Environment & Infrastructure**
- Construct project structure.
- Generate `requirements.txt`.
- Deliver `EnvironmentChecker` and `ProjectScaffolder` scripts.

**PHASE 1: High-Throughput Conformational Sampling**
- Pipeline: SMILES -> Conformer -> GFN2-xTB.
- Diversity Filter.

**PHASE 2: The Quantum Oracle (Delta-Labeling)**
- Automate ORCA .inp generation.
- cclib parsing.

**PHASE 3: Equivariant GNN Development**
- EGNN/SchNet model.
- MC-Dropout.

**PHASE 4: Boltzmann Averaging & Deployment**
- Inference engine.
- Final README with LaTeX equations.
