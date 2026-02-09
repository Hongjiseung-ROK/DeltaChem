import numpy as np
import os
import shutil

class InferenceEngine:
    def __init__(self, model=None, k_b=8.617333262145e-5): # Boltzmann constant in eV/K
        self.model = model
        self.k_b = k_b

    def calculate_boltzmann_weights(self, energies_ev, temperature=298.15):
        """Calculates normalized Boltzmann weights for a set of energies."""
        energies = np.array(energies_ev)
        # Numerical stability: subtract min energy
        energies -= np.min(energies)
        
        exp_factors = np.exp(-energies / (self.k_b * temperature))
        weights = exp_factors / np.sum(exp_factors)
        return weights

    def ensemble_average(self, properties, weights):
        """Calculates weighted average of properties."""
        return np.sum(np.array(properties) * np.array(weights))

class CleanupManager:
    @staticmethod
    def cleanup_orca(work_dir):
        """Removes temporary ORCA files to save space."""
        extensions = [".gbw", ".tmp", ".hess", ".pc", ".prop", ".opt"]
        count = 0
        for root, dirs, files in os.walk(work_dir):
            for file in files:
                if any(file.endswith(ext) for ext in extensions):
                    try:
                        os.remove(os.path.join(root, file))
                        count += 1
                    except Exception as e:
                        print(f"Failed to remove {file}: {e}")
        print(f"Cleanup finished. Removed {count} temporary files.")

if __name__ == "__main__":
    # Test Boltzmann
    inf = InferenceEngine()
    test_energies = [0.0, 0.05, 0.1] # eV
    weights = inf.calculate_boltzmann_weights(test_energies)
    print(f"Boltzmann weights at 298K: {weights}")
    
    test_props = [10.5, 9.2, 8.8]
    avg = inf.ensemble_average(test_props, weights)
    print(f"Ensemble averaged property: {avg}")
