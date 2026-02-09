import cclib
import os

class ResultParser:
    @staticmethod
    def parse_orca(output_path):
        """Parses ORCA output for energy and properties."""
        results = {}
        # Try manual extraction first to avoid cclib's potential decoding issues with ORCA 6
        try:
            with open(output_path, "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()
                if "FINAL SINGLE POINT ENERGY" in content:
                    import re
                    match = re.search(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", content)
                    if match:
                        energy_hartree = float(match.group(1))
                        results["energy_ev"] = energy_hartree * 27.2114
                        print(f"Manual extraction successful: {results['energy_ev']} eV")
                        # Try to get HOMO/LUMO if available via cclib after
        except Exception as e:
            print(f"Manual extraction error: {e}")

        # Try cclib for additional properties if needed, but wrap it
        try:
            data = cclib.io.ccread(output_path)
            if "energy_ev" not in results:
                results["energy_ev"] = data.scfenergies[-1]
            results["homo"] = data.moenergies[0][data.homos[0]]
            results["lumo"] = data.moenergies[0][data.homos[0]+1] if len(data.moenergies[0]) > data.homos[0]+1 else None
            results["dipole"] = data.moments[1] if hasattr(data, "moments") else None
            return results
        except Exception as e:
            print(f"cclib parsing failed: {e}")
            return results if "energy_ev" in results else None

    @staticmethod
    def parse_xtb(work_dir):
        """Parses xTB output files (xtbopt.log or similar)."""
        # xTB energies are usually in Hartree in xtbopt.log
        energy_file = os.path.join(work_dir, "xtbopt.log")
        if os.path.exists(energy_file):
            # Simple extraction for demo purposes
            with open(energy_file, "r") as f:
                content = f.read()
                if "energy:" in content:
                    # Very crude extraction
                    parts = content.split("energy:")
                    energy_hartree = float(parts[-1].split()[0])
                    return energy_hartree * 27.2114 # Convert to eV
        return None
