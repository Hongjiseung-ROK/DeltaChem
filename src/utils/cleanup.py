import os
import glob

class CleanupManager:
    def __init__(self):
        # ORCA temporary file extensions that can be safely deleted after SP calculation
        self.temp_extensions = [
            ".gbw", ".tmp", ".int", ".property", ".opt", 
            ".coords", ".hess", ".pc", ".vcor", ".vwp", ".wfn"
        ]

    def cleanup_orca_temporary_files(self, directory):
        """Removes large ORCA temporary files from the specified directory."""
        print(f"Cleaning up temporary files in {directory}...")
        files_removed = 0
        for ext in self.temp_extensions:
            pattern = os.path.join(directory, f"*{ext}")
            for f in glob.glob(pattern):
                try:
                    os.remove(f)
                    files_removed += 1
                except Exception as e:
                    print(f"Failed to remove {f}: {e}")
        
        if files_removed > 0:
            print(f"Removed {files_removed} temporary files.")
        else:
            print("No temporary files found.")

if __name__ == "__main__":
    manager = CleanupManager()
    manager.cleanup_orca_temporary_files("data/orca_work")
