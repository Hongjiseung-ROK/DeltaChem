import os
from dotenv import load_dotenv

load_dotenv()

class Config:
    # Computational Paths
    ORCA_PATH = os.getenv("ORCA_PATH", r"c:\workspace\222_cc_project\orca_bin\orca.exe")
    XTB_PATH = os.getenv("XTB_PATH", "xtb")
    
    # Directory Structure
    PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    DATA_DIR = os.getenv("DATA_DIR", os.path.join(PROJECT_ROOT, "data"))
    OUTPUT_DIR = os.path.join(DATA_DIR, "production")
    
    # Physics Defaults
    DEFAULT_SOLVENT = os.getenv("DEFAULT_SOLVENT", "water") # GBSA/CPCM
    TEMPERATURE = float(os.getenv("TEMPERATURE", 298.15))
    
    # Model Params
    HIDDEN_DIM = int(os.getenv("HIDDEN_DIM", 64))
    BATCH_SIZE = int(os.getenv("BATCH_SIZE", 32))

    @classmethod
    def get_work_dir(cls, sub_path):
        path = os.path.join(cls.DATA_DIR, sub_path)
        os.makedirs(path, exist_ok=True)
        return path

if __name__ == "__main__":
    print(f"ORCA Path: {Config.ORCA_PATH}")
    print(f"Data Dir: {Config.DATA_DIR}")
