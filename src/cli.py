import argparse
import sys
from src.utils.config import Config
from src.engine.sampling import SamplingEngine
from src.engine.inference import InferenceEngine

def main():
    parser = argparse.ArgumentParser(description="DeltaChem: Enterprise Delta-Learning CLI")
    parser.add_argument("--smiles", type=str, help="SMILES string to process")
    parser.add_argument("--batch", type=str, help="CSV file with SMILES list")
    parser.add_argument("--solvent", type=str, default=Config.DEFAULT_SOLVENT, help="Implicit solvent (e.g., water, methanol)")
    parser.add_argument("--accuracy", action="store_true", help="Run accuracy study instead of production")
    
    args = parser.parse_args()
    
    if not args.smiles and not args.batch and not args.accuracy:
        parser.print_help()
        sys.exit(1)
        
    print(f"DeltaChem CLI Started (Solvent: {args.solvent})")
    
    if args.accuracy:
        from src.analysis.accuracy_study import AccuracyStudy
        study = AccuracyStudy()
        test_set = [("Methane", "C"), ("Water", "O"), ("Methanol", "CO")]
        study.run_study(test_set)
        study.report()
    elif args.smiles:
        print(f"Processing SINGLE molecule: {args.smiles}")
        # Add production logic here
    elif args.batch:
        print(f"Processing BATCH CSV: {args.batch}")

if __name__ == "__main__":
    main()
