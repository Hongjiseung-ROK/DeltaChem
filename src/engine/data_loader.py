class QM9SubsetLoader:
    @staticmethod
    def get_smiles_subset():
        """Returns a representative subset of QM9 (small molecules)."""
        qm9_smiles = [
            ("CH4", "C"),           # Methane
            ("H2O", "O"),           # Water
            ("NH3", "N"),           # Ammonia
            ("C2H6", "CC"),         # Ethane
            ("C2H4", "C=C"),       # Ethene
            ("CH3OH", "CO"),       # Methanol
            ("CH3CN", "CC#N"),     # Acetonitrile
            ("Propane", "CCC"),     # Propane
            ("Acetone", "CC(=O)C"), # Acetone
            ("Benzene", "C1=CC=CC=C1"), # Benzene
            ("Ethanol", "CCO"),      # Ethanol
            ("Cyclopropane", "C1CC1"), # Cyclopropane
            ("Ethyne", "C#C"),       # Ethyne
            ("Formaldehyde", "C=O"), # Formaldehyde
            ("Methylamine", "CN"),    # Methylamine
            ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"), # Aspirin
            ("Nicotine", "CN1CCCC1C2=CN=CC=C2"),      # Nicotine
            ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O") # Ibuprofen
        ]
        return qm9_smiles

if __name__ == "__main__":
    loader = QM9SubsetLoader()
    subset = loader.get_smiles_subset()
    for name, smiles in subset:
        print(f"{name}: {smiles}")
