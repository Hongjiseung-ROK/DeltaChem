import os

def scaffold_project():
    """Constructs the professional Python project structure."""
    dirs = [
        "src",
        "src/utils",
        "data",
        "data/raw",
        "data/processed",
        "tests",
        "notebooks",
        "docs"
    ]
    
    for d in dirs:
        os.makedirs(d, exist_ok=True)
        # Create __init__.py in src and utils
        if "src" in d:
            init_file = os.path.join(d, "__init__.py")
            if not os.path.exists(init_file):
                with open(init_file, "w") as f:
                    pass
    
    print("Project structure scaffolded successfully.")

if __name__ == "__main__":
    scaffold_project()
