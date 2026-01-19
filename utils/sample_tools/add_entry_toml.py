import toml
import sys
from pathlib import Path

def main():
    if len(sys.argv) < 3:
        print("Usage: python add_dataset.py <toml_file> <dataset_name>")
        print("Example: python add_dataset.py datasets.toml JetMET2023A")
        sys.exit(1)
    
    toml_file = sys.argv[1]
    dataset_name = sys.argv[2]

    #load file
    toml_path = Path(toml_file)
    if toml_path.exists():
        with open(toml_path, 'r') as f:
            data = toml.load(f)
    else:
        data = {}
    
    if dataset_name in data:
        print(f"Error: Dataset '{dataset_name}' already exists!")
        sys.exit(1)
    
    # Create empty entry
    data[dataset_name] = {
        "files": [],
        "DAS_dataset_name": "",
        "branch": "Events",
        "is_MC": False,
        "is_UL": False,
        "is_SMS": False
    }
    
    # Write back
    with open(toml_path, 'w') as f:
        toml.dump(data, f)
    
    print(f"Added empty entry '{dataset_name}' to {toml_file}")

if __name__ == "__main__":
    main()