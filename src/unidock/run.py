import os
from typing import List, Callable
from pathlib import Path
import subprocess
import duckdb
import tempfile

VALID_FILE_TYPES = ["smi", "pdb"]

VALID_PARAMS = [
    "receptor",
    "flex",
    "ligand",
    "ligand_index",
    "batch",
    "gpu_batch",
    "scoring",
    "maps",
    "center_x",
    "center_y",
    "center_z",
    "size_x",
    "size_y",
    "size_z",
    "autobox",
    "out",
    "dir",
    "write_maps",
    "cpu",
    "seed",
    "exhaustiveness",
    "max_evals",
    "num_modes",
    "min_rmsd",
    "energy_range",
    "spacing",
    "verbosity",
    "max_step",
    "refine_step",
    "max_gpu_memory",
    "search_mode",
    "config",
    "help",
    "version",
]

def retrieve_smiles(input_path: str) -> List[str]:
    """Return a subset of SMILES from a duckdb database"""

    # Construct the duckdb SQL query
    query = f"SELECT * FROM read_parquet('{input_path}')"

    # Return SMILES as a list of tuples each with a single SMILES string
    results_as_tuples = duckdb.sql(query).fetchall()

    # Converts tuples to single list elements
    results_as_strings = [tup[0] for tup in results_as_tuples]

    return results_as_strings


def smiles_to_smi(smiles_strings: List[str], output_path: str) -> None:
    """Write SMILES strings to smi file"""
     # Save list of SMILES strings to smi file
    with open(output_path, 'w', encoding='utf-8') as file:
        for smiles_string in smiles_strings:
            file.write(f"{smiles_string}\n")


ConvertFn = Callable[[str, str], None]


def context(strategy: ConvertFn, input_path: str, output_path: str) -> None:
    """Converts chemical formats to pdbqt"""
    # Checks if file type is valid
    file_type = Path(input_path).suffix[1:]
    if file_type not in VALID_FILE_TYPES:
        raise ValueError(f"File type not supported: {file_type}")

    # Converts small molecules to pdbqt type
    strategy(input_path, output_path)


def smi_convert(input_path: str, output_path: str) -> None:
    """Convert from .smi format to pdbqt"""
    subprocess.run(
        [
            "obabel",
            "-i",
            "smi",
            input_path,
            "--gen3d",
            "-o",
            "pdbqt",
            "-O",
            output_path,
        ],
        check=False,
    )

class UniDock:
    """Class for running the Uni-Dock application"""

    def __init__(self, config):
        self.config = config

    def run(self):
        """Runs Uni-Dock"""

        # Check that parameters entered by user are valid
        self._check_param_validity()

        # Get all ligand files from the input directory
        if "gpu_batch" in self.config:
            ligand_dir = self.config["gpu_batch"]
            ligand_files = os.listdir(ligand_dir)
            ligand_files = [os.path.join(ligand_dir, file) for file in ligand_files]

        # Construct the command
        command = ["unidock"]
        for param, value in self.config.items():
            # Add the list of filepaths if the parameter is "gpu_batch"
            if param == "gpu_batch":
                command.append("--gpu_batch")
                command.extend(ligand_files)
            else:
                command.append(f"--{param}")
                command.append(str(value))

        # Run the command
        subprocess.run(command, check=False)

    def _check_param_validity(self):
        # Check that all parameters are valid
        for param in self.config:
            if param not in VALID_PARAMS:
                raise ValueError(f"Invalid parameter: {param}")

    def parse_outputs(self):
        """Processes the pdbqt outputs into a list of dictionaries."""

        # Get all ligand files from the output directory
        ligand_files = [
            os.path.join(self.config["dir"], f) for f in os.listdir(self.config["dir"])
        ]
        # Variables to extract
        key_var_names = ["INTER + INTRA", "INTER", "INTRA", "UNBOUND", "NAME"]
        # List of lists of dictionaries containing the results for each ligand model
        ligand_results = []

        # Go through each ligand
        for f in ligand_files:
            with open(f, "r", encoding="utf-8") as file:
                lines = file.readlines()
                # Identify the indices of the lines starting with "MODEL"
                model_line_indices = [
                    i for i, line in enumerate(lines) if line.startswith("MODEL")
                ]
                # List to store the results for each ligand model
                model_results = []

                # Extract key information from each model
                for i in model_line_indices:
                    res = lines[i + 1 : i + 6]
                    # Extracts the values from the key information and whitespace
                    res = [line.split(":")[1].strip() for line in res]
                    # Remove the second and third results for the first item
                    res[0] = res[0].split()[0]
                    # Add the ligand name
                    res.append(lines[i + 6].split("=")[1].strip())
                    # Create a dictionary of the results
                    res_dict = dict(zip(key_var_names, res))
                    # Add the dictionary to the list of results
                    model_results.append(res_dict)

            # Add the results for all models for the given ligand
            ligand_results.append(model_results)

        return ligand_results




def pdb_convert(input_path: str, output_path: str) -> None:
    """Converts .pdb to pdbqt"""
    subprocess.run(
        [
            "obabel",
            "-i",
            "pdb",
            input_path,
            "-o",
            "pdbqt",
            "-O",
            output_path,
        ],
        check=False,
    )

def main():
    # Retrieve small molcule SMILES from database
    smiles = retrieve_smiles("/home/ubuntu/uni-dock/data/inputs/ligands.parquet")

    # Convert small molecule SMILES to .smi file
    smiles_to_smi(smiles, "/home/ubuntu/uni-dock/data/processed/ligands.smi")

    # Convert small molecule smi file to pdbqt
    context(
        smi_convert,
        "/home/ubuntu/uni-dock/data/processed/ligands.smi",
        "/home/ubuntu/uni-dock/data/processed/ligands.pdbqt"
        )

    # Convert target pdb file to pdbqt
    context(
        pdb_convert,
        "/home/ubuntu/uni-dock/data/inputs/target.pdb",
        "/home/ubuntu/uni-dock/data/processed/target.pdbqt"
    )

    # Create Uni-Dock object
    unidock = UniDock(
        {
            "receptor": "data/inputs/mmp13.pdbqt",
            "gpu_batch": "data/inputs/",
            "dir": "data/outputs",
            "center_x": -6.9315,
            "center_y": 26.579,
            "center_z": 54.135999999999996,
            "size_x": 15.341000000000001,
            "size_y": 10.828,
            "size_z": 17.556000000000004,
            "search_mode": "fast",
        }
    )

    # Run UniDock
    unidock.run()

if __name__=="__main__":
    main()
