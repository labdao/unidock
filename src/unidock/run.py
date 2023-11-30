"""Module for implementing Uni-Dock pipeline."""
import os
from typing import List, Callable
from pathlib import Path
import subprocess
from functools import wraps
import tempfile
import re
import pandas as pd
import duckdb


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
    # Create an empty output directory if not present
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Save each SMILES string to individual smi file
    for i, smiles_string in enumerate(smiles_strings):
        output_file = os.path.join(output_path, f"ligand_{i}.smi")
        with open(output_file, "w", encoding="utf-8") as file:
            file.write(f"{smiles_string}\n")


ConvertFn = Callable[[str, str], None]


def context(strategy: ConvertFn, input_path: str, output_path: str) -> None:
    """Converts chemical formats to pdbqt"""
    # Checks to see if input_path is a directory
    if os.path.isdir(input_path):
        # Creates output directory if doesn't exist
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        # Converts each file if input path is directory
        for input_file in os.listdir(input_path):
            # Creates full input file name
            input_file_with_dir = os.path.join(input_path, input_file)
            # Creates output file name with the same stem as the input file
            output_file = os.path.join(output_path, f"{Path(input_file).stem}.pdbqt")
            strategy(input_file_with_dir, output_file)
    else:
        strategy(input_path, output_path)


def check_file_type(func):
    """Decorator to check if file type is valid"""

    @wraps(func)
    def inner(input_path, output_path):
        # Checks if file type is valid
        file_type = Path(input_path).suffix[1:]
        if file_type not in VALID_FILE_TYPES:
            raise ValueError(f"File type not supported: {file_type}")
        func(input_path, output_path)

    return inner


@check_file_type
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


@check_file_type
def pdb_convert(input_path: str, output_path: str) -> None:
    """Converts .pdb to pdbqt"""
    with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
        subprocess.run(["reduce", input_path], stdout=tmp, check=False)
        subprocess.run(["prepare_receptor", "-r", tmp.name, "-o", output_path], check=False)


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

    def save_results(self, output_path: str, best=False) -> None:
        """Processes the pdbqt outputs into a list of dictionaries."""
        # Get all ligand files from the output directory
        ligand_files = [
            os.path.join(self.config["dir"], f) for f in os.listdir(self.config["dir"])
        ]
        # Variables to extract
        key_var_names = ["AFFINITY", "RMSD lower", "RMSD upper", "INTER + INTRA", "INTER", "INTRA", "UNBOUND", "NAME"]
        # List containing results for each ligand model
        ligand_results = []

        # Go through each ligand
        for f in ligand_files:
            with open(f, "r", encoding="utf-8") as file:
                lines = file.readlines()
                # Identify the indices of the lines starting with "MODEL"
                model_line_indices = [
                    i for i, line in enumerate(lines) if line.startswith("MODEL")
                ]

                # Extract key information from each model
                for model_num, i in enumerate(model_line_indices):
                    res_lines = lines[i+1 : i+6]
                    # Get first three values
                    res_vals = self._get_floats(res_lines[0])
                    # Extracts values from lines and add to results list
                    res_vals.extend([float(line.split(":")[1].strip()) for line in res_lines[1:]])
                    # Add the ligand name
                    res_vals.append(lines[i + 6].split("=")[1].strip())
                    # Create a dictionary of the results
                    res_dict = dict(zip(key_var_names, res_vals))
                    # Add a value to the dictionary corresponding to the model number
                    res_dict['MODEL'] = model_num + 1
                    # Add the dictionary to the list of results
                    ligand_results.append(res_dict)

        # Convert the results to a pandas dataframe
        df = pd.DataFrame(ligand_results)

        # Reorder column names
        new_col_order = [ "NAME", "MODEL", "AFFINITY", "RMSD lower", "RMSD upper", "INTER + INTRA", "INTER", "INTRA", "UNBOUND"]
        df = df[new_col_order]

        # Save all poses to csv
        df.to_csv(output_path, index=False)

        # Save the best pose for each ligand if best=True
        if best:
            df = self._get_best_pose(df)
            df.to_csv(output_path.replace(".csv", "_best.csv"), index=False)

    def _get_floats(self, string):
        """Extracts floating point values from a string containing whitespace and text."""
        # Match negative numbers and floating point numbers
        pattern = re.compile(r"-?\d+\.\d+")
        # Find all numbers in the string
        matches = pattern.findall(string)
        # Get numbers and convert them to float
        return [float(num) for num in matches]

    def _get_best_pose(self, all_poses_df: pd.DataFrame):
        """Identifies best configuration for each ligand."""
        # Group by ligand
        grouped_ligands = all_poses_df.groupby("NAME")
        # Get the index of the pose with the highest affinity (most negative score) for each group
        best_poses_indices = grouped_ligands["AFFINITY"].idxmin()
        # Select the poses with the highest affinity for each group
        best_poses_df = all_poses_df.loc[best_poses_indices]
        return best_poses_df

def main():
    # Retrieve small molcule SMILES from database
    smiles = retrieve_smiles("./data/smiles_ligands.parquet")

    # Convert small molecule SMILES to .smi files
    smiles_to_smi(smiles, "./data/processed/smi_ligands")

    # Convert small molecule smi file to pdbqt
    context(
        smi_convert,
        "./data/processed/smi_ligands",
        "./data/processed/pdbqt_ligands",
    )

    # # Convert target pdb file to pdbqt
    # context(
    #     pdb_convert,
    #     "./data/inputs/1adc.pdb",
    #     "./data/processed/target.pdbqt",
    # )

    # Create Uni-Dock object
    unidock = UniDock(
        {
            "receptor": "./data/mmp13_receptor.pdbqt",
            "gpu_batch": "./data/processed/pdbqt_ligands",
            "dir": "./data/outputs",
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

    # Save outputs
    unidock.save_results("./data/results.csv ", best=True)


if __name__ == "__main__":
    main()
