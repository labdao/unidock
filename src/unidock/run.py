"""Module for implementing Uni-Dock pipeline."""
import os
import subprocess
import re
import pandas as pd
from hydra_zen import builds, zen, launch
from ._mol_convert import (
    retrieve_smiles,
    smi_convert,
    pdb_convert,
    smiles_to_smi,
    context,
)


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
        key_var_names = [
            "AFFINITY",
            "RMSD lower",
            "RMSD upper",
            "INTER + INTRA",
            "INTER",
            "INTRA",
            "UNBOUND",
            "NAME",
        ]
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
                    res_lines = lines[i + 1 : i + 6]
                    # Get first three values
                    res_vals = self._get_floats(res_lines[0])
                    # Extracts values from lines and add to results list
                    res_vals.extend(
                        [float(line.split(":")[1].strip()) for line in res_lines[1:]]
                    )
                    # Add the ligand name
                    res_vals.append(lines[i + 6].split("=")[1].strip())
                    # Create a dictionary of the results
                    res_dict = dict(zip(key_var_names, res_vals))
                    # Add a value to the dictionary corresponding to the model number
                    res_dict["MODEL"] = model_num + 1
                    # Add the dictionary to the list of results
                    ligand_results.append(res_dict)

        # Convert the results to a pandas dataframe
        df = pd.DataFrame(ligand_results)

        # Reorder column names
        new_col_order = [
            "NAME",
            "MODEL",
            "AFFINITY",
            "RMSD lower",
            "RMSD upper",
            "INTER + INTRA",
            "INTER",
            "INTRA",
            "UNBOUND",
        ]
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


# Params for running:
# - smiles_file(smiles): str
# - receptor_file(pdb): str
# - bounding box details : dict
# - output_dir: str


def run_unidock_with_smiles(
    smiles_file: str, receptor_file: str, output_dir: str, bounding_box: dict
):
    # Retrieve small molcule SMILES from database
    smiles = retrieve_smiles(smiles_file)

    # Convert small molecule SMILES to .smi files
    smiles_to_smi(smiles, os.path.join(output_dir, "processed/smi_ligands"))

    # Convert small molecule smi file to pdbqt
    context(
        smi_convert,
        os.path.join(output_dir, "processed/smi_ligands"),
        os.path.join(output_dir, "processed/pdbqt_ligands"),
    )

    # Convert target pdb file to pdbqt
    context(
        pdb_convert,
        receptor_file,
        os.path.join(output_dir, "processed/receptor.pdbqt"),
    )

    # Create Uni-Dock object
    unidock = UniDock(
        {
            "receptor": os.path.join(output_dir, "processed/receptor.pdbqt"),
            "gpu_batch": os.path.join(output_dir, "processed/pdbqt_ligands"),
            "dir": output_dir,
            **bounding_box,
            "search_mode": "fast",
        }
    )

    # Run UniDock
    unidock.run()

    # Save outputs
    unidock.save_results(os.path.join(output_dir, "results.csv"), best=True)


def main():
    Config = builds(run_unidock_with_smiles, populate_full_signature=True)
    wrapped_fn = zen(run_unidock_with_smiles)
    job = launch(
        Config,
        wrapped_fn,
        overrides=["smiles_file=", "receptor_file=", "output_dir=", "bounding_box="],
        version_base="1.1",
    )


if __name__ == "__main__":
    main()
