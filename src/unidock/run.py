"""Module for implementing Uni-Dock pipeline."""
import os
import subprocess
import re
import pandas as pd
import hydra
from omegaconf import DictConfig
from ._mol_convert import (
    retrieve_smiles,
    pdb_to_pdbqt,
    smi_to_pdbqt,
    smiles_to_single_smi,
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

        # Create output poses directory
        os.makedirs(self.config["dir"])

        # Construct the command
        command = ["Unidock"]
        for param, value in self.config.items():
            # Add the list of filepaths if the parameter is "gpu_batch"
            if param == "gpu_batch":
                command.append("--gpu_batch")
                command.extend(ligand_files)
            else:
                command.append(f"--{param}")
                command.append(str(value))

        # Run the command
        subprocess.run(command, check=True)

    def _check_param_validity(self):
        # Check that all parameters are valid
        for param in self.config:
            if param not in VALID_PARAMS:
                raise ValueError(f"Invalid parameter: {param}")

    def save_results(self, output_path: str, best=False) -> None:
        """Processes the pdbqt outputs into csv files."""
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
            df = self._get_best_poses(df)
            df.to_csv(output_path.replace(".csv", "_best.csv"), index=False)

    def _get_floats(self, string):
        """Extracts floating point values from a string containing whitespace and text."""
        # Match negative numbers and floating point numbers
        pattern = re.compile(r"-?\d+\.\d+")
        # Find all numbers in the string
        matches = pattern.findall(string)
        # Get numbers and convert them to float
        return [float(num) for num in matches]

    def _get_best_poses(self, all_poses_df: pd.DataFrame):
        """Identifies best configuration for each ligand."""
        # Select all MODEL rows with a value of 1
        best_poses_df = all_poses_df[all_poses_df["MODEL"] == 1]
        return best_poses_df


@hydra.main(version_base=None, config_path="../../data", config_name="config")
def main(cfg: DictConfig) -> None:
    # Convert receptor pdb file to pdbqt
    if cfg.receptor.type == "pdb":
        context(
            pdb_to_pdbqt,
            cfg.ligands.path,
            os.path.join(cfg.output.path, "processed/receptor.pdbqt"),
        )

        cfg.receptor.path = os.path.join(cfg.output.path, "processed/receptor.pdbqt")

    # Convert ligand SMILES to pdbqt files
    if cfg.ligands.type == "smiles":
        # Retrieve small molecule SMILES from database
        smiles = retrieve_smiles(cfg.ligands.path)

        # Convert small molecule SMILES to sdf files
        smiles_to_single_smi(smiles, os.path.join(cfg.output.path, "processed/sdf_ligands"))

        context(
            smi_to_pdbqt,
            os.path.join(cfg.output.path, "processed/smi_ligands"),
            os.path.join(cfg.output.path, "processed/pdbqt_ligands")
        )

        cfg.ligands.path = os.path.join(cfg.output.path, "processed/sdf_ligands")


    # Create Uni-Dock object
    unidock = UniDock(
        {
            "receptor": cfg.receptor.path,
            "gpu_batch": cfg.ligands.path,
            "dir": os.path.join(cfg.output.path, "poses"),
            **cfg.search_space,
            "search_mode": cfg.optional.search_mode,
        }
    )

    # Run UniDock
    unidock.run()

    # Save outputs
    unidock.save_results(os.path.join(cfg.output.path, "results.csv"), best=True)


if __name__ == "__main__":
    main()
