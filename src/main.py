"""Main file for the project."""
import subprocess
import os
from pathlib import Path
import tempfile
import pickle


class MoleculeConvertor:
    """A class to prepare a molecule"""

    ACCEPTED_EXTENSIONS = [".pdb", ".sdf", ".mol2"]

    def __init__(self, input_filepath, output_dir):
        self.input_filepath = Path(input_filepath)
        self.output_dir = Path(output_dir)

    def convert_to_pdbqt(self):
        """Converts either a single file or directory of small molecules to pdbqt"""

        # Convert either single file or directory containing multiple files
        if not self.input_filepath.is_dir():
            self._convert_single_file_to_pdbqt()
        else:
            for filename in os.listdir(self.input_filepath):
                self.input_filepath = os.path.join(self.input_filepath, filename)
                self._convert_single_file_to_pdbqt()

    def _convert_single_file_to_pdbqt(self):
        """Converts single file to pdbqt"""
        file_type = self.input_filepath.suffix
        if file_type == ".pdb":
            self._convert_from_pdb()
        elif file_type == ".sdf":
            self._convert_from_sdf()
        elif file_type == ".mol2":
            self._convert_from_mol2()

    def _convert_from_pdb(self):
        """Converts pdb to pdbqt"""
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            subprocess.run(["reduce", self.input_filepath], stdout=tmp, check=False)
            output_filepath = self.output_dir / (self.input_filepath.stem + ".pdbqt")
            subprocess.run(
                ["prepare_receptor", "-r", tmp.name, "-o", output_filepath], check=False
            )

    def _convert_from_sdf(self):
        """Converts sdf to pdbqt"""
        output_filepath = self.output_dir / (self.input_filepath.stem + ".pdbqt")
        subprocess.run(
            ["obabel", self.input_filepath, "-O", output_filepath], check=False
        )

    def _convert_from_mol2(self):
        """Converts mol2 to pdbqt"""
        output_filepath = self.output_dir / (self.input_filepath.stem + ".pdbqt")
        subprocess.run(
            ["obabel", self.input_filepath, "-O", output_filepath], check=False
        )


class UnidockRunner:
    """Class for running the Uni-Dock application"""

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
            if param not in UnidockRunner.VALID_PARAMS:
                raise ValueError(f"Invalid parameter: {param}")

    def parse_outputs(self):
        """Processes the pdbqt outputs into a list of dictionaries."""

        # Get all ligand files from the output directory
        ligand_files = [
            os.path.join(self.config.dir, f) for f in os.listdir(self.config.dir)
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


def main():
    """Main function for the project."""

    # Instantiate a molecule convertor
    mol_con = MoleculeConvertor("/app/data/inputs/SDF_ideal.sdf", "/app/data/inputs")

    # Covert inputs to pdbqt files
    mol_con.convert_to_pdbqt()

    # Instantiate a Uni-Dock runner
    runner = UnidockRunner(
        {
            "receptor": "data/inputs/mmp13.pdbqt",
            "gpu_batch": "data/inputs/ligands",
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

    # Run the Uni-Dock application
    runner.run()

    # Parse the docked ligands
    results = runner.parse_outputs()

    # Save the results to a pickle file
    with open("data/results.pkl", "wb") as file:
        pickle.dump(results, file)


if __name__ == "__main__":
    main()
