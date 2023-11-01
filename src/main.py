"""Main file for the project."""
import subprocess
import os
from pathlib import Path
import tempfile
import pickle

valid_params = [
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

class Receptor:
    """A class to represent a receptor file"""

    ACCEPTED_EXTENSIONS = [".pdb"]

    def __init__(self, filepath: str):
        self.filepath = Path(filepath)
        self.extension = self.filepath.suffix
        self.receptor_dir = self.filepath.parent
        self.receptor_name = self.filepath.stem
        # Throw an exception if filetype not supported
        if self.extension not in Receptor.ACCEPTED_EXTENSIONS:
            raise ValueError(f"Invalid file type. Acceptable types are {self.ACCEPTED_EXTENSIONS}")
        
    def convert_to_pdbqt(self):
        """Converts file types to pdbqt"""
        if self.extension == ".pdb":
            with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
                subprocess.run(["reduce", self.filepath], stdout=tmp, check=False)
                pdbqt_filepath = str(self.receptor_dir) + str(self.receptor_name) + ".pdbqt"
                subprocess.run(["prepare_receptor", "-r", tmp.name, "-o", pdbqt_filepath], check=False)
                self.filepath = Path(pdbqt_filepath)


def run_unidock(**kwargs):
    """Runs Uni-Dock."""
    
    # Check that all parameters are valid
    for param in kwargs:
        if param not in valid_params:
            raise ValueError(f"Invalid parameter: {param}")

    # Get all ligand files from the input directory
    if "gpu_batch" in kwargs:
        ligand_dir = kwargs["gpu_batch"]
        ligand_files = os.listdir(ligand_dir)
        ligand_files = [os.path.join(ligand_dir, file) for file in ligand_files]

    # Construct the command
    command = ["unidock"]
    for param, value in kwargs.items():
        # Add the list of filepaths if the parameter is "gpu_batch"
        if param == "gpu_batch":
            command.append("--gpu_batch")
            command.extend(ligand_files)
        else:
            command.append(f"--{param}")
            command.append(str(value))

    # Run the command
    subprocess.check_call(command)


def parse_outputs(output_dir):
    """Processes the pdbqt outputs into a list of dictionaries."""
     
    # Get all ligand files from the output directory
    ligand_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir)]
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

    # Prepare filetype
    receptor = Receptor("./data/inputs/7x1i.pdb")
    test = receptor.filepath
    receptor.convert_to_pdbqt()
    print(test)
    print(receptor.filepath)

    # Run Uni-Dock
    run_unidock(receptor="data/inputs/mmp13.pdbqt",
                gpu_batch="data/inputs/ligands",
                dir="data/outputs",
                center_x=-6.9315,
                center_y=26.579,
                center_z=54.135999999999996,
                size_x=15.341000000000001,
                size_y=10.828,
                size_z=17.556000000000004,
                search_mode="fast"
                )

    # Parse the outputs
    results = parse_outputs("data/outputs")

    # Save the results to a pickle file
    with open("data/results.pkl", "wb") as file:
        pickle.dump(results, file)


if __name__ == "__main__":
    main()
