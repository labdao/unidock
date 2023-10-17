"""Main file for the project."""
import subprocess
import os

valid_params = [
    "receptor", "flex", "ligand", "ligand_index", "batch", "gpu_batch", "scoring",
    "maps", "center_x", "center_y", "center_z", "size_x", "size_y", "size_z", "autobox",
    "out", "dir", "write_maps", "cpu", "seed", "exhaustiveness", "max_evals", "num_modes",
    "min_rmsd", "energy_range", "spacing", "verbosity", "max_step", "refine_step", "max_gpu_memory",
    "search_mode", "config", "help", "version"
]

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
        command.append(f"--{param}")
        command.append(str(value))

    # Run the command
    subprocess.check_call(command)


def main():
    """Main function for the project."""
    # Run Uni-Dock
    run_unidock(receptor="data/inputs/mmp13.pdbqt", gpu_batch="data/inputs/ligands", dir="data/outputs",
                center_x=-36.0095, center_y=25.628500000000003, center_z=67.49199999999999,
                size_x=17.201, size_y=14.375000000000004, size_z=12.239999999999995
                )


if __name__ == "__main__":
    main()