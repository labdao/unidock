"""Main file for the project."""
import subprocess
import os
import click

class JobInfo:
    """Class to store information about a job."""
    def __init__(self, test: str):
        self.test = test

    def __str__(self):
        return self.test

@click.command()
@click.argument('receptor_path', type=click.Path(exists=True))
@click.argument('ligand_dir', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path(exists=True))
@click.argument('search_mode', type=click.Choice(['fast', 'balance', 'detail']), default='fast')
def run_unidock(receptor_path: str, ligand_dir: str, output_dir: str, search_mode: str) -> JobInfo:
    """Runs Uni-Dock."""
    # Get all ligand files from the input directory
    ligand_files = os.listdir(ligand_dir)
    ligand_files = [os.path.join(ligand_dir, file) for file in ligand_files]

    # Construct the bash command
    command = [
        "unidock",
        "--receptor", receptor_path,
        "--gpu_batch", *ligand_files[1:300],
        "--dir", output_dir,
        "--center_x", "-36.0095",
        "--center_y", "25.628500000000003",
        "--center_z", "67.49199999999999",
        "--size_x", "17.201",
        "--size_y", "14.375000000000004",
        "--size_z", "12.239999999999995",
        "--search_mode", search_mode
        ]

    # Run the command
    subprocess.check_call(command)

if __name__ == "__main__":
    run_unidock()
