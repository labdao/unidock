"""Convert molecule types"""
import os
import subprocess
import tempfile
from typing import List, Callable
from pathlib import Path
from functools import wraps
from concurrent.futures import ThreadPoolExecutor
import duckdb
from rdkit import Chem
from rdkit.Chem import AllChem


VALID_FILE_TYPES = ["smi", "pdb", "parquet"]


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
    """Convert SMILES strings to smi file"""
    # Create an empty output directory if not present
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Save each SMILES string to individual smi file
    for i, smiles_string in enumerate(smiles_strings):
        output_file = os.path.join(output_path, f"ligand_{i}.smi")
        with open(output_file, "w", encoding="utf-8") as file:
            file.write(f"{smiles_string}\n")

def smiles_to_sdf(smiles_strings: List[str], output_path: str) -> None:
    """Convert SMILES strings to sdf files"""
    # Create an empty output directory if not present
    if not os.path.exists(output_path):
        os.makedirs(output_path)

     # Function to process each SMILES string
    def _process_smiles(i, smiles_string, output_path):
        try:
            mol = Chem.MolFromSmiles(smiles_string)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            file_path = os.path.join(output_path, f'ligand_{i+1}.sdf')
            with Chem.SDWriter(file_path) as writer:
                writer.write(mol)
            return f"Processed: {file_path}"
        except Exception as e:
            return f"Error with {smiles_string}: {str(e)}"

    # Using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=5) as executor:  # Adjust max_workers as needed
        futures = [executor.submit(process_smiles, i, smiles, output_path) for i, smiles in enumerate(smiles_strings)]


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
        check=True,
    )


@check_file_type
def pdb_convert(input_path: str, output_path: str) -> None:
    """Converts .pdb to pdbqt"""
    with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
        subprocess.run(["reduce", input_path], stdout=tmp, check=False)
        subprocess.run(["prepare_receptor", "-r", tmp.name, "-o", output_path], check=False)
