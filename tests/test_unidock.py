from pathlib import Path
import pytest
from unidock.run import retrieve_smiles, smiles_to_smi, smi_convert, context, UniDock

@pytest.fixture(scope="session")
def pdbqt_ligands(tmpdir_factory):
    # Create a temporary directory with the factory
    temp_dir = tmpdir_factory.mktemp("pdbqt")

    # Get the SMILES data
    smiles = retrieve_smiles("data/smiles_ligands.parquet")

    # Create directories for pdbqt and smi within the temporary directory
    pdbqt_dir = temp_dir.mkdir("pdbqt")
    smi_dir = temp_dir.mkdir("smi")

    # Convert SMILES to SMI format and store them in the smi directory
    smiles_to_smi(smiles, smi_dir)

    context(
        smi_convert,
        smi_dir,
        pdbqt_dir,
    )

    return pdbqt_dir


# Check whether a list is returned
def test_retrieve_smiles_returns_list():
    smiles = retrieve_smiles("data/smiles_ligands.parquet")
    assert isinstance(smiles, list)

# Check whether the list contains the correct number of SMILES strings
def test_convert_smiles_to_smi(tmp_path):
    # Arrange
    smiles = retrieve_smiles("data/smiles_ligands.parquet")

    # Act
    smiles_to_smi(smiles, tmp_path)

    # Assert
    assert len(list(tmp_path.iterdir())) == len(smiles)

# Check if number of pdbqt files is equal to number of smi files
@pytest.mark.slow
def test_smiles_to_pdbqt(tmp_path):
     # Arrange
    smiles = retrieve_smiles("data/smiles_ligands.parquet")
    pdbqt = tmp_path / "pdbqt"
    pdbqt.mkdir()
    smi = tmp_path / "smi"
    smi.mkdir()
    smiles_to_smi(smiles, smi)

    # Act
    context(
        smi_convert,
        smi,
        pdbqt,
    )

    # Assert
    assert len(list(pdbqt.iterdir())) == len(list(smi.iterdir()))

# Check to see if pdbqt files are created
def test_unidock_run(tmp_path):
    unidock = UniDock(
        {
            "receptor": "./data/mmp13_receptor.pdbqt",
            "gpu_batch": "./data/pdbqt_ligands",
            "dir":  tmp_path,
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

    assert len(list(tmp_path.iterdir())) == len(list(Path("./data/pdbqt_ligands").iterdir()))



# Check if unidock runs using the generated pdbqt files
@pytest.mark.slow
def test_run_unidock_with_pdbqt_from_smiles(pdbqt_ligands, tmp_path):
    
    unidock = UniDock(
        {
            "receptor": "data/mmp13_receptor.pdbqt",
            "gpu_batch": pdbqt_ligands,
            "dir":  tmp_path,
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

    # Check if the number of files in the output directory is
    # within 5% of the number of files in the input directory
    num_files_tmp_path = float(len(list(tmp_path.iterdir())))
    num_files_pdbqt_ligands = float(len(list(Path("./data/pdbqt_ligands").iterdir())))
    lower_bound = num_files_pdbqt_ligands * 0.95
    upper_bound = num_files_pdbqt_ligands * 1.05

    assert lower_bound <= num_files_tmp_path <= upper_bound
