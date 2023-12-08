from unidock._mol_convert import (
    retrieve_smiles,
    smiles_to_single_smi,
    smiles_to_multiple_smis,
    smiles_to_sdf,
    context,
    smi_to_pdbqt,
)


def test_smiles_to_pdbqt_using_single_smi():
    """Tests running time for converting smiles to pdbqt looping through each file using obabel"""
    smiles = retrieve_smiles("./data/smiles_ligands.parquet")
    smiles_to_single_smi(smiles, "./data/output/smiles_ligand")
    context(smi_to_pdbqt, "./data/output/smiles_ligand", "./data/output/pdbqt_ligands")


def test_smiles_to_pdbqt_using_multiple_smis():
    """Tests running time for converting smiles to pdbqt using a single smi
    file with multiple entries"""
    smiles = retrieve_smiles("./data/smiles_ligands.parquet")
    smiles_to_multiple_smis(smiles, "./data/output/smiles_ligands")
    context(smi_to_pdbqt, "./data/output/smiles_ligands", "./data/output/pdbqt_ligands")


def test_smiles_to_sdf():
    """Tests running time for converting smiles to sdf"""
    smiles = retrieve_smiles("./data/smiles_ligands.parquet")
    smiles_to_sdf(smiles, "./data/output/sdf_ligands")
