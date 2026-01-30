#!/usr/bin/python3
"""
Created on Wed Jul 31 11:55:19 2024

@author: kgr26424
"""
import numpy as np

def calculate_nb_of_probes(
    molecule: str = 'pathtomol_or_smiles',
    format: str = '',
    ratio: float = 0.05,
    box: list = [10,10,10],
    unit: str = 'Ang',
    membrane: bool = False,
    verbosity: bool = False
) -> int:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    if format.lower() == 'pdb':
        mol = Chem.MolFromPDBFile(molecule)
    elif format.lower() == 'mol':
        mol = Chem.MolFromMolFile(molecule)
    elif format.lower() == 'smiles':
        mol = Chem.MolFromSmiles(molecule)
    else:
        raise ValueError("ERROR! Only accept PDB or MOL files, or SMILES string!")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    vol_mol = AllChem.ComputeMolVolume(mol) # Always in Angstrom^3
    if unit == 'Ang':
        converter = 1
    elif unit == 'nm':
        converter = 10
    else:
        raise ValueError("ERROR! Metric unit not valid! Should be 'Ang' or 'nm'")
    assert len(box) == 3
    if membrane:
        vol_box = (box[0]*converter)*(box[1]*converter)*((box[2]*converter)-40)
    else:
        vol_box = (box[0]*converter)*(box[1]*converter)*(box[2]*converter)
    
    n = (ratio*vol_box)/vol_mol
    
    if verbosity:
        print(f'You need {int(n)} probe molecules to reach {ratio*100}% v/v in a {box[0]} {unit} x {box[1]} {unit} x {box[2]} {unit} box')
    
    return int(n)

def calculate_probe_conc(
    n_mol: int = 1,
    box: list = [10,10,10],
    unit: str ='Ang',
    membrane: bool =False,
    verbosity: bool =True
) -> float:
    avogadro=6.02214076e23         # Avogadro constant in mol-1
    if unit == 'nm':
        vol_conv = 1e-24
        memb_height = 4
    elif unit == 'Ang':
        vol_conv = 1e-27
        memb_height = 40
    if membrane:
        conc = (n_mol/(avogadro * (box[0]*box[1]*(box[2]-memb_height))*vol_conv))
    else:
        conc = (n_mol/(avogadro * (box[0]*box[1]*box[2])*vol_conv))
        
    if verbosity:
        print(f'{n_mol} molecules in a {box[0]} {unit} x {box[1]} {unit} x {box[2]} {unit} system correspond to a concentration of: {conc} M.')
    
    return conc