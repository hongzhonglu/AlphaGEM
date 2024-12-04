from Bio import PDB
import numpy as np

def get_plddt_from_pdb(pdb_path):
    p = PDB.PDBParser()
    structure = p.get_structure('input', pdb_path)
    plddt_list = []
    for a in structure.get_residues():
        b = a.get_unpacked_list()
        if len(b) > 0:
            plddt_list.append(b[0].get_bfactor())
    plddt = np.average(np.array(plddt_list))
    plddt = round(plddt, 4)
    return plddt
