
# preamble
import sys
import numpy as np
import MDAnalysis

sys.path.insert(0,'/home/russ/Scripts/git/USalign/pybindings')

import pySOIalign

# function def
def tri2single(resname):
    '''
    CONVERT FROM 3 TO 1 LETTER AA CODES
    Read in an amino acid's three letter code and return the aa's single letter code.
    Input:
        resname: a string of length 3, any case
    Output:
        the single letter code of associated aa; if resname is unexpected, returns 'X' string.
    '''
    return {'ALA': 'A',
            'ARG': 'R',
            'ASN': 'N',
            'ASP': 'D',
            'CYS': 'C',
            'GLN': 'Q',
            'GLU': 'E',
            'GLY': 'G',
            'HIS': 'H',
            'HIE': 'H',
            'HID': 'H',
            'ILE': 'I',
            'LEU': 'L',
            'LYS': 'K',
            'MET': 'M',
            'PHE': 'F',
            'PRO': 'P',
            'SER': 'S',
            'THR': 'T',
            'TRP': 'W',
            'TYR': 'Y',
            'VAL': 'V'}.get(resname.upper(), 'X')

# prep universe and input parameters
u = MDAnalysis.Universe('/home/russ/Scripts/git/USalign/tests/data/2pel_chainA.pdb')
u_cas = u.select_atoms('name CA')
u_len = u_cas.n_atoms
u_cas_coords = u_cas.positions.flatten()
u_seq = [tri2single(a.resname) for a in u_cas]
u_sec = pySOIalign.make_sec_py(u_cas_coords, u_len)
u_neighbor_list = pySOIalign.getCloseK_py(u_cas_coords, u_len, 5)
u_sec_bonds = pySOIalign.assign_sec_bond_py(u_sec, u_len)
u_alnStruct = pySOIalign.alnStruct(u_cas_coords, u_neighbor_list, u_sec_bonds, u_seq, u_sec, u_len)

v = MDAnalysis.Universe('/home/russ/Scripts/git/USalign/tests/data/3cna_chainA.pdb')
v_cas = v.select_atoms('name CA')
v_len = v_cas.n_atoms
v_cas_coords = v_cas.positions.flatten()
v_seq = [tri2single(a.resname) for a in v_cas]
v_sec = pySOIalign.make_sec_py(v_cas_coords, v_len)
v_neighbor_list = pySOIalign.getCloseK_py(v_cas_coords, v_len, 5)
v_sec_bonds = pySOIalign.assign_sec_bond_py(v_sec, v_len)
v_alnStruct = pySOIalign.alnStruct(v_cas_coords, v_neighbor_list, v_sec_bonds, v_seq, v_sec, v_len)

uv_alnParameters = pySOIalign.alnParameters(5, -2, 6, 0, False, 0, False, 0, False)

alnResults = pySOIalign.runSOIalign(u_alnStruct, v_alnStruct, uv_alnParameters)

