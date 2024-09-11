
#
# this is just a quick script to calculate the input features and alignment
# between two structures in `data/`. Useful to test that the bindings are
# working. won't be used to run tests enmasse.
#
# !!! NOTE: remove before pushing to main !!!
#

import sys
sys.path.insert(0,'/home/russ/Scripts/git/USalign/pybindings') 

import pySOIalign
import numpy as np
import MDAnalysis

mm_opt = 6
closeK_opt = 0

def tri2single(resname):
    '''
    CONVERT FROM 3 TO 1 LETTER AA CODES
    Read in an amino acid's three letter code and return the aa's single 
    letter code.

    Input:
        resname: a string of length 3, any case
    Output:
        the single letter code of associated aa; if resname is unexpected, 
        returns 'X' string.
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

u = MDAnalysis.Universe(
        '/home/russ/Scripts/git/USalign/tests/data/2pel_chainA.pdb')
u_CAs = u.select_atoms('name CA')
u_len = u_CAs.n_atoms
u_CAs_coords = np.array(u_CAs.positions.flatten(),dtype=np.float64)
u_seq = ''.join([tri2single(a.resname) for a in u_CAs])
u_sec = pySOIalign.make_sec_py(u_CAs_coords, u_len)
u_sec_bonds = pySOIalign.assign_sec_bond_py(u_sec, u_len)
u_neighbor_list = pySOIalign.getCloseK_py(u_CAs_coords, u_len, 5)
u_alnStruct = pySOIalign.alnStruct(u_CAs_coords, 
                                   u_neighbor_list, 
                                   u_sec_bonds, 
                                   u_seq, 
                                   u_sec, 
                                   u_len)

v = MDAnalysis.Universe(
        '/home/russ/Scripts/git/USalign/tests/data/3cna_chainA.pdb')
v_CAs = v.select_atoms('name CA')
v_CAs_coords = np.array(v_CAs.positions.flatten(),dtype=np.float64)
v_len = v_CAs.n_atoms
v_seq = ''.join([tri2single(a.resname) for a in v_CAs])
v_sec = pySOIalign.make_sec_py(v_CAs_coords, v_len)
v_sec_bonds = pySOIalign.assign_sec_bond_py(v_sec, v_len)
v_neighbor_list = pySOIalign.getCloseK_py(v_CAs_coords, v_len, 5)
v_alnStruct = pySOIalign.alnStruct(v_CAs_coords, 
                                   v_neighbor_list, 
                                   v_sec_bonds, 
                                   v_seq, 
                                   v_sec, 
                                   v_len)

uv_alnParameters = pySOIalign.alnParameters(0, -2, 6, 0, False, 0.0, False, 0.0, False)

results = pySOIalign.runSOIalign(u_alnStruct, v_alnStruct, uv_alnParameters)

print(results.rmsd0)
print(results.L_ali)
print(results.n_ali8)
print(results.TM1)
print(results.TM2)
print(results.translation_vector)
print(results.rotation_matrix.reshape((3,3)))

