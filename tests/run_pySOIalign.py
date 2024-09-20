
#
# this is just a quick script to calculate the input features and alignment
# between two structures in `data/`. Useful to test that the bindings are
# working. won't be used to run tests enmasse.
#
# !!! NOTE: remove before pushing to main !!!
#

import sys
sys.path.insert(0,'/home/russ/Scripts/git/USalign/pybindings') 
import enum
import pySOIalign
import numpy as np
import MDAnalysis

closeK_opt = 0
mm_opt = 6


def readPDB(filename):
    """
    function to parse a PDB file to gather CA atom coords and resnames

    :param filename: string, path to a PDB file that is to be parsed
    :return coords: numpy array of float64 values, shape (nAtoms x 3)
    :return seq: string, sequence string
    """
    with open(filename,'r') as test_data:
        ca_atom_lines = [line for line in test_data.readlines() if line[:4] == 'ATOM' and '  CA  ' in line]
    
    coords = [[line[30:38].strip(),line[38:46].strip(),line[46:54].strip()] for line in ca_atom_lines]
    coords = np.array(coords, dtype = np.float64)
    seq = ''.join([AAtri2single(line[16:20].strip().upper()) for line in ca_atom_lines])

    return coords, seq


class AAResnameMapping(enum.Enum):
    """
    only a subset of amino acid residue names understood by MDAnalysis:
    https://userguide.mdanalysis.org/stable/standard_selections.html#protein-selection

    !!! denotes special/ambiguous cases
    """
    # amino acids
    ALA = 'A'
    ARG = 'R'
    ASN = 'N'
    ASP = 'D'
    ASX = 'B'   # !!!
    CYS = 'C'
    GLN = 'Q'
    GLU = 'E'
    GLX = 'Z'   # !!!
    GLY = 'G'
    HIS = 'H'   # !!!
    HIE = 'H'   # !!!
    HID = 'H'   # !!!
    HIP = 'H'   # !!!
    ILE = 'I'
    LEU = 'L'
    LYS = 'K'
    MET = 'M'
    PHE = 'F'
    PRO = 'P'
    SER = 'S'
    THR = 'T'
    TRP = 'W'
    TYR = 'Y'
    VAL = 'V'
    UNK = 'X'   # !!!
    

class NAResnameMapping(enum.Enum):
    """
    only a subset of nucleic acid residue names understood by MDAnalysis:
    https://userguide.mdanalysis.org/stable/standard_selections.html#nucleic-acids

    !!! denotes special/ambiguous cases
    """
    # nucleic acids
    ADE = 'A'
    CYT = 'C'
    GUA = 'G'
    THY = 'T'
    URA = 'U'
    UNK = 'X'   # !!!


def AAtri2single(resname):
    """
    CONVERT FROM 3 TO 1 LETTER AA CODES
    Read in an amino acid's three letter code and return the aa's single 
    letter code.

    Input:
        :param resname: a string.upper(), associated with a name in the 
                        AAResnameMapping Enum object
    Output:
        :return: the single letter code for the associated aa. 
                 If resname is unexpected, returns 'X' string to denote an 
                 unknown resname. Or an empty string, if the Enum fails for an
                 unexpected error.
    """
    try: 
        return AAResnameMapping[resname].value
    except KeyError:
        print(f"{resname} not in the AAResnameMapping. Returning 'X' to denote an unknown residue type.")
        return 'X'
    except Exception as e:
        print(f"AAtri2single({resname}) failed unexpectedly. Returning '' to remove this {resname} from the analysis.")
        return ''


def AAsingle2tri(resname):
    """
    CONVERT FROM 1 TO 3 LETTER AA CODES
    Read in an amino acid's single letter code and return the aa's tri letter
    code.

    !!! NOTE: There are ambiguous values in the AAResnameMapping Enum object.
              Due to the organization of the Enum, the first entry with the 
              value will be used.

    Input:
        :param resname: a string.upper(), associated with a value in the 
                        AAResnameMapping Enum object
    Output:
        :return: the tri letter code for the associated aa. 
                 If resname is unexpected, returns 'UNK' string to denote an 
                 unknown resname. Or an empty string, if the Enum fails for an
                 unexpected error.
    """
    try:
        return AAResnameMapping(resname).name
    except ValueError:
        print(f"{resname} not in the AAResnameMapping. Returning 'UNK' to denote an unknown residue type.")
        return 'UNK'
    except Exception as e:
        print(f"AAsingle2tri({resname}) failed unexpectedly. Returning '' to remove this {resname} from the analysis.")
        return ''

#
#def tri2single(resname):
#    """
#    no real error handling was done here. just giving 'X' back if resname was 
#    unexpected.
#    """
#    return {'ALA': 'A',
#            'ARG': 'R',
#            'ASN': 'N',
#            'ASP': 'D',
#            'CYS': 'C',
#            'GLN': 'Q',
#            'GLU': 'E',
#            'GLY': 'G',
#            'HIS': 'H',
#            'HIE': 'H',
#            'HID': 'H',
#            'ILE': 'I',
#            'LEU': 'L',
#            'LYS': 'K',
#            'MET': 'M',
#            'PHE': 'F',
#            'PRO': 'P',
#            'SER': 'S',
#            'THR': 'T',
#            'TRP': 'W',
#            'TYR': 'Y',
#            'VAL': 'V'}.get(resname.upper(), 'X')
#

u = MDAnalysis.Universe(
        '/home/russ/Scripts/git/USalign/tests/data/2pel_chainA.pdb')
u_CAs = u.select_atoms('name CA')
u_len = u_CAs.n_atoms
u_CAs_coords = np.array(u_CAs.positions.flatten(),dtype=np.float64)
u_seq = ''.join([AAtri2single(a.resname.upper()) for a in u_CAs])
#u_CAs_coords, u_seq = readPDB('/home/russ/Scripts/git/USalign/tests/data/2pel_chainA.pdb')
#u_len = u_CAs_coords.shape[0]
u_sec = pySOIalign.wrap_make_sec(u_CAs_coords, u_len)
u_sec_bonds = pySOIalign.wrap_assign_sec_bond(u_sec, u_len)
u_neighbor_list = pySOIalign.wrap_getCloseK(u_CAs_coords, u_len, closeK_opt)
u_alnStruct = pySOIalign.alnStruct(u_CAs_coords, 
                                   u_neighbor_list, 
                                   u_sec_bonds, 
                                   u_seq, 
                                   u_sec, 
                                   u_len)

print(u_len, u_seq, u_len == len(u_seq), u_CAs_coords.shape)

v = MDAnalysis.Universe(
        '/home/russ/Scripts/git/USalign/tests/data/3cna_chainA.pdb')
v_CAs = v.select_atoms('name CA')
v_CAs_coords = np.array(v_CAs.positions.flatten(),dtype=np.float64)
v_len = v_CAs.n_atoms
v_seq = ''.join([AAtri2single(a.resname.upper()) for a in v_CAs])
#v_CAs_coords, v_seq = readPDB('/home/russ/Scripts/git/USalign/tests/data/3cna_chainA.pdb')
#v_len = v_CAs_coords.shape[0]
v_sec = pySOIalign.wrap_make_sec(v_CAs_coords, v_len)
v_sec_bonds = pySOIalign.wrap_assign_sec_bond(v_sec, v_len)
v_neighbor_list = pySOIalign.wrap_getCloseK(v_CAs_coords, v_len, closeK_opt)
v_alnStruct = pySOIalign.alnStruct(v_CAs_coords, 
                                   v_neighbor_list, 
                                   v_sec_bonds, 
                                   v_seq, 
                                   v_sec, 
                                   v_len)

print(v_len, v_seq, v_len == len(v_seq), v_CAs_coords.shape)

uv_alnParameters = pySOIalign.alnParameters(closeK_opt, -1*(u_len+v_len), mm_opt)
#uv_alnParameters = pySOIalign.alnParameters(0, -1*(u_len+v_len), 6, 0, False, 0.0, False, 0.0, False)

results = pySOIalign.runSOIalign(u_alnStruct, v_alnStruct, uv_alnParameters)

print(results.rmsd0)
print(results.L_ali)
print(results.n_ali8)
print(results.TM1)
print(results.TM2)
print(results.translation_vector)
print(results.rotation_matrix.reshape((3,3)))

