#!/usr/bin/env python3
""" Task manager for running the dask pipeline for running structure alignment, 
    as defined in the input list file.
    USAGE: 
        python3 dask_taskmgr.py [-h] --tskmgr-log-file TSKMGR.log \
                                     --alignments-pdb-list-file ALN_PATH \
                                     --scheduler-file SCHEDULER_FILE \
                                     --nIterations 1
                                     
    INPUT: 
        -h, --help      show this help message and exit
        --tskmgr-log-file TSKMGR.log, -log TSKMGR.log 
                        path for a file within which logging output for the 
                        workflow will be written
        --alignments-pdb-list-file ALN_PATH, -inp ALN_PATH
                        path to a list file that contains the paths to the pair
                        of structure models to be aligned; third element in each
                        line is the alignment method to be used for the calc,
                        mirroring the options for USalign -mm parameter. 5 and
                        6 for fNS and sNS, respectively.
        --usalign-path /path/to/USalign, -path /path/to/USalign 
                        path to the USalign executable
        --scheduler-file SCHEDULER_FILE, -s SCHEDULER_FILE
                        dask scheduler file; optional, default = None
        --nIterations, 1, -niters 1 
                        integer, number of iterations of each alignment to 
                        perform, default = 1
"""

import sys
import time
import io
import ast
import re
import json
import argparse
import platform
import logging
import itertools
import numpy as np
import csv
from uuid import uuid4
from pathlib import Path

import subprocess
from subprocess import CalledProcessError

import dask
import dask.config
from distributed import Client, Worker, as_completed, get_worker

#######################################
### LOGGING FUNCTIONS
#######################################

def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""
    formatter = logging.Formatter('%(asctime)s    %(levelname)s       %(message)s')
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def clean_logger(logger):
    """To cleanup the logger instances once we are done with them"""
    for handle in logger.handlers:
        handle.flush()
        handle.close()
        logger.removeHandler(handle)


#######################################
### PARSING FUNCTION
#######################################

def tri2single(resname):
    '''
    CONVERT FROM 3 TO 1 LETTER AA CODES
    Read in an amino acid's three letter code and return the aa's single 
    letter code.
    
    Input:
        :param resname: a string of length 3, any case
    
    Output:
        :return: the single letter code of associated aa; if resname is 
                 unexpected, returns 'X' string.
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


def parse_usalign_file(file_object):
    """ Parse a USalign alignment file

    #BUG: regex to get struct1 and struct2 expect local or global paths that 
          include some slashes... so this propagates to the alignment workflow 
          where the list files need to include the paths... grumble

    INPUT:
    :param file_object: file object that can be read. many methods of creating 
                        such an object; ex: io.StringIO or "with open()" 
    RETURNS:
    :return results: dictionary of quantitative results associated with the
                     alignment. 
            Keys:
                'struct1': string, path to the first structure in the alignment;
                           the mobile or query structure
                'struct2': string, path to the second structure in the 
                           alignment; the target structure
                'struct1_chainID': string, chain ID string associated with the 
                                   query structure
                'struct2_chainID': string, chain ID string associated with the 
                                   target structure
                'TMscore1': float, TMscore value normalized by the query 
                            structure's length
                'TMscore2': float, TMscore value normalized by the target 
                            structure's length
                'RMSD': float, root mean square deviation between the aligned 
                        atoms
                'SeqIDAli': float, sequence identity of aligned residues
                'Len1': float, the query structure's length
                'Len2': float, the target structure's length
                'LenAligned': float, the number of residues used in the 
                'd0_1': float, the normalized distance for query structure; used 
                        in TMscore metric calculation
                'd0_2': float, the normalized distance for target structure; used 
                        in TMscore metric calculation
                'map_1_to_2': optional, dictionary, keys are query structure's 
                              ('1') residue indices (string) that map to a 
                              tuple with query residue name, target structure's
                              ('2') residue index (string), target residue name,
                              and alignment distance (float, only available for 
                              SNS and FNS alignments); if alignment_type is not 
                              in ['CP', 'SNS', 'FNS'], this mapping will not be 
                              created
                'dt': float, reported time for alignment; units: seconds
                'trans_vector': list, shape = (3), cart coords to translate 
                                mobile onto target.
                'rot_vector': list, shape = (3,3), rot matrix to align 
                              mobile to target; always apply after translation
    """
    start_time = time.time()
    results = {}
    # file_object is a TextIOBase text stream created by open() or io.StringIO 
    # or some other way; here we just collect all lines
    lines = file_object.readlines()

    ### gather relevant lines common across all USalign outputs
    # query structure lines
    struct1_lines = [line.strip() for line in lines if 'Structure_1:' in line]
    # finds file name; finds the first instance of a match; takes the stem of
    # that path. 
    results['struct1'] = Path(re.search(r'[\w-]+\.', struct1_lines[0])[0]).stem
    #results['struct1'] = re.search(r'((?<!\w)(\.{1,2})?(?<!\/)(\/((\\\b)|[^ \b%\|:\n\"\\\/])+)+\/?)', struct1_lines[0])[0]    # finds absolute or relative paths; finds the first instance of a match. 
    # finds the chainID that was used in the alignment
    results['struct1_chainID'] = re.search(r'(?<=[:])\w+',struct1_lines[0])[0]
    # gathering quantitative metrics associated with TM-score, Length of query 
    # structure, and d0 value for the alignment
    results['TMscore1'], results['Len1'], results['d0_1'] = [float(elem) for elem in re.findall(r'(?<=[=\s](?=\d))\d*[\.]?\d*',struct1_lines[2])]
    
    # target structure lines
    struct2_lines = [line.strip() for line in lines if 'Structure_2:' in line]
    # finds file name; finds the first instance of a match; takes the stem of
    # that path. 
    results['struct2'] = Path(re.search(r'[\w-]+\.', struct2_lines[0])[0]).stem 
    #results['struct2'] = re.search(r'((?<!\w)(\.{1,2})?(?<!\/)(\/((\\\b)|[^ \b%\|:\n\"\\\/])+)+\/?)', struct2_lines[0])[0]    # finds absolute or relative paths; finds the first instance of a match. 
    # find the chainID that was used in the alignment
    results['struct2_chainID'] = re.search(r'(?<=[:])\w+',struct2_lines[0])[0] 
    # gathering quantitative metrics associated with TM-score, Length of target
    # structure, and d0 value for the alignment
    results['TMscore2'], results['Len2'], results['d0_2'] = [float(elem) for elem in re.findall(r'(?<=[=\s](?=\d))\d*[\.]?\d*',struct2_lines[2])] 
    
    # Alignment overview line
    aln_lines = [line.strip() for line in lines if 'Aligned length=' in line]
    # gathering quantitative metrics associated with Length of Alignment, RMSD,
    # and sequence ID for the aligned residues
    results['LenAligned'], results['RMSD'], results['SeqIDAli'] = [float(elem) for elem in re.findall(r'(?<=[=\s](?=\d))\d*[\.]?\d*',aln_lines[0])]

    # dictionary of tuples; keys are mobile residue index with values being a tuple of (target residue index, mobile resname, target resname, distance)
    results['map_1_to_2'] = {}
    
    # gather the alignment mapping, strip white space on both sides
    map_lines    = [line.strip() for line in lines if ' CA ' == line[:4]]
    for line in map_lines:
        # line format: 'CA  LEU A 111 \t CA  GLN A  97 \t    1.568'
        # awkward white spaces:      ^^               ^^
        temp = [elem.strip() for elem in line.split('\t')]
        
        if len(temp) == 2:
            # CP-align does not print distance btw atoms
            dist = 0
        elif len(temp) == 3:
            # SOIalign methods print distance btw atoms
            dist = float(temp[2])
        else:
            # something funky going on, skip the map_lines
            break
        
        # from the above example: {'111': ('97','LEU','GLN',1.568)}, where 
        # dict key is mobile resid (type int) that maps to a value (type: 
        # tup). The tup is filled with target resid, mobile resname, target
        # resname, and distance btw the two atoms.
        results['map_1_to_2'].update(
            {int(temp[0][-4:].strip()): (int(temp[1][-4:].strip()),
                                         temp[0][4:7],
                                         temp[1][4:7],
                                         float(temp[2]))})

    # collect USalign reported timing
    #time_lines = [line.strip() for line in lines if 'Total CPU time' in line]
    #results['dt'] = float(re.findall(r'\d*[\.]?\d+',time_lines[0])[0])

    # collect translation and rotation lines if present; these lines are not
    # written in standard USalign output. Gotta `cat` the -m file...
    trans_rot_array = np.array([line.split()[1:] for line in lines if line[0] in ['0','1','2']],dtype=float)
    if trans_rot_array.size > 0:
        results['trans_vector'] = trans_rot_array[:,0].tolist()
        results['rot_matrix'] = trans_rot_array[:,1:].tolist()

    results['dt'] = [time.time() - start_time]
    return results


#######################################
### NATIVE USALIGN CALLS
#######################################

def submit_preprocessing_pipeline(structure):
    """
    preprocessing task function for the dask workflow
    INPUT: 
        :param structure: tuple, a string and an int
            elem[0]: structure_file_path, string, path to the protein structure 
                     file; assumed to be a .pdb file since USalign has limited 
                     parsing for file-input formats.
            elem[1]: closeK_opt, int, value used to set how many nearest 
                     neighbors will be identified for each residue.
    RETURNS:
        :return: results_dict, dictionary of results. Keys:
             'dt': total time of task, units: seconds
             'first_coords': list of xyz coords as floats
             'last_coords':  list of xyz coords as floats
             'seq': string, structure amino acid sequence
             'sec': string, structure 2ndary structure string
             'sec_bounds': list of lists of shape (len x 2)
             'k_neighbors': list of lists of shape (len*closeK_opt x 3)
    """
    start_time = time.time()
    structure_file_path = structure[0]
    closeK_opt = structure[1]
    results_dict = {'task_type': 'preprocess',
                    'structure': structure_file_path,
                    'closeK_opt': closeK_opt}
    # run the preprocessing executable that recreates the steps performed in
    # USalign's SOIalign function to prep a structure for alignment
    try:
        completed_process = subprocess.run(
            f'./preprocessing_tests {structure_file_path} {closeK_opt}',
            shell=True,
            capture_output=True,
            check=True)
    # if the alignment script errors out, skip the parsing steps
    except Exception as e:
        print(f'submit_preprocessing, calc, {structure_file_path}, {closeK_opt}', e, file=sys.stderr, flush=True)
        return results_dict

    # parse the preprocessing_tests output
    try:
        # creates a file-like-object text stream
        stdout_string = completed_process.stdout.decode() 
        stdout_file = io.StringIO(stdout_string)
        lines = stdout_file.readlines()
        # parse the text stream
        results_dict['first_coords'] = ast.literal_eval(lines[0][26:-1]) 
        results_dict['last_coords']  = ast.literal_eval(lines[1][26:-1])
        results_dict['seq'] = lines[2][10:-1]
        results_dict['sec'] = lines[3][25:-1]
        results_dict['sec_bounds']  = ast.literal_eval(lines[5])
        if closeK_opt >= 3:
            results_dict['k_neighbors'] = ast.literal_eval(lines[7])
    
    except Exception as e:
        print(f'submit_preprocessing, proc, {structure_file_path}, {closeK_opt}', e, file=sys.stderr, flush=True)

    # this time encompasses the submit_preprocessing subprocess call 
    # as well as processing of results
    results_dict['dt'] = time.time() - start_time
    
    return results_dict


def submit_subprocess_pipeline(alignments, 
                               nIterations, 
                               USalign_executable_path):
    """
    alignment and parsing function for the dask workflow
    INPUT:
        :param alignments: tuple of two strings and two ints, where the strings 
                  are paths to the two proteins to be aligned nIterations times.
                  elem[0] is the path string pointing at the query protein
                  elem[1] is the path string associated with target protein,
                          being aligned to.
                  elem[2] is the USalign -mm parameter value to be used. 5 for 
                          fNS and 6 for sNS. 
                  elem[3] is the USalign -closeK_opt parameter value to be used. 
                          fNS and 6 for sNS. 
        :param nIterations: int, number of iterations to perform.
        :param USalign_executable_path: str, path to USalign executable.
    RETURNS:
        :return: results_dict, dictionary of results. Keys:
             'avg_dt': average time of task, units: seconds
             'std_dt': standard deviation time of task, units: seconds
             See the parse_usalign_file() function documentation for full list
             of keys from the USalign results.
    """
    results_dict = {'task_type': 'alignment',
                    'mm_opt': int(alignments[2]),
                    'closeK_opt': int(alignments[3]),
                    'mobile': alignments[0],
                    'target': alignments[1]}

    # do many iterations
    for _iter in range(nIterations):
        start_time = time.time()
        # run the alignment script
        try:
            # 1) run USalign between the two structure files, using alignment 
            #    method, outputting alignment results to stdout, and saving the
            #    translocation/rotation matrix to file
            # 2) cat the trans/rot matrix file to stdout
            # 3) delete the trans/rot matrix file so we don't get a ton of files 
            temp_string = str(uuid4())
            completed_process = subprocess.run(
                    f'{USalign_executable_path} {alignments[0]} {alignments[1]} -mm {alignments[2]} -closeK {alignments[3]}  -outfmt 0 -m {temp_string}.dat; cat {temp_string}.dat; rm {temp_string}.dat',
                    shell=True,
                    capture_output=True,
                    check=True)
        # if the alignment script errors out, skip the parsing steps
        except Exception as e:
            print(f'submit_pipeline, aln, {_iter}, {alignments[0]}, {alignments[1]}, {alignments[2]}, {alignments[3]}', e, file=sys.stderr, flush=True)
            continue

        # parse the alignment output
        try:
            # creates a file-like-object text stream
            stdout_string = completed_process.stdout.decode() 
            stdout_file = io.StringIO(stdout_string)
            # put this text stream through the parser
            # hard coding the aln_algo as '' because we are not interested in
            # the residue mapping at the present
            parsed_results = parse_usalign_file(stdout_file)
           
            # this dt time encompasses the USalign subprocess call as well as 
            # post-processing of the alignment results

            # fill the results_dict with the alignment results
            if _iter == 0:
                results_dict.update(parsed_results)
                results_dict['dt'] = [time.time() - start_time]
            
            # otherwise check values against others and add to dt
            else:
                results_dict['dt'].append(time.time() - start_time)
                # these quant metrics should demonstrate that each iteration is
                # recovering the observed values from iter 0; USalign 
                # alignments should be deterministic
                assert results_dict['TMscore1'] == parsed_results['TMscore1']
                assert results_dict['TMscore2'] == parsed_results['TMscore2']
                assert results_dict['RMSD'] == parsed_results['RMSD']
                assert results_dict['SeqIDAli'] == parsed_results['SeqIDAli']
        
        # if parsing or asserts fail, return an error message and move to the 
        # next iteration
        except Exception as e:
            print(f'submit_pipeline, parsing, {_iter}, {alignments[0]}, {alignments[1]}, {alignments[2]}, {alignments[3]}', e, file=sys.stderr, flush=True)
            continue
    
    # calculate the average and stdev dt for the set of iterations
    if nIterations > 2:
        results_dict['avgdt'] = np.mean(results_dict['dt'])
        results_dict['stddt'] = np.std( results_dict['dt'])
    else:
        results_dict['avgdt'] = 0
        results_dict['stddt'] = 0
    
    return results_dict


#######################################
### PYBINDING FUNCTIONS
#######################################

def pySOIalign_prep_structure(struct_params):
    """
    function to create the alnStruct object that will be used as data entry for
    the pySOIalign.runSOIalign() function
    INPUT:
        :param struct_params: tuple, a path string and an int...
    OUTPUT:
        :return: u_alnStruct, a pySOIalign.alnStruct object filled with the 
                 relevant data for the structure. 
        :return: struct_file_path, string used as input; return only so that
                 the results can be identified
        :return: closeK_opt, int used as input; return only so that the results 
                 can be identified
        :return: dt, float value for elapsed time while preparing the 
                 u_alnStruct object; units: seconds
    """
    struct_file_path = struct_params[0]
    closeK_opt = struct_params[1]

    start_time = time.time()

    # preprocess the structure via an MDAnalysis Universe object; there are so
    # many different ways to do these steps. MDA might not be the most 
    # efficient; BUT, YOU ONLY HAVE TO DO IT ONCE!!!
    u = MDAnalysis.Universe(struct_file_path)
    u_CAs = u.select_atoms('name CA')
    u_len = u_CAs.n_atoms
    u_CAs_coords = np.array(u_CAs.positions.flatten(),dtype=np.float64)
    u_seq = ''.join([tri2single(a.resname) for a in u_CAs])
    
    # run the python bindings of the preprocessing calculations from USalign
    u_sec = pySOIalign.make_sec_py(u_CAs_coords, u_len)
    u_sec_bonds = pySOIalign.assign_sec_bond_py(u_sec, u_len)
    if closeK_opt >= 3:
        u_neighbor_list = pySOIalign.getCloseK_py(u_cas_coords, 
                                                  u_len, 
                                                  closeK_opt)
    # otherwise, the neighbor list will go unused... still need to create an
    # object that the alnStruct class expects though.
    else: 
        u_neighbor_list = np.array([[0,0,0]],dtype=np.float64)

    # create the alnStruct object
    u_alnStruct = pySOIalign.alnStruct(u_CAs_coords, 
                                       u_neighbor_list, 
                                       u_sec_bonds, 
                                       u_seq, 
                                       u_sec, 
                                       u_len)
    
    return u_alnStruct, struct_file_path, closeK_opt, time.time() - start_time


def run_alignment(alnTuple, nIterations):
    """
    function to run the pySOIalign.runSOIalign() function
    INPUT:
        :param alnObjects: tuple/list of three alignment objects
            elem[0]: mobile structure's alnStruct object
            elem[1]: target structure's alnStruct object
            elem[2]: the alnParameters object containing all other parameters
            elem[3]: mobile structure identifier
            elem[4]: target structure identifier
            elem[5]: closeK_opt value
            elem[6]: mm_opt value
    OUTPUT:
        :returns: results dictionary filled with alignment results. 
            keys:
                'translation_vector', np.array of shape (3)
                'rotation_matrix', np.array of shape (9)
                'seqM', mapping string between mobile and target seqs
                'seqA_mobile', seqA_target', aligned sequences for both structs
                'TM1', 'TM2', TMscore floats normalized by mobile and target
                              respectively
                'TM3', 'TM4', 'TM5', only used if a_opt, u_opt, or d_opt
                'd0a', 'd0u', 'd0_scale', only used if above are used
                'rmsd0', final Kabsch RMSD value 
                'd0_out', float for final normalizing distance value
                'Liden', int, number of identical restypes in aligned seq
                'n_ali8', int, number of residues w/in d0_out distance
                'L_ali', int, number of aligned residues in aligned seq
                'd0_mobile', 'd0_target', floats for normalizing distances
                                          for the two structures...
                'dts', list of elapsed times for each iteration
                'avgdt', 'stddt', average and stdev for dts
    """
    identifier = f"aln_{alnTuple[3]}_{alnTuple[4]}_{alnTuple[6]}_{alnTuple[5]}"
    results = {}
    dts = [] 
    for _iter in range(nIterations):
        start_time = time.time()
        results_object = pySOIalign.runSOIalign(u_alnStruct, v_alnStruct, uv_alnParameters)
        dts.append(time.time() - start_time)
        if _iter == 0:
            results['translation_vector'] = results_object.translation_vector
            results['rotation_matrix'] = results_object.rotation_matrix.reshape((3,3))
            results['TMscore1'] = results_object.TM1
            results['TMscore2'] = results_object.TM2
            results['RMSD'] = results_object.rmsd0
            #results['SeqIDAli'] = results_object.
            #results['Len1'] = results_object.
            #results['Len2'] = results_object.
            results['LenAligned'] = results_object.L_ali
            # !!! NOTE, not passing almost all of the important results... need
            # to implement some output formatting function. Think about this 
            # some more
        else:
            assert results['translation_vector'] == results_object.translation_vector
            assert results['rotation_matrix'] == results_object.rotation_matrix
            assert results['TMscore1'] == results_object.TM1
            assert results['TMscore2'] == results_object.TM2
            assert results['RMSD'] == results_object.rmsd0

    # calculate the average and stdev dt for the set of iterations
    if nIterations > 2:
        results['avgdt'] = np.mean(dts)
        results['stddt'] = np.std( dts)
    else:
        results_dict['avgdt'] = 0
        results_dict['stddt'] = 0
    
    return results, identifier


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Structural alignment task manager')
    # essential parameters
    ##parser.add_argument('--pySOIalign-path', 
    ##                    '-pySOIalign', 
    ##                    required=True, 
    ##                    help='global path string pointing to the directory of the pySOIalign .so file',
    ##                    type=str)
    parser.add_argument('--tskmgr-log-file', 
                        '-log', 
                        required=True, 
                        help='string that will be used to store logging info for this run',
                        type=str)
    parser.add_argument('--alignments-list-file', 
                        '-inp', 
                        required=True, 
                        help='list file where each line consists of the paths for the two protein models that are to be aligned',
                        type=str)
    parser.add_argument('--usalign-path', 
                        '-path', 
                        required=True,
                        help='path to the USalign executable',
                        type=str)
    # optional parameters
    parser.add_argument('--scheduler-file', 
                        '-s', 
                        required=False,
                        default=None, 
                        help='dask scheduler file')
    parser.add_argument('--nIterations', 
                        '-niters', 
                        required=False,
                        default=1, 
                        help='number of iterations of each alignment to perform',
                        type=int)
    args = parser.parse_args()


    #########################################
    ##### IMPORTING SOIALIGN MODULE
    #########################################
    ##try:
    ##    sys.path.insert(0, pySOIalign_path)
    ##    import pySOIalign
    ##except Exception as e:
    ##    print("Failed to load pySOIalign shared file as a module. Check compilation and/or path of the .so file.", e, file=sys.stderr, flush=True)
    ##    sys.exit(1)

    #######################################
    ### 
    #######################################
    # set up the main logger file and list all relevant parameters.
    main_logger = setup_logger('tskmgr_logger',args.tskmgr_log_file)
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Path to the USalign executable: {args.usalign_path}')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Alignments are listed in {args.alignments_list_file}.')
    main_logger.info(f'Each alignment will be performed {args.nIterations} times.')

    ## unimportant for testing suite
    #dask_parameter_string = ''
    #for key, value in dask.config.config.items():
    #    dask_parameter_string += f"'{key}': '{value}'\n"
    #dask_parameter_string += '################################################################################'
    #main_logger.info(f'\n################################################################################\nDask parameters:\n{dask_parameter_string}')

    # start dask client.
    # if a scheduler file is input, then the CLI is being used, where the user 
    # has defined a scheduler and number of workers. Likely has also defined 
    # what resources each worker gets.
    if args.scheduler_file:
        client = Client(scheduler_file=args.scheduler_file,
                        timeout=120,
                        name='AlignmentTaskMgr')
    # no scheduler file, so assume that the default Client behavior is expected
    # (spin up a scheduler and worker within the Client call); since tasks in 
    # this workflow are known to only need 1 thread, we hard wire the client to
    # spin up the appropriate resourced workers.
    else:
        client = Client(timeout=120,
                        name='AlignmentTaskMgr',
                        threads_per_worker=1)
   
    # parsing the alignments_list_file
    # alignment_list represents the list of structure pairs that will be aligned 
    main_logger.info(f'Reading the alignments file.')
    with open(args.alignments_list_file,'r') as structures_file:
        # get the specific alignments to be performed
        alignments_list = [line.split() for line in structures_file.readlines() if line[0] != '#']
        # get the structures that need to be pre-processed
        structure_closeK_list = set([(aln[i],int(aln[3])) for aln in alignments_list for i in range(2)])
    
    main_logger.info(f'A total of {len(structure_closeK_list):,} structures will be preprocessed.')
    main_logger.info(f'A total of {len(alignments_list)*args.nIterations:,} alignments will be performed.')
    
    #######################################
    ### Native USalign Functions, run via subprocess
    #######################################
    
    # create the preprocess futures
    preprocess_futures = client.map(submit_preprocessing_pipeline,
                                    structure_closeK_list) 

    # create the alignment futures
    aln_futures = client.map(submit_subprocess_pipeline, 
                             alignments_list, 
                             nIterations = args.nIterations,
                             USalign_executable_path = args.usalign_path,
                             pure=False)

    futures = preprocess_futures + aln_futures

    # subproc_results_dict will be a dict of dicts, with keys being a str to
    # identify the process and values being the respective results subdict
    subproc_results_dict = {}
    completed = as_completed(futures)
    # loop over completed tasks
    for finished_task in completed:
        # gather return values from the task
        results = finished_task.result()
        if results['task_type'] == 'preprocess':
            struct = Path(results['structure']).stem
            identifier = f"pre_{struct}_{results['closeK_opt']}"
            main_logger.info(f"Preprocessing for {results['structure']} finished; took {results['dt']} seconds.")
            
        elif results['task_type'] == 'alignment':
            mobile = Path(results['mobile']).stem
            target = Path(results['target']).stem
            identifier = f"aln_{mobile}_{target}_{results['mm_opt']}_{results['closeK_opt']}"
            main_logger.info(f"-mm {results['mm_opt']} -closeK_opt {results['closeK_opt']} alignment between {mobile} and {target} has finished.")
        
        subproc_results_dict[identifier] = results
    
    # all native USalign tasks have finished, write out results to json file
    with open('USalign_subprocess_results.json','w') as json_out:
        json.dump(subproc_results_dict, json_out, indent = 4)

    with open(f'USalign_subprocess_results.dat','w') as out_file:
        out_file.write('Query_Target_Method,TMscore1,TMscore2,RMSD,SeqIDAli,Len1,Len2,LenAligned),Avg Time (s),Stdev Time (s)\n')
        sorted_keys = [key for key in subproc_results_dict.keys() if key[:3] == 'aln']
        for key in sorted_keys:
            results  = subproc_results_dict[key]
            TMscore1 = results['TMscore1']
            TMscore2 = results['TMscore2']
            RMSD     = results['RMSD']
            SeqIDAli = results['SeqIDAli']
            Len1     = results['Len1']
            Len2     = results['Len2']
            LenAlign = results['LenAligned']
            avg_time = results['avgdt']
            std_time = results['stddt']
            out_file.write(f'{key},{TMscore1},{TMscore2},{RMSD},{SeqIDAli},{Len1},{Len2},{LenAlign},{avg_time:.3f},{std_time:.3f}\n')


   ## #######################################
   ## ### pySOIalign Functions
   ## #######################################

   ## # run the preprocessing calculations before creating the alnment futures
   ## # so that the alnStruct objects can be handed off efficiently
   ## preprocess_futures = client.map(pySOIalign_prep_structure,
   ##                                 structure_closeK_list)

   ## # prep the results_dict collection for the pySOIalign results
   ## pySOIalign_results_dict = {}
   ## completed = as_completed(preprocess_futures)
   ## # loop over completed tasks
   ## for finished_task in completed:
   ##     # gather return values from the task
   ##     alnStruct_obj, struct_file, closeK_opt, dt = finished_task.result()
   ##     struct = Path(struct_file).stem
   ##     identifier = f"pre_{struct}_{closeK_opt}"
   ##     main_logger.info(f"Python-based preprocessing for {struct_file} finished; took {dt} seconds.")
   ##     pySOIalign_results_dict[identifier] = {'alnStruct': alnStruct_obj,
   ##                                            'dt': dt}
   ## 
   ## # creating the collection of input data for the pySOIalign.runSOIalign() 
   ## # function; need the two alnStruct objects and an alnParameters object
   ## alignment_tuples = []
   ## for alignment in alignments_list:
   ##     # get the alnStruct object associated with the mobile structure
   ##     mobile_identifier = f"pre_{Path(alignment[0]).stem}_{alignment[3]}"
   ##     mobile_alnStruct = pySOIalign_results_dict[mobile_identifier]['alnStruct']
   ##     
   ##     # get the alnStruct object associated with the target structure
   ##     target_identifier = f"pre_{Path(alignment[1]).stem}_{alignment[3]}"
   ##     target_alnStruct = pySOIalign_results_dict[target_identifier]['alnStruct']
   ## 
   ##     # create the alnParameters object to be fed into the alignment tasks
   ##     alnParameters = pySOIalign.alnParameters(int(alignment[3]), # closeK_opt
   ##                                              -2,    # molec_type
   ##                                              int(alignment[2]), # mm_opt
   ##                                              0,     # a_opt
   ##                                              False, # u_opt
   ##                                              0.0,   # Lnorm_ass
   ##                                              False, # d_opt
   ##                                              0.0,   # d0_scale
   ##                                              False) # fast_opt
   ##     # fill the list with tuples of the grouped input objects
   ##     alignment_tuples.append((mobile_alnStruct,
   ##                              target_alnStruct,
   ##                              alnParameters,
   ##                              Path(alignment[0]).stem,
   ##                              Path(alignment[1]).stem,
   ##                              alignment[3],
   ##                              alignment[2]))
   ## 
   ## # create the future objects for the pySOIalign tasks
   ## aln_futures = client.map(pySOIalign_pipeline,
   ##                          alignment_tuples,
   ##                          nIterations = args.nIterations,)

   ## completed = as_completed(aln_futures)
   ## # loop over completed tasks
   ## for finished_task in completed:
   ##     # gather return values from the task
   ## 







    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    clean_logger(main_logger)

