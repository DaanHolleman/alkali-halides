"""
Import the displace function into a python terminal.
Then call e.g. displace('-m step') to generate a step displacement.
"""

#%% IMPORT
import numpy as np
import os, sys, pickle
import json
from pymatgen.core import Structure, Lattice
from tabulate import tabulate
from ..toolbox.helper import leading_zeros
from ..toolbox.files import select
from .scripts_displace.methods import line_cell, line_cart, mag_cell, mag_cart, zero, plane_cell, volume_cell, step_cell, step_cartesian, shell_cell, shell_cartesian, shell_cartesian_oct, shell_cell_oct
import argparse, shlex, glob

CARDNAMES = ['&CONTROL', '&SYSTEM', '&ELECTRONS', '&IONS', '&CELL', '&FCP', '&RISM', 
             'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS', 'ADDITIONAL_K_POINTS', 'CELL_PARAMETERS', 
             'CONSTRAINTS', 'OCCUPATIONS', 'ATOMIC_VELOCITIES', 'ATOMIC_FORCES', 'SOLVENTS', 'HUBBARD']

#%%

def get_method_keys():
    method_keys = ['line','line-cart','mag','mag-cart',
                   'zero',
                   'plane','volume',
                   'step','step-cart',
                   'shell','shell-cart','shell-cart-oct','shell-oct']
    
    return method_keys

def parse_method(method:str):
    """
    Ask for user input and return relative quantities.
    Returns primitive cell, species, coords, and displacement 
    """
    # Create method dictionary
    method_keys = get_method_keys()
    routines = [line_cell, line_cart, mag_cell, mag_cart,
                zero,
                plane_cell, volume_cell, 
                step_cell, step_cartesian, 
                shell_cell, shell_cartesian, shell_cartesian_oct, shell_cell_oct]
    method_dict = { key:func for key, func in zip(method_keys, routines) }
    
    # Default option is first in method_keys (line)
    if method is None:
        method = method_keys[0]
    
    # Clean up user input
    method = method.lower().strip()
    
    # Raise error if not an option
    if method not in method_keys:
        raise ValueError(f'Option {method} is not a valid option. Please choose from:\n\t{method_keys}')
    
    # Choose routine
    method_routine = method_dict[method]
    
    # Gather input data and return
    QE_data = ask_QE_atom()
    input_data = method_routine(QE_data)
    
    return input_data

def load_input():
    """
    Load a previously used input. Binary file must be present in current working directory.
    """
    with open(cf.SAVEFILE, 'rb') as f:
        data = pickle.load(f)
    return data

def save_input(data):
    """
    Save inputs to binary file.
    """
    with open(cf.SAVEFILE, 'wb') as f:
        pickle.dump(data, f)

def argv_value(argv:list, keys:list):
    [short, long] = keys
    lookfor = short if short in argv else long
    return (argv + [''])[argv.index(lookfor) + 1]


def parse_argv():
    parser = argparse.ArgumentParser(
        prog = 'py_displace',
        description = 'Displaces atoms in various ways',
        epilog = 'The following methods are available:\n' + get_help_string()
    )
    parser.add_argument('-n','--nodir', action='store_true',
                        help='Do not create subdirectories for each file')
    
    parser.add_argument('-l','--load', action='store_true', 
                        help='Load a previous result from displace.bin')
    
    parser.add_argument('-s','--save', action='store_true',
                        help='Save results to displace.bin')
    
    parser.add_argument('-m','--method', default=get_method_keys()[0], choices=get_method_keys(), metavar = '',
                        help='Method used in displacing')
    
    parser.add_argument('-f','--find', action='store_true',
                        help='Automatically find Quantum Espresso input files')

    parser.add_argument('--SAVEFILE', action='store_const', default='./displace.bin', const='./displace.bin')
    
    cf = parser.parse_args()
    
    return cf
    

#%%

def ask_QE_atom():
    # Find QE input file
    if cf.find:
        fn = select( glob.glob('*.in') )
    else:
        user = input('Enter Quantum Espresso file.\n>>> ')
        fn = user.strip()
    
    rprim, species, coords = read_QE(fn)
    coupled = [[label, *list(abc)] for label, abc in zip(species, coords)]
    headers = ['Atom','A','B','C']
    table = tabulate(coupled, tablefmt='plain', headers=headers, showindex = True)
    print(table)
    user = input('Choose atom index to displace.\n>>> ')
    move_index = int(user)
    
    return fn, move_index, rprim, species, coords

def get_card(lines, cardname):
    """
    Retrieve cards from a Quantum Espresso file.
    lines: list of strings.
    cardname: string with the name of the card.
    """
    contents = []
    walker = iter(lines)
    for line in walker:
        if cardname in line:
            break
    else:
        return contents
    for line in walker:
        line = line.strip()
        if len(line) == 0 or any([ line.index(cardname) == 0 for cardname in CARDNAMES if cardname in line.upper()]):
            break
        contents += [line]
    return contents

def read_QE(filename):
    """
    Read relevant cards from a Quantum Espresso file.
    filename: string with the filename.
    """
    lines = open(filename,'r').readlines()
    # CARDS
    cell_pars = get_card(lines, 'CELL_PARAMETERS')
    atomic_pos = get_card(lines, 'ATOMIC_POSITIONS')
    
    species = [ line.split()[0] for line in atomic_pos ]
    coords = np.array([ line.split()[1:4] for line in atomic_pos ], float)
    rprim = np.array([ line.split() for line in cell_pars ], float)
    
    return rprim, species, coords

def create_json(fn, rprim, species, coords):
    """
    Create a json file containing the pymatgen structure.
    fn : str
        filename of the Quantum Espresso file used to make the structure.
    rprim : array (3,3)
        the real space primitive cell.
    species : list (N)
        which atoms are present.
    coords : array (N,3)
        the coordinates of the atoms in the cell corresponding to species.
    """
    structure = Structure(
        lattice = Lattice(rprim),
        species = species,
        coords = coords
    )
    fn = fn[:-3] + '.json' if fn[-3:] == '.in' else fn
    with open(fn,'w') as file:
        json.dump(structure.as_dict(), file)
    return fn

def loop_displacements(fn, move_index, rprim, species, pos_abc, dis_abc):
    """
    Loop through displacements and write to json files.
    """
    # Leading zeros in dir/file names
    N10 = leading_zeros(dis_abc)
    
    ## CREATE JASONS
    npos_abc = pos_abc.copy()
    json_files = ['' for ii in range(len(dis_abc))]
    for ii, dis in enumerate(dis_abc):
        # adjust position
        npos_abc[move_index] = pos_abc[move_index] + dis
        # location of file
        new_dir = f'D{ii:0{N10}d}'
        if not cf.nodir:
            new_fn = new_dir + '/' + fn
            if not os.path.exists(new_dir): 
                os.makedirs(new_dir)
        else:
            new_fn = new_dir + '-' + fn
        json_fn = create_json(new_fn, rprim, species, npos_abc)
        json_files[ii] = json_fn
    
    return json_files

def read_coords(fn, move_index):
    file = json.load(open(fn))
    atom = file['sites'][move_index]
    coord = np.array(atom['xyz'], float)
    return coord

#%%

def stdout(input_data, json_files):
    """
    Print file of user inputs.
    """
    fn, move_index, rprim, species, pos_abc, dis_abc = input_data
    
    out  = 'USER INPUT\n'
    out += f'Input file:\n\t{fn}\n'
    out += f'Moved atom:\n\t[{move_index}] {species[move_index]}\n'
    out += f'Atom positions [abc]:\n{pos_abc}\n'
    out += f'Steps:\n\t{len(dis_abc)}\n'
    out += f'Displacement [abc]:\n{dis_abc}\n'
    out += 'Displacement [xyz]:\n'
    coords = np.array([ read_coords(json_fn, move_index) for json_fn in json_files ])
    out += f'{coords}\n'
    with open('displace.out','w') as file:
        file.write(out)
    print('Written user input to displace.out')

def get_help_string():
    return '\t'.join(get_method_keys())

#%%

def displace():
    """
    Main method for creating displaced structure files in json format.
    """
    ## ARGV
    global cf
    cf = parse_argv()
    
    ## Data handling
    if cf.load:
        input_data = load_input() 
    else:
        input_data = parse_method(cf.method)
    if cf.save: 
        save_input(input_data)
    
    # Create json files
    json_files = loop_displacements(*input_data)
    stdout(input_data, json_files)

def main():
    displace()
