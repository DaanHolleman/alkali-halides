# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 13:37:49 2024

@author: dholl
"""
import numpy as np
from io import StringIO
from .crystal import Crystal
from .attrdict import AttrDict
import os, sys

def load_database(filename):
    filename = os.path.join( os.path.dirname(__file__), 'data', filename ) # TODO : FIND SOMETHING BETTER FOR THIS UGLY PATH
    with open(filename) as file:
        header = file.readline().strip().split(';')
        dtypes = [ dtype for dtype in file.readline().strip().split(';') if len(dtype) != 0 ]
        remainder = ''.join(file.readlines())
    file_io = StringIO(remainder)
    database = np.genfromtxt(file_io, dtypes, delimiter = ';')
    return header, database

def convert_to_dictionary(header, database, dict_database = {}):
    for entry in database:
        key_crystal = entry[0]
        if key_crystal == '': continue
        dict_crystal = dict_database.get(key_crystal, {})
        runner = zip(header, entry)
        next(runner)
        for key, value in runner:
            dict_crystal[key] = value
        dict_database[key_crystal] = dict_crystal
    return dict_database

def construct_database():
    dict_database = {}
    for fn in FILENAMES:
        header, database = load_database(fn)
        dict_database = convert_to_dictionary(header, database, dict_database)
    return dict_database

def make_species(crystal_key):
    crystal = DICT_DATABASE[crystal_key]
    alkali = str(crystal['alkali'])
    halide = str(crystal['halide'])
    return [alkali, halide]

def get_crystal(crystal_key):
    """
    Build the crystal that is defined by the crystal_key (e.g. 'LiF')

    Parameters
    ----------
    crystal_key : str
        Label of the alkali halide.

    Returns a crystal object.

    """
    return Crystal(make_species(crystal_key), **DICT_DATABASE[crystal_key])

def get_all_crystals():
    """
    Builds all 20 alkali halides.

    Returns a list of these.
    """
    crystals = AttrDict()
    for specie, [crystal, kwargs] in zip(species, DICT_DATABASE.items()):
        crystals[crystal] = Crystal(specie, **kwargs)
    return crystals
    
# Define filenames
FILENAME_CALCULATED = 'calculated.csv'
FILENAME_SETTINGS   = 'settings.csv'
FILENAME_LITERATURE = 'literature.csv'
FILENAMES = [FILENAME_CALCULATED, FILENAME_SETTINGS, FILENAME_LITERATURE]

# Construct database from these
DICT_DATABASE = construct_database()
species = [ make_species(crystal_key) for crystal_key in DICT_DATABASE.keys() ]
crystals = get_all_crystals()
