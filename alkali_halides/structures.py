# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 13:34:18 2024

@author: dholl
"""

import numpy as np
from .attrdict import AttrDict

structures = AttrDict()
structures.fcc = AttrDict()
structures.fcc.high_symmetry = AttrDict(
    G = [0.000, 0.000, 0.000],
    K = [0.375, 0.375, 0.750],
    L = [0.500, 0.500, 0.500],
    U = [0.625, 0.250, 0.625],
    W = [0.500, 0.250, 0.750],
    W2= [0.750, 0.250, 0.500],
    X = [0.500, 0.000, 0.500],
)
structures.fcc.rprim = np.array([[0,1,1],[1,0,1],[1,1,0]], float)
structures.fcc.coordinates = np.array([[0,0,0],[.5, .5, .5]], float)
structures.fcc.basic_to_primitive = 0.5

# TODO : Add the high symmetry points of the BCC brillouin zone. These are from FCC.
structures.bcc = AttrDict()
structures.bcc.high_symmetry = AttrDict(
    G = [0.000, 0.000, 0.000],
    # K = [0.375, 0.375, 0.750],
    # L = [0.500, 0.500, 0.500],
    # U = [0.625, 0.250, 0.625],
    # W = [0.500, 0.250, 0.750],
    # W2= [0.750, 0.250, 0.500],
    # X = [0.500, 0.000, 0.500],
)
structures.bcc.rprim = np.array([[0,1,1],[1,0,1],[1,1,0]], float)
structures.bcc.coordinates = np.array([[0,0,0],[.5, .5, .5]], float)
structures.bcc.basic_to_primitive = 0.5
