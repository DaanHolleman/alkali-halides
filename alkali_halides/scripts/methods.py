# -*- coding: utf-8 -*-

import numpy as np
from scipy.linalg import inv


#%% GRID CREATION FUNCTIONS

def line(QE_data:tuple, cart = False, mag = False):
    """
    Generate structures with displacements along a single vector.
    """
    # Unpack
    fn, move_index, rprim, species, pos_abc = QE_data
    
    # Vector of displacement
    if mag:
        msg = f'Enter displacement vector [{"xyz m" if cart else "abc m"}].\n>>> '
        user = input(msg)
        dis_vec = np.array(user.strip().split(),float)
        u_vec = dis_vec[0:3] / np.sqrt( np.sum(dis_vec[0:3]**2) )
        mag = dis_vec[3]
        dis_vec = mag * u_vec
    else:
        msg = f'Enter displacement vector [{"xyz" if cart else "abc"}].\n>>> '
        user = input(msg)
        dis_vec = np.array(user.strip().split(),float)
    
    # Displacement range
    user = input('Enter displacement range. (Start Stop Number)\n>>> ')
    [start, stop, number] = np.array(user.strip().split(),float)
    number = int(number)
    
    # Create displacement coordinates
    mult = np.linspace(start, stop, number)
    dis_abc = dis_vec * mult[:,None]
    
    return fn, move_index, rprim, species, pos_abc, dis_abc

def line_cart(QE_data:tuple):
    """
    Generate structures with displacements along a single vector.
    """
    fn, move_index, rprim, species, pos_abc, dis_xyz = line(QE_data)
    
    # Project onto cell basis
    inv_rprim = inv(rprim)
    dis_abc = dis_xyz @ inv_rprim

    return fn, move_index, rprim, species, pos_abc, dis_abc
def line_cell(QE_data:tuple):
    """
    Generate structures with displacements along a single vector.
    """
    return line(QE_data)
def mag_cart(QE_data:tuple):
    """
    Generate structures by specifying a direction in cartesian coordinates and magnitude.
    """
    fn, move_index, rprim, species, pos_abc, dis_xyz = line(QE_data, cart = True, mag = True)
    # Project onto cell basis
    inv_rprim = inv(rprim)
    dis_abc = dis_xyz @ inv_rprim
    return fn, move_index, rprim, species, pos_abc, dis_abc
def mag_cell(QE_data:tuple):
    """
    Generate structures by specifying a direction in cell coordinates and magnitude.
    """
    return line(QE_data, mag = True)


def zero(QE_data:tuple):
    """
    Zero displacement. Effectively converts QE input to JSON.
    """
    # Unpack
    fn, move_index, rprim, species, pos_abc = QE_data
    
    # No displacement
    dis_abc = np.zeros((1,3), float)

    return fn, move_index, rprim, species, pos_abc, dis_abc


def plane_cell(QE_data:tuple):
    """
    Generate structures with displacements distributed according to a plane made from two cell basis 
    vectors.
    """
    make_vec = lambda user : np.array(user.strip().split(),float)
    
    fn, move_index, rprim, species, pos_abc = QE_data
    
    # Vectors of displacement
    dis_vecs = np.zeros((2,3))
    user = input('Enter displacement vectors.\nV1\n>>> ')
    dis_vecs[0] = make_vec(user)
    user = input('V2\n>>> ')
    dis_vecs[1] = make_vec(user)
    
    # Displacement ranges
    ranges = np.zeros((2,3), object)
    user = input('Enter displacement range. (Start Stop Number)\nV1\n>>> ')
    ranges[0] = make_vec(user)
    prev_user = user
    user = input('V2 [previous]\n>>> ')
    if len(user.strip()) == 0: 
        user = prev_user
    ranges[1] = make_vec(user)
    ranges[:,2] = ranges[:,2].astype(int)
    mults = [ np.linspace(*ssn) for ssn in ranges ]
    
    dis_abcs = [ mult[:,None] * dis_vec for mult, dis_vec in zip(mults, dis_vecs) ]
    dis_abc = dis_abcs[0]
    for next_dis_abc in dis_abcs[1:]:
        dis_abc = add_convolve(dis_abc, next_dis_abc)
    
    return fn, move_index, rprim, species, pos_abc, dis_abc
    
    

def volume_cell(QE_data:tuple):
    """
    Generate structures with displacements distributed according to a volume made up of the cell basis
    vectors.
    """
    make_vec = lambda user : np.array(user.strip().split(),float)
    
    fn, move_index, rprim, species, pos_abc = QE_data
    
    # Vectors of displacement
    dis_vecs = np.zeros((3,3))
    user = input('Enter displacement vectors.\nV1\n>>> ')
    dis_vecs[0] = make_vec(user)
    user = input('V2\n>>> ')
    dis_vecs[1] = make_vec(user)
    user = input('V3\n>>> ')
    dis_vecs[2] = make_vec(user)
    
    # Displacement ranges
    ranges = np.zeros((3,3), object)
    user = input('Enter displacement range. (Start Stop Number)\nV1\n>>> ')
    ranges[0] = make_vec(user)
    prev_user = user
    user = input('V2 [previous]\n>>> ')
    if len(user.strip()) == 0: 
        user = prev_user
    ranges[1] = make_vec(user)
    user = input('V3 [previous]\n>>> ')
    if len(user.strip()) == 0: 
        user = prev_user
    ranges[2] = make_vec(user)
    ranges[:,2] = ranges[:,2].astype(int)
    mults = [ np.linspace(*ssn) for ssn in ranges ]
    
    dis_abcs = [ mult[:,None] * dis_vec for mult, dis_vec in zip(mults, dis_vecs) ]
    dis_abc = dis_abcs[0]
    for next_dis_abc in dis_abcs[1:]:
        dis_abc = add_convolve(dis_abc, next_dis_abc)
    
    return fn, move_index, rprim, species, pos_abc, dis_abc
    

def step_cartesian(QE_data:tuple):
    """
    Generate 4 structures with displacement along cartesian axes. 
    The first structure has no displacement.
    """
    
    fn, move_index, rprim, species, pos_abc = QE_data
    
    # Displacement vectors are the inverse of the rprim matrix
    dis_vecs = inv(rprim)
    
    # Displacement ranges
    user = input('Step size in Angstrom:\n>>> ')
    step = float( user.strip() )
    
    dis_abc = np.array([ dis_vec * step for dis_vec in dis_vecs ], float)
    dis_abc = np.append([[0,0,0]], dis_abc, axis=0)
    
    return fn, move_index, rprim, species, pos_abc, dis_abc

def step_cell(QE_data:tuple):
    """
    Generate 4 structures with displacement along cartesian axes. 
    The first structure has no displacement.
    """
    fn, move_index, rprim, species, pos_abc = QE_data
    
    # Displacement vectors are the inverse of the rprim matrix
    dis_vecs = np.identity(3)
    
    # Displacement ranges
    user = input('Step size in alat:\n>>> ')
    step = float( user.strip() )
    
    dis_abc = np.array([ dis_vec * step for dis_vec in dis_vecs ], float)
    dis_abc = np.append([[0,0,0]], dis_abc, axis=0)
    
    return fn, move_index, rprim, species, pos_abc, dis_abc

def shell_cartesian(QE_data:tuple):
    """
    Create shell by projecting a cube onto a cartesian sphere of a given radius
    """
    
    fn, move_index, rprim, species, pos_abc, dis_abc = shell(QE_data, 'Angstrom')
    
    # Project onto cell basis
    inv_rprim = inv(rprim)
    dis_abc = dis_abc @ inv_rprim # move to cell basis
    
    return fn, move_index, rprim, species, pos_abc, dis_abc

def shell_cartesian_oct(QE_data:tuple):
    """
    Create shell by projecting a cube onto a cartesian sphere of a given radius
    """
    
    fn, move_index, rprim, species, pos_abc, dis_abc = shell(QE_data, 'Angstrom', True)
    
    # Project onto cell basis
    inv_rprim = inv(rprim)
    dis_abc = dis_abc @ inv_rprim # move to cell basis
    
    return fn, move_index, rprim, species, pos_abc, dis_abc

def shell_cell(QE_data:tuple):
    """
    Create shell by projecting a cube onto a cell-basis sphere of a given radius
    """
    fn, move_index, rprim, species, pos_abc, dis_abc = shell(QE_data, 'alat')
    return fn, move_index, rprim, species, pos_abc, dis_abc

def shell_cell_oct(QE_data:tuple):
    """
    Create shell by projecting a cube onto a cell-basis sphere of a given radius
    """
    fn, move_index, rprim, species, pos_abc, dis_abc = shell(QE_data, 'alat', True)
    return fn, move_index, rprim, species, pos_abc, dis_abc

def shell(QE_data:tuple, units, octant:bool = False):
    """
    Create shell by projecting a cube onto a sphere of a given radius
    """    
    fn, move_index, rprim, species, pos_abc = QE_data
    
    dis_vecs = np.identity(3)
    
    # Displacement ranges
    user = input(f'Shell size in {units}:\n>>> ')
    shell = float( user.strip() )
    
    # Number of points on one axis
    user = input('Number of points per axis:\n>>> ')
    N = int( user.strip() )
    
    # Generate cube
    if octant:
        mult = np.linspace(0,1,N)
    else:
        mult = np.linspace(0,1,N)
    mults = [mult] * 3
    
    dis_abcs = [ mult[:,None] * dis_vec for mult, dis_vec in zip(mults, dis_vecs) ]
    xy = dis_abcs[0]
    xy = add_convolve(xy, dis_abcs[1])
    xz = dis_abcs[0]
    xz = add_convolve(xz, dis_abcs[2])
    yz = dis_abcs[1]
    yz = add_convolve(yz, dis_abcs[2])
    
    if octant:
        dis_abc = np.concatenate((xy+dis_vecs[2],
                                  xz+dis_vecs[1],
                                  yz+dis_vecs[0]))
    else:
        dis_abc = np.concatenate((xy-dis_vecs[2], xy+dis_vecs[2],
                                  xz-dis_vecs[1], xz+dis_vecs[1],
                                  yz-dis_vecs[0], yz+dis_vecs[0]))
    
    norm = np.sqrt(np.sum(dis_abc**2,axis=1))
    dis_abc = dis_abc / norm[:,None] * shell
    
    return fn, move_index, rprim, species, pos_abc, dis_abc
