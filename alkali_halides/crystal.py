# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 13:38:34 2024

@author: dholl
"""
from .attrdict import AttrDict
from .structures import structures

from pymatgen.core import Structure, Lattice
from pymatgen.transformations.standard_transformations import PerturbStructureTransformation, SupercellTransformation
import numpy as np

class Crystal(object):
    def __init__(self, species, *args, **kwargs):
        # TODO : Use @property @setter for name protection
        self.species = species
        self.crystal = ''.join(species)
        self.alkali, self.halide = species
        self.valence = kwargs.get('valence')
        
        # Set properties from literature
        lit = AttrDict()
        lit.structure = kwargs.get('lit_structure')
        lit.a0 = kwargs.get('lit_a0')
        lit.Eg = kwargs.get('lit_Eg')
        lit.E1s = kwargs.get('lit_E1s')
        lit.eps0 = kwargs.get('lit_eps0')
        lit.epsinf = kwargs.get('lit_epsinf')
        self.lit = lit
        
        # Set properties calculated by self
        calc = AttrDict()
        calc.a0 = kwargs.get('calc_a0')
        calc.eps0 = kwargs.get('calc_eps0')
        calc.epsinf = kwargs.get('calc_epsinf')
        calc.Eg = kwargs.get('calc_Eg')
        calc.E1s = kwargs.get('calc_E1s')
        self.calc = calc
        
        # Set properties found from convergence
        conv = AttrDict()
        conv.pressure = kwargs.get('conv_pressure')
        conv.eps = kwargs.get('conv_eps')
        conv.total_energy = kwargs.get('conv_total_energy')
        self.conv = conv
        
        # BGW / QE Settings
        settings = AttrDict()
        settings.structure = kwargs.get('set_structure')
        settings.nbnd = kwargs.get('set_nbnd')
        settings.ecutwfc = kwargs.get('set_ecutwfc')
        settings.ngkpt_scf = np.array([ kwargs.get('set_kpoints_scf') ] * 3)
        settings.ngkpt_co  = np.array([ kwargs.get('set_kpoints_co')  ] * 3)
        settings.ngkpt_fi  = np.array([ kwargs.get('set_kpoints_fi')  ] * 3)
        settings.fft = np.array([ kwargs.get('set_fft', 0) ] * 3)
        settings.ecuteps = kwargs.get('set_ecuteps')
        settings.screened_cutoff = kwargs.get('set_screened_cutoff')
        settings.ecutsig = settings.screened_cutoff
        self.settings = settings
        
        # Duplicates
        self.settings.valence = self.valence
        self.total_energy = self.conv.total_energy
        self.prefix = self.crystal
        self.nbnd = self.settings.nbnd
        
        # Structure is fcc, bcc, etc.
        self.set_structure( kwargs.get('use_literature_structure', False) )
    
    def set_structure(self, use_literature:bool = False):
        code = self.lit.structure if use_literature else self.settings.structure
        if code is None:
            raise ValueError(f'Structure of {self} is not defined.')
        self.structure = structures[ code ]
    
    @property
    def bgwpy_kwargs(self): # TODO: return bgwpy ready kwargs
        bgwpy = dict(
        )
        return bgwpy
    
    def __repr__(self):
        return f'<Crystal({self.crystal})>'
    
    @property
    def pseudos(self):
        ext = '.upf'
        return [self.alkali + ext, self.halide + ext]
    
    def build_structure(self, supercell:int = None, perturbed:float = None, round_to_em8:bool = True):
        """
        Generate a crystal structure to build using pymatgen.

        Parameters
        ----------
        supercell : int(3), optional
            The supercell specifications. The default of None means no supercell. 
            If only 1 integer N is provided, the result will be a (N,N,N) supercell.
            If 3 integers X, Y, Z are provided, the results will be a (X,Y,Z) supercell.
        perturbed : float(2), optional
            Perturbation applied to the crystal (after any supercell operations). 
            The default is no perturbation.
            See PerturbStructureTransformation for more information.
        round_to_em8 : bool (false), optional
            Whether to round to 8 decimal places after applying any transformations.
        
        Returns
        -------
        structure : Structure
            pymatgen Structure after applying all appropriate transformations.

        """
        # Create the initial crystal structure
        
        species = self.species
        a0 = self.calc.a0 * self.structure.basic_to_primitive
        
        rprim = Lattice( a0 * self.structure.rprim )
        coords = self.structure.coordinates
        structure = Structure(rprim, species, coords)
        
        # Apply supercell transformation
        if supercell is not None:
            if type(supercell) is int:
                supercell = (supercell, supercell, supercell)
            trans_supercell = SupercellTransformation(supercell)
            structure = trans_supercell.apply_transformation(structure)
        
        # Apply perturbations
        if perturbed is not None:
            if type(perturbed) is tuple:
                max_dist, min_dist = perturbed
            else:
                max_dist, min_dist = perturbed, None
            trans_perturb = PerturbStructureTransformation(max_dist, min_dist)
            structure = trans_perturb.apply_transformation(structure)
        
        # Round coordinates
        if round_to_em8:
            # Round fractional coordinates to 8 decimal places
            from numpy import around
            rounded_coords = around(structure.frac_coords, decimals=8)
            # Create new structure
            structure = Structure(
                lattice=structure.lattice,  # Use the same lattice
                species=structure.species,  # Use the same species
                labels=structure.labels,    # Use the same labels
                coords=rounded_coords,      # Use the rounded coordinates
                coords_are_cartesian=False  # Indicate that the coords are cell coordinates
            )
        
        return structure
