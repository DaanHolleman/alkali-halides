# alkali-halides
*A Python package that implements an alkali halide database for use with BGWpy*

For ease of use with BGWpy, a python package called `alkali_halides` was written that
implements twenty alkali halide crystals as objects that hold the information needed to
generate input scripts for calculations, as well as some additional information about convergence, and properties from Song & Williams.  

The crystal structures of CsCl, CsBr, and CsI that were used are face centered cubic. However,
the ground state structures of these crystals are body centered cubic. For the sake of being able
to compare all alkali halides, it was decided to keep these as face centered cubic.  

Each crystal of the alkali halides is given four attribute dictionaries called `settings`, `lit`,
`calc`, and `conv`. All settings needed for the input script are stored in the settings attribute
dictionary. The values from literature are stored in lit and their calculated counterparts in `calc`. The values found from convergence tests are stored in `conv`.  

The implementation of the attribute dictionary allows the user to obtain the values from the
dictionary, not just from keywords, but also from attribute collection. For example, the value for the wavefunction cutoff, stored in the settings dictionary, can be obtained via `settings.ecutwfc` in addition to the regular method `settings[’ecutwfc’]`.  

Linear displacements of single atoms can be done by invoking the script `AH_displace` from the terminal. See --help for more information on the options that are available. The default option for displacements is line, which displaces the atom in a line, the direction of which is decided from the lattice vectors (abc). If line-cart is used, the direction of the path is decided from the cartesian coordinates (xyz). The options mag and mag-cart allows you to specify the magnitude of the vector.  
For example, mag-cart 1 1 1 0.66 will use the vector $0.66 \cdot (1,1,1) / \sqrt{1+1+1}$.