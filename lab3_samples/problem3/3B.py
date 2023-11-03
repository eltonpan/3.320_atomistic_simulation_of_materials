import sys
import os
sys.path.append("/home/gridsan/{}/".format(os.environ['USER']))
from labutil.src.plugins.pwscf import *
from ase.spacegroup import crystal
from ase.io import write
# from ase.build import *
from ase.build import bulk
from ase import Atoms
import numpy
import matplotlib.pyplot as plt


def make_struc(a, c):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    lattice = [a,a,c] * numpy.identity(3)
    print(lattice)
    symbols = ['Cu', 'Au']
    sc_pos = [[0,0,0], [0.5,0.5,0.5]]
    cuaucell = Atoms(symbols=symbols, scaled_positions=sc_pos, cell=lattice)
    write('CuAu.cif', cuaucell)
    # fecell.set_atomic_numbers([26, 27])
    structure = Struc(ase2struc(cuaucell))
    print(structure.species)
    return structure


if __name__ == '__main__':
    make_struc(a = 3.5, c = 4.0)
