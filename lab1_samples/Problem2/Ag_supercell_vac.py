# Problem 2
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append("/home/gridsan/{}/".format(os.environ['USER']))
from ase.spacegroup import crystal
from ase.build import make_supercell
from ase.io import write, read
from labutil.src.plugins.lammps import *
os.environ['LAMMPS_COMMAND']    = '/home/gridsan/{}/3320_atomistic_shared/lammps/src/lmp_serial'.format(os.environ['USER'])
os.environ['LAMMPS_POTENTIALS'] = '/home/gridsan/{}/3320_atomistic_shared/lammps/potentials'.format(os.environ['USER'])
os.environ['WORKDIR']           = os.getcwd()

# 2A

# Generate slabs of different sizes
unit_cell_dims = np.arange(1,15)
for bool in [True, False]:
    for unit_cell_dim in unit_cell_dims:
        # Build supercell
        unitcell = crystal('Ag', [(0, 0, 0)], spacegroup=225, cellpar=[4.090 , 4.090 ,4.090 , 90, 90, 90])
        multiplier = np.identity(3) * unit_cell_dim
        ase_supercell = make_supercell(unitcell, multiplier)

        # Remove 1 atom
        with_vac = bool
        
        if with_vac:
            ase_supercell.pop(0)
            write(f'Ag_slab_sc{unit_cell_dim}{unit_cell_dim}{unit_cell_dim}_with_vac.cif', ase_supercell)
        else:
            write(f'Ag_slab_sc{unit_cell_dim}{unit_cell_dim}{unit_cell_dim}.cif', ase_supercell)

# Compute cohesive energies
# Without optimization
input_template = """
# ---------- 1. Initialize simulation ---------------------
units metal
atom_style atomic
dimension  3
boundary   p p p
read_data $DATAINPUT

# ---------- 2. Specify interatomic potential ---------------------
pair_style eam
pair_coeff * * $POTENTIAL

# pair_style lj/cut 4.5
# pair_coeff 1 1 0.3450 2.6244 4.5

# ---------- 3. Run single point calculation  ---------------------
thermo_style custom step pe lx ly lz press pxx pyy pzz
run 0

# -- include optimization of the unit cell parameter
fix 1 all box/relax iso 0.0 vmax 0.001

# -- enable optimization of atomic positions (and the cell) 
min_style cg
minimize 1e-10 1e-10 10000 100000

# ---- 4. Define and print useful variables -------------
variable natoms equal "count(all)"
variable totenergy equal "pe"
variable length equal "lx"

print "Total energy (eV) = ${totenergy}"
print "Number of atoms = ${natoms}"
print "Lattice constant (Angstoms) = ${length}"
"""

def compute_energy_from_cif(fname, template, element, potpath):
    '''
    Computes energy from cif file input
    '''
    struc = read(fname) # read cif file
    struc = Struc(ase2struc(struc)) # important to convert from ASE to struc
    potpath = os.path.join(os.environ['LAMMPS_POTENTIALS'],potpath)
    potential = ClassicalPotential(path=potpath, ptype='eam', element=[element])
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab1", fname.split('.')[0]))
    output_file = lammps_run(struc=struc, runpath=runpath, potential=potential, intemplate=template, inparam={})
    energy, lattice = get_lammps_energy(outfile=output_file)

    n_atoms = struc.n_atoms
    return energy, lattice, n_atoms

if __name__ == '__main__':
    energy_list, size_list = [], []
    for d in unit_cell_dims:
        print(d)
        energy, _, n_atoms = compute_energy_from_cif(fname = f'Ag_slab_sc{d}{d}{d}.cif', template = input_template, element = 'Ag', potpath = 'Ag_u3.eam')
        E_wo_vac = (n_atoms-1)/n_atoms*energy
        print('E_cohesive: ', energy/n_atoms)

        energy, _, n_atoms = compute_energy_from_cif(fname = f'Ag_slab_sc{d}{d}{d}_with_vac.cif', template = input_template, element = 'Ag', potpath = 'Ag_u3.eam')
        E_w_vac = energy
        E_vac = E_w_vac - E_wo_vac

        energy_list.append(E_vac)
        size_list.append(d)

    plt.scatter(size_list, energy_list)
    plt.xlabel("Supercell size")
    plt.ylabel("Vacancy formation energy (eV)")
    plt.savefig('Ag_supercell_vac_w_opt_eam.png')
    plt.show()


