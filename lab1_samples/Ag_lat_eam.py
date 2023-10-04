import sys
import os
sys.path.append("/home/gridsan/{}/".format(os.environ['USER']))
from labutil.src.plugins.lammps import *
from ase.build import make_supercell
from ase.spacegroup import crystal
import numpy
import matplotlib.pyplot as plt

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
# fix 1 all box/relax iso 0.0 vmax 0.001

# -- enable optimization of atomic positions (and the cell) 
# min_style cg
# minimize 1e-10 1e-10 1000 10000

# ---- 4. Define and print useful variables -------------
variable natoms equal "count(all)"
variable totenergy equal "pe"
variable length equal "lx"

print "Total energy (eV) = ${totenergy}"
print "Number of atoms = ${natoms}"
print "Lattice constant (Angstoms) = ${length}"
       """

# With optimization
input_template_for_opt = """
# ---------- 1. Initialize simulation --------------------- 
units metal
atom_style atomic
dimension 3
boundary   p p p
read_data $DATAINPUT

# ---------- 2. Specify interatomic potential ------------- 
# pair_style lj/cut 4.5
# pair_coeff 1 1 0.3450 2.6244 4.5
pair_style eam
pair_coeff * * $POTENTIAL

# ---------- 3. Run the calculation ----------------
# -- perform a single-point energy calculation only 
run 0

# # -- include optimization of the unit cell parameter
fix 1 all box/relax iso 0.0 vmax 0.001

# # -- enable optimization of atomic positions (and the cell) 
min_style cg
minimize 1e-10 1e-10 10000 100000 # need to increase maxiter and maxeval for optimizer to converge

# ---- 4. Define and print useful variables ------------- 
variable natoms equal "count(all)"
variable totenergy equal "pe"
variable length equal "lx"
print "Total energy (eV) = ${totenergy}"
print "Number of atoms = ${natoms}"
print "Lattice constant (Angstoms) = ${length}"
"""

def make_struc(alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    unitcell = crystal('Ag', [(0, 0, 0)], spacegroup=225, cellpar=[alat, alat, alat, 90, 90, 90])
    
    # For supercell
    multiplier = numpy.identity(3) * 2
    ase_supercell = make_supercell(unitcell, multiplier)
    structure = Struc(ase2struc(ase_supercell))
    
    # # For unit cell
    # structure = Struc(ase2struc(unitcell))
    return structure


def compute_energy(alat, template):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potpath = os.path.join(os.environ['LAMMPS_POTENTIALS'],'Ag_u3.eam')
    potential = ClassicalPotential(path=potpath, ptype='eam', element=["Ag"])
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab1", str(alat)))
    struc = make_struc(alat=alat)
    output_file = lammps_run(struc=struc, runpath=runpath, potential=potential, intemplate=template, inparam={})
    energy, lattice = get_lammps_energy(outfile=output_file)
    return energy, lattice


def lattice_scan():
    alat_list = numpy.linspace(3.8, 4.4, 10)
    energy_list = [compute_energy(alat=a, template=input_template)[0] for a in alat_list]
    print('---Without optimization---')
    print('lattice params:', alat_list)
    print('energies:', energy_list)
    plt.scatter(alat_list, energy_list)
    plt.xlabel('Lattice constant ($\AA$)')
    plt.ylabel('Energy (eV)')
    plt.savefig('Ag_lat_eam.png', dpi = 100)
    plt.show()


def lattice_scan_with_optimization():
    alat_list = numpy.linspace(3.8, 4.4, 10)
    results = [compute_energy(alat=a, template=input_template_for_opt) for a in alat_list]
    optimized_energy_list = [a[0] for a in results] # optimized energies
    optimized_lattice_list = [a[1] for a in results] # optimized lattice params
    print('---With optimization---')
    print('lattice params:', alat_list)
    print('optimized energies:', optimized_energy_list)
    print('optimized lattice:', optimized_lattice_list)


if __name__ == '__main__':
    # put here the function that you actually want to run
    os.environ['LAMMPS_COMMAND']    = '/home/gridsan/{}/3320_atomistic_shared/lammps/src/lmp_serial'.format(os.environ['USER'])
    os.environ['LAMMPS_POTENTIALS'] = '/home/gridsan/{}/3320_atomistic_shared/lammps/potentials'.format(os.environ['USER'])
    os.environ['WORKDIR']           = os.getcwd()
    lattice_scan() # without optimization
    print()
    lattice_scan_with_optimization() # with optimization
