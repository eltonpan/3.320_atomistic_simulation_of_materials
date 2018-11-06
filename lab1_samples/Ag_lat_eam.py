from labutil.src.plugins.lammps import *
from ase.spacegroup import crystal
from ase.build import *
import numpy
import matplotlib.pyplot as plt


input_template = """
# ---------- 1. Initialize simulation ---------------------
units metal
atom_style atomic
dimension  3
boundary   p p p
read_data $DATAINPUT

# ---------- 2. Specify interatomic potential ---------------------
pair_style eam/alloy
pair_coeff * * $POTENTIAL  Ag

# pair_style lj/cut 4.5
# pair_coeff 1 1 0.3450 2.6244 4.5

# ---------- 3. Run single point calculation  ---------------------
thermo_style custom step pe lx ly lz press pxx pyy pzz
run 0

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
    #multiplier = numpy.identity(3) * 2
    #ase_supercell = make_supercell(unitcell, multiplier)
    structure = Struc(ase2struc(unitcell))
    return structure


def compute_energy(alat, template):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potpath = os.path.join(os.environ['LAMMPS_POTENTIALS'],'Ag_u3.eam')
    potential = ClassicalPotential(path=potpath, ptype='eam', element=["Al"])
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab1", str(alat)))
    struc = make_struc(alat=alat)
    output_file = lammps_run(struc=struc, runpath=runpath, potential=potential, intemplate=template, inparam={})
    energy, lattice = get_lammps_energy(outfile=output_file)
    return energy, lattice


def lattice_scan():
    alat_list = numpy.linspace(3.8, 4.3, 7)
    energy_list = [compute_energy(alat=a, template=input_template)[0] for a in alat_list]
    print(energy_list)
    plt.plot(alat_list, energy_list)
    plt.show()


if __name__ == '__main__':
    # put here the function that you actually want to run
    lattice_scan()
