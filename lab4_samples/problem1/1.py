import sys
import os
sys.path.append("/home/gridsan/{}/".format(os.environ['USER']))
from labutil.src.plugins.lammps import *
from ase.spacegroup import crystal
# from ase.build import *
from ase.build import make_supercell
import matplotlib.pyplot as plt
import numpy as np
import pdb

def make_struc(size):
    """
    Creates the crystal structure using ASE.
    :param size: supercell multiplier
    :return: structure object converted from ase
    """
    alat = 4.10
    unitcell = crystal('Ag', [(0, 0, 0)], spacegroup=225, cellpar=[alat, alat, alat, 90, 90, 90])
    multiplier = numpy.identity(3) * size
    supercell = make_supercell(unitcell, multiplier)
    structure = Struc(ase2struc(supercell))
    return structure


def compute_dynamics(size, timestep, nsteps, temperature):
    """
    Make an input template and select potential and structure, and input parameters.
    Return a pair of output file and RDF file written to the runpath directory.
    """
    intemplate = """
    # ---------- Initialize simulation ---------------------
    units metal
    atom_style atomic
    dimension  3
    boundary   p p p
    read_data $DATAINPUT

    pair_style eam
    pair_coeff * * $POTENTIAL

    velocity  all create $TEMPERATURE 87287 dist gaussian

    # ---------- Describe computed properties------------------
    compute msdall all msd
    thermo_style custom step pe ke etotal temp press density c_msdall[4]
    thermo $TOUTPUT

    # ---------- Specify ensemble  ---------------------
    fix  1 all nve
    #fix  1 all nvt temp $TEMPERATURE $TEMPERATURE $TDAMP

    # --------- Compute RDF ---------------
    compute rdfall all rdf 100 1 1
    fix 2 all ave/time 1 $RDFFRAME $RDFFRAME c_rdfall[*] file $RDFFILE mode vector

    # --------- Run -------------
    timestep $TIMESTEP
    run $NSTEPS
    """

    potential = ClassicalPotential(ptype='eam', element='Ag', name='Ag_u3.eam')
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab4/Problem1", "size_" + str(size)))
    struc = make_struc(size=size)
    inparam = {
        'TEMPERATURE': temperature,
        'NSTEPS': nsteps,
        'TIMESTEP': timestep,
        'TOUTPUT': 100,                 # how often to write thermo output
        'TDAMP': 50 * timestep,       # thermostat damping time scale
        'RDFFRAME': int(nsteps / 4),   # frames for radial distribution function
        'RDFFILE': "/home/gridsan/epan1/labutil/lab4_samples/lammps.rdf"
    }
    outfile = lammps_run(struc=struc, runpath=runpath, potential=potential,
                                  intemplate=intemplate, inparam=inparam)
    output = parse_lammps_thermo(outfile=outfile)
    rdffile = get_rdf(runpath=Dir(path = "/home/gridsan/epan1/labutil/lab4_samples"))
    rdfs = parse_lammps_rdf(rdffile=rdffile)
    return output, rdfs


def md_run(size=3, timestep=0.001, nsteps=10000, temperature=300):
    output, rdfs = compute_dynamics(size=size, timestep=timestep, nsteps=nsteps, temperature=temperature)
    [simtime, pe, ke, energy, temp, press, dens, msd] = output
    simtime = [int(s) for s in simtime]
    energy = [float(e) for e in energy]
    temp = [float(t) for t in temp]
    simtime = np.array(simtime)
    cutoff_idx = np.argwhere(simtime >= 2000.)[0][0] # idx of cutoff step to start averaging energy

    avg_e = np.mean(energy[cutoff_idx:])
    print('avg_e:', avg_e)
    ## ------- plot output properties
    #plt.plot(simtime, temp)
    #plt.show()
    plt.plot(simtime, press)
    plt.show()

    # plt.figure()
    # plt.plot(simtime, energy)
    # plt.xlabel('steps')
    # plt.ylabel('energy (eV)')
    # plt.savefig('tot_e_vs_t.png')
    return avg_e, simtime, temp


    # ----- plot radial distribution functions
    for rdf in rdfs:
        plt.plot(rdf[0], rdf[1])

    plt.show()

def one_A():
    timesteps = np.linspace(0.001, 0.02, 10)
    avg_e, _, _ = [md_run(timestep=t) for t in timesteps]

    plt.figure()
    plt.plot(timesteps, avg_es)
    plt.xlabel('Timestep size (ps)')
    plt.ylabel('Time-averaged energy (eV)')
    plt.savefig('1A.png')

def one_B():

    _, simtime, temp = md_run(timestep=0.001)

    plt.figure(figsize = (7,5))
    plt.plot(simtime, temp)
    plt.xlabel('Time (ps)')
    plt.ylabel('Temperature (K)')
    plt.savefig('1B.png')

def one_C():

    _, simtime, temp = md_run(timestep=0.001, size=8)

    plt.figure(figsize = (7,5))
    plt.plot(simtime, temp)
    plt.xlabel('Time (ps)')
    plt.ylabel('Temperature (K)')
    plt.savefig('1C.png')


if __name__ == '__main__':
    # put here the function that you actually want to run
    os.environ['LAMMPS_COMMAND']    = '/home/gridsan/{}/3320_atomistic_shared/lammps/src/lmp_serial'.format(os.environ['USER'])
    os.environ['LAMMPS_POTENTIALS'] = '/home/gridsan/{}/3320_atomistic_shared/lammps/potentials'.format(os.environ['USER'])
    os.environ['WORKDIR']           = os.getcwd()

    # one_A()
    # one_B()
    one_C()
    
