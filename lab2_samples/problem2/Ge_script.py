import sys
import os
sys.path.append("/home/gridsan/{}/".format(os.environ['USER']))
from labutil.src.plugins.pwscf import *
from ase.spacegroup import crystal
from ase.io import write
import matplotlib.pyplot as plt
import numpy as np
import time


def make_struc(alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    # set primitive_cell=False if you want to create a simple cubic unit cell with 8 atoms
    gecell = crystal('Ge', [(0, 0, 0)], spacegroup=227, cellpar=[alat, alat, alat, 90, 90, 90], primitive_cell=True)
    # check how your cell looks like
    # write('s.cif', gecell)
    structure = Struc(ase2struc(gecell))
    return structure


def compute_energy(alat, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    # potname = 'Ge.pz-bhs.UPF'
    # pseudopath = os.environ['ESPRESSO_PSEUDO']

    potname = 'Ge.pz-dn-rrkjus_psl.0.2.2.UPF'
    pseudopath = '../'

    potpath = os.path.join(pseudopath, potname)

    pseudopots = {'Ge': PseudoPotential(name=potname, path=potpath, ptype='uspp', element='Ge', functional='LDA')}
    struc = make_struc(alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab2/Problem1", str(alat)))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'scf',
            'pseudo_dir': pseudopath,
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none',
        },
        'SYSTEM': {
            'ecutwfc': ecut,
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.5,
            'conv_thr': 1e-7,
        },
        'IONS': {},
        'CELL': {},
        })

    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts)
    output = parse_qe_pwscf_output(outfile=output_file)
    return output


def lattice_scan(nk = 3, ecut = 30, alat = 5.0):
    output = compute_energy(alat=alat, ecut=ecut, nk=nk)
    print(output)
    energy = output['energy']/13.6057039763/2 # convert from Rydberg to eV/atom (2 atoms per primitive cell)
    n_unique_k_points = output['n_unique_k_points']

    return energy, n_unique_k_points


if __name__ == '__main__':
    os.environ['WORKDIR']       = os.getcwd()
    os.environ['PWSCF_COMMAND'] = "~/3320_atomistic_shared/qe-7.2/bin/pw.x"
    # put here the function that you actually want to run

    energies, times = [], []
    nks = [int(x) for x in np.arange(1,10)]
    nks_unique = []
    for nk in nks:
        start = time.time()
        e, nk_unique = lattice_scan(nk = nk)
        print(nk, e)
        end = time.time()

        energies.append(e)
        times.append(end-start)
        nks_unique.append(nk_unique)
    
    print('energies', energies)
    print('times', times)

    nks = np.array(nks)**3

    plt.figure()
    plt.scatter(nks_unique, energies)
    plt.xlabel('Number of k-points')
    plt.ylabel('Energy (eV)')
    plt.savefig('2A.png', dpi = 100)
    plt.show()

    plt.figure()
    plt.scatter(nks_unique, times)
    plt.xlabel('Number of k-points')
    plt.ylabel('Wall time (s)')
    plt.savefig('2B.png', dpi = 100)
    plt.show()

    
