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

    # Displace a Ge atom +0.05 in the z direction (fractional coordinates)
    frac_coords = gecell.get_scaled_positions()
    frac_coords[-1][-1] += 0.05 # add to the z fractional coordinate of the 2nd atom
    gecell.set_scaled_positions(frac_coords)
    
    structure = Struc(ase2struc(gecell))
    return structure


def compute_force(alat, nk, ecut):
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


def lattice_scan(nk = 4, ecut = 30, alat = 5.0):
    output = compute_force(alat=alat, ecut=ecut, nk=nk)
    print(output)
    force = output['force']
    # *25.71104309541616 # convert from Ryd/Bohr to eV/Angstrom
    
    return force


if __name__ == '__main__':
    os.environ['WORKDIR']       = os.getcwd()
    os.environ['PWSCF_COMMAND'] = "~/3320_atomistic_shared/qe-7.2/bin/pw.x"
    # put here the function that you actually want to run

    forces, times = [], []
    ecuts = [int(x) for x in np.arange(5,80,5)]
    for ecut in ecuts:
        start = time.time()
        f = lattice_scan(ecut = ecut)
        end = time.time()

        print(ecut, f)
        forces.append(f)
        times.append(end-start)
    
    print('forces', forces)
    print('times', times)

    plt.figure()
    plt.scatter(ecuts, forces)
    plt.xlabel('Plane-wave cutoff energy (Ry)')
    plt.ylabel('Force (eV/Angstrom)')
    plt.savefig('3.png', dpi = 100)
    plt.show()

    # plt.figure()
    # plt.scatter(ecuts, times)
    # plt.xlabel('Plane-wave cutoff energy (Ry)')
    # plt.ylabel('Wall time (s)')
    # plt.savefig('1B.png', dpi = 100)
    # plt.show()

    
