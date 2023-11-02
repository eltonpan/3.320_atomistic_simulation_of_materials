import sys
import os
sys.path.append("/home/gridsan/{}/".format(os.environ['USER']))
from labutil.src.plugins.pwscf import *
from ase.spacegroup import crystal
from ase.io import write
# from ase.build import *
from ase.build import bulk
import numpy as np
import matplotlib.pyplot as plt



def make_struc(alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    fecell = bulk('Fe', 'bcc', a=alat)
    # check how your cell looks like
    #write('s.cif', gecell)
    print(fecell, fecell.get_atomic_numbers())
    structure = Struc(ase2struc(fecell))
    print(structure.species)
    return structure


def compute_energy(alat, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potname = 'Fe.pbe-nd-rrkjus.UPF'
    potpath = os.path.join(os.environ['ESPRESSO_PSEUDO'], potname)
    pseudopots = {'Fe': PseudoPotential(path=potpath, ptype='uspp', element='Fe',
                                        functional='GGA', name=potname),
                #   'Co': PseudoPotential(path=potpath, ptype='uspp', element='Fe',
                #                         functional='GGA', name=potname)
                  }
    struc = make_struc(alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    dirname = 'Fe_bcc_a_{}_ecut_{}_nk_{}'.format(alat, ecut, nk)
    # runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem1", dirname))
    runpath = Dir(path=os.path.join(os.getcwd(), "Lab3", dirname))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'scf',
            'pseudo_dir': '/home/gridsan/{}/labutil/lab3_samples/pseudopotentials/'.format(os.environ['USER']),
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none',
        },
        'SYSTEM': {
            'ecutwfc': ecut,
            'ecutrho': ecut * 8,
            'nspin': 2,
            'starting_magnetization(1)': 0.7,
            'occupations': 'smearing',
            'smearing': 'mp',
            'degauss': 0.02
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
                               params=input_params, kpoints=kpts, ncpu=4)
    output = parse_qe_pwscf_output(outfile=output_file, return_vol=True)
    return output



def lattice_scan():
    ecut = 30
    nk = 10

    alat_list = [
                2.2, 
                2.5, 
                2.8, 
                3.1, 
                3.4, 
                3.7, 
                4.0
                ]
    energies = []
    volumes =  []
    for alat in alat_list:
        output = compute_energy(alat=alat, ecut=ecut, nk=nk)
        print(output)
        e, f, p, n_unique_k_points, v = output['energy']/13.6057039763, output['force'], output['pressure'], output['n_unique_k_points'], output['vol'] # no dividing since bcc primitive cell only has 1 atom
        energies.append(e)
        volumes.append(v)
    
    print('energies:', energies)
    print('volumes:', volumes)
    
    plt.figure()
    plt.scatter(volumes, energies)
    plt.xlabel('Volume (a.u.^3)')
    plt.ylabel('Energy per atom (eV)')
    plt.savefig('1B_bcc.png', dpi = 100)
    plt.show()

    

if __name__ == '__main__':
    os.environ["PWSCF_COMMAND"] = '/home/gridsan/epan1/3320_atomistic_shared/qe-7.2/bin/pw.x'
    # put here the function that you actually want to run
    lattice_scan()
