import sys
import os
sys.path.append("/home/gridsan/{}/".format(os.environ['USER']))
from labutil.src.plugins.pwscf import *
from ase.spacegroup import crystal
from ase.io import write
# from ase.build import *
from ase.build import bulk
import numpy
import matplotlib.pyplot as plt



def make_struc(alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    fecell = bulk('Au', 'fcc', a=alat)
    # check how your cell looks like
    #write('s.cif', gecell)
    print(fecell, fecell.get_atomic_numbers())
    # fecell.set_atomic_numbers([26, 27])
    structure = Struc(ase2struc(fecell))
    print(structure.species)
    return structure


def compute_energy(alat, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potname = 'Au.pz-d-rrkjus.UPF'
    potpath = os.path.join(os.environ['ESPRESSO_PSEUDO'], potname)
    pseudopots = {'Au': PseudoPotential(path=potpath, ptype='uspp', element='Au',
                                        functional='LDA', name=potname),
                #   'Co': PseudoPotential(path=potpath, ptype='uspp', element='Fe',
                #                         functional='GGA', name=potname)
                  }
    struc = make_struc(alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    dirname = 'Au_a_{}_ecut_{}_nk_{}'.format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.getcwd(), "Lab3", dirname))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'vc-relax',
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
                               params=input_params, kpoints=kpts, ncpu=8)
    output = parse_qe_pwscf_output(outfile=output_file)
    return output



def lattice_scan():
    ecut = 40
    alat = 4.0
    
    nk_list = numpy.arange(1,16)
    energy_list = []
    nk_unique_list = []
    for nk in nk_list:
        output = compute_energy(alat=alat, ecut=ecut, nk=nk)
        e, nk_unique = output['energy'], output['n_unique_k_points']
        energy_list.append(e)
        nk_unique_list.append(nk_unique)
    print(output)

    plt.figure()
    plt.scatter(nk_unique_list, energy_list)
    plt.xlabel('Num. k-points')
    plt.ylabel('Energy (eV)')
    plt.savefig('3A_Au.png', dpi = 100)
    plt.show()

if __name__ == '__main__':
    os.environ["PWSCF_COMMAND"] = '/home/gridsan/epan1/3320_atomistic_shared/qe-7.2/bin/pw.x'
    # put here the function that you actually want to run
    lattice_scan()
