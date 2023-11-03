import sys
import os
sys.path.append("/home/gridsan/{}/".format(os.environ['USER']))
from labutil.src.plugins.pwscf import *
from ase.io import write
from ase import Atoms
import matplotlib.pyplot as plt



def make_struc(alat, displacement=0):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    lattice = alat * numpy.identity(3)
    symbols = ['Pb', 'Ti', 'O', 'O', 'O']
    sc_pos = [[0,0,0], [0.5,0.5,0.5 + displacement], [0,0.5,0.5], [0.5,0,0.5], [0.5,0.5,0]]
    perov = Atoms(symbols=symbols, scaled_positions=sc_pos, cell=lattice)
    # check how your cell looks like
    # write('s.cif', perov)
    structure = Struc(ase2struc(perov))
    return structure



def compute_energy(alat, nk, ecut, displ=0):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    # pseudopots = {'Pb': PseudoPotential(ptype='uspp', element='Pb', functional='LDA', name='Pb.pz-d-van.UPF'),
    #               'Ti': PseudoPotential(ptype='uspp', element='Ti', functional='LDA', name='Ti.pz-sp-van_ak.UPF'),
    #               'O': PseudoPotential(ptype='uspp', element='O', functional='LDA', name='O.pz-rrkjus.UPF')}
    pseudopots = {'Pb': PseudoPotential(path=os.path.join(os.environ['ESPRESSO_PSEUDO'], 'Pb.pz-d-van.UPF'),     ptype='uspp', element='Pb', functional='LDA', name='Pb.pz-d-van.UPF'),
                  'Ti': PseudoPotential(path=os.path.join(os.environ['ESPRESSO_PSEUDO'], 'Ti.pz-sp-van_ak.UPF'), ptype='uspp', element='Ti', functional='LDA', name='Ti.pz-sp-van_ak.UPF'),
                   'O': PseudoPotential(path=os.path.join(os.environ['ESPRESSO_PSEUDO'], 'O.pz-rrkjus.UPF'),     ptype='uspp', element='O',  functional='LDA', name='O.pz-rrkjus.UPF')}
    struc = make_struc(alat=alat, displacement=displ)
    # fix the Pb and Ti atoms in place during relaxation
    constraint = Constraint(atoms={'0': [0,0,0], '1': [0,0,0]})
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=True)
    dirname = 'PbTiO3_a_{}_ecut_{}_nk_{}_displ_{}'.format(alat, ecut, nk, displ)
    # runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab3/Problem2", dirname))
    runpath = Dir(path=os.path.join(os.getcwd(), "Lab3", dirname))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'relax',
            'pseudo_dir': '/home/gridsan/{}/labutil/lab3_samples/pseudopotentials/'.format(os.environ['USER']),
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none'
        },
        'SYSTEM': {
            'ecutwfc': ecut,
            'ecutrho': ecut * 8,
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.7,
            'conv_thr': 1e-7,
        },
        'IONS': {
            'ion_dynamics': 'bfgs'
        },
        'CELL': {},
        })

    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts, constraint=constraint, ncpu=2)
    output = parse_qe_pwscf_output(outfile=output_file)
    return output



def lattice_scan():
    nk = 4
    ecut = 30
    alat = 3.88
    displ_list = numpy.linspace(-0.1, 0.1, 20)
    print(displ_list)
    energy_list = []
    for displ in displ_list:
        output = compute_energy(alat=alat, ecut=ecut, nk=nk, displ = displ)
        energy_list.append(output['energy']/13.6057039763)
        print(output)
    print(displ_list)
    print(energy_list)
    plt.plot(displ_list, energy_list)
    plt.xlabel('Fraction displacement')
    plt.ylabel('Energy (eV)')
    plt.savefig('2B.png', dpi = 100)
    plt.show()


if __name__ == '__main__':
    os.environ["PWSCF_COMMAND"] = '/home/gridsan/epan1/3320_atomistic_shared/qe-7.2/bin/pw.x'
    # put here the function that you actually want to run
    lattice_scan()
