import sys
import os
sys.path.append("/home/gridsan/{}/".format(os.environ['USER']))
from labutil.src.plugins.ising import *
import matplotlib.pyplot as plt
import pdb

# Variable initialization
def one_A(T, n_eq, n_mc):
    Ns = np.arange(5, 16)
    eq_idxs = []
    for N in Ns:
        #### RUN MC CODE ####
        E, M, E_eq, M_eq, imagelist = mc_run(N,n_eq,n_mc,T)

        x_eq = np.arange(n_eq*N*N)
        x_mc = np.arange(n_mc*N*N)+(n_eq*N*N)

        fig, (ax_energy, ax_mag) = plt.subplots(2,1, sharex=True, sharey=False)

        E_eq, M_eq = np.array(E_eq), np.array(M_eq)

        # pdb.set_trace()
        eq_idx = np.argwhere(np.abs(M_eq)>0.95)[0][0] 
        eq_idxs.append(eq_idx)
        print(N, eq_idx)

    # ax_energy.plot(x_eq,E_eq,label='Equilibration')
    # ax_energy.plot(x_mc,E,label='Production')
    # ax_energy.axvline(n_eq*(N**2),color = 'black')
    # ax_energy.legend()
    # ax_energy.set_ylabel('Energy per site/ J')

    # ax_mag.plot(x_eq,M_eq,label='Equilibration')
    # ax_mag.plot(x_mc,M,label='Production')
    # ax_mag.axvline(n_eq*(N**2),color = 'black')
    # ax_mag.legend()
    # ax_mag.set_ylabel('Magnetization per site')
    # ax_mag.set_xlabel('Number of flip attempts')

    # animator(imagelist)
    # plt.show()

    plt.figure(figsize = (10,7))
    plt.plot(Ns, eq_idxs)
    plt.xlabel('N')
    plt.ylabel('Num. steps needed for equilibration')
    plt.savefig('1A.png')



def one_B():

    N = 12
    n_eq = 750
    n_mc = 750

    avg_Ms = []
    Ts = np.arange(1.7, 3.5, 0.1)
    for T in Ts:
        #### RUN MC CODE ####
        E, M, E_eq, M_eq, imagelist = mc_run(N,n_eq,n_mc,T)

        x_eq = np.arange(n_eq*N*N)
        x_mc = np.arange(n_mc*N*N)+(n_eq*N*N)

        fig, (ax_energy, ax_mag) = plt.subplots(2,1, sharex=True, sharey=False)

        E_eq, M_eq = np.array(E_eq), np.array(M_eq)
        E, M = np.array(E), np.array(M)

        # pdb.set_trace()
        avg_M = np.mean(np.abs(M))
        avg_Ms.append(avg_M)
        print(T, avg_M)

    plt.figure(figsize = (7,5))
    plt.plot(Ts, avg_Ms)
    plt.xlabel('T')
    plt.ylabel('Average magnetization')
    plt.savefig('1B.png')

if __name__ == "__main__":
    # one_A(T = 1, n_eq = 750, n_mc = 750)
    one_B()

