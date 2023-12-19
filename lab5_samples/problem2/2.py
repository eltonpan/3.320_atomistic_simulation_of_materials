import sys
import os
sys.path.append("/home/gridsan/{}/".format(os.environ['USER']))
from labutil.src.plugins.ising import *
import matplotlib.pyplot as plt
import pdb

# Variable initialization
def two_A():

    N = 12
    n_eq = 750
    n_mc = 2000

    avg_Es = []
    Ts = np.arange(1.7, 3.5, 0.1)
    for T in Ts:
        #### RUN MC CODE ####
        E, M, E_eq, M_eq, imagelist = mc_run(N,n_eq,n_mc,T)

        x_eq = np.arange(n_eq*N*N)
        x_mc = np.arange(n_mc*N*N)+(n_eq*N*N)

        E_eq, M_eq = np.array(E_eq), np.array(M_eq)
        E, M = np.array(E), np.array(M)

        # pdb.set_trace()
        avg_E = np.mean(E)
        avg_Es.append(avg_E)
        print(T, avg_E)


    # pdb.set_trace()
    grads = []
    for idx in range(len(Ts)-1):
        grad = (avg_Es[idx+1]-avg_Es[idx])/(Ts[idx+1]-Ts[idx])
        grads.append(grad)

    plt.figure(figsize = (7,5))
    plt.plot(Ts[:-1], grads)
    plt.xlabel('T')
    plt.ylabel('Cv')
    plt.savefig('2A.png')



if __name__ == "__main__":
    two_A()

