import argparse
from qutip import *
from scipy import *
from matplotlib import *
import numpy as np
import loadgmt,kevent
from pylab import *


def LoadInteractionsAdv(fName):
    inFile = open(fName, 'r')
    data   = loadtxt(fName)
    (nSz, nSx, nSzSz) = tuple(data[0,:])
    (nSz, nSx, nSzSz) = (int(nSz), int(nSx), int(nSzSz))
    Sz = []
    for i in range(1,1+nSz):
        Sz.append([int(data[i,0]),data[i,2]])
    Sx = []
    for i in range(1+nSz,1+nSz+nSx):
        Sx.append([int(data[i,0]),data[i,2]])
    SzSz = []
    for i in range(1+nSz+nSx,1+nSz+nSx+nSzSz):
        SzSz.append([int(data[i,0]),int(data[i,1]),data[i,2]])

    return Sz, Sx, SzSz

def LoadInteractions(fName):
    inFile = open(fName, 'r')
    data   = loadtxt(fName)
    (nSz, nSx, nSzSz) = tuple(data[0,:])
    (nSz, nSx, nSzSz) = (int(nSz), int(nSx), int(nSzSz))
    Sz = []
    for i in range(1,1+nSz):
        Sz.append(data[i,2])
    Sx = []
    for i in range(1+nSz,1+nSz+nSx):
        Sx.append(data[i,2])
    SzSz = []
    for i in range(1+nSz+nSx,1+nSz+nSx+nSzSz):
        SzSz.append(data[i,2])

    return Sz, Sx, SzSz

def main(): 
    parser = argparse.ArgumentParser(description='Finite temperature exact diagonalization of TFIM')
    parser.add_argument('--interaction','-i', help='File with interactions', type=str)
    args = parser.parse_args()
    fName = args.interaction

    Sz, Sx, SzSz = LoadInteractions(fName)
    print 'Sz: ', Sz
    print 'Sx: ', Sx
    print 'SzSz: ', SzSz
    N     = len(Sx) 


    # figure = figure(1)
    #ax = subplot(111)
    #connect('key_press_event',kevent.press)
    #colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]

    I  = qeye(2)

    for j,beta in enumerate([0.1, 1.0, 2.0, 10.0]):
        Es     = []
        # Add one site interactions 
        Ht = tensor(Sx[0]*sigmax(), tensor([I]*(N-1)))
        for i in range(1,N-1):
            pref = tensor([I]*i)
            suff = tensor([I]*(N-i-1))
            Ht   += tensor(pref, Sx[i]*sigmax(), suff)
        Ht += tensor(tensor([I]*(N-1)), Sx[N-1]*sigmax())

        H = Ht
        # Add two sites bonds
        if N>2: H += tensor(SzSz[0]*tensor(sigmaz(),sigmaz()) ,tensor([I]*(N-2)))
        else  : H += tensor(SzSz[0]*sigmaz(),sigmaz())
 
        
        for i in range(1,N-2):
            pref = tensor([I] * i)
            suff = tensor([I] * (N-i-2))
            H   += tensor(pref, SzSz[i]*tensor(sigmaz(),sigmaz()) ,suff)
        if N>2: H += tensor(tensor([I]*(N-2)),SzSz[N-2]*tensor(sigmaz(),sigmaz()))

        # Add periodic boundary conditions
        if N>2:
            bulk = tensor([I] * (N-2))
            H += tensor(sigmaz(),bulk,SzSz[N-1]*sigmaz()) 
        else:
            H += tensor(sigmaz(),SzSz[N-1]*sigmaz()) 

        rho  = (-H*beta).expm()
        Z    = rho.tr()
        rho /= Z
        E    = (rho*H).tr()/float(N)
        Sxi   = tensor(sigmax(), tensor([I]*(N-1)))
        M    = (Sxi*rho).tr()#/float(N*N)

        #Es     += [np.sqrt(M)]
        Es     += [E]
        print ' N = %02d' %N,' beta = %0.4f' %beta, ' Energy : %0.7f' %real(E), ' Magentization : %0.7f' %real(M) 

    #xlabel(r'$\mathrm{h/J}$')
    #ylabel(r'$\mathrm{E/N}$')
    #xlim([-2.05,2.05])
    #lg = legend(frameon=False) 
    # lg.draggable(state=True)
    #show()


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
