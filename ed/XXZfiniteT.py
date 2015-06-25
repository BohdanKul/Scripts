from qutip import *
from scipy import *
from matplotlib import *
import numpy as np
import loadgmt,kevent
from pylab import *

N     = 8    # Number of spins in a chain
T     = 1    # Temperature
A     = 4 
Alist = range(A)
Alist = [0,1,4,5]
#delta = 0


figure = figure(1)
ax = subplot(111)
title(r'$L=%d, A=%d$' %(N,A)) 
connect('key_press_event',kevent.press)
colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]

Sz = 0.5*sigmaz()
Sx = 0.5*sigmax()
Sy = 0.5*sigmay()
Sp = sigmap()
Sm = sigmam()
I  = qeye(2)

# Define model's bond

for j,T in enumerate([0.01]):
    Es     = []
    deltas = []
    for delta in [0.5,1.0]:
        xxz_bond = tensor(Sx,Sx) + tensor(Sy,Sy) + delta*tensor(Sz,Sz)
        #print xxz_bond
        # Construct Hamiltonian matrix
        if N>2: H = tensor(xxz_bond,tensor([I]*(N-2)))
        else  : H = xxz_bond 
        #print H
        for i in range(1,N-2):
            pref = tensor([I] * i)
            suff = qeye(1)
            suff = tensor([I] * (N-i-2))
            H   += tensor(pref,xxz_bond,suff)
        if N>2: H += tensor(tensor([I]*(N-2)),xxz_bond)
        #print tensor(tensor([I]*(N-2)),xxz_bond)
        # Add periodic boundary conditions
        if N>2:
            bulk = tensor([I] * (N-2))
            H += tensor(Sx,bulk,Sx) 
            H += tensor(Sy,bulk,Sy) 
            H += tensor(Sz,bulk,delta*Sz) 
        else:
            H += tensor(Sx,Sx) 
            H += tensor(Sy,Sy) 
            H += tensor(Sz,delta*Sz) 

        #H = tensor(Sz,Sz)
        #print H
        rho  = (-H/T).expm()
        Z    = rho.tr()
        rho /= Z
        E    = (rho*H).tr()/float(N)
        if  (A==0): rhoA = rho
        else:       rhoA = rho.ptrace(Alist)
        rhoA /= rhoA.tr()
        ##print rhoA
        S2   = -log((rhoA*rhoA).tr())

        deltas += [delta]
        #Es     += [S2]
        Es     += [E]
        print 'Jz = %0.3f' %delta, ' N = %02d' %N,' A = %02d' %A,' T = %0.4f' %T, ' delta = %0.2f' %delta, ' Energy : %0.7f' %E, ' Entropy: %0.7f' %S2
xlabel(r'$\mathrm{\Delta}$')
#ylabel(r'$\mathrm{S_2}$')
ylabel(r'$\mathrm{E}$')
xlim([-2.05,2.05])
legend(frameon=False) 
show()



