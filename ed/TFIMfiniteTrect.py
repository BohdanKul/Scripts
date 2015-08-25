from qutip import *
from scipy import *
from matplotlib import *
import numpy as np
import loadgmt,kevent
from pylab import *





figure = figure(1)
ax = subplot(111)
connect('key_press_event',kevent.press)
colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]

N     = 2    # Number of spins in a chain
J     = -0.0
delta = 0.0

I  = qeye(2)
#P = tensor(projection(2,0,0), tensor([I]*(N-Nv)))
ket = tensor(basis(2,1), basis(2,0))
P2 =  ket*ket.dag()


clamped = [0,1]
if  len(clamped)>0:
    ket = basis(2,clamped[0])
    for qbit in clamped[1:]:
        ket = tensor(ket, basis(2,qbit))
    P = ket*ket.dag()

if N > len(clamped):
   if len(clamped)>0: P = tensor(P, tensor([qeye(2)]*(N-len(clamped))))
   else:              P = tensor(tensor([qeye(2)]*N))
else:
    print "Length of clamped bits exceed the total system size"


oper = tensor(sigmaz(),tensor([qeye(2)]*(N-1)))
#print P
print oper
PBC = True 
PBC = False 
for j,beta in enumerate([1.0]):
    Es     = []
    hs = []
    #for h in np.linspace(0,1,10):
    for h in [1.0]:
        ZZbond = 1.0*J*tensor(sigmaz(),sigmaz())
        dBond  = 1.0*delta*sigmax()
        hBond  = 1.0*h*sigmaz()
        
        # Add one site interactions 
        if delta!=0:
            Ht = tensor(dBond, tensor([I]*(N-1)))
            for i in range(1,N-1):
                pref = tensor([I]*i)
                suff = tensor([I]*(N-i-1))
                Ht   += tensor(pref, dBond, suff)
            Ht += tensor(tensor([I]*(N-1)), dBond)
            H = Ht
        else:
            H = tensor([I]*(N)) - tensor([I]*(N))
        # Add one site interactions 
        if h!=0.0:
            Ht = tensor(hBond, tensor([I]*(N-1)))
            for i in range(1,N-1):
                pref = tensor([I]*i)
                suff = tensor([I]*(N-i-1))
                Ht   += tensor(pref, hBond, suff)
            Ht += tensor(tensor([I]*(N-1)), hBond)
            H += Ht

        # Add two sites bonds
        if J!=0:
            if N>2: H += tensor(ZZbond,tensor([I]*(N-2)))
            else  : H += ZZbond 
            
            for i in range(1,N-2):
                pref = tensor([I] * i)
                suff = tensor([I] * (N-i-2))
                H   += tensor(pref,ZZbond,suff)
            if N>2: H += tensor(tensor([I]*(N-2)),ZZbond)
        # Add periodic boundary conditions
        if PBC and J!=0:
            if N>2:
                bulk = tensor([I] * (N-2))
                H += tensor(sigmaz(),bulk,J*sigmaz()) 
            else:
                H += tensor(sigmaz(),J*sigmaz()) 
#        print H

        rho  = (-H*beta).expm()
        Z    = rho.tr()
        rho /= Z
#        print rho
        #E = (rho*H*P).tr()/float(N)
        E = 0
        O    = (rho*oper*P).tr()/float(N)
        #M = (Ht*Ht*rho).tr()/float(N*N)
        M = 0
        hs += [h]
        #Es     += [np.sqrt(M)]
        Es     += [E]
        print 'PBC = ',PBC, ' N = %02d' %N,' beta = %0.4f' %beta, ' J = %0.2f' %J,' h = %0.2f' %h, ' delta = %0.2f' %delta, ' Energy : %0.7f' %real(E), ' Magentization : %0.7f' %np.sqrt(real(M)), ' Operator: %0.7f' % O

xlabel(r'$\mathrm{h/J}$')
ylabel(r'$\mathrm{E/N}$')
#xlim([-2.05,2.05])
#lg = legend(frameon=False) 
#if lg: lg.draggable(state=True)
show()
