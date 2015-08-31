from qutip import *
from scipy import *
from scipy.integrate import simps
from matplotlib import *
import numpy as np
import loadgmt,kevent
from pylab import *





#figure = figure(1)
#ax = subplot(111)
#connect('key_press_event',kevent.press)
#colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]

N     = 4    # Number of spins in a chain
J     = 1.0
delta = 1.0

I  = qeye(2)
#P = tensor(projection(2,0,0), tensor([I]*(N-Nv)))
ket = tensor(basis(2,1), basis(2,0))
P2 =  ket*ket.dag()


clamped = []
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


PBC = True 
#PBC = False 
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
        
        rho  = (-H*beta).expm()
        rho *= P
        Z    = rho.tr()
        rho /= Z
        
        E = 0
        E = (rho*H).tr()/float(N)
        
        
        #M = (Ht*Ht*rho).tr()/float(N*N)
        M = 0
        hs += [h]
        #Es     += [np.sqrt(M)]
        print 'PBC = ',PBC, ' N = %02d' %N,' beta = %0.4f' %beta, ' J = %0.2f' %J,' h = %0.2f' %h, ' delta = %0.2f' %delta, ' Energy : %0.7f' %real(E), ' Magnetization : %0.7f' %np.sqrt(real(M))
        print "  sigma_z: ",
        for i in range(N):
            if   i==0:     oper = tensor(sigmaz(),tensor([qeye(2)]*(N-1)))
            elif i==(N-1): oper = tensor(tensor([qeye(2)]*(N-1)), sigmaz())
            else:          oper = tensor(tensor([qeye(2)]*i), sigmaz(), tensor([qeye(2)]*(N-1-i)))
            O    = (rho*oper).tr()
            print  ' %d = %+0.7f' %(i,real(O)),


        print 
        print "  sigma_x: ",
        for i in range(N):
            if   i==0:     oper = tensor(sigmax(),tensor([qeye(2)]*(N-1)))
            elif i==(N-1): oper = tensor(tensor([qeye(2)]*(N-1)), sigmax())
            else:          oper = tensor(tensor([qeye(2)]*i), sigmax(), tensor([qeye(2)]*(N-1-i)))
            O    = (rho*oper).tr()
            print  ' %d = %+0.7f' %(i,real(O)),
        print

        print "  sigma_z x sigma_z: ",
        if N>1:
            if N==2:
               oper = tensor(sigmaz(), sigmaz())
               O    = (rho*oper).tr()
               print  '0-%d = %+0.7f' %(i,real(O)),
            else:    
                for i in range(N-1):
                    if   i==0:      oper = tensor(sigmaz(),                sigmaz(), tensor([qeye(2)]*(N-2)))
                    elif i==(N-2):  oper = tensor(tensor([qeye(2)]*(N-2)), sigmaz(), sigmaz())
                    else:           oper = tensor(tensor([qeye(2)]*i),     sigmaz(), sigmaz(), tensor([qeye(2)]*(N-i-2)))
                    O    = (rho*oper).tr()
                    print  '%d-%d = %+0.7f' %(i,i+1,real(O)),
                if PBC:          
                   oper = tensor(sigmaz(), tensor([qeye(2)]*(N-2)), sigmaz())
                   O    = (rho*oper).tr()
                   print  '%d-%d = %+0.7f' %(0,N-1,real(O)),
        print

        print "--- Time-evolved ---"
        print "  sigma_x: ",
        for i in range(N):
            if   i==0:     oper = tensor(sigmax(),tensor([qeye(2)]*(N-1)))
            elif i==(N-1): oper = tensor(tensor([qeye(2)]*(N-1)), sigmax())
            else:          oper = tensor(tensor([qeye(2)]*i), sigmax(), tensor([qeye(2)]*(N-1-i)))
        
            Nsamples = 100
            xs = []
            ys = []
            for step in range(Nsamples+1):
                tau = 1.0*step*beta/(1.0*Nsamples)
                xs  += [tau]
                Ub   = ( H*tau).expm()
                Uf   = (-H*tau).expm()
                rho  = (-H*beta).expm()
                rho *= P
                Z    = rho.tr()
                rho /= Z
                ys  += [(Uf*oper*Ub*rho).tr()]
            O = simps(np.array(ys), np.array(xs))
            
            plot(real(xs), real(ys))
            print  ' %d = %+0.7f' %(i,real(O)),

        print

        print "  sigma_z: ",
        for i in range(N):
            if   i==0:     oper = tensor(sigmaz(),tensor([qeye(2)]*(N-1)))
            elif i==(N-1): oper = tensor(tensor([qeye(2)]*(N-1)), sigmaz())
            else:          oper = tensor(tensor([qeye(2)]*i), sigmaz(), tensor([qeye(2)]*(N-1-i)))
        
            Nsamples = 100
            xs = []
            ys = []
            for step in range(Nsamples+1):
                tau = 1.0*step*beta/(1.0*Nsamples)
                xs  += [tau]
                Ub   = ( H*tau).expm()
                Uf   = (-H*tau).expm()
                rho  = (-H*beta).expm()
                rho *= P
                Z    = rho.tr()
                rho /= Z
                ys  += [(Uf*oper*Ub*rho).tr()]
            O = simps(np.array(ys), np.array(xs))
            
            print  ' %d = %+0.7f' %(i,real(O)),

        print

        #print "  sigma_z x sigma_z: ",
        #for i in range(N-1+int(PBC)):
        #    if N==2: oper = tensor(sigmaz(), sigmaz())
        #    else:  
        #        if i < N-1:
        #            if   i==0:      oper = tensor(sigmaz(),                sigmaz(), tensor([qeye(2)]*(N-2)))
        #            elif i==(N-2):  oper = tensor(tensor([qeye(2)]*(N-2)), sigmaz(), sigmaz())
        #            else:           oper = tensor(tensor([qeye(2)]*i),     sigmaz(), sigmaz(), tensor([qeye(2)]*(N-i-2)))
        #        else:        
        #           oper = tensor(sigmaz(), tensor([qeye(2)]*(N-2)), sigmaz())
        #
        #    Nsamples = 100
        #    xs = []
        #    ys = []
        #    for step in range(Nsamples+1):
        #        tau = 1.0*step*beta/(1.0*Nsamples)
        #        xs  += [tau]
        #        Ub   = ( H*tau).expm()
        #        Uf   = (-H*tau).expm()
        #        rho  = (-H*beta).expm()
        #        rho *= P
        #        Z    = rho.tr()
        #        rho /= Z
        #        ys  += [(Uf*oper*Ub*rho).tr()]
        #    O = simps(np.array(ys), np.array(xs))
        #    
        #    #plot(real(xs), real(ys))
        #    if i<N-1: print  '%d-%d = %+0.7f' %(i  , i+1, real(O)),
        #    else:     print  '%d-0 = %+0.7f'  %(N-1, real(O)),

        #print
        print

           
#xlabel(r'$\mathrm{h/J}$')
#ylabel(r'$\mathrm{E/N}$')
#xlim([-2.05,2.05])
#lg = legend(frameon=False) 
#if lg: lg.draggable(state=True)
show()
