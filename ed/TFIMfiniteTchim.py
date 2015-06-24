from qutip import *
from scipy import *
from matplotlib import *
import numpy as np
import loadgmt,kevent
from pylab import *

x  = 1
y  = 2 
ux = 2
uy = 2
N = x*y*ux*uy

figure = figure(1)
ax = subplot(111)
connect('key_press_event',kevent.press)
colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]

#sigmaz = sigmaz()
#sigmax = sigmax()
I  = qeye(2)

J = 1.0
h = 1.0
# Define model's bond

#for j,beta in enumerate([0.1, 1.0, 2.0, 10.0]):
for j,beta in enumerate([1.0]):
    Es     = []
    hs = []
    #for h in np.linspace(0,1,10):
    for h in [0,0.25, 0.5, 0.75, 1]:
        hBond  = -1.0*h*sigmax()
        
        # Add one site interactions 
        Ht = tensor(hBond, tensor([I]*(N-1)))
        for i in range(1,N-1):
            pref = tensor([I]*i)
            suff = tensor([I]*(N-i-1))
            Ht   += tensor(pref, hBond, suff)
        Ht += tensor(tensor([I]*(N-1)), hBond)

        H = Ht
        # Add two sites bonds within all unit cells
        for l in range(y):
            for k in range(x):
                if l!=0   or k!=0:   bpref = tensor([I]*(ux*uy*(l+y*k)))
                if l!=y-1 or k!=x-1: bsuff = tensor([I]*(ux*uy*(x*y-l-y*k-1))) 
                
                uH = tensor([I]*ux*uy)*0
                for i in range(uy):
                    for j in range(uy):
                        if i    !=   0 : pref = tensor([I]*i)
                        if uy+j != i+1 : bulk = tensor([I]*(uy-i-1+j))
                        if j    != uy-1: suff = tensor([I]*(uy-j-1))
                        
                        if   uy==1:                inc = tensor(     sigmaz(),       -1.0*J*sigmaz())
                        elif (i==0) and (j==uy-1): inc = tensor(     sigmaz(), bulk, -1.0*J*sigmaz())
                        elif (i==0):               inc = tensor(     sigmaz(), bulk, -1.0*J*sigmaz(), suff)
                        elif (j==uy-1):            inc = tensor(pref,sigmaz(), bulk, -1.0*J*sigmaz())
                        elif (i==uy-1) and (j==0): inc = tensor(pref,sigmaz(),       -1.0*J*sigmaz(), suff)
                        else:                      inc = tensor(pref,sigmaz(), bulk, -1.0*J*sigmaz(), suff)
                        #print " Hamiltonian: ", H
                        #print " Increament: ", inc
                        uH += inc
                if   x==1   and y==1:   H+=uH
                elif l==0   and k==0:   H+=tensor(       uH, bsuff)  
                elif l==y-1 and k==x-1: H+=tensor(bpref, uH       )
                else:                   H+=tensor(bpref, uH, bsuff)

        # Add two sites vertical intra-cell bonds
        if  y>1:
            for l in range(y-1):
                for k in range(x):
                    for i in range(uy):
                        if not(l==0 and k==0 and i==0): pref = tensor([I]*(ux*uy*(l+y*k)+i))
                        bulk = tensor([I]*(ux*uy-1))
                        suff = tensor([I] *(ux*uy*(x*y-l-y*k-2)+ux*uy-i-1))
                       
                        if (l==0 and k==0 and i==0): H+=tensor(      sigmaz(), bulk, -1.0*J*sigmaz(), suff)
                        else:                        H+=tensor(pref, sigmaz(), bulk, -1.0*J*sigmaz(), suff)
       
        # Add two sites horizontal intra-cell bonds
        if  x>1:
            for l in range(y):
                for k in range(x-1):
                    for i in range(uy):
                        pref = tensor([I]*(ux*uy*(l+y*k)+uy+i))
                        bulk = tensor([I]*(ux*uy*y-1))
                        if not(l==y-1 and k==x-2 and i==uy-1): suff = tensor([I]*(ux*uy*(x*y-l-y*k-y-1)+uy-i-1))
                      
                        if (l==y-1 and k==x-2 and i==uy-1): H+=tensor(pref, sigmaz(), bulk, -1.0*J*sigmaz()      )
                        else:                               H+=tensor(pref, sigmaz(), bulk, -1.0*J*sigmaz(), suff)

        rho  = (-H*beta).expm()
        Z    = rho.tr()
        rho /= Z
        E    = (rho*H).tr()/float(N)
        M    = (Ht*Ht*rho).tr()/float(N*N)

        hs += [h]
        #Es     += [np.sqrt(M)]
        Es     += [E]
        print ' N = %02d' %N, ' x = %02d'%x, ' y = %02d'%y,' ux = %02d'%ux,' uy = %02d'%uy, ' beta = %0.4f' %beta, ' J = %0.2f' %J,' h = %0.2f' %h, ' Energy : %0.7f' %real(E), ' Magentization : %0.7f' %np.sqrt(real(M)) 
    if  beta == 10.0:
        ax.errorbar([0.00, 0.25, 0.50, 0.75, 1.00],
                    [-0.99997896, -1.01568132, -1.06360697, -1.14672478, -1.27850453],
                    [3.50260383e-05,   3.15993799e-05,   5.49630521e-05,   6.10981861e-05, 8.66366407e-05],
                    mec = colors[j], ls = '', marker = 's', mfc = 'white')
    ax.plot(hs, Es, 'b',ls='-', color = colors[j], marker = '', label=r'$\beta=%0.2f$' %beta)
    print Es

xlabel(r'$\mathrm{h/J}$')
ylabel(r'$\mathrm{E/N}$')
#xlim([-2.05,2.05])
lg = legend(frameon=False) 
lg.draggable(state=True)
show()
