from qutip import *
from scipy import *
from matplotlib import *
import numpy as np
import loadgmt,kevent
from pylab import *

N     = 8    # Number of spins in a chain
T     = 1    # Temperature
#delta = 0


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

for j,beta in enumerate([0.1, 1.0, 2.0, 10.0]):
    Es     = []
    hs = []
    #for h in np.linspace(0,1,10):
    for h in [0,0.25, 0.5, 0.75, 1]:
        ZZbond = -1.0*J*tensor(sigmaz(),sigmaz())
        hBond  = -1.0*h*sigmax()
        
        # Add one site interactions 
        Ht = tensor(hBond, tensor([I]*(N-1)))
        for i in range(1,N-1):
            pref = tensor([I]*i)
            suff = tensor([I]*(N-i-1))
            Ht   += tensor(pref, hBond, suff)
        Ht += tensor(tensor([I]*(N-1)), hBond)

        H = Ht
        # Add two sites bonds
        if N>2: H += tensor(ZZbond,tensor([I]*(N-2)))
        else  : H += ZZbond 
        
        for i in range(1,N-2):
            pref = tensor([I] * i)
            suff = tensor([I] * (N-i-2))
            H   += tensor(pref,ZZbond,suff)
        if N>2: H += tensor(tensor([I]*(N-2)),ZZbond)

        # Add periodic boundary conditions
        if N>2:
            bulk = tensor([I] * (N-2))
            H += tensor(sigmaz(),bulk,-1.0*J*sigmaz()) 
        else:
            H += tensor(sigmaz(),-1.0*J*sigmaz()) 

        rho  = (-H*beta).expm()
        Z    = rho.tr()
        rho /= Z
        E    = (rho*H).tr()/float(N)
        M    = (Ht*Ht*rho).tr()/float(N*N)

        hs += [h]
        #Es     += [np.sqrt(M)]
        Es     += [E]
        print ' N = %02d' %N,' beta = %0.4f' %beta, ' J = %0.2f' %J,' h = %0.2f' %h, ' Energy : %0.7f' %real(E), ' Magentization : %0.7f' %np.sqrt(real(M)) 
    if  beta == 0.1:
        ax.errorbar([0.00, 0.25, 0.50, 0.75, 1.00],
                    [-0.09871421, -0.10625461, -0.12401147, -0.15468257, -0.19793142],
                    [ 0.0006086,   0.00068187,  0.00074412,  0.00083378,  0.0009409 ],
                    mec = colors[j], ls = '', marker = 's', mfc = 'white')

    if  beta == 1.0:
        ax.errorbar([0.00, 0.25, 0.50, 0.75, 1.00],
                    [-0.81752932, -0.83560363, -0.89074228, -0.9884731,  -1.13415408],
                    [ 0.00030496,  0.00033013,  0.00031177,  0.00033864,  0.00039329],
                    mec = colors[j], ls = '', marker = 's', mfc = 'white')
    if  beta == 2.0:
        print Es
        ax.errorbar([0.00, 0.25, 0.50, 0.75, 1.00],
                    [-0.99524406, -1.00980205, -1.05466637, -1.13287314, -1.25568052],
                    [ 0.00013192,  0.00012351,  0.00018119,  0.00021636,  0.00022715],
                    mec = colors[j], ls = '',  marker = 's', mfc = 'white')
    if  beta == 10.0:
        ax.errorbar([0.00, 0.25, 0.50, 0.75, 1.00],
                    [-0.99997896, -1.01568132, -1.06360697, -1.14672478, -1.27850453],
                    [3.50260383e-05,   3.15993799e-05,   5.49630521e-05,   6.10981861e-05, 8.66366407e-05],
                    mec = colors[j], ls = '', marker = 's', mfc = 'white')
    ax.plot(hs, Es, 'b',ls='-', color = colors[j], marker = '', label=r'$\beta=%0.2f$' %beta)


xlabel(r'$\mathrm{h/J}$')
ylabel(r'$\mathrm{E/N}$')
#xlim([-2.05,2.05])
lg = legend(frameon=False) 
lg.draggable(state=True)
show()
