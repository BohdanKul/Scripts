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
        #print 'Jz = %0.3f' %delta,  ' Energy = %0.7f' %E, ' Entropy = %0.7f' %S2
    #ax.plot(deltas, Es, 'b',ls='-', color = colors[j], marker = '.', label=r'$ED\,T=%0.2f$' %T)
    #deltas = np.arange(1.01,2.05,0.01)
    #ax.plot(deltas,np.sqrt(0.5/np.pi-0.5*np.arccos(deltas)/(np.pi**2)), 'y',ls='-', marker = '.')

xlabel(r'$\mathrm{\Delta}$')
#ylabel(r'$\mathrm{S_2}$')
ylabel(r'$\mathrm{E}$')
xlim([-2.05,2.05])
legend(frameon=False) 
show()


#ax.errorbar(np.arange(0,1.2,0.1), [-0.12388921991361787, -0.12821130350306337, -0.13389134934374727, -0.14101739987268161, -0.14889584691971477, -0.15710154978034604, -0.16669939975951903, -0.17713359968615464, -0.18907219961716926, -0.20236244953773375,-0.21616369957117534, -0.23044809334836355], [0.00013667796255575787, 0.00034572820985246701, 0.00036745606712788078, 0.00033380153932379921, 0.00033760634965878007, 0.00031836417172304503, 0.00034599036536937577, 0.00031850614364982183, 0.00032181081791655367, 0.0003193244130496919,0.00030886062451879659, 0.00034093836706992762], color='r', ls='', marker='s',mfc='None',mec = 'r', label = r"$MC$")
#ax.errorbar(2.0, -0.412937529942215, 0.00040939009499613692, color='r', ls='', marker='s',mfc='None',mec = 'r')

#ax.plot([0,0.25,0.5,0.75,1.0],
#        [0.68655103964,0.669437326132,0.658374050747, 0.652260692226, 0.650362166049],
#         marker = '.', ls='--',color=colors[4],  label = r'$ED\,T=0$')
#
#
#ax.errorbar([0,0.25,0.5,0.75,1.0],
#             [0.6857888921524403, 0.6688304172489111, 0.6567963678655567, 0.650764371410471, 0.651991734103279],
#             [0.0014768972604316717, 0.0010943507096154192, 0.0009784606815572308, 0.0008589020028244075, 0.0017278815509074313],
#             color=colors[0], ls='', marker='s',mfc='None',mec = colors[0], label = r'$MC\,T=0.01$')
#
#ax.errorbar([0,0.25,0.5,0.75,1.0],
#            [0.7121281957230885, 0.6874751454605332, 0.6679465993729068, 0.6597669728899697, 0.6550564079283808],
#            [0.0013802139169074978, 0.0027246818806691124, 0.0032315535081278357, 0.0027270357318958864, 0.0007336368485037216],
#             color=colors[1], ls='', marker='s',mfc='None',mec = colors[1], label=r'$MC\,T=0.1$')
#
#ax.errorbar([0,0.25,0.5,0.75,1.0],
#             [1.2814853663777406, 1.2560053721611615, 1.236301598802633, 1.2083679558130265, 1.1844370698040998],
#             [0.0022972750903136694, 0.002420015969385813, 0.0021509712033470534, 0.0024032878790228467, 0.0009184628185913591],
#             color=colors[2], ls='', marker='s',mfc='None',mec = colors[2], label=r'$MC\,T=1$')
#
#ax.errorbar([0,0.25,0.5,0.75,1.0],
#            [1.3863446131336763, 1.3825799076416605, 1.3799621706607206, 1.384404728859941, 1.3843063886983913],
#            [0.0019133408741475707, 0.0016603134862451091, 0.0017267660653117766, 0.0017560127043574926, 7.773619257301776e-05],
#             color=colors[3], ls='', marker='s',mfc='None',mec = colors[3], label=r'$MC\,T=10$')


