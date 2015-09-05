from qutip import *
from Hbuilder import *
from bmachine import *

import pyutils
import loadgmt,kevent
import random as rm
from pylab import *
import argparse
import numpy as np

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def bitfield(n):
    return np.array([1 if digit=='1' else 0 for digit in bin(n)[2:]])

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--data',       help='Data file',         type=str)
    parser.add_argument('--inter','-I', help='Interactions file', type=str)

    parser.add_argument('--beta', '-b', help='Inverse temperature ', type=float)
    
    parser.add_argument('-X', help='Lattice width ', type=int)
    parser.add_argument('-Y', help='Lattice height (set to 1 for a 1-d system)', type=int)
    parser.add_argument('--OBC', help='Set open boundary conditions ', action='store_true', default=False)
    parser.add_argument('-J', help='Ising interaction',  type=float)
    parser.add_argument('-D', help='Transverse field',   type=float)
    parser.add_argument('-H', help='Longitudinal field', type=float)

    args = vars(parser.parse_args())

    #  Error handling ---------------------------------------------------------
    if not((args['inter'] is None) or (('X' in args) and 
                                       ('Y' in args) and
                                       ('J' in args) and
                                       ('D' in args) and
                                       ('H' in args)
                                      )):
       print 'Define interactions'
       return 0

    print args['inter']
    if (args['inter'] is not None) and (('X' in args) or 
                                        ('Y' in args) or
                                        ('J' in args) or
                                        ('D' in args) or
                                        ('H' in args)
                                       ):
       print 'Conflicting arguments '
       return 0
    
    
    # Initialize couplings and fields -----------------------------------------
    N = 0 
    beta = args['beta']
    if (args['inter'] is not None):
        Z, X, ZZ, bonds = LoadInters(args.inter)
        N = len(Z)
    else:
        N = int(args['X'])*int(args['Y'])
        Z = np.ones(N)*float(args['H'])
        X = np.ones(N)*float(args['D'])
        if args['OBC']: bonds, ZZ  = OBC_Rect(int(args['X']), int(args['Y']), float(args['J']))
        else:           bonds,  ZZ = PBC_Rect(int(args['X']), int(args['Y']), float(args['J']))
    Nbonds = len(ZZ)
   
    #print X
    #print Z
    #print ZZ
    
    print '---------Data acquisition----------'
    # Load or generate training set -------------------------------------------
    Nclamped = N 
    Ndata = 50 
    if (args['data'] is not None):
        data = np.loadtxt(args['data'])
        Nclamped = data.shape[1]
        if Nclamped>N: 
            print 'Training set vectors exceed the graph size'
            return 0
        Ndata = data.shape[0]
    else:
        if Nclamped>N: 
            print 'Training set vectors exceed the graph size'
            return 0
        BM = BoltzmannMachine(N, bonds, Z, X, ZZ , beta)
        probTable = np.zeros(2**Nclamped)
        for i in range(2**Nclamped):
            cbits = bitfield(i)                              # convert i to a list of bits
            cbits = [0]*(Nclamped-len(cbits))+cbits.tolist() # keep the list length constant   
            BM.setProjector(cbits)
            if i==0: probTable[i] = real(BM.evaluateProjector())
            else:    probTable[i] = probTable[i-1] + real(BM.evaluateProjector())
       
        rm.seed(0)
        data = []
        index = 1
        for i in range(Ndata):
            RN = rm.random()
            index = searchsorted(probTable, RN)
            #index = 0
            cbits = bitfield(index)
            cbits = [0]*(Nclamped-len(cbits))+cbits.tolist() # keep the list length constant 
            data += [cbits]
        del BM

    data = [[1, 0], [1,0], [1,0], [1,0], [0,1], [0,1], [0,1]]
    Ndata = 7
    print 'Data: ', data
    print 
    print '---------Gradient ascent----------'
    # Gradient descent --------------------------------------------------------
    
    # Compute the log-likelihood of the generating Hamiltonian
    BM = BoltzmannMachine(N, bonds, Z, X, ZZ, beta)
    gLL = 0
    for cbits in data:
        BM.setProjector(cbits)
        gLL -= np.log(real(BM.evaluateProjector()))/(1.0 * Ndata)

    # General an initial guess for the Hamiltonian
    gdZ  = np.random.randn(N)*max([max(Z),0.5])+Z
    gdX  = np.random.randn(N)*max([max(Z),0.5])+X
    #gdX  = X
    gdZZ = np.random.randn(Nbonds)*max([max(Z),0.5])+ZZ
   
    gdZ  = np.array([-0.03040965, -0.04240908])
    gdX  = np.array([-1.0,-1.0])
    gdZZ = np.array([-0.004542226])
    print "Initial guess for Z: ",  gdZ
    print "Initial guess for X: ",  gdX
    print "Initial guess for ZZ: ", gdZZ

    LL = .5/beta
    RL = 0.05
    nsteps = 100
    b    = (1.0*nsteps*RL)/(LL - RL)
    a    = LL*b
    step = 0
    Norm = 100.0
    while ((step < nsteps) and (Norm > 0.01)):
        print '----------------step ', step
        # Initialize
        Zavers  = np.zeros((Ndata+1, N))
        Xavers  = np.zeros((Ndata+1, N))
        ZZavers = np.zeros((Ndata+1, Nbonds))
        BM = BoltzmannMachine(N, bonds, gdZ, gdX, gdZZ, beta)
        
        # Compute the log-likelihood
        LL = 0
        for cbits in data:
            BM.setProjector(cbits)
            LL -= np.log(real(BM.evaluateProjector()))/(1.0 * Ndata)
        print '   LL = %0.4f vs generating H LL = %0.4f' %(LL, gLL)
        
        # Accumulate averages
        for i, cbits in enumerate(data):
            BM.setProjector(cbits)
            Zavers[i, :], Xavers[i,:], ZZavers[i,:] = BM.computeLocalAverages()
            #Zavers[i, :], ZZavers[i, :] = -1.0*(np.array(cbits)*2-1)*beta, (cbits[0]*2-1)*(cbits[1]*2-1)*beta
        BM.setProjector([])
        Zavers[Ndata,:], Xavers[Ndata,:], ZZavers[Ndata,:] = BM.computeLocalAverages()

        #print Zavers[:Ndata]
        #print Zavers[Ndata]
        #print Zavers[:Ndata] - Zavers[Ndata]
        # Compute derivatives and the norm
        dZ   = np.sum((Zavers[:Ndata]  - beta*Zavers[Ndata]),  axis=0)/(1.0*Ndata) 
        #dX   = np.sum((Xavers[:Ndata]  - beta*Xavers[Ndata]),  axis=0)/(1.0*Ndata) 
        dX   = np.zeros_like(Xavers[0])
        dZZ  = np.sum((ZZavers[:Ndata] - beta*ZZavers[Ndata]), axis=0)/(1.0*Ndata) 
        Norm = np.sqrt(np.sum(dZ*dZ) + np.sum(dX*dX) + np.sum(dZZ*dZZ))
        #print '   Norm = %0.3f' %Norm  
        #print '  Unclamped <Z>:          ', Zavers[Ndata] 
        #print '  Unclamped <X>:          ', Xavers[Ndata] 
        #print '  Unclamped <ZZ>:         ', ZZavers[Ndata] 
        #print '  Clamped <Z>:            ', Zavers[0], ' for ', data[0] 
        #print '  Clamped <X>:            ', Xavers[0]
        #print '  Clamped <ZZ>:           ', ZZavers[3],' for ', data[3]
        #print '  Average clamped <Z>:    ', np.sum(Zavers[:Ndata] , axis = 0)/(1.0*Ndata)#*beta)
        #print '  Average clamped <X>:    ', np.sum(Xavers[:Ndata] , axis = 0)/(1.0*Ndata)#*beta)
        #print '  Average clamped <ZZ>:   ', np.sum(ZZavers[:Ndata], axis = 0)/(1.0*Ndata)#*beta)
        
        #eta   = a/(b+1.0*step)
        eta = 0.1
        #print '  eta                     ', eta 
        #print '  dZ:                     ', dZ
        #print '  dX:                     ', dX
        #print '  dZZ:                    ', dZZ

        # Follow the negative gradient 
        gdZ   += -eta*dZ
        gdX   += -eta*dX
        gdZZ  += -eta*dZZ
        #print '  Z field:                ', gdZ
        #print '  X field:                ', gdX
        #print '  ZZ field:               ', gdZZ
        # Clear the memory required to store ED solver
        del BM
        step += 1 


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()



