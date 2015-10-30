from qutip import *
from Hbuilder import *
from bmachine import *

import Hfile
import pyutils
import loadgmt,kevent
import numpy.random as rm
from pylab import *
import argparse
import numpy as np
import time
import collections
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def bitfield(n):
    return np.array([1 if digit=='1' else 0 for digit in bin(n)[2:]])

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--data',       help='Data file',                  type=str)
    parser.add_argument('--seed',       help='Seed used to generate data', type=int, default=0)
    parser.add_argument('--classic',    help='Perform classical training', action='store_true', default=False)
    parser.add_argument('--time',       help='Benchmark execution time ',  action='store_true', default=False)
    parser.add_argument('--verbose',    help='Verbose ',                   action='store_true', default=False)
    parser.add_argument('--inter','-I', help='Interactions file',          type=str)
    parser.add_argument('--beta', '-b', help='Inverse temperature ',       type=float)
    parser.add_argument('-X',           help='Lattice width ',             type=int)
    parser.add_argument('-Y', help='Lattice height (set to 1 for a 1-d system)', type=int)
    parser.add_argument('--OBC', help='Set open boundary conditions ', action='store_true', default=False)
    parser.add_argument('-J', help='Ising interaction',  type=float)
    parser.add_argument('-D', help='Transverse field',   type=float)
    parser.add_argument('-H', help='Longitudinal field', type=float)

    args = vars(parser.parse_args())
    print '--------------', args['seed']
    rm.seed(args['seed'])

    #  Error handling ---------------------------------------------------------
    if not((args['inter'] is None) or (('X' in args) and 
                                       ('Y' in args) and
                                       ('J' in args) and
                                       ('D' in args) and
                                       ('H' in args)
                                      )):
       print 'Define interactions'
       return 0

    if (args['inter'] is not None) and ((args['X'] is not None) or 
                                        (args['Y'] is not None) or
                                        (args['J'] is not None) or
                                        (args['D'] is not None) or
                                        (args['H'] is not None)
                                       ):
       print 'Conflicting arguments '
       return 0
    
    
    # Initialize couplings and fields -----------------------------------------
    N = 0 
    beta = args['beta']
    if (args['inter'] is not None):
        Z, X, ZZ, bonds = Hfile.LoadInters(args['inter'])
        N = len(Z)
    else:
        N = int(args['X'])*int(args['Y'])
        Z = np.ones(N)*float(args['H'])
        X = np.ones(N)*float(args['D'])
        if args['OBC']: bonds, ZZ  = OBC_Rect(int(args['X']), int(args['Y']), float(args['J']))
        else:           bonds,  ZZ = PBC_Rect(int(args['X']), int(args['Y']), float(args['J']))
    Nbonds = len(ZZ)
   
    
    print '---------Data acquisition----------'
    # Load or generate training set -------------------------------------------
    #Nclamped = N 
    #Ndata = 50000 
    #if (args['data'] is not None):
    #    data = np.loadtxt(args['data'], dtype='int32')
    #    Nclamped = data.shape[1]
    #    if Nclamped>N: 
    #        print 'Training set vectors exceed the graph size'
    #        return 0
    #    Ndata = data.shape[0]
    #    data  = data.tolist()
    #else:
    #    if Nclamped>N: 
    #        print 'Training set vectors exceed the graph size'
    #        return 0
    #    BM = BoltzmannMachine(N, bonds, Z, X, ZZ , beta)
    #    probTable = np.zeros(2**Nclamped)
    #    for i in range(2**Nclamped):
    #        cbits = bitfield(i)                              # convert i to a list of bits
    #        cbits = [0]*(Nclamped-len(cbits))+cbits.tolist() # keep the list length constant   
    #        BM.setProjector(cbits)
    #        if i==0: probTable[i] = real(BM.evaluateProjector())
    #        else:    probTable[i] = probTable[i-1] + real(BM.evaluateProjector())
    #   
    #    data = []
    #    index = 1
    #    for i in range(Ndata):
    #        RN = rm.random()
    #        index = searchsorted(probTable, RN)
    #        cbits = bitfield(index)
    #        cbits = [0]*(Nclamped-len(cbits))+cbits.tolist() # keep the list length constant 
    #        data += [cbits]
    #    del BM

    data = [[1, 0], [1,0], [1,0], [1,0], [0,1], [0,1], [0,1]]
    Ndata = 7
    
    
    # find unique states and count them
    udata = []
    cdata = collections.OrderedDict()
    for i,d in enumerate(data):
        if not(d in udata): udata += [d]; cdata[repr(d)]  = 1
        else:                             cdata[repr(d)] += 1
    weights = np.array(cdata.values())/float(Ndata)
    data    = udata
    Ndata   = len(data)
    print '--- Data (%d):    ' %len(data), data
    print '--- Weights (%4.2f) : '%np.sum(weights) , weights 
    weights = weights.reshape((Ndata,1))


    print '---------Gradient ascent----------'
    # Gradient descent --------------------------------------------------------
    
    # Compute the log-likelihood of the generating Hamiltonian
    #BM = BoltzmannMachine(N, bonds, Z, X, ZZ, beta)
    gLL = 0
    #for i, cbits in enumerate(data):
    #    BM.setProjector(cbits)
    #    gLL -= np.log(real(BM.evaluateProjector()))*weights[i]

    # General an initial guess for the Hamiltonian
    gdZ  = np.random.randn(N)*max([max(Z),0.5])+Z
    if args['classic']:  gdX  = np.zeros_like(Z)
    else:                gdX  = np.random.randn(N)*max([max(X),0.5])+X
    gdZZ = np.random.randn(Nbonds)*max([max(Z),0.5])+ZZ
   
    gdZ  = np.array([[0, -0.03040965],
                     [1, -0.04240908]])
    gdX  = np.array([[0, -1.0],
                     [1 ,-1.0]])
    gdZZ = np.array([[0,1, -0.004542226]])
    print "Initial guess for Z: ",  gdZ
    print "Initial guess for X: ",  gdX
    print "Initial guess for ZZ: ", gdZZ
    kwargs = {'X': gdX, 'Z1': gdZ, 'Z2': gdZZ}
    BM = BoltzmannMachine(N, beta, **kwargs)

    # Initialize proportionality factor schedule
    LL = 1.8/beta    # left limit
    RL = 0.4/beta   # right limit
    nsteps = 200    # number of steps
    b    = (1.0*nsteps*RL)/(LL - RL)
    a    = LL*b
    step = 0

    Ts   = []
    Norm = 100.0 # Norm of the gradient 
    while ((step < nsteps) and (Norm > 0.00001)):
        
        if args['verbose']: print '----------------step ', step, '-------------------' 
        # Averages
        Zavers  = np.zeros((Ndata+1, N))
        Xavers  = np.zeros((Ndata+1, N))
        ZZavers = np.zeros((Ndata+1, Nbonds))
        
        
        
        if args['time']: t0 = time.time()
        
        kwargs = {'X': gdX, 'Z1': gdZ, 'Z2': gdZZ}
        BM = BoltzmannMachine(N, beta, **kwargs)
        
        if args['time']: Ts += [time.time() - t0]
        


        # Compute the log-likelihood
        if args['time']: t0 = time.time()
        
        LL = 0
        for i,cbits in enumerate(data):
            BM.setProjector(cbits)
            LL -= np.log(real(BM.evaluateProjector())) * weights[i]
        if args['verbose']: print '--- LL = %0.4f vs %0.4f' %(LL, gLL)
        else:               print '--- step=%03d LL = %0.4f vs %0.4f' %(step, LL, gLL)
        
        if args['time']: Ts += [time.time() - t0]
        

        # Accumulate averages
        if args['time']: t0 = time.time()
        
        for i, cbits in enumerate(data):
            BM.setProjector(cbits)
            t = []
            if args['classic']: 
                Zavers[i, :] = -1.0*(np.array(cbits)*2-1)*beta
                for j,bond in enumerate(bonds): 
                    ZZavers[i, j] = (cbits[bond[0]]*2-1)*(cbits[bond[1]]*2-1)*beta
                    t += [(cbits[bond[0]]*2-1)*(cbits[bond[1]]*2-1)*beta]
            else: 
                (Zavers[i, :], Xavers[i,:], ZZavers[i,:]) = np.split(BM.computeLocalAverages(), [N, N+N])
        BM.setProjector([]) 
        (Zavers[Ndata,:], Xavers[Ndata,:], ZZavers[Ndata,:]) = np.split(BM.computeLocalAverages(), [N, N+N])
        
        if args['time']: Ts += [time.time() - t0]
        
        
        # Compute derivatives and the norm
        if args['classic']: dX = np.zeros(N)
        else:               dX = np.sum((Xavers[:Ndata,:]*weights),  axis=0) - beta*Xavers[Ndata]
        dZ   = np.sum((Zavers[:Ndata,:]*weights),  axis=0) - beta*Zavers[Ndata]
        dZZ  = np.sum((ZZavers[:Ndata,:]*weights), axis=0) - beta*ZZavers[Ndata]
        Norm = np.sqrt(np.sum(dZ*dZ) + np.sum(dX*dX) + np.sum(dZZ*dZZ))
        #eta   = a/(b+1.0*step)
        eta = 0.1
       
        # Follow the negative gradient 
        #if not args['classic']: gdX   += -eta*dX
        gdZ[:,-1]+= -eta*dZ
        gdZZ[:,-1]  += -eta*dZZ
       
       
        if  args['time']:
            print '--- times:  BM=%0.3f LL=%0.3f E=%0.3f ' %(Ts[0], Ts[1], Ts[2])
        if  args['verbose']:
            print '--- data vs model averages: '

            print '    Z:  ',
            for i in range(N): 
                zdata  = np.sum(Zavers[:Ndata,:]*weights , axis = 0)[i]
                zmodel = Zavers[Ndata,i]
                print '(%+0.3f vs %+0.3f) ' %(zdata, zmodel),

            print '\n    X:  ',
            for i in range(N):
                xdata  = np.sum(Xavers[:Ndata,:]*weights , axis = 0)[i]
                xmodel = Xavers[Ndata,i]
                print '(%+0.3f vs %+0.3f) ' %(xdata, xmodel),

            print '\n    ZZ: ',
            for i in range(Nbonds): 
                zzdata  = np.sum(ZZavers[:Ndata,:]*weights , axis = 0)[i]
                zzmodel = ZZavers[Ndata,i]
                print '(%+0.3f vs %+0.3f) ' %(zzdata, zzmodel),
            
            print '\n--- algotirhm: '
            print '    norm = %0.3f' %Norm  
            print '    eta  = %0.3f' %eta  

        
            print '--- couplings: '
            print '    Z:  ',
            for i in range(N): print '%+0.3f ' %gdZ[i],
            
            print '\n    X:  ',
            for i in range(N): print '%+0.3f ' %gdX[i],
            
            print '\n    ZZ: ',
            for i in range(Nbonds): print '%+0.3f ' %gdZZ[i],
            print
            print
        # Clear the memory required to store ED solver
        del BM
        step += 1 


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()



