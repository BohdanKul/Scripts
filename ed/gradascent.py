from qutip import *
from Hbuilder import *
from bmachine import *

import pyutils
import loadgmt,kevent
from pylab import *
import argparse
import numpy as np

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('data',         help='Data file', type=str)
    parser.add_argument('--inter','-I', help='Interactions file', type=str)

    parser.add_argument('--beta', '-b', help='Inverse temperature ', type=float)
    
    parser.add_argument('-X', help='Lattice width ', type=int)
    parser.add_argument('-Y', help='Lattice height (set to 1 for a 1-d system)', type=int)
    parser.add_argument('--OBC', help='Set open boundary conditions ', action='store_true', default=False)
    parser.add_argument('-J', help='Ising interaction',   type=int)
    parser.add_argument('-D', help='Transverse field',   type=int)
    parser.add_argument('-H', help='Longitudinal field', type=int)

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
    if (args['inter'] is not None):
        Zfield, Xfield, Inter = LoadInters(args.inter)
        N = len(Zfield)
    else:
        N = int(args['X'])*int(args['Y'])
        Zfield = np.ones(N)*float(args['H'])
        Xfield = np.ones(N)*float(args['D'])
        if args['OBC']: Inter = OBC_Rect(int(args['X']), int(args['Y']), float(args['J']))
        else:           Inter = PBC_Rect(int(args['X']), int(args['Y']), float(args['J']))

   
    #print Xfield
    #print Zfield
    #print Inter
    
    # Load training set ------------- -----------------------------------------
    data = np.loadtxt(args['data'])
    if data.shape[1]>N: 
        print 'Training set vectors exceed the graph size'
        return 0
    
    
    # Initialize the Hamiltonian ----------------------------------------------
    beta = args['beta']
    BM = BoltzmannMachine(N, Zfield, Zfield, Inter, beta)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()



