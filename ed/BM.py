import Hfile, Hbuilder, bmachine
import argparse, collections, time
import numpy.random as rm
import numpy as np
from scipy.optimize import fmin_l_bfgs_b as LBFGS

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def bitfield(n):
    return np.array([1 if digit=='1' else 0 for digit in bin(n)[2:]])
        
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def Bits2Int(bitlist):
    out = 0
    for bit in bitlist:
        out = (out << 1) | bit

    return out


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--clamped',    help='File with a clamped vector', type=str)
    parser.add_argument('--inter','-I', help='Interactions file',          type=str)
    parser.add_argument('--beta', '-b', help='Inverse temperature ',       type=float, default=1.0)

    args = vars(parser.parse_args())

   
    #  Error handling ---------------------------------------------------------
    if (args['inter'] is None):
       print 'Define interactions'
       return 0

    # Initialize couplings and fields -----------------------------------------
    beta = args['beta']
    Hs, Ds, Js, bonds = Hfile.LoadInters(args['inter'])
    Ns = len(Hs)
    Nb = len(Js)

    # Initialize BM 
    kwargs = {'X': Ds, 'Z1': Hs, 'Z2': Js}
    BM = bmachine.BoltzmannMachine(Ns, beta, **kwargs)
   

    # Load projective vector --------------------------------------------------
    if (args['clamped'] is not None):
        cbits = np.loadtxt(args['clamped'], dtype='int32')
        cbits = Hfile.MCtoED(cbits)
    else: 
        cbits = []

    # Pass it to BM and compute correlations of interest ----------------------
    BM.setProjector(cbits)
    aves = BM.computeLocalAverages()
    del BM
    
    # Record results to a file
    fname = 'ED_N-%02d_P-%04d.dat' %(Ns, Bits2Int())
    f = open(fname, 'w')
    
    aves = np.hstack((aves[Ns:Ns+Ns), aves[:Ns], aves[Ns+Ns:]))
    for i in range(Ds.shape[0]):   header += '%20s' %('<X'+str(i)+'>') 
    for i in range(Hs.shape[0]):   header += '%20s' %('<Z'+str(i)+'>') 
    for bond in Js[:,:2].tolist(): header += '%20s' %('<ZZ(%d, %d)>' %(bond[0], bond[1]))
    header += '\n'
    f.write(header)

    for meas in aves[]:  st+= '    %16.8E' %meas
    st += '\n'
        
    f.write(st)
    f.close()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()



