import argparse
import bmachine, Hbuilder

import numpy  as np
import numpy.random as rm

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def bitfield(n):
    return np.array([1 if digit=='1' else 0 for digit in bin(n)[2:]])


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--seed',        help='RN seed', default=0,    type=int)
    parser.add_argument('--beta',  '-b', help='Inverse temperature ',  type=int)
    parser.add_argument('--Ndata', '-M', help='The size of data set ', type=int)
    parser.add_argument('--Nsites', '-N', help='Number of spins ',      type=int)

    args = vars(parser.parse_args())
    
    Nsites = args['Nsites']
    Ndata  = args['Ndata']
    beta   = args['beta']

    bonds  = Hbuilder.FullyConnected(Nsites)
    Nbonds = len(bonds)
    
    rm.seed(args['seed'])
    Js     = rm.random_integers( 0, 1, Nbonds)*2-1
    Hs     = rm.random_integers(-1, 1, Nsites)
    Ds     = np.zeros(Nsites)

    BM = bmachine.BoltzmannMachine(Nsites, bonds, Hs, Ds, Js, beta)

    probTable = np.zeros(2**Nsites)
    for i in range(2**Nsites):
        cbits = bitfield(i)                              # convert i to a list of bits
        cbits = [0]*(Nsites-len(cbits))+cbits.tolist() # keep the list length constant   
        BM.setProjector(cbits)
        if i==0: probTable[i] = np.real(BM.evaluateProjector())
        else:    probTable[i] = probTable[i-1] + np.real(BM.evaluateProjector())

    data = []
    index = 1
    for i in range(Ndata):
        RN = rm.random()
        index = np.searchsorted(probTable, RN)
        cbits = bitfield(index)
        cbits = [0]*(Nsites-len(cbits))+cbits.tolist() # keep the list length constant 
        data += [cbits]
    del BM

    fname = 'data_N-%02d_b-%05.2f_%05d.dat' %(Nsites, beta, args['seed'])
    np.savetxt(fname, np.array(data), fmt='%1d')
    
    head  = ' Sz Sx SzSz\n %d  %d  %d' %(Nsites, Nsites, Nbonds)
    fname = 'ham_FC_N-%02d_%05d.dat' %(Nsites, args['seed'])
    fJs   = np.hstack((np.array(bonds), np.expand_dims(Js, (1))))
    fHs   = np.transpose(np.vstack((np.arange(Nsites), np.arange(Nsites), Hs))) 
    fDs   = np.transpose(np.vstack((np.arange(Nsites), np.arange(Nsites), Ds))) 
    np.savetxt(fname, np.vstack((fHs, fDs, fJs)), fmt=' %1d', header=head)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
