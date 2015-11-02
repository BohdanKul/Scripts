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
    parser.add_argument('--mode', choices=['pair', 'triple', 'mean'], default='pair',    type=str)
    parser.add_argument('--alpha',     help='Amplitude of triple bond', default=0,    type=float)
    parser.add_argument('--seed',      help='RN seed', default=0,    type=int)
    parser.add_argument('--beta','-b', help='Inverse temperature ',  type=float)
    parser.add_argument('--mS',        help='Minimum  data entropy', type=float, default =0)
    parser.add_argument('--Nd',  '-M', help='The size of data set ', type=int)
    parser.add_argument('--Ns',  '-N', help='Number of spins ',      type=int)

    args = vars(parser.parse_args())
   
    Ns   = args['Ns']
    Nd   = args['Nd']
    beta = args['beta']
    
    rm.seed(args['seed'])
   
    sites = np.arange(Ns)
    sites = sites.reshape((Ns,1))
    Ds    = np.zeros(Ns)
    Ds    = Ds.reshape((Ns,1))
    Ds    = np.hstack((sites, Ds))
    

    Hs    = rm.random_integers(-1, 1, Ns)
    Hs    = Hs.reshape((Ns,1))
    Hs    = np.hstack((sites, Hs))

    bonds = Hbuilder.pairFullyConnected(Ns)
    Js  = rm.random_integers( 0, 1, len(bonds))*2-1
    Js  = Js.reshape((len(bonds),1))
    Js  = np.hstack((np.array(bonds), Js ))
    kwargs = {}

    kwargs['X']  = Ds 
    kwargs['Z1'] = Hs
    kwargs['Z2'] = Js 
   
    JJs = []
    if args['mode'] == 'triple':
       bonds3 = Hbuilder.tripleFullyConnected(Ns)
       JJs    = (rm.random_integers(0, 1, len(bonds3))*2-1)*args['alpha']
       JJs    = JJs.reshape((len(bonds3),1))
       JJs    = np.hstack((np.array(bonds3), JJs)) 
       kwargs['Z3'] = JJs

    BM = bmachine.BoltzmannMachine(Ns, beta, **kwargs)

    probTable = np.zeros(2**Ns)
    for i in range(2**Ns):
        cbits = bitfield(i)                        # convert i to a list of bits
        cbits = [0]*(Ns-len(cbits))+cbits.tolist() # keep the list length constant   
        BM.setProjector(cbits)
        if i==0: probTable[i] = np.real(BM.evaluateProjector())
        else:    probTable[i] = probTable[i-1] + np.real(BM.evaluateProjector())

    data = []
    index = 1
    for i in range(Nd):
        RN = rm.random()
        index = np.searchsorted(probTable, RN)
        cbits = bitfield(index)
        cbits = [0]*(Ns-len(cbits))+cbits.tolist() # keep the list length constant 
        data += [cbits]
    del BM

    fname = 'data_N-%02d_b-%05.2f_alpha-%04.2f_%05d.dat' %(Ns, beta, args['alpha'], args['seed'])
    np.savetxt(fname, np.array(data), fmt='%1d')
   
    head  = ' Sz Sx SzSz SzSzSz\n %d  %d    %d      %d' %(Ns, Ns, len(Js), len(JJs))
    Hfname = 'gH_mode-%s_N-%02d_alpha-%04.2f_%05d.dat' %(args['mode'], Ns, args['alpha'], args['seed'])
    if  args['mode'] == 'pair':
        Hs   = np.hstack((np.arange(Ns).reshape((Ns,1)), Hs))
        Ds   = np.hstack((np.arange(Ns).reshape((Ns,1)), Ds)) 
        np.savetxt(Hfname, np.vstack((Hs, Ds, Js)), fmt=' %1d', header=head)
    if  args['mode'] == 'triple':
        Hs   = np.hstack((np.arange(Ns).reshape((Ns,1)), np.arange(Ns).reshape((Ns,1)), Hs))
        Ds   = np.hstack((np.arange(Ns).reshape((Ns,1)), np.arange(Ns).reshape((Ns,1)), Ds)) 
        Js   = np.hstack((np.arange(Js.shape[0]).reshape((Js.shape[0],1)), Js)) 
        np.savetxt(Hfname, np.vstack((Hs, Ds, Js, JJs)), fmt=' %1d', header=head)
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
