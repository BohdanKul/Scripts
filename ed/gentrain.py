import bmachine, Hbuilder
import argparse, collections

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
    parser.add_argument('--minE',   help='Minimum desired dataset entropy', default=0,    type=float)
    parser.add_argument('--maxI',      help='Maximum number of iterations', default=1,    type=int)
    parser.add_argument('--maxD',      help='Maximum number of vectors in the dataset ', default=100,    type=int)
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
   

    Entropy  = 0
    bentropy = 0
    bsize    = 1000
    nI = 0
    weights = []
    while (((not (Entropy > args['minE'])) or (len(weights)>args['maxD'])) and (nI < args['maxI'])):
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
           _bonds3 = Hbuilder.tripleFullyConnected(Ns)
           bonds3 = []
           used   = []
           for i in range(int(np.ceil(len(_bonds3)*0.1))):
               index = rm.randint(0, len(_bonds3))
               if not(index in used): 
                   bonds3 += [_bonds3[index]]
                   used   += [index]
           
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


        # Find unique states and count them
        udata = []
        cdata = collections.OrderedDict()
        for i,d in enumerate(data):
            if not(d in udata): udata += [d]; cdata[repr(d)]  = 1
            else:                             cdata[repr(d)] += 1
        global weights
        weights = np.array(cdata.values())/float(Nd)

        Entropy = -1.0*np.sum(weights*np.log(weights))
        
        if (Entropy>bentropy) and (len(weights)<args['maxD']):
           bentropy = Entropy 
           bdata = data
           print "Entropy: ", Entropy, " iteration: ", nI, " size: ", len(weights)," weights: ", weights
        

        nI += 1

    if (nI == args['maxI']): 
        print "Reached the maximum number of iterations"
    data = bdata

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
