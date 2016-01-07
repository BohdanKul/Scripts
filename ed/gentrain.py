import bmachine, Hbuilder, kevent
import argparse, collections
import pylab as pl

import numpy        as np
import numpy.random as rm
import scipy.stats  as stats

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def bitfield(n):
    return np.array([1 if digit=='1' else 0 for digit in bin(n)[2:]])

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def getCouplings(Ns, mode, alpha):
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
    if mode == 'triple':
       _bonds3 = Hbuilder.tripleFullyConnected(Ns)
       bonds3 = []
       used   = []
       for i in range(int(np.ceil(len(_bonds3)*0.1))):
           index = rm.randint(0, len(_bonds3))
           if not(index in used): 
               bonds3 += [_bonds3[index]]
               used   += [index]
       
       JJs    = (rm.random_integers(0, 1, len(bonds3))*2-1)*alpha
       JJs    = JJs.reshape((len(bonds3),1))
       JJs    = np.hstack((np.array(bonds3), JJs)) 
       kwargs['Z3'] = JJs

    return kwargs

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def probBoltzmann(v, BM):
    BM.setProjector(v)
    return BM.evaluateProjector()


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def probProduct(modes, iweight, v):
    v = np.array(v)
    if len(modes.shape)>1: P = np.prod((1+(modes*2-1)*(v*2-1)*iweight)/2.0, axis=1)
    else:                  P = np.prod((1+(modes*2-1)*(v*2-1)*iweight)/2.0)

    return P

#------------------------------------------------------------------------------
#----------------- -------------------------------------------------------------
def probDiscreteGauss(mean, Ns, sigma=1.0):
    N = 2**Ns
    w  = N//2 - 0.5
    (LB, RB) = (-w, w)
    Es = np.linspace(LB, RB, N+1)
    #X = stats.truncnorm(LB/sigma, RB/sigma)
    Ps = np.diff(stats.truncnorm.cdf(Es, LB/sigma, RB/sigma, scale=sigma))
    
    return np.roll(Ps, -N//2+mean)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def probMMDiscreteGauss(v, modes, Ns, sigma=1.0):
    Nmodes = len(modes)
    
    modePs = np.zeros((Nmodes, 2**Ns))
    for i,mode in enumerate(modes):
        modePs[i,:] = probDiscreteGauss(mode, Ns, sigma)
    
    P = 1.0/(1.0*Nmodes)*np.sum(modePs, axis=0)[v]

    return P 
 


#------------------------------------------------------------------------------
#----------------- -------------------------------------------------------------
def probMM(v, modes, iweight = 1.0):

    if len(modes.shape)>1: Nmodes = modes.shape[0]
    else:                  Nmodes = 1
    
    P = 1.0/(1.0*Nmodes)*np.sum(probProduct(modes, iweight, v))

    return P 
    
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def probrMM(v, modes, Ns, rNs, iweight = 1.0):
    if len(modes.shape)>1: Nmodes = modes.shape[0]
    else:                  Nmodes = 1
    
    rmodes = np.ones((Nmodes, Ns))
    rmodes[:, :rNs] = modes

    sigmas = np.ones((1, Ns))*0.5
    sigmas[:rNs] *= iweight

    P = 1.0/(1.0*Nmodes)*np.sum(probProduct(rmodes, sigmas, v))

    return P



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--mode',      type=str,   default='pair', choices=['pair', 'triple', 'mean', 'MM', 'dMM', 'rMM'])
    parser.add_argument('--alpha',     type=float, default=0,   help='Amplitude of triple bond')
    parser.add_argument('--modes',     type=int,   default=1,   help='Number of modes in MM distribution')
    parser.add_argument('--rsize',     type=int,   default=1,   help='Number of clamped qubits in rMM')
    parser.add_argument('--iweight',   type=float, default=0,   help='Inverse weight of multi-modal distribution')
    parser.add_argument('--minE',      type=float, default=0,   help='Minimum desired dataset entropy')
    parser.add_argument('--maxI',      type=int,   default=100, help='Maximum number of iterations')
    parser.add_argument('--maxD',      type=int,   default=100, help='Maximum number of vectors in the dataset ')
    parser.add_argument('--seed',      type=int,   default=0,   help='RN seed')
    parser.add_argument('--batch',    type=int,    default=1,   help='Batch size')
    parser.add_argument('--visual',    type=int,   default=0,   help='Visualize probability distribution')
    parser.add_argument('--beta','-b', type=float, default=1.0, help='Inverse temperature ')
    parser.add_argument('--Nd',  '-M', type=int,   help='The size of data set ')
    parser.add_argument('--Ns',  '-N', type=int,   help='Number of spins ')

    args = vars(parser.parse_args())
   
    Ns      = args['Ns']
    Nd      = args['Nd']
    beta    = args['beta']
    Nmodes  = args['modes'] 
    iweight = args['iweight'] 
    
    for i in range(args['batch']):
        print "\n\n---------- Batch ", i, " ------------"
        seed = args['seed'] + i 
        rm.seed(seed)
        Entropy  = 0
        bentropy = 0
        bsize    = 1000
        nI = 0
        weights = []
        
        while (((not (Entropy > args['minE'])) or (len(weights)>args['maxD'])) and (nI < args['maxI'])):

            params = {}
            if   args['mode'] == 'MM':
                modes = np.random.randint(2, size = Nmodes*Ns).reshape((Nmodes, Ns))
            elif args['mode'] == 'dMM' :
                modes = np.random.randint(2**Ns, size = Nmodes)
            elif args['mode'] == 'rMM' :
                modes = np.random.randint(2, size = Nmodes*args['rsize']).reshape((Nmodes, args['rsize']))
            else:
                kwargs = getCouplings(Ns, args['mode'], args['alpha'])
                BM = bmachine.BoltzmannMachine(Ns, beta, **kwargs)

            cdf = np.zeros(2**Ns)
            pdf = np.zeros(2**Ns)
            for i in range(2**Ns):
                cbits = bitfield(i)                        # convert i to a list of bits
                cbits = [0]*(Ns-len(cbits))+cbits.tolist() # keep the list length constant   
                if   args['mode']=='MM':    P = probMM(cbits,  modes, iweight)
                elif args['mode']=='dMM':   P = probMMDiscreteGauss(i, modes, Ns, sigma=iweight)
                elif args['mode']=='rMM':   P = probrMM(cbits, modes, Ns, args['rsize'], iweight)
                else:                       P = probBoltzmann(cbits, BM)
                
                if i==0: cdf[i] = P 
                else:    cdf[i] = cdf[i-1] + P 
                pdf[i] = P
            data = []
            index = 1
            for i in range(Nd):
                RN = rm.random()
                index = np.searchsorted(cdf, RN)
                cbits = bitfield(index)
                cbits = [0]*(Ns-len(cbits))+cbits.tolist() # keep the list length constant 
                data += [cbits]

            if not(args['mode'] in ['MM', 'dMM', 'rMM']): del BM


            # Find unique states and count them
            udata = []
            cdata = collections.OrderedDict()
            for i,d in enumerate(data):
                if not(d in udata): udata += [d]; cdata[repr(d)]  = 1
                else:                             cdata[repr(d)] += 1
            global weights
            weights = np.array(cdata.values())/float(Nd)

            DS = -1.0*np.sum(weights*np.log(weights))
            nonzero = np.nonzero(pdf)
            MS = -1.0*np.sum(pdf[nonzero]*np.log(pdf[nonzero]))
           
            if (DS>bentropy) and (len(weights)<args['maxD']):
               bentropy = DS 
               bdata = data
               print "Model entropy: ", MS, "Data entropy: ", DS, " iteration: ", nI, " size: ", len(weights)," weights: ", weights
               if (args['visual']==1): 
                   fig = pl.figure(1, figsize=(13,6))
                   pl.connect('key_press_event',kevent.press)
                   ax  = pl.subplot(111)
                   ax.plot(pdf)
                   pl.show()    
            

            nI += 1

        if (nI == args['maxI']): 
            print "Reached the maximum number of iterations"
        data = bdata

        if args['mode'] == 'MM':
            fname = 'data_MM_N-%02d_modes-%02d_iw-%04.2f_%05d.dat' %(Ns, Nmodes, iweight, seed)

        elif args['mode'] == 'dMM':
            fname = 'data_dMM_N-%02d_modes-%02d_iw-%04.2f_%05d.dat' %(Ns, Nmodes, iweight, seed)
        elif args['mode'] == 'rMM':
            fname = 'data_rMM_N-%02d_Nr-%02d_modes-%02d_iw-%04.2f_%05d.dat' %(Ns, args['rsize'],Nmodes, iweight, seed)
        else:
            fname = 'data_Hamiltonian_N-%02d_b-%05.2f_alpha-%04.2f_%05d.dat' %(Ns, beta, args['alpha'], seed)
        
        np.savetxt(fname, np.array(data), fmt='%1d')
   
    Z2 = kwargs['Z2'].shape[0]
    
    if 'Z3' in kwargs.keys(): Z3 = kwargs['Z3'].shape[0]
    else:                     Z3 = 0
    head  = ' Sz Sx SzSz SzSzSz\n %d  %d    %d      %d' %(Ns, Ns, Z2, Z3)
    Hfname = 'gH_mode-%s_N-%02d_alpha-%04.2f_%05d.dat' %(args['mode'], Ns, args['alpha'], seed)
    if  args['mode'] == 'pair':
        Hs   = np.hstack((np.arange(Ns).reshape((Ns,1)), kwargs['Z1']))
        Ds   = np.hstack((np.arange(Ns).reshape((Ns,1)), kwargs['X'])) 
        np.savetxt(Hfname, np.vstack((Hs, Ds, kwargs['Z2'])), fmt=' %1d', header=head)
    if  args['mode'] == 'triple':
        Hs   = np.hstack((np.arange(Ns).reshape((Ns,1)), np.arange(Ns).reshape((Ns,1)), Hs))
        Ds   = np.hstack((np.arange(Ns).reshape((Ns,1)), np.arange(Ns).reshape((Ns,1)), Ds)) 
        Js   = np.hstack((np.arange(Js.shape[0]).reshape((Js.shape[0],1)), Js)) 
        np.savetxt(Hfname, np.vstack((Hs, Ds, Js, JJs)), fmt=' %1d', header=head)
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
