import Hfile, Hbuilder, bmachine
import argparse, collections
import numpy.random as rm
import numpy as np
from scipy.optimize import fmin_l_bfgs_b as LBFGS

(Ns, Nb, counter, beta, fname, mode, delta) = (0,0,0,0,'','class', 0)
(data, weights, bonds) = ([],[],[])

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def getLL(Hs, Ds, Js):
    kwargs = {'X': Ds, 'Z1': Hs, 'Z2': Js}
    BM = bmachine.BoltzmannMachine(Ns, beta, **kwargs)
    gLL = 0
    for i, cbits in enumerate(data):
        BM.setProjector(cbits)
        gLL -= np.log(np.real(BM.evaluateProjector()))*weights[i]
    del BM

    return gLL

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def progtracker(inter):
    global Ns, counter, beta, fname, mode
    global data, weights, bonds

    Nb = len(bonds)
    counter += 1
    Hs, Ds, Js = unpackInter(inter, bonds, Ns, mode)
    
    # Compute the log-likelihood of the generating Hamiltonian
    kwargs = {'X': Ds, 'Z1': Hs, 'Z2': Js}
    BM = bmachine.BoltzmannMachine(Ns, beta, **kwargs)
    LL = 0
    for i, cbits in enumerate(data):
        BM.setProjector(cbits)
        LL -= np.log(np.real(BM.evaluateProjector()))*weights[i]
    del BM

    # store the values
    f = open(fname, 'a')
    str = '    %8.4f' %LL
    for meas in np.hstack([Hs[:,-1], Ds[:,-1], Js[:,-1]]):
        str+= '    %8.4f' %meas
    str += '\n'
    f.write(str)
    f.close()

    return 0

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def LLgrad(inter, *args):
    ''' Train all couplings and fields'''
    (Ns, beta, bonds, data, weights, cmode) = args   
    Nd = len(data)
    Nb = len(bonds) 
    
    # initialize Boltzmann machine 
    Hs, Ds, Js = unpackInter(inter, bonds, Ns, cmode)
    kwargs = {'X': Ds, 'Z1': Hs, 'Z2': Js}
    BM = bmachine.BoltzmannMachine(Ns, beta, **kwargs)
        
    # compute the log-likelihood
    vLL = 0
    for i,cbits in enumerate(data):
        BM.setProjector(cbits)
        vLL -= np.log(np.real(BM.evaluateProjector())) * weights[i]
     
    # compute the gradient
    aves  = np.zeros((Nd, Ns+Ns+Nb))
    for i, cbits in enumerate(data):
        BM.setProjector(cbits)
        if (cmode=='class') and (i!=Nd-1): 
           aves[i, :Ns] = -1.0*(np.array(cbits)*2-1)*beta
           for j,bond in enumerate(bonds): 
               aves[i, 2*Ns+j] = (cbits[bond[0]]*2-1)*(cbits[bond[1]]*2-1)*beta
        else: 
            aves[i,:] = BM.computeLocalAverages()
    gLL = np.sum((aves*weights),  axis=0)
    gLL = packInter(gLL, Ns, cmode)
    
    return vLL, gLL


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def unpackInter(inter, bonds, Ns, mode):
    if  mode=="class":     
        (Hs, Js)     = np.split(inter, [Ns])
        Ds = np.zeros(Ns)
    elif mode=="quant":     
        (Hs, D, Js)  = np.split(inter, [Ns, Ns+1])
        Ds = np.ones(Ns)*D
    elif mode=="quant_all": 
        (Hs, Ds, Js) = np.split(inter, [Ns, Ns+Ns])
    elif mode=="quant_stat": 
        (Hs, Js) = np.split(inter, [Ns])
        global delta
        Ds = np.ones(Ns)*delta

    
    Nb = len(bonds) 
    Ds = np.transpose(np.vstack((np.arange(Ns), Ds)))
    Hs = np.transpose(np.vstack((np.arange(Ns), Hs)))
    Js = np.hstack((bonds, Js.reshape(Nb,1)))

    return Hs, Ds, Js

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def packInter(gLL, Ns, mode):
    (gHs, gDs, gJs) = np.split(gLL, [Ns, Ns+Ns])
    if   mode=="class":      return np.hstack([gHs, gJs])
    elif mode=="quant":      return np.hstack([gHs, np.sum(gDs), gJs])
    elif mode=="quant_stat": return np.hstack([gHs, gJs]) 
    elif mode=="quant_all":  return gLL
    
    
        

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def bitfield(n):
    return np.array([1 if digit=='1' else 0 for digit in bin(n)[2:]])

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main(): 

    modes  =  ['class', 'quant', 'quant_all', 'quant_stat'] 
    
    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--seed',       help='Seed used to generate data', type=int, default=0)
    parser.add_argument('--data',       help='Data file',                  type=str)
    parser.add_argument('--inter','-I', help='Interactions file',          type=str)
    parser.add_argument('--delta',      help='Transverse field',           type=float)
    parser.add_argument('--beta', '-b', help='Inverse temperature ',       type=float)
    parser.add_argument('--mode',       help='Training mode', choices=modes, default='class')

    args = vars(parser.parse_args())
    rm.seed(args['seed'])

    #  Error handling ---------------------------------------------------------
    if (args['inter'] is None):
       print 'Define interactions'
       return 0

    if ((args['mode'] in ['quant', 'quant_all', 'quant_stat']) and (args['delta'] is None)):
       print 'Define the transverse field of the quantum model'
       return 0
    global mode 
    mode = args['mode']

    # Initialize couplings and fields -----------------------------------------
    Ns = 0 
    global beta
    beta = args['beta']
    if (args['inter'] is not None):
        global bonds
        Hs, Ds, Js, bonds = Hfile.LoadInters(args['inter'])
        if (args['mode'] in ['quant', 'quant_all']): 
            Ds[:,1] = np.ones(len(Ds))*args['delta']
        global Ns
        Ns = len(Hs)
    else:
        Ns = int(args['X'])*int(args['Y'])
        Hs = np.ones(N)*float(args['H'])
        Ds = np.ones(N)*float(args['D'])
        if args['OBC']: bonds, Js = OBC_Rect(int(args['X']), int(args['Y']), float(args['J']))
        else:           bonds, Js = PBC_Rect(int(args['X']), int(args['Y']), float(args['J']))
    Nb = len(Js)
   
   

    # Load or generate training set -------------------------------------------
    print '---------Data acquisition----------'
    Nclamped = Ns 
    Nd = 50000 
    if (args['data'] is not None):
        global data
        data = np.loadtxt(args['data'], dtype='int32')
        Nclamped = data.shape[1]
        if Nclamped>Ns: 
            print 'Training set vectors exceed the graph size'
            return 0
        Nd = data.shape[0]
        data  = data.tolist()
    else:
        if Nclamped>N: 
            print 'Training set vectors exceed the graph size'
            return 0
        kwargs = {'X': Ds, 'Z1': Hs, 'Z2': Js}
        BM = bmachine.BoltzmannMachine(Ns, beta, **kwargs)
        probTable = np.zeros(2**Nclamped)
        for i in range(2**Nclamped):
            cbits = bitfield(i)                              # convert i to a list of bits
            cbits = [0]*(Nclamped-len(cbits))+cbits.tolist() # keep the list length constant   
            BM.setProjector(cbits)
            if i==0: probTable[i] = np.real(BM.evaluateProjector())
            else:    probTable[i] = probTable[i-1] + np.real(BM.evaluateProjector())
       
        data = []
        index = 1
        for i in range(Nd):
            RN = rm.random()
            index = searchsorted(probTable, RN)
            cbits = bitfield(index)
            cbits = [0]*(Nclamped-len(cbits))+cbits.tolist() # keep the list length constant 
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
    data    = udata
    data += [[]]          # add a vector with no bits clamped
    entropy = -1.0*np.sum(weights*np.log(weights))
    print type(weights), entropy
    weights = np.append(weights, -beta) # and its weight
    Nd   = len(data)
    print '--- Data (%d):    '     %len(data), data
    print '--- Weights (%4.2f) : ' %np.sum(weights) , weights 
    print '--- Entropy : %4.2f'    %entropy
    weights = weights.reshape((Nd,1))

    # General an initial guess for the Hamiltonian -----------------------------
    Hs[:,-1] = np.random.randn(Ns)*np.abs(Hs[:,-1]).max() 
    Js[:,-1] = np.random.randn(Nb)*np.abs(Js[:,-1]).max()
    
    if   args['mode'] == 'class':      iparams = np.hstack([Hs[:,-1],                Js[:,-1]])
    elif args['mode'] == 'quant':      iparams = np.hstack([Hs[:,-1], args['delta'], Js[:,-1]])
    elif args['mode'] == 'quant_all':  iparams = np.hstack([Hs[:,-1], Ds[:,-1],      Js[:,-1]])
    elif args['mode'] == 'quant_stat': iparams = np.hstack([Hs[:,-1],                Js[:,-1]])

 
    # impose limits on the trained values
    #lims    = []
    #f = 5.0
    #for i,par in enumerate(iparams):
    #    apar = np.abs(par)*f
    #    lims.append((-apar, apar))

    # Compute the log-likelihood of the generating Hamiltonian
    gLL = getLL(Hs, Ds, Js)

    # Creat a log file 
    global delta
    if args['delta']!= None: delta=args['delta']
    global fname 
    fname = 'train_delta-%04.2f_%s' %(args['delta'], args['data'][5:])
    
    f = open(fname, 'w')
    header = '#%11s' %('LL')
    for i in range(Hs.shape[0]):   header += '%12s' %('H'+str(i)) 
    for i in range(Ds.shape[0]):   header += '%12s' %('D'+str(i)) 
    for bond in Js[:,:2].tolist(): header += '%12s' %('J(%d, %d)' %(bond[0], bond[1]))
    header +='\n'
    f.write(header)

    st = '    %8.4f' %gLL
    for meas in np.hstack([Hs[:,-1], Ds[:,-1], Js[:,-1]]):
        st+= '    %8.4f' %meas
    st += '\n'
    f.write(st)
    f.close()


    print '---------Gradient ascent----------'
    # Gradient descent --------------------------------------------------------
    
    
    counter = 0
    fx, fLL, info = LBFGS(LLgrad, x0=iparams, args = (Ns, beta, bonds, data, weights, args['mode']), iprint = 1, pgtol=0.001, factr=1e13/2.2, callback=progtracker)  

    print fx
    print info
    #print lims

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()



