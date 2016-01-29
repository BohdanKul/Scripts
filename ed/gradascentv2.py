import Hfile, Hbuilder, bmachine
import argparse, collections, time
import numpy.random as rm
import numpy as np
from scipy.optimize import fmin_l_bfgs_b as LBFGS

(Ns, Nb, counter, beta, fname, mode, delta, gradEval) = (0,0,0,0,'','class', 0, True)
(data, weights, bonds) = ([],[],[])

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def getLLandAves(Hs, Ds, Js):
    kwargs = {'X': Ds, 'Z1': Hs, 'Z2': Js}
    BM = bmachine.BoltzmannMachine(Ns, beta, **kwargs)
    gLL = 0
    for i, cbits in enumerate(data[:-1]):
        BM.setProjector(cbits)
        gLL -= np.log(np.real(BM.evaluateProjector()))*weights[i]

    BM.setProjector([])
    aves = BM.computeLocalAverages()
    del BM

    return gLL, aves


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def progtracker(inter):
    global gradEval 
    gradEval = True
    
    return 0


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def LL(inter, *args):
    (Ns, beta, bonds, data, weights, cmode) = args   
    Nd = len(data)
    Nb = len(bonds) 
    
    # initialize Boltzmann machine 
    Hs, Ds, Js = unpackInter(inter, bonds, Ns, cmode)
    kwargs = {'X': Ds, 'Z1': Hs, 'Z2': Js}
    BM = bmachine.BoltzmannMachine(Ns, beta, **kwargs)
        
    # compute the log-likelihood
    vLL = 0
    for i,cbits in enumerate(data[:-1]):
        BM.setProjector(cbits)
        vLL -= np.log(np.real(BM.evaluateProjector())) * weights[i]
 
    del BM
    return vLL

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def grad(inter, *args):
    (Ns, beta, bonds, data, weights, cmode) = args   
    Nd = len(data)
    Nb = len(bonds) 
    
    # initialize Boltzmann machine 
    Hs, Ds, Js = unpackInter(inter, bonds, Ns, cmode)
    kwargs = {'X': Ds, 'Z1': Hs, 'Z2': Js}
    BM = bmachine.BoltzmannMachine(Ns, beta, **kwargs)
        
    # compute the gradient
    aves  = np.zeros((Nd, Ns+Ns+Nb))
    
    for i, cbits in enumerate(data):
        BM.setProjector(cbits)
        if (cmode in ['class', 'aquant']) and (i!=Nd-1): 
           aves[i, :Ns] = -1.0*(np.array(cbits)*2-1)*beta
           for j,bond in enumerate(bonds): 
               aves[i, 2*Ns+j] = (cbits[bond[0]]*2-1)*(cbits[bond[1]]*2-1)*beta
        else: 
            aves[i,:] = BM.computeLocalAverages()
    
    del BM
    gLL = np.sum((aves*weights),  axis=0)
    
    gLL = packInter(gLL, Ns, cmode)
    return gLL
    
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def LLgrad(inter, *args):
    ''' Train all couplings and fields'''
    (Ns, beta, bonds, data, weights, cmode, totime, noise) = args   
    Nd = len(data)
    Nb = len(bonds) 
    
    # initialize Boltzmann machine 
    if totime: t0 = time.time()
    Hs, Ds, Js = unpackInter(inter, bonds, Ns, cmode)
    kwargs = {'X': Ds, 'Z1': Hs, 'Z2': Js}
    BM = bmachine.BoltzmannMachine(Ns, beta, **kwargs)
    if totime: print "--- BM init: %0.4f"  %(time.time() - t0)
        
    # compute the log-likelihood
    if totime: t0 = time.time()
    vLL = 0
    for i,cbits in enumerate(data[:-1]):
        BM.setProjector(cbits)
        vLL -= np.log(np.real(BM.evaluateProjector())) * weights[i]
    if totime: print "--- LL comp: %0.4f" %(time.time() - t0)
     
    # compute the gradient
    if totime: 
        t0 = time.time()
        ts = []
    aves  = np.zeros((Nd, Ns+Ns+Nb))
    for i, cbits in enumerate(data):
        if totime: t0i = time.time()
        BM.setProjector(cbits)
        if (cmode in ['class', 'aquant']) and (i!=Nd-1): 
           aves[i, :Ns] = -1.0*(np.array(cbits)*2-1)*beta
           for j,bond in enumerate(bonds): 
               aves[i, 2*Ns+j] = (cbits[bond[0]]*2-1)*(cbits[bond[1]]*2-1)*beta
        else: 
            aves[i,:] = BM.computeLocalAverages()
        if totime: ts += [time.time() -t0i]
    if totime: 
        print "--- GR eval: %0.4f" %(time.time()-t0)
        print "---          ",
        for T in ts:
            print "%0.4f, " %T,
    
    gLL = np.sum((aves*weights),  axis=0)
    if noise !=0: gLL = gLL + rm.randn((gLL.shape[0]))*np.abs(gLL)*noise

    # record the progress into a file
    global gradEval
    if gradEval:
        st = '    %16.8E' %vLL
        for meas in np.hstack([Hs[:,-1], Ds[:,-1], Js[:,-1]]): st+= '    %16.8E' %meas
        for meas in aves[-1,:]:                                st+= '    %16.8E' %meas
        st += '\n'
        
        global fname
        f = open(fname, 'a')
        f.write(st)
        f.close()

    gradEval = False
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
    elif mode in ["quant_stat", "aquant"]: 
        (Hs, Js) = np.split(inter, [Ns])
        global delta
        Ds = np.ones(Ns)*delta
    elif mode=="quant_free":
        (Hs, D)  = np.split(inter, [Ns, Ns+1])
        Ds = np.ones(Ns)*D
    
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
    elif mode=="quant_free":       return np.hstack([gHs, np.sum(gDs)])
    elif mode in ["quant_stat", 'aquant']: return np.hstack([gHs, gJs]) 
    elif mode=="quant_all":  return gLL

    
    
        

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def bitfield(n):
    return np.array([1 if digit=='1' else 0 for digit in bin(n)[2:]])

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main(): 

    modes  =  ['class', 'quant', 'quant_all', 'quant_stat','quant_free', 'aquant'] 
    
    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--seed',       help='Seed used to generate data', type=int, default=0)
    parser.add_argument('--data',       help='Data file',                  type=str)
    parser.add_argument('--inter','-I', help='Interactions file',          type=str)
    parser.add_argument('--delta',      help='Transverse field',           type=float)
    parser.add_argument('--beta', '-b', help='Inverse temperature ',       type=float)
    parser.add_argument('--mode',       help='Training mode', choices=modes, default='class')
    parser.add_argument('--time',       help='Benchmark execution time ',  action='store_true', default=False)
    parser.add_argument('--noise',      help='Add noise to the gradient',  type=float, default=0)

    args = vars(parser.parse_args())
    rm.seed(args['seed'])

   
    #  Error handling ---------------------------------------------------------
    if (args['inter'] is None):
       print 'Define interactions'
       return 0

    if ((args['mode'] in ['quant', 'quant_all', 'quant_stat', 'aquant']) and (args['delta'] is None)):
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
    print '--- Entropy : %6.4f'    %entropy
    weights = weights.reshape((Nd,1))

    # General an initial guess for the Hamiltonian -----------------------------
    Hs[:,-1] = (rm.rand(Ns)-0.5)*0.1
    Js[:,-1] = (rm.rand(Nb)-0.5)*0.1
    #Hs[:,-1] = rm.randn(Ns)*np.abs(Hs[:,-1]).max() 
    #Js[:,-1] = rm.randn(Nb)*np.abs(Js[:,-1]).max()
    
    if   args['mode'] == 'class':      iparams = np.hstack([Hs[:,-1],                Js[:,-1]])
    elif args['mode'] == 'quant':      iparams = np.hstack([Hs[:,-1], args['delta'], Js[:,-1]])
    elif args['mode'] == 'quant_all':  iparams = np.hstack([Hs[:,-1], Ds[:,-1],      Js[:,-1]])
    elif args['mode'] in ['quant_stat', 'aquant']: 
                                       iparams = np.hstack([Hs[:,-1],                Js[:,-1]])
    elif args['mode'] == 'quant_free': iparams = np.hstack([Hs[:,-1], args['delta']])

 
    # impose limits on the trained values
    #lims    = []
    #f = 5.0
    #for i,par in enumerate(iparams):
    #    apar = np.abs(par)*f
    #    lims.append((-apar, apar))


    # Creat log file 
    global delta
    if args['delta']!= None: delta=args['delta']
    
    global fname 
    fname = 'train_'
    if args['mode'] == 'class': fname += 'mode-%s_' %args['mode']
    else:                       fname += 'mode-%s_delta-%04.2f_' %(args['mode'], delta)
    fname += 'noise-%04.2f_' %(args['noise'])
    fname += args['data'][5:]
    
    
    # Record the header
    f = open(fname, 'w')
    header = '#%19s' %('LL')
    for i in range(Hs.shape[0]):   header += '%20s' %('H'+str(i)) 
    for i in range(Ds.shape[0]):   header += '%20s' %('D'+str(i)) 
    for bond in Js[:,:2].tolist(): header += '%20s' %('J(%d, %d)' %(bond[0], bond[1]))
    
    for i in range(Hs.shape[0]):   header += '%20s' %('<Z'+str(i)+'>') 
    for i in range(Ds.shape[0]):   header += '%20s' %('<X'+str(i)+'>') 
    for bond in Js[:,:2].tolist(): header += '%20s' %('<ZZ(%d, %d)>' %(bond[0], bond[1]))
    
    header +='\n'
    f.write(header)

    # Compute the log-likelihood of the generating Hamiltonian
    gLL, aves = getLLandAves(Hs, Ds, Js)
    st = '    %16.8E' %gLL
    for meas in np.hstack([Hs[:,-1], Ds[:,-1], Js[:,-1]]): st+= '    %16.8E' %meas
    for meas in aves:                                      st+= '    %16.8E' %meas
    st += '\n'
    f.write(st)
    f.close()


    print '---------Gradient ascent----------'
    # Gradient descent --------------------------------------------------------
    
    
    counter = 0
    fx, fLL, info = LBFGS(LLgrad, x0=iparams, args = (Ns, beta, bonds, data, weights, args['mode'], args['time'], args['noise']), iprint = 1, pgtol=1e-05, factr=1e13/2.2, maxiter=50, callback=progtracker)  
    #fx, fLL, info = LBFGS(LL, x0=iparams, fprime=grad, args = (Ns, beta, bonds, data, weights, args['mode']), iprint = 1, pgtol=1e-08, factr=1e13/2.2/10000.0, maxiter=500, callback=progtracker)  

    print fx
    print info
    #print lims

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()



