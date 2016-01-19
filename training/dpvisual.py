# Visually compare learnt probability distributions between 
# classical, quantum  BMs and data

import Hfile, bmachine, kevent
import argparse, collections
import matplotlib as mpl
import numpy as np
import pylab as pl
from matplotlib import pyplot as plt

# ----------------------------------------------------------------------
def HammingDistance1(v1, v2 = [1,0,1,1,0,1,1,1,0,1]):
    v1 = v1[1]
    d = 0
    ones1 = 0
    ones2 = 0
    sign  = 0
    for c1, c2 in zip(v1, v2):
        if c1!=c2: 
           d += 1
           if (sign==0) and (c1==1): sign = +1
           if (sign==0) and (c2==1): sign = -1

    return d*sign

# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
def HammingDistance(v1, v2 = [1,0,1,1,0,1,1,1,0,1]):
    d = 0
    ones1 = 0
    ones2 = 0
    sign  = 0
    for c1, c2 in zip(v1, v2):
        if c1!=c2: 
           d += 1
           if (sign==0) and (c1==1): sign = +1
           if (sign==0) and (c2==1): sign = -1

    return d*sign

# ----------------------------------------------------------------------
#def sortMultiModal(weights, data, Nmodes):
#
#    wind = np.argsort(weights)
#    for mode in data[wind][-Nmodes,:]
#        
#    for i,v in enumerate(np.array(data)[wind][-10:]):
#        print v, HammingDistance((i,v)), data.index(v.tolist())

# ----------------------------------------------------------------------
def main(): 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--data', '-D', help='Data file',         type=str)
    parser.add_argument('--inter','-I', help='Interactions file', type=str)
    parser.add_argument('--quant','-Q', help='Quantum training file',     type=str)
    parser.add_argument('--class','-C', help='Classical training file',     type=str)

    args = vars(parser.parse_args())

    # ----------------------------------------------------------------------
    data = np.loadtxt(args['data'], dtype='int32')
    udata = []
    cdata = collections.OrderedDict()
    for i,a in enumerate(data):
        d = a.tolist()
        if not(d in udata): udata += [d]; cdata[repr(d)]  = 1
        else:                             cdata[repr(d)] += 1
    weights = np.array(cdata.values())/float(data.shape[0])
    data    = udata
    
    #sortind = [i[0] for i in sorted(enumerate(data), key=HammingDistance1)]
    #sortind  = np.argsort(weights)
    #weights = weights[sortind]
    #data    = (np.array(data)[sortind]).tolist()
    #for v in data:
    #    print v, HammingDistance(v)
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]
    fig = pl.figure(1, figsize=(13,6))
    pl.connect('key_press_event',kevent.press)
    ax  = pl.subplot(111)
    #ax.plot( np.sum(-(data*2-1), axis=0)/(1.0*data.shape[0]))
    
    # ----------------------------------------------------------------------
    Hs, Ds, Js, bonds = Hfile.LoadInters(args['inter'])
    (Ns, Nb) = (len(Hs), len(bonds)) 
    
    # ----------------------------------------------------------------------
    headers = Hfile.getHeaders(args['quant'])
    tdata = Hfile.cleanData(np.loadtxt(args['quant']))
    hs = tdata[-1, headers.index('H0')    :headers.index('D0')]
    ds = tdata[-1, headers.index('D0')    :headers.index('J(0,1)')]     
    js = tdata[-1, headers.index('J(0,1)'):headers.index('<Z0>')]     
    ds = np.transpose(np.vstack((np.arange(Ns), ds)))
    hs = np.transpose(np.vstack((np.arange(Ns), hs)))
    js = np.hstack((bonds, js.reshape(Nb,1)))

    # ----------------------------------------------------------------------

    qPs = []
    kwargs = {'X': ds, 'Z1': hs, 'Z2': js}
    BM = bmachine.BoltzmannMachine(Ns, 1.0, **kwargs)
    for cbits in data:
        BM.setProjector(cbits)
        qPs += [BM.evaluateProjector()]
    del BM
    qPs = np.array(qPs)
    #ax.plot((-np.log(qPs)+np.log(weights)), color=colors[0], lw=1, ls='-', label = r'$quantum$')
    ax.plot((-weights*(np.log(qPs)+np.log(weights))), color=colors[0], lw=1, ls='-', label = r'$quantum$')

    # ----------------------------------------------------------------------
    headers = Hfile.getHeaders(args['class'])
    tdata = Hfile.cleanData(np.loadtxt(args['class']))
    hs = tdata[-1, headers.index('H0')    :headers.index('D0')]
    ds = tdata[-1, headers.index('D0')    :headers.index('J(0,1)')]     
    js = tdata[-1, headers.index('J(0,1)'):headers.index('<Z0>')]     
    ds = np.transpose(np.vstack((np.arange(Ns), ds)))
    hs = np.transpose(np.vstack((np.arange(Ns), hs)))
    js = np.hstack((bonds, js.reshape(Nb,1)))

    cPs = []
    kwargs = {'X': ds, 'Z1': hs, 'Z2': js}
    BM = bmachine.BoltzmannMachine(Ns, 1.0, **kwargs)
    for cbits in data:
        BM.setProjector(cbits)
        cPs += [BM.evaluateProjector()]
    del BM
    cPs = np.array(cPs)
    #ax.plot((-np.log(cPs)+np.log(weights)), color=colors[1], lw=1, ls='-', label = r'$classical$')
    ax.plot((-weights*(np.log(cPs)+np.log(weights))), color=colors[1], lw=1, ls='-', label = r'$classical$')

    ax.plot([0, len(weights)],[0,0], color='black', lw=2, ls='-')
    #ax.plot(-weights*np.log(weights), color='black', lw=2, ls='-', label = r'$data$')
    pl.xlabel('Vector')
    pl.ylabel('P')
    lgd = pl.legend(loc = 'best')
    lgd.draggable(state=True)
    pl.tight_layout()
    pl.show()
# ------------------------------------------------------------ ----------
# ------------------------------------------------------------ ----------
if __name__ == "__main__": 
    main()

