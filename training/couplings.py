# Extract and plot one and two body correlations from 
# a data file and quant & class training logs. 
# The correlations of interest for the data file are 
# extracted from the inter file. 

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

    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]
    fig = pl.figure(1, figsize=(13,6))
    pl.connect('key_press_event',kevent.press)
    ax  = pl.subplot(111)

    # ----------------------------------------------------------------------
    Hs, Ds, Js, bonds = Hfile.LoadInters(args['inter'])
    (Ns, Nb) = (len(Hs), len(bonds)) 
    
    # ----------------------------------------------------------------------
    data = np.loadtxt(args['data'], dtype='int32')
    udata = []
    cdata = collections.OrderedDict()
    for i,a in enumerate(data):
        d = a.tolist()
        if not(d in udata): udata += [d]; cdata[repr(d)]  = 1
        else:                             cdata[repr(d)] += 1
    weights = np.array(cdata.values())/float(data.shape[0])
    print weights.shape
    data    = udata
    Nd      = len(data)
    weights = weights[np.newaxis]
    weights = weights.T 
    print weights.shape

    # ----------------------------------------------------------------------
    aves  = np.zeros((Nd, Ns+Ns+Nb))
    for i,cbits in enumerate(data):
        aves[i, :Ns] = -1.0*(np.array(cbits)*2-1)
        for j,bond in enumerate(bonds): 
            aves[i, 2*Ns+j] = (cbits[bond[0]]*2-1)*(cbits[bond[1]]*2-1)
    
    aves = np.sum(aves*weights, axis=0)

    ax.plot(aves, color=colors[0], lw=2, ls='-', label=r'$data$')
   
    
    # ----------------------------------------------------------------------
    #headers = Hfile.getHeaders(args['quant'])
    #tdata   = Hfile.cleanData(np.loadtxt(args['quant']))
    #qaves   = tdata[-1, headers.index('<Z0>'):]     
    #
    #ax.plot(qaves, color=colors[1], lw=2, ls='-', label=r'$quantum$')

    ## ----------------------------------------------------------------------
    #headers = Hfile.getHeaders(args['class'])
    #tdata   = Hfile.cleanData(np.loadtxt(args['class']))
    #caves   = tdata[-1, headers.index('<Z0>'):]     
    #
    #ax.plot(caves, color=colors[2], lw=2, ls='-', label=r'$classical$')

    #lheaders = []
    #for head in headers:
    #    lhead    = r'$%s$' %head
    #    lheaders += [lhead]
    #pl.xticks(range(len(caves)), lheaders[headers.index('<Z0>'):], rotation='vertical')
    pl.xlabel(r'$Averages$')
    pl.ylabel(r'$Strength$')
    lgd = pl.legend(loc = 'best')
    lgd.draggable(state=True)
    pl.tight_layout()
    pl.show()
# ------------------------------------------------------------ ----------
# ------------------------------------------------------------ ----------
if __name__ == "__main__": 
    main()

