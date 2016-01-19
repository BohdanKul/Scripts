# visualize BM training in terms of the convergence of the 
# associated probability distribution to the data distribution 

import Hfile, bmachine, kevent
import argparse, collections
import matplotlib as mpl
import numpy as np
import pylab as pl
from matplotlib import pyplot as plt

def toDecimal(cbits):
    D = 0
    for i, bit in enumerate(cbits):
        if bit==1:
            D += 2**(len(cbits)-i-1)

    return D

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
    parser.add_argument('--train','-T', help='Training file',     type=str)

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
    
    eweights = np.zeros(2**len(udata[0]))
    for i,cbits in enumerate(data):
        d = toDecimal(cbits)
        eweights[d] = weights[i]

    #sortind = [i[0] for i in sorted(enumerate(data), key=HammingDistance1)]
    #sortind  = np.argsort(weights)
    #weights = weights[sortind]
    #data    = (np.array(data)[sortind]).tolist()
    #for v in data:
    #    print v, HammingDistance(v)
    #colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]
    fig = pl.figure(1, figsize=(13,6))
    pl.connect('key_press_event',kevent.press)
    ax  = pl.subplot(111)
    #ax.plot( np.sum(-(data*2-1), axis=0)/(1.0*data.shape[0]))
   
    if ((args['inter'] != None) and (args['train'] != None)):
        # ----------------------------------------------------------------------
        Hs, Ds, Js, bonds = Hfile.LoadInters(args['inter'])
        
        # ----------------------------------------------------------------------
        headers = Hfile.getHeaders(args['train'])
        tdata = Hfile.cleanData(np.loadtxt(args['train']))
        Hs = tdata[:, headers.index('H0')    :headers.index('D0')]     
        Ds = tdata[:, headers.index('D0')    :headers.index('J(0,1)')]     
        Js = tdata[:, headers.index('J(0,1)'):headers.index('<Z0>')]     

        # ----------------------------------------------------------------------
        (Ns, Nb) = (Hs.shape[1], len(bonds)) 

        norm = mpl.colors.Normalize(vmin=0, vmax=tdata.shape[0])
        c_m = mpl.cm.cool
        s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        s_m.set_array([])

        for i in range(tdata.shape[0]):
            ds = np.transpose(np.vstack((np.arange(Ns), Ds[i,:])))
            hs = np.transpose(np.vstack((np.arange(Ns), Hs[i,:])))
            js = np.hstack((bonds, Js[i,:].reshape(Nb,1)))

            #Ps = []
            #kwargs = {'X': ds, 'Z1': hs, 'Z2': js}
            #BM = bmachine.BoltzmannMachine(Ns, 1.0, **kwargs)
            #for cbits in data:
            #    BM.setProjector(cbits)
            #    Ps += [BM.evaluateProjector()]
            #del BM
            #ax.plot(-weights*np.log(Ps), color=s_m.to_rgba(i), lw=1, ls='-')#            Ps = []
    
            kwargs = {'X': ds, 'Z1': hs, 'Z2': js}
            BM = bmachine.BoltzmannMachine(Ns, 1.0, **kwargs)
            Ps = np.zeros(2**len(udata[0]))
            for cbits in data:
                BM.setProjector(cbits)
                d = toDecimal(cbits)
                Ps[d] = BM.evaluateProjector()
            del BM
            ax.plot(-eweights*np.log(Ps), color=s_m.to_rgba(i), lw=1, ls='-')#, label = r'$epoch \, %d$' %i)

        plt.colorbar(s_m)
    ax.plot(-eweights*np.log(eweights), color='black', lw=2, ls='-', label = r'$data$')
    #ax.plot(-weights*np.log(weights), color='black', lw=2, ls='-', label = r'$data$')
    pl.xlabel(r'$Data$')
    pl.ylabel(r'$P_{data}\log P_{model}$')
    lgd = pl.legend(loc = 'best')
    lgd.draggable(state=True)
    pl.tight_layout()
    pl.show()
# ------------------------------------------------------------ ----------
# ------------------------------------------------------------ ----------
if __name__ == "__main__": 
    main()

