# Visualize best trained KL divergence between classical and quantum models. 
# KL-divergence is extracted from the corresponding training files.

import Hfile, bmachine, kevent, ssexyhelp
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
def main(): 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--quant','-Q', help='Quantum training file',   nargs='+')
    parser.add_argument('--class','-C', help='Classical training file', nargs='+')

    args = vars(parser.parse_args())

    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]
    fig = pl.figure(1, figsize=(10,5))
    pl.connect('key_press_event',kevent.press)
    ax  = pl.subplot(111)

    
    # ----------------------------------------------------------------------
    if 'class' in args.keys():
        cdata = np.zeros((5, len(args['class'])//5))
        ssizes = []
        i = -1
        for filename in args['class']:
            data = Hfile.cleanData(np.loadtxt(filename))
            
            cKL             = np.amin(data[:,1])
            fparams         = ssexyhelp.getReduceParamMap(filename)
            (bsize, seed)   = (int(fparams['ssise']), int(filename[:-4].split('_')[-1]))
            
            if not(bsize in ssizes): 
                i += 1
                ssizes += [bsize]
            cdata[seed, i] = cKL  
        
        cKLs = np.sum(cdata, axis = 0)
        ax.plot(ssizes, cKLs, color=colors[1], lw=1, ls='-', label=r'$classical$')

    # ----------------------------------------------------------------------
    if 'quant' in args.keys():
        qdata = np.zeros((5, len(args['quant'])//5))
        ssizes = []
        i = -1
        for filename in args['quant']:
            data = cleanData(np.loadtxt(filename))
            
            qKL             = np.amin(data[:,1])
            fparams         = ssexyhelp.getReduceParamMap(filename)
            (bsize, seed)   = (int(fparams['ssise']), int(filename[:-4].split('_')[-1]))
            
            if not(bsize in ssizes): 
                i += 1
                ssizes += [bsize]
            qdata[seed, i] = qKL  

        qKLs = np.sum(qdata, axis = 0)
        ax.plot(ssizes, qKLs, color=colors[2], lw=1, ls='-', label=r'$quantum$')


    
    pl.xlabel(r'$Training \, sample \, size$')
    pl.ylabel(r'$KL-divergence$')
    lgd = pl.legend(loc = 'best')
    lgd.draggable(state=True)
    lgd.draw_frame(False)
    pl.tight_layout()
    pl.show()
   
    #diff = np.array(cLL) - np.array(qLL)
    #plt.hist(diff, 10, normed=1, alpha=0.75)
    #plt.show()
   
# ------------------------------------------------------------ ----------
# ------------------------------------------------------------ ----------
if __name__ == "__main__": 
    main()

