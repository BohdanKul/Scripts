import kevent
import argparse
from   matplotlib import pyplot as plt
import numpy as np
import pylab as pl

def main(): 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--quant', '-Q', help='Quantum training file',   nargs='+')
    parser.add_argument('--class', '-C', help='Classical training file', nargs='+')

    args = vars(parser.parse_args())

    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]
    #fig = pl.figure(1, figsize=(10,5))
    f, (ax1, ax2)  = pl.subplots(2)
    pl.connect('key_press_event',kevent.press)

    # ----------------------------------------------------------------------
    cdata = {}
    for filename in args['class']:
        data = np.loadtxt(filename)
        LL      = np.amin(data[1:,0])
        seed = int(filename[:-4].split('_')[-1])
        cdata[seed] = LL

    # ----------------------------------------------------------------------
    qdata = {}
    for filename in args['quant']:
        data = np.loadtxt(filename)
        LL      = np.amin(data[1:,0])
        seed = int(filename[:-4].split('_')[-1])
        qdata[seed] = LL
    
    # ----------------------------------------------------------------------
    cLL   = []
    qLL   = []
    for cseed in cdata.keys():
        if cseed in qdata.keys():
            cLL += [cdata[cseed]]
            qLL += [qdata[cseed]]
    cLL = np.array(cLL)
    qLL = np.array(qLL)
    ax1.scatter(cLL, qLL)
    ax1.set_xlabel(r'$LL_{class}$')
    ax1.set_ylabel(r'$LL_{quant}$')
    
    ax2.scatter(cLL, cLL-qLL)
    ax2.set_xlabel(r'$Classical \, LL$')
    ax2.set_ylabel(r'$LL_{class} - LL_{quant}$')
    
    #lgd = pl.legend(loc = 'best')
    #lgd.draggable(state=True)
    #lgd.draw_frame(False)
    pl.tight_layout()
    pl.show()
 
# ------------------------------------------------------------ ----------
# ------------------------------------------------------------ ----------
if __name__ == "__main__": 
    main()
