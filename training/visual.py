import argparse 
import kevent, ssexyhelp
import numpy     as np
import pylab     as pl

def findKbV(dic, svalue):
    for key, value in dic.iteritems():
        if value==svalue:
            return key
            break

    return None

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('filenames', help='Training output files', nargs='+')

    args = parser.parse_args()
    filenames = args.filenames
    pl.connect('key_press_event',kevent.press)

    ddata = {}
    for filename in filenames:
        rdata    = np.loadtxt(filename) 
        if len(rdata.shape)> 1: pindex = (-1,0)
        else:                   pindex = (0)
        
        #(nr, nc) = rdata.shape 

        fparams = ssexyhelp.getReduceParamMap(filename)
        (alpha, Ns, beta, delta) = (fparams['alpha'], fparams['N'] , fparams['b'], fparams['delta'])
        seed = int(findKbV(fparams, ''))
        if not(delta in ddata.keys()): ddata[delta] = np.zeros(30) 
        else:                          ddata[delta][seed] = rdata[pindex]

    fig = pl.figure(1, figsize=(13,6))
    ax  = pl.subplot(111)
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']

    i = 0
    for delta, data in ddata.iteritems():
        ax.plot(range(len(data)), data,  
                color=colors[i], lw=3, ls='-',
                label = r'$\Delta = %0.2f$' %delta)
        i += 1

    pl.xlabel('Instance')
    pl.ylabel('LL')
    pl.legend(loc = 'best')
    pl.title(r'$\beta=%0.2f \, \alpha=%0.2f $' %(beta, alpha), fontsize = 20)
    pl.tight_layout()
    pl.show()
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()


