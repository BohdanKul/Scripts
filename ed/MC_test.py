import ssexyhelp, MCstat, kevent, Hfile
import argparse
import numpy as np
import pylab as pl

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--MC','-Q', help='Quantum training file',   nargs='+')
    parser.add_argument('--ED','-C', help='Classical training file', nargs='+')
    parser.add_argument('--clamped', help='Clamped simulation', default=False, action='store_true')
    args = vars(parser.parse_args())


    if 'MC' in args.keys():
        skip = 3 # skip estimators we don't care about 
        for filename in args['MC']:
            fparams = ssexyhelp.getReduceParamMap(filename)
            data    = np.loadtxt(filename)
            headers = Hfile.getHeaders(filename)
            if not(args['clamped']): (first, last) = (headers.index('X0'),  headers.index('sX0'))
            else:                    (first, last) = (headers.index('sX0'), len(headers))
            aves    = np.zeros(last-first)
            stds    = np.zeros(last-first)
            for i, c in enumerate(range(first, last)):
                cdata   = data[~np.isnan(data[:,c]),c]
                aves[i] = np.mean(cdata)
                stds[i] = np.amax(MCstat.bin(cdata[np.newaxis].T))
            
    if 'ED' in args.keys():
        for filename in args['ED']:
            fparams = ssexyhelp.getReduceParamMap(filename)
            ED    = np.loadtxt(filename)
    
    print ED
    print aves
    print stds
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]


    fig = pl.figure(1, figsize=(10,5))
    pl.connect('key_press_event',kevent.press)
    ax  = pl.subplot(111)
    ax.plot((ED-aves)/stds, color=colors[0])#, lw=2, m='o', ls='')#, label=r'$data$')
    
    pl.ylabel(r'$(ED-MC)/\Delta_{MC}$')
    pl.xlabel(r'$Averages$')
    lgd = pl.legend(loc = 'best')
    
    lheaders = []
    for head in Hfile.getHeaders(args['ED'][0]):
        lheaders += [r'$%s$' %head]
    pl.xticks(range(len(lheaders)), lheaders, rotation='vertical')

    #lgd.draggable(state=True)
    #lgd.draw_frame(False)
    pl.tight_layout()
    pl.show()
# ---------------------------------------------------------- ----------
# ------------------------------------------------------------ ----------
if __name__ == "__main__": 
    main()
