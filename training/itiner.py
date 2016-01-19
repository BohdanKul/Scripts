# Visualize the training path on a classical vs quantum energy plot

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
    
        inFile = open(filename,'r');
        inLine = inFile.readlines();
        headers = inLine[0].split('   ')
        headers.pop(0)
        header = []
        for head in headers:
            if head!='': header += [head]
        inFile.close()

        split1 = header.index(' <Z0>')
        split2 = header.index(' <X0>')
        split3 = header.index(' <ZZ(0, 1)>')

        fparams = ssexyhelp.getReduceParamMap(filename)
        #(alpha, Ns, beta, mode) = (fparams['alpha'], fparams['N'] , fparams['b'], fparams['mode'])
        (mode) = (fparams['mode']+str(fparams['modes'])+str(fparams['iw']))
       
        pindex = len(rdata[:, 0])
        if np.isnan(rdata[-1,0]) or rdata[-2,0]<rdata[-1,0]: 
           print "Weird LL: ", filename
           pindex += -1
        
        if not (mode in ddata.keys()): 
           ddata[mode] = {} 

        #seed = int(findKbV(fparams, ''))
        seed = int(fparams['hseed'])
        #seed = int(findKbV(fparams, ''))
        #print rdata[:pindex, 0]
        #print ddata[mode][seed, 0, :]
        #print 'pindex', pindex,rdata[: pindex, 1:split1]
        Es = rdata[: pindex, 1:split1]*rdata[:pindex, split1:]
        split2 -= split1
        split3 -= split1
        split1  = 0
        ddata[mode][seed] = {}
        ddata[mode][seed][0] = rdata[:pindex,0]
        ddata[mode][seed][1] = np.sum(Es[:, split1:split2], axis=1)\
                        + np.sum(Es[:, split3:], axis=1)
        ddata[mode][seed][2] = np.sum(Es[:, split2:split3], axis=1) 

        
        #if not(delta in ddata.keys()): ddata[delta] = np.zeros(30) 
        #else:                          ddata[delta][seed] = rdata[pindex]

    fig = pl.figure(1, figsize=(13,6))
    ax  = pl.subplot(111)
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']

    i = 0
    print ddata
    for label, data in ddata.iteritems():
        if label == 'quant': col = colors[0] 
        else:                col = colors[1] 
        for seed,sdata in data.iteritems():
            ax.plot(sdata[1][0], sdata[2][0], marker='.')
            ax.plot(sdata[1], sdata[2],
                    ls = ['-','--',':'][(seed+2)%3])

                    #color=col, lw=3, ls='-')
                    #, label = r'$mode = %s$' %label)

    pl.xlabel(r'$E_C$')
    pl.ylabel(r'$E_Q$')
    pl.legend(loc = 'best')
    #pl.title(r'$\beta=%0.2f \, \alpha=%0.2f $' %(beta, alpha), fontsize = 20)
    pl.tight_layout()
    pl.show()
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()


