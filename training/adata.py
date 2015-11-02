import argparse, collections 
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


    for filename in filenames:
        data = np.loadtxt(filename, dtype='int32')
        Ndata = data.shape[0]
        data  = data.tolist()
        udata = []
        cdata = collections.OrderedDict()
        for i,d in enumerate(data):
            if not(d in udata): udata += [d]; cdata[repr(d)]  = 1
            else:                             cdata[repr(d)] += 1

        weights = np.array(cdata.values())/float(Ndata)
        data    = udata
        fstr    = ('%0.3f '*len(data))
        lsw     = tuple((np.sort(weights)[::-1]/np.amax(weights)).tolist())
        print 'Weights: ', fstr % lsw  #, ' data: ', data
        #print weights 
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()


