import numpy    as np
import random   as rm
import bmachine as bm
from pylab    import real
import collections

#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
def bitfield(n):                                                                   
    return np.array([1 if digit=='1' else 0 for digit in bin(n)[2:]]) 

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def LoadInters(fname):  
    header = open(fname, 'r').readlines()[1]
    (nSz, nSx, nSzSz) = header[1:].split()[0:3]
    (nSz, nSx, nSzSz) = (int(nSz), int(nSx), int(nSzSz))

    data = np.loadtxt(fname, skiprows=2)
    Zfield = data[:nSz,-2:]
    Xfield = data[nSz:nSz+nSx, -2:]
    bonds  = data[nSz+nSx:, 0:2]
    bonds  = bonds.astype(int)
    Inter  = data[nSz+nSx:, :]

    return Zfield, Xfield, Inter,  bonds


#------------------------------------------------------------------------------
# Acquire data
#------------------------------------------------------------------------------
def GetData(datafile, N, seed, Nsamples = 0, bonds = [], Z = [], X = [], ZZ = [], beta = 1.0): 
    # Either from a pre-generated set
    if (datafile is not None):
        data = np.loadtxt(args['data'])
        Nclamped = data.shape[1]
        if Nclamped>N: 
            print 'Training set vectors exceed the graph size'
            return 0
        Nsamples = data.shape[0]

    # Or generate it from a  Hamiltonian at a given beta with help of ED
    else:
        Nclamped = N
        if Nclamped>N: 
            print 'Training set vectors exceed the graph size'
            return 0
        BM = bm.BoltzmannMachine(N, bonds, Z, X, ZZ , beta)
        probTable = np.zeros(2**Nclamped)
        for i in range(2**Nclamped):
            cbits = bitfield(i)                              # convert i to a list of bits
            cbits = [0]*(Nclamped-len(cbits))+cbits.tolist() # keep the list length constant   
            BM.setProjector(cbits)
            if i==0: probTable[i] = real(BM.evaluateProjector())
            else:    probTable[i] = probTable[i-1] + real(BM.evaluateProjector())
       
        rm.seed(seed)
        data = []
        for i in range(Nsamples):
            RN = rm.random()
            index = np.searchsorted(probTable, RN)
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
    weights = np.array(cdata.values())/float(Nsamples)
    data    = udata


    return data, weights


