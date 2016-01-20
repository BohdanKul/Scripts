import numpy    as np
import random   as rm
import bmachine as bm
from pylab    import real
import collections


#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
def EDtoMC(bits):
    ''' 
        Convert between ED and MC representation of Ising degrees of freedom 
        0 ->  1
        1 -> -1
    '''
    return (-np.array(bits)*2+1).tolist()

#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
def MCtoED(bits):
    ''' 
        Convert between MC and ED representation of Ising degrees of freedom 
        1 ->  0
       -1 ->  1
    '''
    return ((-(np.array(bits)-1))//2).tolist()

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
    bonds  = data[nSz+nSx:nSz+nSx+nSzSz, -3:-1]
    bonds  = bonds.astype(int)
    Inter  = data[nSz+nSx:nSz+nSz+nSzSz, -3:]

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

# ----------------------------------------------------------------------
def getHeaders(filename):
    inFile = open(filename,'r');
    inLine = inFile.readlines();
    inLine[0] = inLine[0].rstrip()
    oheaders = inLine[0].split('   ')
    oheaders.pop(0)
    nheaders = []
    for head in oheaders:
        if head!='': nheaders += [head.replace(' ','')]
    inFile.close()

    return nheaders

# ----------------------------------------------------------------------
def cleanData(view):
    if np.isnan(view[-1, 0]) or (view[-2, 0] < view[-1, 0]):
        view = view[:-1,:]

    return view


