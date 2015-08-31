from qutip import *
import numpy as np

N = 8
Bonds  = []
Lfield = []
Tfield = []
Inter  = []

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def EmbedSiteOper(site_oper, field, site, N):
    I = qeye(2)
    site_oper *= field

    if N == 1:      return site_oper
    elif site==0:   return tensor(site_oper, tensor([I]*(N-1)))
    elif site==N-1: return tensor(tensor([I]*(N-1)), site_oper)
    else:           return tensor(tensor([I]*site), site_oper, tensor([I]*(N-site-1))) 

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def Embed2SiteOper(site_oper, strength, siteA, siteB, N):
    I = qeye(2)
    if siteB < siteA: 
       siteA, siteB = siteB, siteA
    
    operA = site_oper
    operB = site_oper*strength
    
    pre  = EmbedSiteOper(operA, 1.0, siteA, siteA+1)
    suff = EmbedSiteOper(operB, strength, siteB-siteA-1, N-siteA-1)

    return tensor(pre, suff)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def LoadInters(fname):  
    header = open(fname, 'r').readlines()[1]
    (nSz, nSx, nSzSz) = header[1:].split()
    (nSz, nSx, nSzSz) = (int(nSz), int(nSx), int(nSzSz))

    data = np.loadtxt(fname, skiprows=2)
    Zfield = data[:nSz,2]
    Xfield = data[nSz:nSz+nSx, 2]
    Inter  = data[nSz+nSx:, :]
    
    return Zfield, Xfield, Inter

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def OBC_Line(N, strength):
    inter = np.zeros((N-1, 3))

    for i in range(N-1):
        inter[i] = [i, i+1, strength]
    
    return inter

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def PBC_Line(N, strength):
    inter = np.zeros((N, 3))
    inter[:N-1, :] = OBC_Line(N, strength)
    inter[N-1] = [N-1, 0, strength]

    return inter

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def OBC_Rect(x, y, strength):
    if y==1: return OBC_Line(x, strength)

    inter = np.zeros( ((x-1)*(y-1)*2+x-1+y-1, 3) )
    i = 0
    for Y in range(y):
        for X in range(x):
            if (X!=x-1) and (Y!=y-1):
                inter[i] = [x*Y+X, x*Y+X+1,   strength]
                inter[i+1] = [x*Y+X, x*(Y+1)+X, strength]
                i += 2
            elif (X==x-1) and (Y!=y-1):
                inter[i] = [x*Y+X, x*(Y+1)+X, strength]
                i += 1
            elif (X!=x-1) and (Y==y-1):
                inter[i] = [x*Y+X, x*Y+X+1,   strength]
                i += 1
    return inter


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def PBC_Rect(x, y, strength):
    if y==1: return PBC_Line(x, strength)

    inter = np.zeros( (x*y*2, 3) ) 
    
    # bulk connections
    i = (x-1)*(y-1)*2+x-1+y-1
    inter[: i, : ] = OBC_Rect(x, y, strength)
    
    # vertical PBC
    for X in range(x):
        inter[i] = [x*(y-1)+X, X,   strength]
        i += 1
    
    # horizontal PBC
    for Y in range(y):
        inter[i] = [x*(Y+1)-1, x*Y,   strength]
        i += 1
    
    return inter



