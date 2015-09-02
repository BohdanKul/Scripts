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
    operB = site_oper
    
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
    bonds  = data[nSz+nSx:, 0:1]
    Inter  = data[nSz+nSx:, 2]
    
    return Zfield, Xfield, Inter, bonds

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def OBC_Line(N, strength):
    inter = np.ones(N-1)*strength

    bonds = []
    for i in range(N-1):
        bonds   += [[i, i+1]] 
    return bonds, inter 

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def PBC_Line(N, strength):
    bonds, inter = OBC_Line(N, strength)
    bonds += [[N-1, 0]]
    inter = np.append(inter, strength)

    return bonds, inter

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def OBC_Rect(x, y, strength):
    if y==1: return OBC_Line(x, strength)

    inter = np.ones((x-1)*(y-1)*2+x-1+y-1)*strength
    bonds = []
    for Y in range(y):
        for X in range(x):
            if (X!=x-1) and (Y!=y-1):
                bonds += [[x*Y+X, x*Y+X+1]]
                bonds += [[x*Y+X, x*(Y+1)+X]]
            elif (X==x-1) and (Y!=y-1): bonds += [[x*Y+X, x*(Y+1)+X]]
            elif (X!=x-1) and (Y==y-1): bonds += [[x*Y+X, x*Y+X+1]]
    return bonds, inter


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def PBC_Rect(x, y, strength):
    if y==1: return PBC_Line(x, strength)

    inter = np.ones(x*y*2) * strength
    
    # bulk connections
    #i = (x-1)*(y-1)*2+x-1+y-1
    bonds, t = OBC_Rect(x, y, strength)

    # vertical PBC
    for X in range(x):
        bonds += [[x*(y-1)+X, X]]
    
    # horizontal PBC
    for Y in range(y):
        bonds += [[x*(Y+1)-1, x*Y]]
    
    return bonds, inter



