import numpy        as np
import scipy.sparse as sparse 

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def Embed1SiteOper(site_oper, field, N, site):
    site_oper = site_oper.copy()*field
    site = int(site)
    if N == 1:      return site_oper
    elif site==0:   return sparse.kron(site_oper, sparse.eye(2**(N-1)))
    elif site==N-1: return sparse.kron(sparse.eye(2**(N-1)), site_oper)
    else:           return sparse.kron(sparse.kron(sparse.eye(2**site), site_oper), sparse.eye(2**(N-site-1))) 

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def Embed2SiteOper(site_oper, strength, N, siteA, siteB):
    I = sparse.eye(2)
    siteA, siteB = int(siteA), int(siteB)
    if siteB < siteA: 
       siteA, siteB = siteB, siteA
    
    operA = site_oper.copy()
    operB = site_oper.copy()
    pre  = Embed1SiteOper(operA, 1.0, siteA+1, siteA)
    suff = Embed1SiteOper(operB, strength, N-siteA-1, siteB-siteA-1)
    return sparse.kron(pre, suff)


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def Embed3SiteOper(site_oper, strength, N, siteA, siteB, siteC):
    I = sparse.eye(2)
    siteA, siteB, siteC = sorted([int(siteA), int(siteB), int(siteC)]) 
    operA = site_oper.copy()
    operC = site_oper.copy()
    
    pre =  Embed2SiteOper(operA, 1.0, siteB+1, siteA, siteB)
    suff = Embed1SiteOper(operC, strength, N-siteB-1, siteC-siteB-1)

    return sparse.kron(pre, suff)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def EmbedOper(site_oper, N, bond):
    strength, bond = float(bond[-1]), bond[:-1]
    if len(bond) == 1: return Embed1SiteOper(site_oper, strength, N, *bond)
    if len(bond) == 2: return Embed2SiteOper(site_oper, strength, N, *bond)
    if len(bond) == 3: return Embed3SiteOper(site_oper, strength, N, *bond)
    

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


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def pairFullyConnected(N):

    bonds = []
    for i in range(N-1):
        for j in range(i+1, N):
            bonds += [[i,j]]

    return bonds

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def tripleFullyConnected(N):

    bonds = []
    for i in range(N-2):
        for j in range(i+1, N-1):
            for k in range(j+1, N):
                bonds += [[i, j, k]]

    return bonds


