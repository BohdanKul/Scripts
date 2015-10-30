import time
import numpy               as np
import numpy.random        as rm
import scipy.linalg as sl
import mexp as m
import scipy.sparse as sparse
import Hbuilder
#import scipy.sparse.linalg as sl


def bsexp(M):
    e = np.ceil(np.log(sl.norm(M,np.inf)))
    r = int(max(0, e+1))
    M = M/(2.0**r)
    FM = np.eye(M.shape[0])
    F  = FM + M
    (l, n, incF) = (1,1,1)
    while incF>1e-16:
        n += 1
        l += 1
        FM = np.dot(FM, M)
        F  = F + np.dot(FM, M/float(l))
        incF = np.sum(FM)

    for k in range(r):
        F = np.dot(F,F)

    return F



rm.seed(0)
Ns    = 2 
Niter = 1 
times1 = np.zeros(Niter)
times2 = np.zeros(Niter)
X = sparse.csc_matrix(np.array([[0, 1], [1, 0]]))
Z = sparse.csc_matrix(np.array([[1, 0], [0, -1]]))
for i in range(Niter):
    H = np.matrix((rm.random_integers(0,1, 2**Ns * 2**Ns)*2 - 1).reshape((2**Ns, 2**Ns)))
    
    #kwargs = {}
    #H = sparse.csc_matrix((2**Ns,2**Ns))
    #print H.shape   
    #sites = np.arange(Ns)
    #sites = sites.reshape((Ns,1))
    #Ds    = np.ones(Ns)
    #Ds    = Ds.reshape((Ns,1))
    #Ds    = np.hstack((sites, Ds))
    #

    #Hs    = rm.random_integers(-1, 1, Ns)
    #Hs    = Hs.reshape((Ns,1))
    #Hs    = np.hstack((sites, Hs))

    #bonds = Hbuilder.pairFullyConnected(Ns)
    #Js  = rm.random_integers( 0, 1, len(bonds))*2-1
    #Js  = Js.reshape((len(bonds),1))
    #Js  = np.hstack((np.array(bonds), Js ))
    #
    #kwargs['X']  = Ds 
    #kwargs['Z1'] = Hs
    #kwargs['Z2'] = Js 

    #for itype, bonds in kwargs.iteritems():
    #    if itype == 'X': oper = X 
    #    else:            oper = Z
    #    for bond in bonds:
    #        H = H + Hbuilder.EmbedOper(oper, Ns, bond)  
 
    #H = H.todense()
    
    t0   = time.time()
    rho1 = m.expm(H)
    times1[i] = time.time() - t0
    print 'new : %0.4f ' %times1[i]

    t0   = time.time()
    rho2 = sl.expm(H)
    times2[i] = time.time() - t0
    print 'old: %0.4f ' %times2[i]
  
    #print rho1.todense()
    #print
    #print rho2.todense()
    #rho2 = rho2.todense()
    #rho1 = rho1.todense()
    print rho2
    print rho1
    print 'diff: ', np.amax(np.abs(rho2-rho1))/np.amax(np.abs(rho2))

print 'New mean: ', np.mean(times1), '+/-', np.std(times1)
print 'Old mean: ', np.mean(times2), '+/-', np.std(times2)



