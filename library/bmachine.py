import numpy         as np
import scipy.sparse  as sparse
import scipy.linalg  as slin
import Hbuilder
import mexp                
from   scipy.integrate import simps, trapz

class BoltzmannMachine:
    def __init__(self, N, beta, **kwargs):
        self.Ns = N
        self.Nb = len(kwargs['Z2'])
        I = np.eye(2)
        X = np.array([[0, 1], [1, 0]])
        Z = np.array([[1, 0], [0, -1]])
        self.H = np.zeros((2**N, 2**N))
        #I = sparse.eye(2, format="csc")
        #X = sparse.csc_matrix(np.array([[0, 1], [1, 0]]))
        #Z = sparse.csc_matrix(np.array([[1, 0], [0, -1]]))
        #self.H = sparse.csc_matrix((2**N, 2**N))
       
        for itype, bonds in kwargs.iteritems():
            if itype == 'X': oper = X 
            else:            oper = Z
            for bond in bonds:
                self.H = self.H + Hbuilder.EmbedOper(oper, N, bond)  
        
        self.beta = beta
        
        #self.rho  = np.matrix(slin.expm(-self.H*beta), copy=False)
        #self.H = np.matrix(self.H)
        #self.rho  = mexp.expm(-self.H*beta)
        
        self.rho  = np.matrix(mexp.expm(-self.H*beta), copy=False)
        
        self.rhoN = self.rho/np.sum(self.rho.diagonal())
        #print 'X1: ', kwargs['X'][:,-1]
        #print 'Z1: ', kwargs['Z1'][:,-1]
        #print 'Z2: ', kwargs['Z2'][:,-1]
        #print "max H: ",   np.amax(np.abs(np.array([self.H.max(), self.H.min()]))),\
        #      "max rho: ", np.amax(np.abs(np.array([self.rho.max(), self.rho.min()]))),\
        #      "Z: ",       np.sum(self.rho.diagonal())
        self.operZs = []
        self.operXs = []
        for site in range(N):
            self.operZs += [Hbuilder.EmbedOper(Z, N, [site, 1.0], True)]
            self.operXs += [Hbuilder.EmbedOper(X, N, [site, 1.0], True)]
            
        self.operZZs = []
        for bond in kwargs['Z2'].copy():
            bond[-1] = 1.0 
            self.operZZs += [Hbuilder.EmbedOper(Z, N, bond, True)]
        
        #print
        #print self.H
        #print
        #print self.rho
    def setProjector(self, cbits):
        self.clamped = False
        if  len(cbits)>0:
            self.clamped = True
            ket = sparse.csc_matrix(np.array((1-cbits[0], cbits[0])))
            for qbit in cbits[1:]:
                ket = sparse.kron(ket, sparse.csc_matrix(np.array((1-qbit, qbit))))
            self.P = sparse.kron(ket.transpose(), ket)

        if (len(cbits)>0) and (len(cbits)!=self.Ns): self.P = sparse.kron(self.P, sparse.qye(2**(self.N-len(cbits))))
        elif (len(cbits)==0):                        self.P = sparse.eye(2**self.Ns, format="csc")

    def evaluateProjector(self):
        #print self.rhoN.diagonal()
        #print self.rhoN*self.P.todense()
        #print np.sum((self.rhoN*self.P).diagonal())
        #print 
        return np.sum((self.rhoN*self.P).diagonal())
    
    def computeLocalAverages(self, test = False):
        U  = self.rho*self.P
        U /= np.sum(U.diagonal())
        if  ((not self.clamped) or (test)):
            aves = np.zeros(self.Ns+self.Ns+self.Nb) 
            for site in range(self.Ns):
                aves[site] = np.real(np.sum(
                                            (U*self.operZs[site]).diagonal()
                                           )
                                    )
            
                aves[self.Ns+site] += np.real(np.sum(
                                                    (U*self.operXs[site]).diagonal()
                                                    )
                                             )

            for i in range(self.Nb):
                aves[2*self.Ns+i] += np.real(np.sum(
                                                   (U*self.operZZs[i]).diagonal()
                                                   )
                                            )
            return aves
        else:
            Nsamples = 8 
            eZ  = np.zeros((Nsamples+1,self.Ns))
            eX  = np.zeros((Nsamples+1,self.Ns)) 
            eZZ = np.zeros((Nsamples+1,self.Nb))
            # Evaluate time-evolved operators on a discrete grid
            tau = 1.0*self.beta/(2.0*Nsamples)
            #Ub0 = np.matrix(slin.expm( self.H*tau), copy=False)
            #Uf0 = np.matrix(slin.expm(-self.H*tau), copy=False)
            Ub0 = np.matrix(mexp.expm( self.H*tau), copy=False)
            Uf0 = np.matrix(mexp.expm(-self.H*tau), copy=False)
            for step in range(Nsamples+1):
                if  step>0: U = Ub0*U*Uf0

                for site in range(self.Ns):
                    eZ[step][site] = np.real(np.sum(
                                                    (self.operZs[site]*U).diagonal()
                                                   )
                                            )
                    eX[step][site] = np.real(np.sum(
                                                    (self.operXs[site]*U).diagonal()
                                                    )
                                            )
                
                for i in range(self.Nb):
                    eZZ[step][i] = np.real(np.sum(
                                                  (self.operZZs[i]*U).diagonal()
                                                 )
                                          )
        
            # Numerically integrate 
            UE = np.zeros(self.Ns+self.Ns+self.Nb) 
            xs = np.linspace(0.0, Nsamples+1-1, Nsamples+1)*tau
            for site in range(self.Ns):
                UE[site]            = trapz(eZ[:,site], xs)*2.0
                UE[self.Ns + site] += trapz(eX[:,site], xs)*2.0
            
            for bond in range(self.Nb):
                UE[2*self.Ns+bond] += [trapz(eZZ[:,bond], xs)*2.0]
         
            return UE 
        
