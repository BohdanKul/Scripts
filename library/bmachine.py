import numpy               as np
import scipy.sparse        as sparse
import scipy.sparse.linalg as sl
import Hbuilder

from   scipy.integrate import simps, trapz

class BoltzmannMachine:
    def __init__(self, N, beta, **kwargs):
        self.Ns = N
        self.Nb = len(kwargs['Z2'])
        I = sparse.eye(2)
        X = sparse.csc_matrix(np.array([[0, 1], [1, 0]]))
        Z = sparse.csc_matrix(np.array([[1, 0], [0, -1]]))
        self.H = sparse.eye(2**N) - sparse.eye(2**N)
       
        for itype, bonds in kwargs.iteritems():
            if itype == 'X': oper = X 
            else:            oper = Z
            for bond in bonds:
                #print itype, bond
                #print Hbuilder.EmbedOper(oper, N, bond).todense()  
                self.H = self.H + Hbuilder.EmbedOper(oper, N, bond)  
        
        self.beta = beta

        self.rho  = sl.expm(-self.H*beta)

        self.operZs = []
        self.operXs = []
        for site in range(N):
            self.operZs += [Hbuilder.EmbedOper(Z, N, [site, 1.0])]
            self.operXs += [Hbuilder.EmbedOper(X, N, [site, 1.0])]
            
        self.operZZs = []
        for bond in kwargs['Z2'].copy():
            bond[2] = 1.0 
            self.operZZs += [Hbuilder.EmbedOper(Z, N, bond)]
        #print 'Hamiltonian:'
        #print kwargs
        #print self.H.todense()

    def setProjector(self, cbits):
        self.clamped = False
        if  len(cbits)>0:
            self.clamped = True
            ket = sparse.csc_matrix(np.array((1-cbits[0], cbits[0])))
            for qbit in cbits[1:]:
                ket = sparse.kron(ket, sparse.csc_matrix(np.array((1-qbit, qbit))))
            self.P = sparse.kron(ket.transpose(), ket)

        if (len(cbits)>0) and (len(cbits)!=self.Ns): self.P = sparse.kron(self.P, sparse.qye(2**(self.N-len(cbits))))
        elif (len(cbits)==0):                        self.P = sparse.eye(2**self.Ns)

    def evaluateProjector(self):
        return np.sum((self.rho.multiply(self.P)).diagonal())/np.sum(self.rho.diagonal())
    
    def computeLocalAverages(self, test = False):
        U  = self.rho.multiply(self.P)
        Z  = np.sum(U.diagonal())
        U /= Z
        if  ((not self.clamped) or (test)):
            aves = np.zeros(self.Ns+self.Ns+self.Nb) 
            for site in range(self.Ns):
                aves[site] = np.real(np.sum(
                                            (U.multiply(self.operZs[site])).diagonal()
                                           )
                                    )
            
                aves[self.Ns+site] += np.real(np.sum(
                                                    (U.multiply(self.operXs[site])).diagonal()
                                                    )
                                             )

            for i in range(self.Nb):
                aves[2*self.Ns+i] += np.real(np.sum(
                                                   (U.multiply(self.operZZs[i])).diagonal()
                                                   )
                                            )
            return aves
        else:
            Nsamples = 10 
            eZ  = np.zeros((Nsamples+1,self.Ns))
            eX  = np.zeros((Nsamples+1,self.Ns)) 
            eZZ = np.zeros((Nsamples+1,self.Nb))
            # Evaluate time-evolved operators on a discrete grid
            tau = 1.0*self.beta/2.0/(Nsamples)
            Ub0 = sl.expm( self.H*tau)
            Uf0 = sl.expm(-self.H*tau)
            for step in range(Nsamples+1):
                if  step>0: U = Ub0.multiply(U.multiply(Uf0))

                for site in range(self.Ns):
                    eZ[step][site] = np.real(np.sum(
                                                    (self.operZs[site].multiply(U)).diagonal()
                                                   )
                                            )
                    #eX[step][site] = np.real((self.operXs[site]*U).tr())
                
                for i in range(self.Nb):
                    eZZ[step][i] = np.real(np.sum(
                                                  (self.operZZs[i].multiply(U)).diagonal()
                                                 )
                                          )
        
            # Numerically integrate 
            UE = np.zeros(self.Ns+self.Ns+self.Nb) 
            xs  = np.linspace(0.0, Nsamples+1-1, Nsamples+1)*tau
            for site in range(self.Ns):
                UE[site]            = trapz(eZ[:,site], xs)*2.0
                #UE[self.Ns + site] += trapz(eX[:,site], xs)*2.0
            
            for bond in range(self.Nb):
                UE[2*self.Ns+bond] += [trapz(eZZ[:,bond], xs)*2.0]
         
            return UE 
        
