import numpy         as np
import scipy.sparse  as sparse
import scipy.linalg  as slin
import Hbuilder
import mexp
from   scipy.integrate import simps, trapz

def Bits2Int(bitlist):
    out = 0
    for bit in bitlist:
        out = (out << 1) | bit

    return out

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

        for otype in ['X', 'Z1', 'Z2']:
            if otype == 'X': oper = X
            else:            oper = Z
            for bond in kwargs[otype]:
                self.H = self.H + Hbuilder.EmbedOper(oper, N, bond)

        self.beta = beta

        #self.rho  = np.matrix(slin.expm(-self.H*beta), copy=False)
        #self.H = np.matrix(self.H)
        #self.rho  = mexp.expm(-self.H*beta)

        self.rho  = np.matrix(mexp.expm(-self.H*beta), copy=False)

        #self.rhoN = self.rho/np.sum(self.rho.diagonal())
        self.rhoN = np.diagonal(self.rho).copy()
        self.rhoN /= np.sum(self.rhoN)

        #print 'X1: ', kwargs['X'][:,-1]
        #print 'Z1: ', kwargs['Z1'][:,-1]
        #print 'Z2: ', kwargs['Z2'][:,-1]
        #print "max H: ",   np.amax(np.abs(np.array([self.H.max(), self.H.min()]))),\
        #      "max rho: ", np.amax(np.abs(np.array([self.rho.max(), self.rho.min()]))),\
        #      "Z: ",       np.sum(self.rho.diagonal())
        self.oper1Zs = []
        self.operXs = []
        for site in range(N):
            self.oper1Zs += [Hbuilder.EmbedOper(Z, N, [site, 1.0], True).diagonal()]
            self.operXs += [sparse.csr_matrix(
                                    Hbuilder.EmbedOper(X, N, [site, 1.0], True)
                                ).nonzero()[1]]
            #self.oper1Zs += [Hbuilder.EmbedOper(Z, N, [site, 1.0], True)]
            #self.operXs += [Hbuilder.EmbedOper(X, N, [site, 1.0], True)]

        self.oper2Zs = []
        for bond in kwargs['Z2'].copy():
            bond[-1] = 1.0
            self.oper2Zs += [Hbuilder.EmbedOper(Z, N, bond, True).diagonal()]
            #self.oper2Zs += [Hbuilder.EmbedOper(Z, N, bond, True)]

        self.oper3Zs = []
        if 'tribody' in kwargs.keys():
            for bond in list(kwargs['tribody']):
                bond += [1.0]
                self.oper3Zs += [Hbuilder.EmbedOper(Z, N, bond, True).diagonal()]
        self.N3b = len(self.oper3Zs)

        self.oper4Zs = []
        if 'quadbody' in kwargs.keys():
            for bond in list(kwargs['quadbody']):
                bond += [1.0]
                self.oper4Zs += [Hbuilder.EmbedOper(Z, N, bond, True).diagonal()]
        self.N4b = len(self.oper4Zs)

    def setProjector(self, cbits):
        self.bindex = Bits2Int(cbits)
        self.clamped = False
        if  len(cbits)>0:
            self.clamped = True
        #    ket = sparse.csc_matrix(np.array((1-cbits[0], cbits[0])))
        #    for qbit in cbits[1:]:
        #        ket = sparse.kron(ket, sparse.csc_matrix(np.array((1-qbit, qbit))))
        #    self.P = sparse.kron(ket.transpose(), ket)

        #if (len(cbits)>0) and (len(cbits)!=self.Ns): self.P = sparse.kron(self.P, sparse.qye(2**(self.N-len(cbits))))
        #elif (len(cbits)==0):                        self.P = sparse.eye(2**self.Ns, format="csc")

    def evaluateProjector(self):
        #return np.sum((self.rhoN*self.P).diagonal())
        return self.rhoN[self.bindex]

    def computeLocalAverages(self, test = False):
        #U  = self.rho*self.P
        #U /= np.sum(U.diagonal())

        if  ((not self.clamped) or (test)):
            U = self.rho/np.sum(self.rho.diagonal())
            aves = np.zeros(self.Ns + self.Ns + self.Nb + self.N3b + self.N4b)

            # Z1 operator
            for site in range(self.Ns): aves[site] = np.sum(self.oper1Zs[site]*U.diagonal().T)
                #aves[site] = np.real(np.sum((U*self.oper1Zs[site]).diagonal()))

            # X1 operator
            shift = self.Ns
            for site in range(self.Ns):
                for r,c in enumerate(self.operXs[site]):
                    aves[shift+site] += U[c, r]
                #aves[self.Ns+site] += np.real(np.sum((U*self.operXs[site]).diagonal()))

            # Z2 operator
            shift += self.Ns
            for i in range(self.Nb): aves[shift+i] = np.sum(self.oper2Zs[i]*U.diagonal().T)
                #aves[2*self.Ns+i] += np.real(np.sum((U*self.oper2Zs[i]).diagonal()))

            # Z3 operator
            shift += self.Nb
            for i in range(self.N3b): aves[shift+i] = np.sum(self.oper3Zs[i]*U.diagonal().T)

            # Z4 operator
            shift +=self.N3b
            for i in range(self.N4b): aves[shift+i] = np.sum(self.oper4Zs[i]*U.diagonal().T)

            return aves
        else:
            U  = np.matrix(np.zeros((2**self.Ns, 2**self.Ns)), copy=False)
            Prob = self.rho[self.bindex, self.bindex]
            U[:, self.bindex]  = self.rho[:,  self.bindex]/Prob

            Nsamples = 12
	    eO   = np.zeros((Nsamples+1, self.Ns + self.Ns + self.Nb + self.N3b + self.N4b))

            # Evaluate time-evolved operators on a discrete grid
            tau = 1.0*self.beta/(2.0*Nsamples)
            #Ub0 = np.matrix(slin.expm( self.H*tau), copy=False)
            #Uf0 = np.matrix(slin.expm(-self.H*tau), copy=False)
            Ub0 = np.matrix(mexp.expm( self.H*tau), copy=False)
            Uf0 = np.matrix(mexp.expm(-self.H*tau), copy=False)
            for step in range(Nsamples+1):
                if  step>0: U = Ub0*U*Uf0

                shift = 0
                for site in range(self.Ns):
                    eO[step][site] = np.sum(self.oper1Zs[site]*U.diagonal().T)
                    #eZ[step][site] = np.real(np.sum((self.oper1Zs[site]*U).diagonal()))

                shift += self.Ns
                for site in range(self.Ns):
                    for r,c in enumerate(self.operXs[site]):
                        eO[step][shift+site] += U[c, r]
                    #eX[step][site] = np.real(np.sum((self.operXs[site]*U).diagonal()))

                shift += self.Ns
                for i in range(self.Nb):  eO[step][shift+i] = np.sum(self.oper2Zs[i]*U.diagonal().T)

                shift += self.Nb
                for i in range(self.N3b): eO[step][shift+i] = np.sum(self.oper3Zs[i]*U.diagonal().T)

                shift += self.N3b
                for i in range(self.N4b): eO[step][shift+i] = np.sum(self.oper4Zs[i]*U.diagonal().T)

            # Numerically integrate
            UE = np.zeros(self.Ns + self.Ns + self.Nb + self.N3b + self.N4b)
            xs = np.linspace(0.0, Nsamples+1-1, Nsamples+1)*tau
            shift = self.Ns
            for site in range(self.Ns):
                UE[site]       = trapz(eO[:,site], xs)*2.0
                UE[shift+site] = trapz(eO[:,shift+site], xs)*2.0

            shift += self.Ns
            for bond in range(self.Nb):  UE[shift+bond] = trapz(eO[:,shift+bond], xs)*2.0

            shift += self.Nb
            for bond in range(self.N3b): UE[shift+bond] = trapz(eO[:,shift+bond], xs)*2.0

            shift += self.N3b
            for bond in range(self.N4b): UE[shift+bond] = trapz(eO[:,shift+bond], xs)*2.0

            return UE

