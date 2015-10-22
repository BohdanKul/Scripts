from qutip    import *
import Hbuilder

from scipy.integrate import simps, trapz
import numpy as np

class BoltzmannMachine:
    def __init__(self, N, beta, **kwargs):
        self.Ns = N
        self.Nb = len(kwargs['Z2'])
        I = qeye(2)
        self.H = tensor([I]*N) - tensor([I]*N)
       
        for itype, bonds in kwargs.iteritems():
            if itype == 'X': oper = sigmax()
            else:            oper = sigmaz()
            for bond in bonds:
                self.H += Hbuilder.EmbedOper(oper, N, bond)  
        
        self.beta = beta

        self.rho  = (-self.H*beta).expm()

        self.operZs = []
        self.operXs = []
        for site in range(N):
            self.operZs += [Hbuilder.EmbedOper(sigmaz(), N, [site, 1.0])]
            self.operXs += [Hbuilder.EmbedOper(sigmax(), N, [site, 1.0])]
            
        self.operZZs = []
        for bond in kwargs['Z2'].copy():
            bond[2] = 1.0 
            self.operZZs += [Hbuilder.EmbedOper(sigmaz(), N, bond)]


    def setProjector(self, cbits):
        self.clamped = False
        if  len(cbits)>0:
            self.clamped = True
            ket = basis(2,cbits[0])
            for qbit in cbits[1:]:
                ket = tensor(ket, basis(2,qbit))
            self.P = ket*ket.dag()

        if (len(cbits)>0) and (len(cbits)!=self.Ns): self.P = tensor(self.P, tensor([qeye(2)]*(self.N-len(cbits))))
        elif (len(cbits)==0):                       self.P = tensor(tensor([qeye(2)]*self.Ns))

    def evaluateProjector(self):
        return (self.rho*self.P).tr()/self.rho.tr()
    
    def computeLocalAverages(self, test = False):
        U  = self.rho*self.P
        Z  = U.tr()
        U /= Z
        if  ((not self.clamped) or (test)):
            aves = np.zeros(self.Ns+self.Ns+self.Nb) 
            for site in range(self.Ns):
                aves[site] = np.real((U*self.operZs[site]).tr())
            
                aves[self.Ns+site] += np.real((U*self.operXs[site]).tr())

            for i in range(self.Nb):
                aves[2*self.Ns+i] += np.real((U*self.operZZs[i]).tr())
            
            return aves
        else:
            Nsamples = 10 
            eZ  = np.zeros((Nsamples+1,self.Ns))
            eX  = np.zeros((Nsamples+1,self.Ns)) 
            eZZ = np.zeros((Nsamples+1,self.Nb))
            # Evaluate time-evolved operators on a discrete grid
            tau = 1.0*self.beta/2.0/(Nsamples)
            Ub0 = ( self.H*tau).expm()
            Uf0 = (-self.H*tau).expm()
            for step in range(Nsamples+1):
                if  step>0: U = Ub0*U*Uf0   

                for site in range(self.Ns):
                    eZ[step][site] = np.real((self.operZs[site]*U).tr())
                    #eX[step][site] = np.real((self.operXs[site]*U).tr())
                
                for i in range(self.Nb):
                    eZZ[step][i] = np.real((self.operZZs[i]*U).tr())
        
            # Numerically integrate 
            UE = np.zeros(self.Ns+self.Ns+self.Nb) 
            xs  = np.linspace(0.0, Nsamples+1-1, Nsamples+1)*tau
            for site in range(self.Ns):
                UE[site]            = trapz(eZ[:,site], xs)*2.0
                #UE[self.Ns + site] += trapz(eX[:,site], xs)*2.0
            
            for bond in range(self.Nb):
                UE[2*self.Ns+bond] += [trapz(eZZ[:,bond], xs)*2.0]
         
            return UE 
        
