from qutip    import *
from Hbuilder import *

from scipy.integrate import simps
import numpy as np


class BoltzmannMachine:
    def __init__(self, N, Zfield, Xfield, Inter, beta):
        self.N = N
        I = qeye(2)
        self.H = tensor([I]*N)

        for site, zfield in enumerate(Zfield):
            if  zfield!=0:
                self.H += EmbedSiteOper(sigmaz(), zfield, site, N)

        for site, xfield in enumerate(Xfield):
            if  xfield!=0:
                self.H += EmbedSiteOper(sigmax(), xfield, site, N)

        self.bonds = []
        for inter in Inter:
            (siteA, siteB, inter) = inter
            bonds += [ [siteA,siteB] ]
            if  inter!=0:
                self.H += Embed2SiteOper(sigmaz(), inter, siteA, siteB, N)
        
        self.beta = beta
        self.rho  = (-self.H*beta).expm()


    def setProjector(cbits):
        cbits = []
        self.clamped = False
        if  len(cbits)>0:
            self.clamped = True
            ket = basis(2,cbits[0])
            for qbit in clamped[1:]:
                ket = tensor(ket, basis(2,qbit))
            self.P = ket*ket.dag()

        if len(cbits)>0: self.P = tensor(self.P, tensor([qeye(2)]*(self.N-len(cbits))))
        else:            self.P = tensor(tensor([qeye(2)]*self.N))


    def computeLocalAverages():
        rho *= self.rho*self.P
        Z    = rho.tr()
        rho /= Z
        
        if  not self.clamped:
            aveZ = []
            aveX = []
            for site in range(self.N):
                oper = EmbedSiteOper(sigmaz(), 1.0, site, N) 
                aveZ += [real((rho*oper).tr())]
            
                oper = EmbedSiteOper(sigmax(), 1.0, site, N) 
                aveX += [real((rho*oper).tr())]

            aveZZ = []
            for bond in self.bonds:
                (siteA, siteB) = (bond[0], bond[1])
                oper = Embed2SiteOper(sigmaz(), 1.0, siteA, siteB, N)
                aveZZ += [real((rho*oper).tr())]
            
            return aveZ, aveX, aveZZ
        else:
            Nsamples = 100
            eZ  = [np.zeros(Nsampes)]*self.N
            eX  = [np.zeros(Nsampes)]*self.N
            eZZ = [np.zeros(Nsampes)]*len(self.bonds)
            xs = []
            # Evaluate time-evolved operators on a discrete grid
            for step in range(Nsamples+1):
                tau = 1.0*step*self.beta/(1.0*Nsamples)
                xs  += [tau]
                Ub   = ( H*tau).expm()
                Uf   = (-H*tau).expm()
                
                for site in range(N):
                    oper = EmbedSiteOper(sigmaz(), 1.0, site, self.N) 
                    eZ[site][step] = real((Uf*oper*Ub*rho).tr())
                
                for site in range(N):
                    oper = EmbedSiteOper(sigmax(), 1.0, site, self.N) 
                    eX[site][step] = real((Uf*oper*Ub*rho).tr())
                
                for i,bond in enumerate(self.bonds):
                    (siteA, siteB) = (bond[0], bond[1])
                    oper = Embed2SiteOper(sigmaz(), 1.0, siteA, siteB, self.N)
                    eZZ[i][step] = real((Uf*oper*Ub*rho).tr())
            
            # Numerically integrate 
            (UaveZ, UaveX, UaveZZ) = ([], [], [])
            for site in range(N):
                UaveZ += [sims(np.array(xs), eZ[site])]
                UaveX += [sims(np.array(xs), eX[site])]
            
            for bond in range(len(self.bonds)):
                UaveZZ += [sims(np.array(xs), eZZ[bond])]
            
            return UaveZ, UaveX, UaveZZ
        
