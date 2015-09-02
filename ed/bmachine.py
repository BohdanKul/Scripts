from qutip    import *
from Hbuilder import *

from scipy.integrate import simps
import numpy as np
from pylab import real 


class BoltzmannMachine:
    def __init__(self, N, bonds, Zfield, Xfield, Inter, beta):
        self.N = N
        I = qeye(2)
        self.H = tensor([I]*N) - tensor([I]*N)
       
        for site, zfield in enumerate(Zfield):
            if  zfield!=0:
                self.H += EmbedSiteOper(sigmaz(), zfield, site, N)

        for site, xfield in enumerate(Xfield):
            if  xfield!=0:
                self.H += EmbedSiteOper(sigmax(), xfield, site, N)

        self.bonds = bonds
        for i,inter in enumerate(Inter):
            siteA, siteB = bonds[i]
            if  inter!=0:
                self.H += Embed2SiteOper(sigmaz(), inter, siteA, siteB, N)
       
        self.beta = beta
        self.rho  = (-self.H*beta).expm()


    def setProjector(self, cbits):
        self.clamped = False
        if  len(cbits)>0:
            self.clamped = True
            ket = basis(2,cbits[0])
            for qbit in cbits[1:]:
                ket = tensor(ket, basis(2,qbit))
            self.P = ket*ket.dag()

        if (len(cbits)>0) and (len(cbits)!=self.N): self.P = tensor(self.P, tensor([qeye(2)]*(self.N-len(cbits))))
        elif (len(cbits)==0):                       self.P = tensor(tensor([qeye(2)]*self.N))


    def evaluateProjector(self):
        return (self.rho*self.P).tr()/self.rho.tr()
    
    def computeLocalAverages(self):
        rho  = self.rho*self.P
        Z    = rho.tr()
        rho /= Z
        if  not self.clamped:
            aveZ = []
            aveX = []
            for site in range(self.N):
                oper = EmbedSiteOper(sigmaz(), 1.0, site, self.N) 
                aveZ += [real((rho*oper).tr())]
            
                oper = EmbedSiteOper(sigmax(), 1.0, site, self.N) 
                aveX += [real((rho*oper).tr())]

            aveZZ = []
            for bond in self.bonds:
                (siteA, siteB) = (bond[0], bond[1])
                oper = Embed2SiteOper(sigmaz(), 1.0, siteA, siteB, self.N)
                aveZZ += [real((rho*oper).tr())]
            
            return aveZ, aveX, aveZZ
        else:
            Nsamples = 100
            eZ  = [np.zeros(Nsamples+1)]*self.N
            eX  = [np.zeros(Nsamples+1)]*self.N
            eZZ = [np.zeros(Nsamples+1)]*len(self.bonds)
            xs = []
            # Evaluate time-evolved operators on a discrete grid
            for step in range(Nsamples+1):
                tau = 1.0*step*self.beta/(1.0*Nsamples)
                xs  += [tau]
                Ub   = ( self.H*tau).expm()
                Uf   = (-self.H*tau).expm()
                
                for site in range(self.N):
                    oper = EmbedSiteOper(sigmaz(), 1.0, site, self.N) 
                    eZ[site][step] = real((Uf*oper*Ub*rho).tr())
                
                for site in range(self.N):
                    oper = EmbedSiteOper(sigmax(), 1.0, site, self.N) 
                    eX[site][step] = real((Uf*oper*Ub*rho).tr())
                
                for i,bond in enumerate(self.bonds):
                    (siteA, siteB) = (bond[0], bond[1])
                    oper = Embed2SiteOper(sigmaz(), 1.0, siteA, siteB, self.N)
                    eZZ[i][step] = real((Uf*oper*Ub*rho).tr())
            
            # Numerically integrate 
            (UaveZ, UaveX, UaveZZ) = ([], [], [])
            for site in range(self.N):
                UaveZ += [simps(eZ[site], np.array(xs))]
                UaveX += [simps(eX[site], np.array(xs))]
            
            for bond in range(len(self.bonds)):
                UaveZZ += [simps(eZZ[bond], np.array(xs))]
           
            return UaveZ, UaveX, UaveZZ
        
