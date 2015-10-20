from qutip    import *
from Hbuilder import *

from scipy.integrate import simps
from scipy.integrate import trapz 
import numpy as np
from pylab import real 
import time

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

        self.operZs = []
        self.operXs = []
        for site in range(self.N):
            self.operZs += [EmbedSiteOper(sigmaz(), 1.0, site, self.N)]
            self.operXs += [EmbedSiteOper(sigmax(), 1.0, site, self.N)]
            
        self.operZZs = []
        for i,bond in enumerate(self.bonds):
            (siteA, siteB) = (bond[0], bond[1])
            self.operZZs += [Embed2SiteOper(sigmaz(), 1.0, siteA, siteB, self.N)]


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
    
    def computeLocalAverages(self, test = False):
        U  = self.rho*self.P
        Z  = U.tr()
        U /= Z
        if  ((not self.clamped) or (test)):
            aveZ = []
            aveX = []
            for site in range(self.N):
                aveZ += [real((U*self.operZs[site]).tr())]
            
                oper = EmbedSiteOper(sigmax(), 1.0, site, self.N) 
                aveX += [real((U*self.operXs[site]).tr())]

            aveZZ = []
            for i in range(len(self.bonds)):
                aveZZ += [real((U*self.operZZs[i]).tr())]
            
            return aveZ, aveX, aveZZ
        else:
            Nsamples = 10 
            eZ  = np.zeros((Nsamples+1,self.N))
            eX  = np.zeros((Nsamples+1,self.N)) 
            eZZ = np.zeros((Nsamples+1,len(self.bonds)))
            # Evaluate time-evolved operators on a discrete grid
            tau = 1.0*self.beta/2.0/(Nsamples)
            Ub0 = ( self.H*tau).expm()
            Uf0 = (-self.H*tau).expm()
            for step in range(Nsamples+1):
                if  step>0: U = Ub0*U*Uf0   

                for site in range(self.N):
                    #t = self.rho*self.P
                    #t /= t.tr()
                    #eZ[step][site] = real((self.operZs[site]*t).tr())
                    eZ[step][site] = real((self.operZs[site]*U).tr())
                    eX[step][site] = real((self.operXs[site]*U).tr())
                    #eX[step][site] = real((self.operXs[site]*U).tr())
                
                for i in range(len(self.bonds)):
                    #print 'bond ', i, ' oper: ', self.operZZs[i]
                    eZZ[step][i] = real((self.operZZs[i]*U).tr())
        
            # Numerically integrate 
            (UaveZ, UaveX, UaveZZ) = ([], [], [])
            xs  = np.linspace(0.0, Nsamples+1-1, Nsamples+1)*tau
            #print eZ[:,0]
            for site in range(self.N):
                UaveZ += [trapz(eZ[:,site], xs)*2.0]
                UaveX += [trapz(eX[:,site], xs)*2.0]
                #UaveZ += [simps(eZ[:,site], xs)]
                #UaveX += [simps(eX[:,site], xs)]
            
            for bond in range(len(self.bonds)):
                UaveZZ += [trapz(eZZ[:,bond], xs)*2.0]
           
            return UaveZ, UaveX, UaveZZ
        
