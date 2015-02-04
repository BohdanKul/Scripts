import os,sys,glob
import loadgmt,kevent
import ssexyhelp
import MCstat
from optparse import OptionParser
import argparse
from scipy import integrate
from pylab import *
import mplrc
import numpy as np
from scipy import interpolate
from sets import Set
from simpserror import simps_error
from uncertainties   import unumpy
from utils import *
from collections import OrderedDict

def GetDicKey(dic, val=''):

    keys = []
    for key,value in dic.items():
            if  value==val: 
                keys.append(key)

    return keys

def GetZ_Simpson(Beta,E,dE, offset):        
    
    Sleng = len(Beta) - offset 
    Z  = np.zeros(Sleng)
    dZ = np.zeros(Sleng)
    Z[0]  = 0
    dZ[0] = 0

    Beta = Beta[offset:]
    E    = E[offset:]
    dE   = dE[offset:]

    for i in range(1,len(Beta)):
        #print i, E[i], Beta[i]
        
        Z[i]  = integrate.simps((-1)*E[:i],Beta[:i], even='first')
        if  i<3:  dZ[i-offset] = 0    
        else:     dZ[i-offset] = simps_error(dE[:i], Beta[:i], even='first')
    return unumpy.uarray(Z,dZ),Beta

def GetZ_Spline(ospline,Beta,offset,nbetas):        
    
    if nbetas == 0: nbetas = 200
    minB   = 0 # Beta[offset]
    maxB   = 8 #Beta[-1]

    #Beta   = np.linspace(minB,maxB,nbetas)
    Z      = np.zeros_like(Beta)
    dZ     = np.zeros_like(Beta)

    ispline = ospline.antiderivative()
    a = ispline(Beta[0])
    for i,b in enumerate(Beta):
        Z[i] = ispline(b) - a
    
    return -Z,Beta


# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main():


    parser = argparse.ArgumentParser(description='Test integration techniques based on a spline model of energy data')
    parser.add_argument('fileNames',   help='Reduce estimator files', nargs='+')
    args = parser.parse_args() 

    rcParams.update(mplrc.aps['params'])
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
    figure(1,(7,6))
    connect('key_press_event',kevent.press)
    
    ax1  = subplot(3,1,1)
    ax2  = subplot(3,1,2)
    ax3  = subplot(3,1,3)
    ax1.set_ylabel(r'$\mathrm{lnZ}$')
    ax3.set_ylabel(r'$\mathrm{E}$')
    ax3.set_xlabel(r'$\mathrm{\beta}$')
       


    
    Lxs = {}
    offset = 0

    for fileName in args.fileNames:
        scalarhelp = ssexyhelp.ScalarReduce(fileName)
        parmap = scalarhelp.getParmap()
        Lx = int(parmap['Lx'])
        if  not Lx in Lxs:
            Lxs[Lx] = {}

        geom = GetDicKey(parmap)[0]
        r    = int(parmap['r'])
        Lxs[Lx][(geom,r)] = scalarhelp

   
    data = {}
    i = 0
    j = 0
    equalA = False
    shift = 0
    for Lx, dic in Lxs.items():
        
        # ====================== Modify original data =====================
        LxData = OrderedDict()
        print "Lx = %2.0f" %Lx
        for geom, rscalar in reversed(dic.items()):
            
            rscalar.loadData()
            rb, Beta = rscalar.getrParams()
            E, dE    = rscalar.getAverages('ET')
          
            LxData[geom] = {}
            LxData[geom]['Beta'] = Beta
            LxData[geom]['E']    = E
            LxData[geom]['dE']   = dE
      
        


# ====================== Process data =============================

        for geom, data in reversed(dic.items()):
            
            # Get the raw data
            Beta = LxData[geom]['Beta']  
            E    = LxData[geom]['E']*float(Lx)**2     
            dE   = LxData[geom]['dE']*float(Lx)**2     
            dE[0] = dE[1]
            
            
            # Interpolate the raw data on nBeta
            nBeta  = np.linspace(min(Beta),max(Beta),len(Beta))
            nE,ndE = SplineGenerate(Beta,E,dE,nBeta)
            
            # Integrate it
            sZ = SplineIntegrate(nBeta,nE,ndE)/float(Lx) 
            
            # Perform a bootstrap analysis on the spline integration
            bsZ, bsdZ  = Bootstrap(SplineIntegrate,nBeta,nE,ndE,250)
            (bsZ,bsdZ) = (bsZ/float(Lx),bsdZ/float(Lx))
            
            
            
            # Display results
            print geom
            ax3.errorbar(Beta, E, dE, 
                         ls = '', marker = '', color = colors[-2], 
                         label=r"$\mathrm{Raw \, data}$")
            
            FancyErrorbar(ax3,nBeta, unumpy.uarray(nE,ndE), 
                          colors[-1],r"$\mathrm{Interpolation}$")
            
            ax2.plot(nBeta, (sZ-bsZ)/bsdZ,
                     ls = '-', color = colors[i%len(colors)], 
                     label=r"$\mathrm{(Z_{bs} - Z_{raw})/dZ_{bs}}$")

            #ax2.plot(nBeta, bsdZ,
            #         ls = '--', color = colors[i%len(colors)], 
            #         label=r"$\mathrm{dZ_{bs}}$")

            ax1.plot(nBeta, sZ,
                     ls = '-', color = colors[i%len(colors)], 
                     label=r"$\mathrm{Raw \, data}$")
            
            ax1.plot(nBeta, bsZ,
                     ls = '--', color = colors[i%len(colors)-1], 
                     label=r"$\mathrm{Bootstrap}$")
            
            j += 1
           
        i += 1

    ax1.legend(loc='best',frameon=False)
    ax2.legend(loc='best',frameon=False)
    ax3.legend(loc='best',frameon=False)
    tight_layout()
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

