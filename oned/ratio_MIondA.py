import os,sys,glob
import loadgmt,kevent
import ssexyhelp
#Study of the ratio trick dependence upon the size of region A increament


import MCstat
import argparse
from pylab import *
import mplrc
import numpy as np
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import AutoMinorLocator
from numpy import array
from uncertainties       import unumpy
from uncertainties       import ufloat
from utils import *
import numpy as np
from scipy.optimize import curve_fit

def EE1d(l):
    # Entanglement entropy for a single region l with PBC
    # l is a fractional size of region A

    global L
    c  = 1.0 # conformal charge
    S2 = c/4.0*np.log( (L/np.pi) * np.sin(np.pi*l) )
    #S2 = np.sin(np.pi*l) 
    
    return S2

def fDeltaJz(j_z):
    return 0.5*(1.0-np.arccos(1.0*j_z)/np.pi)

def EETcorrection(l):

    global L
    g = 2
    n = 2
    global j_z
    Delta = fDeltaJz(j_z)
    A     = g/(1.0-n)*( (1.0/(n**(2.0*Delta-1.0)))*
                        ((np.sin(np.pi*l/L))**(2.0*Delta))/
                        ((np.sin(np.pi*l/(n*L)))**(2.0*Delta))
                        -n
                      ) 
    #A = 4.0*(1.0 - np.cos(np.pi*l/2.0/L))
    
    return A

def EE1dFit(l,C):
    return EE1d(l)+C

def EETcorrFit(l,A,C):
    print A,C
    global Beta
    global L
    return (EETcorrection(l)+A)*np.exp(-np.pi*Beta/L) + C


def EETfullFitNew(l,A,B):
    global Beta
    global L
    S2 = EE1d(l) #*np.exp(-np.pi*Beta/L)

def EETfullFit(l,A,B):
    global Beta
    global L
    S2 = EE1d(l) + A + (EETcorrection(l*L)*B)#*np.exp(-np.pi*Beta/L)

    return S2 

def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('fileNames', help='Reduce scalar estimator', nargs='+')
    args = parser.parse_args() 

   
    rcParams.update(mplrc.aps['params'])
    connect('key_press_event',kevent.press)
    
 
    figure(1)
    ax  = subplot(1,1,1)
    ax.set_xlabel(r'$l/L$')#, labelpad=-1.5)
    ax.set_ylabel(r'$S_2$')#, labelpad=-1.0)
    #ax.set_ylim([yl,yh]) 
    ax.set_xlim([-0.05,1.05]) 
    minorLocator   = AutoMinorLocator(3)                                           
    ax.xaxis.set_minor_locator(minorLocator)                                    
    minorLocator   = AutoMinorLocator(5)                                        
    ax.yaxis.set_minor_locator(minorLocator)
    #subplots_adjust(top=0.94,bottom=0.15, right=0.995, left=0.139)
    
    Xs     = []
    Values = []
    Errors = []
    used = []
    global L
    for j,fileName in enumerate(args.fileNames):
        scalarhelp = ssexyhelp.ScalarReduce(fileName)
    
        # Get parameters of a simulation
        fileparams = scalarhelp.getParmap()
        #dA = 3
        dA         = fileparams['dA']
        Lx         = fileparams['Lx']
        Ly         = fileparams['Ly']
        global Beta
        Beta = float(fileparams['b'])
        global L
        L = float(Lx)
        global j_z
        j_z = float(fileparams['delta'])

        #if  j==0: m = j
        #else:     m = j-1
        #Load the files
        scalarhelp.loadData()
    
        #Get the range of region-A sizes
        rp, As = scalarhelp.getrParams()

        #Get the ratios from the advanced loop column
        ratios, dratios  = scalarhelp.getAverages('ALRatio')
        ratios = unumpy.uarray(ratios,dratios)

        #Get the products from the const increament data
        if len(ratios.shape)==1:
           for i in range(1,len(ratios)):
               ratios[i] = ratios[i-1]*ratios[i]

        #Plot the results
        ratios = -unumpy.log(ratios)#/(float(Lx))
        
        colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]
        for choice in [0]:
            if  len(As.shape) == 1:
                if  choice==1:
                    ratios = ratios[:] + ratios[-1::-1]
                    ratios = ratios - ratios[-1]
                    ratios /=2.0
                else:
                    ratios = np.insert(ratios,0,0)
                x = (As+(As[1]-As[0]))
                x /=float(x[-1]) 
                x = np.insert(x,0,0)
               
                if choice==1:
                    ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  
                            color = colors[j+choice], ls = ['--','','-'][choice], marker='s',
                            mec=colors[j+choice],mfc='None',ms=3,capsize=3)
                else:
                    # Subtract off ground state entropy
                    #S2       = EE1dFit(x,0.5486)
                    #maxS2    = EE1dFit(x,0.5486+0.0015)
                    #S2ground = unumpy.uarray(S2,maxS2-S2)
                    #ratios  -= S2ground
                    
                    # Ignore the thermal entropy point
                    x      = x[:-1]
                    ratios = ratios[:-1]#*np.exp(np.pi*float(Beta)/L)
                    mratios = ratios - EETfullFit(x,0,0)

                    # Plot data
                    ax.errorbar(x,  unumpy.nominal_values(mratios),  unumpy.std_devs(mratios),  
                            color = colors[j+choice], ls = ['','--','-'][choice], marker='s',
                            mec=colors[j+choice],mfc='None',ms=3,capsize=3,
                            label=r'$L=%0.0f \, \beta = %0.0f \, J_z = %0.1f$' %(L,Beta,j_z)) 
                    

                    ##
                    (A,B) = (1,100)
                    last = -5
                    coeff, var_matrix = curve_fit(EETfullFit,x[1:last],unumpy.nominal_values(ratios)[1:last],p0=(A,B))
                    (nCoef,eCoef) = (coeff,np.sqrt(var_matrix.diagonal()))
                    #Plot the fit
                    l = np.linspace(0,x[-1]*0.9991,1000)
                    S2T  = EETfullFit(l,*nCoef)
                    S2T -= EETfullFit(l,0,.0)
                    eBeta = -1.0*unumpy.log(ufloat(nCoef[-1], eCoef[-1]))/np.pi*L
                    print 'Effective beta: %3.2f +/- %0.2f' %(eBeta.n, eBeta.s), " Real beta: %0.2f " %Beta
                    ax.plot(l,S2T,ls='-',color = colors[j+choice],alpha=0.5)
                    
    
    #for p,text in enumerate(lgd.get_texts()):
    #    if  p==0: m = p
    #    else:     m = p+#1
    #    text.set_color(colors[m])
    lgd = legend(loc='best',ncol=1,frameon=False)
    lgd.draggable(state=True)
    tight_layout() 
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

                    #minS2 = EE1d(l,nC-eC)
                    #maxS2 = EE1d(l,nC+eC)
                    #S2    = EE1d(l,nC)
                    #ax.plot(l,S2,ls='-',color = colors[j+choice],alpha=0.5)
                    #ax.fill_between(l, minS2,maxS2, facecolor = 'r', alpha =0.25, edgecolor='None')
