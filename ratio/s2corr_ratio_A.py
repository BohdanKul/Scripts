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
from matplotlib.ticker import MultipleLocator
from scipy import interpolate
from sets import Set
from simpserror import simps_error
from uncertainties   import unumpy
from utils import *
from scipy.optimize import curve_fit
from collections import OrderedDict

# -----------------------------------------------------------------------------
# Renyi finite-beta, finite-L correction  
# -----------------------------------------------------------------------------
def RenyiCorrection_C(l, C):
    c     = 1.0
    g     = 2.0
    n     = 2.0
    Delta = 0.25
    global beta
    global L
    B = 2.0*np.pi*Delta/L
    A     = g/(1.0-n)*( (1.0/(n**(2.0*Delta-1.0)))*
                        ((np.sin(np.pi*l/L))**(2.0*Delta))/
                        ((np.sin(np.pi*l/(n*L)))**(2.0*Delta))
                        -n
                      )

    return A*np.exp(-B*beta) +C

def RenyiCorrection2_AC(l, A,C):
    global beta
    c     = 1.0
    n     = 2.0
    h     = (c/24.0)*(1.0-1.0/(n**2)) 
    Delta = 2*h
    L     = 64.0;
    B = 2.0*np.pi*Delta/L
    return A*np.exp(-B*beta)-C

def RenyiCorrection_DeltaC(l, Delta,C):
    c     = 1.0
    n     = 2.0
    g     = 2.0
    global beta
    global L
    B = 2.0*np.pi*Delta/L
    #A     = g/(1.0-n)*(1.0/pow(n,2.0*Delta-1.0)*(pow(np.sin(np.pi*l/L),2.0*Delta)/pow(np.sin(np.pi*l/(n*L)),2.0*Delta))-n)
    A     = g/(1.0-n)*( (1.0/(n**(2.0*Delta-1.0)))*
                        ((np.sin(np.pi*l/L))**(2.0*Delta))/
                        ((np.sin(np.pi*l/(n*L)))**(2.0*Delta))
                        -n
                      )
    return A*np.exp(-B*beta)-C 

def RenyiCorrection_DeltaBC(l, Delta,B,C):
    c     = 1.0
    n     = 2.0
    g     = 2.0
    global beta
    global L
    A     = g/(1.0-n)*( (1.0/(n**(2.0*Delta-1.0)))*
                        ((np.sin(np.pi*l/L))**(2.0*Delta))/
                        ((np.sin(np.pi*l/(n*L)))**(2.0*Delta))
                        -n
                      )
    return A*B-C 

def RenyiZero(l):
    global L
    n = 2.0
    c = 1.0
    S2_zero = c/(1.0-n)*(1.0-n*n)/(6.0*n)*np.log(L/np.pi*np.sin(np.pi*l/L))
    
    return S2_zero

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main():

    # define the mapping between short names and label names 
    parMap = {'x': r'L_x',
              'y': r'L_y',
              'b': r'\beta',
              'T': r'T',
              'r': r'r'}


    parser = argparse.ArgumentParser(description='Plot Raw MC Equilibration Data for Scalar Estimators.')
    parser.add_argument('fileNames', help='Scalar estimator files', nargs='+')
    args = parser.parse_args() 

    rcParams.update(mplrc.aps['params'])
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
    figure(1,(7,6))
    connect('key_press_event',kevent.press)
    ax1 = subplot(1,1,1,)

    global L 
    L = 64.0
    #L = 256.0
    global beta 
    beta = 150.0
    geom = 'ful'
    for fileName in args.fileNames:
        scalarhelp = ssexyhelp.ScalarReduce(fileName)
        parmap = scalarhelp.getParmap()
        Lx = int(parmap['Lx'])
           
        scalarhelp.loadData()
        rb, As = scalarhelp.getrParams()
        As  = As + As[1]-As[0]
        As = np.insert(As,0,0)
        nAs = np.linspace(As[1],As[-2],100)
        Zr, dZr  = scalarhelp.getAverages('ALRatio')
        Zr   =  unumpy.uarray(Zr,dZr)
        if len(Zr.shape)==1:
           for i in range(1,len(Zr)):
               Zr[i] = Zr[i-1]*Zr[i]
        S2A = -unumpy.log(Zr)
        S2A = np.insert(S2A,0,0)
        S2A = S2A[:] + S2A[-1::-1]
        S2A = S2A - S2A[-1]
        S2A /=2.0

        As  = As[1:-1]
        S2A = S2A[1:-1] 
        S2A = S2A-RenyiZero(As) 
        S2n = unumpy.nominal_values(S2A)
        S2d = unumpy.std_devs(S2A)
        ax1.errorbar(As,S2n,S2d,
              color = colors[i%len(colors)],label = ((r'$\mathrm{L=%2.0d}$' %Lx) +'\n'+ (r'$\mathrm{S_{2}^{\beta=%3.2f}}$' %(beta) )))
        #ax1.plot(nAs,RenyiZero(nAs),color='black',label = r'$S_2(T=0)$')
            
        (A,B,C) = (0.4,1,0.4)
       
        #coeff, var_matrix = curve_fit(RenyiCorrection_C,As,S2n,p0=(C))
        #(C) = coeff
        #errs = np.sqrt(var_matrix.diagonal())
        #S2pred  = RenyiCorrection_C(nAs,C)
        #ax1.plot(nAs, S2pred, linewidth = 2, label = r"$\mathrm{A(\Delta=\frac{1}{4})e^{-\beta*B(\Delta=\frac{1}{4})}-C}$" )

        #(A,B,C) = (0.4,1,0.4)
        #coeff, var_matrix = curve_fit(RenyiCorrection_DeltaC,As,S2n,p0=(B,C))
        #(Delta,C) = coeff
        #errs = np.sqrt(var_matrix.diagonal())
        #S2pred  = RenyiCorrection_DeltaC(nAs,Delta,C)
        #ax1.plot(nAs, S2pred, linewidth = 2, label = r"$\mathrm{A(\Delta) e^{-\beta*B(\Delta)}-C, \, \Delta=%0.3f(%0.3f)}$" %(Delta,errs[0]))

        #(A,B,C) = (0.4,1,0.4)
        #coeff, var_matrix = curve_fit(RenyiCorrection_DeltaBC,As,S2n,p0=(A,B,C))
        #(Delta,B,C) = coeff
        #errs = np.sqrt(var_matrix.diagonal())
        #
        #S2pred  = RenyiCorrection_DeltaBC(nAs,Delta,B,C)
        #ax1.plot(nAs, S2pred, linewidth = 2, label = r"$\mathrm{A(\Delta) e^{-\beta*B(\Delta)}-C, \, \Delta=%0.3f(%0.3f)}$" %(Delta,errs[0]))


        i += 1

    ax1.legend(loc='best',frameon=False)
    tight_layout()
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

