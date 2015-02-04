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
def RenyiCorrection_C(betas, C):
    c     = 1.0
    g     = 2.0
    n     = 2.0
    Delta = 0.25
    global l
    global L
    B = 2.0*np.pi*Delta/L
    A     = g/(1.0-n)*( (1.0/(n**(2.0*Delta-1.0)))*
                        ((np.sin(np.pi*l/L))**(2.0*Delta))/
                        ((np.sin(np.pi*l/(n*L)))**(2.0*Delta))
                        -n
                      )

    return A*np.exp(-B*betas) -C

def RenyiCorrection_DeltaC(betas, Delta,C):
    c     = 1.0
    n     = 2.0
    g     = 2.0
    global l
    global L
    B = 2.0*np.pi*Delta/L
    #A     = g/(1.0-n)*(1.0/pow(n,2.0*Delta-1.0)*(pow(np.sin(np.pi*l/L),2.0*Delta)/pow(np.sin(np.pi*l/(n*L)),2.0*Delta))-n)
    A     = g/(1.0-n)*( (1.0/(n**(2.0*Delta-1.0)))*
                        ((np.sin(np.pi*l/L))**(2.0*Delta))/
                        ((np.sin(np.pi*l/(n*L)))**(2.0*Delta))
                        -n
                      )
    return A*np.exp(-B*betas)-C 

def RenyiCorrection_DeltaBC(betas, Delta,B,C):
    c     = 1.0
    n     = 2.0
    g     = 2.0
    global l
    global L
    print l, L
    A     = g/(1.0-n)*( (1.0/(n**(2.0*Delta-1.0)))*
                        ((np.sin(np.pi*l/L))**(2.0*Delta))/
                        ((np.sin(np.pi*l/(n*L)))**(2.0*Delta))
                        -n
                      )
    return A*np.exp(-B*betas)-C 

def RenyiZero():
    global l
    global L
    n = 2.0
    c = 1.0
    eta = 1.0
    S2_zero = c/(1.0-n)*(1.0-n*n)/(6.0*n)*np.log(L/np.pi*np.sin(np.pi*l/L))
    S2_zero = 2.0*c/(1.0-n)*(1.0-n*n)/(6.0*n)*np.log(L/np.pi*np.sin(np.pi*l/L))
    S2_zero = c*(1.0+n)/(3.0*n*eta)*np.log(eta*L/np.pi*np.sin(np.pi*l/L))
    
    return S2_zero

def RenyiCorrection_ADeltaC(betas, A,Delta,C):
    global L
    B = 2.0*np.pi*Delta/L
    A = A
    return A*np.exp(-B*betas)-C

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
    ax1.set_xlabel(r'$\beta$')
    ax1.set_ylabel(r'$S_2^l$')
    global L 
    global l 
    i = 0
    geom = 'ful'
    for fileName in args.fileNames:
        scalarhelp = ssexyhelp.ScalarReduce(fileName)
        parmap = scalarhelp.getParmap()
        dA = int(parmap['dA'])
        Lx = int(parmap['Lx'])
        L = Lx
        l = dA
        
        scalarhelp.loadData()
        rb, Beta = scalarhelp.getrParams()
        Zr, dZr  = scalarhelp.getAverages('ALRatio')
        Z   =  unumpy.uarray(Zr,dZr)
        S2A = -unumpy.log(Z)

        
        S2n = unumpy.nominal_values(S2A)
        S2d = unumpy.std_devs(S2A)
        #ax1.plot([Beta[0],Beta[-1]],[RenyiZero(),RenyiZero()],color='black',label = r'$S_2(T=0)$')
        mins = S2n - S2d                                       
        maxs = S2n + S2d                                       
        fill_between(Beta, mins,maxs, facecolor = colors[i], alpha =0.25, edgecolor='None')
        ax1.plot(Beta[0],S2n[0],
                color = colors[i],linewidth=4)#,label = r"$Data \, L=%2d$" %Lx) 
        (A,B,C) = (0.4,1,0.4)
        #flshift = np.max(np.where(Beta<1.14)) 
        minB = 40
        if (Beta[0] > minB): flshift = np.max(np.where(Beta<65.14)) 
        else:                 flshift = np.max(np.where(Beta<minB))
        fhshift = np.min(np.where(Beta>395.91)) #len(Beta)# 
        ZBeta = np.linspace(Beta[flshift],Beta[fhshift-1],1000) 
        #coeff, var_matrix = curve_fit(RenyiCorrection_ADeltaC,np.array(Beta)[flshift:fhshift],S2n[flshift:fhshift],p0=(A,B,C))
        #(A,B,C) = coeff


        #(A,B,C) = (0.4,0.01,0.4)
        #coeff, var_matrix = curve_fit(RenyiCorrection_DeltaBC,np.array(Beta)[flshift:fhshift],S2n[flshift:fhshift],p0=(A,B,C))
        #(Delta,B,C) = coeff
        #errs = np.sqrt(var_matrix.diagonal())
        #S2pred  = RenyiCorrection_DeltaBC(ZBeta,Delta,B,C)
        #ax1.plot(ZBeta, S2pred, linewidth = 1.5, label = r"$\mathrm{A_{\Delta} e^{-\beta*B}-C, \, \Delta=%0.3f(%0.3f)}$" %(Delta,errs[1]))
         
        #coeff, var_matrix = curve_fit(RenyiCorrection_ADeltaC,np.array(Beta)[flshift:fhshift],S2n[flshift:fhshift],p0=(A,B,C))
        #(A,Delta,C) = coeff
        #errs = np.sqrt(var_matrix.diagonal())
        #S2pred  = RenyiCorrection_ADeltaC(ZBeta,A,Delta,C)
        #ax1.plot(ZBeta, S2pred, linewidth = 1.5, label = r"$\mathrm{A e^{-\beta*B_{\Delta}}-C, \, \Delta=%0.3f(%0.3f)}$" %(Delta,errs[1]))

        coeff, var_matrix = curve_fit(RenyiCorrection_DeltaC,np.array(Beta)[flshift:fhshift],S2n[flshift:fhshift],p0=(0.5,C))
        (Delta,C) = coeff
        errs = np.sqrt(var_matrix.diagonal())
        S2pred  = RenyiCorrection_DeltaC(ZBeta,Delta,C)
        #ax1.plot(ZBeta, S2pred, linewidth = 1.5, label = r"$\mathrm{A_{\Delta} e^{-\beta*B_{\Delta}}-C, \, \Delta=%0.3f(%0.3f)}$" %(Delta,errs[0]))
        ax1.plot(ZBeta, S2pred, linewidth = 1.5, label = r"$\mathrm{L=%2d \,  \Delta=%0.3f(%0.3f)}$" %(Lx,Delta,errs[0]))
        
        #coeff, var_matrix = curve_fit(RenyiCorrection_C,np.array(Beta)[flshift:fhshift],S2n[flshift:fhshift],p0=(C))
        #(C) = coeff
        #errs = np.sqrt(var_matrix.diagonal())
        #S2pred  = RenyiCorrection_C(ZBeta,C)
        #ax1.plot(ZBeta, S2pred, linewidth = 1.5, label = r"$\mathrm{A_{\Delta=\frac{1}{4}}e^{-\beta*B_{\Delta=\frac{1}{4}}}-C}$" )

        i += 1

    ax1.set_title(r'$A_{\Delta} e^{-\beta*B_{\Delta}}-C$')
    lg = ax1.legend(loc='best',ncol = 2,frameon=False)
    lg.draggable(state=True)
    tight_layout()
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

