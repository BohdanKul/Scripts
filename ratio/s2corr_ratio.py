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
from uncertainties   import ufloat
from utils import *
from scipy.optimize import curve_fit
from collections import OrderedDict
from scipy import special

# -----------------------------------------------------------------------------
# Renyi finite-beta, finite-L correction  
# -----------------------------------------------------------------------------
def RenyiCorrection_C(betas, C):
    c     = 1.0
    g     = 2.0
    n     = 2.0
    Delta = fDeltaJz(j_z)
    global l
    global L
    print Delta
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


def fDeltaJz(j_z):
    return 0.5*(1.0-np.arccos(1.0*j_z)/np.pi)
    #return np.pi*(0.5/np.pi-0.5/(np.pi**2)*np.arccos(1.0*j_z))

def CutOff(j_z):       
    global L
    #if not(j_z<0):  minB =  L*(1.25-j_z*0.5)
    if not(j_z<0):  minB =  L*(1.5)
    elif j_z>-0.39: minB =  50
    else:           minB = 50 + 10*(0.45-j_z) #Lx*(3.6-j_z*1.75)

    return minB
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
    parser.add_argument('--cutoff',  action='store_true', default=False,  help='Plot Delta as function of the cut-off')
    args = parser.parse_args() 

    #rcParams.update(mplrc.aps['params'])
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
    if  args.cutoff:
        figure(2,(7,6))
        ax2 = subplot(211)
        connect('key_press_event',kevent.press)
        #ax2.set_xlabel(r'$\beta_{\mathrm{cutoff}}$')
        xticks([],[])
        ax2.set_ylabel(r'$\Delta$')
        ax4 = subplot(212)
        ax4.set_xlabel(r'$\beta_{\mathrm{cutoff}}$')
        ax4.set_ylabel(r'$\Delta$')
        subplots_adjust(hspace=0)


    figure(1,(7,6))
    connect('key_press_event',kevent.press)
    ax1  = subplot2grid((2,2), (0,0),colspan=2)
    ax1.set_xlabel(r'$\beta$')
    ax1.set_ylabel(r'$S_2^l$')
    global L 
    global l 
    global j_z
    i = 0
    geom = 'ful'
    LDelta  = []
    LDeltaE = []
    LJz     = []
    doubleJz = 1.0

    for fileName in args.fileNames:
        #print fileName
        scalarhelp = ssexyhelp.ScalarReduce(fileName)
        parmap = scalarhelp.getParmap()
        dA = int(parmap['dA'])
        Lx = int(parmap['Lx'])
        if 'delta' in parmap.keys(): 
            j_z = float(parmap['delta'])*doubleJz
            delta_fit = fDeltaJz(j_z)
        else: 
            delta_fit = .5
        LJz += [j_z]
        L = Lx
        l = dA
        
        scalarhelp.loadData()
        rb, Beta = scalarhelp.getrParams()
        if 'ALRatio' in scalarhelp.headers:
            Zr, dZr  = scalarhelp.getAverages('ALRatio')
            Zr       =  unumpy.uarray(Zr,dZr)
        elif 'nAred' in scalarhelp.headers:
            Ar, dAr  = scalarhelp.getAverages('nAred')
            Ar       = unumpy.uarray(Ar,dAr)
            Ae, dAe  = scalarhelp.getAverages('nAext')
            Ae       = unumpy.uarray(Ae,dAe)
            Zr       = Ae/Ar
        S2A = -unumpy.log(Zr)
        
        S2n = unumpy.nominal_values(S2A)
        S2d = unumpy.std_devs(S2A)
        #ax1.plot([Beta[0],Beta[-1]],[RenyiZero(),RenyiZero()],color='black',label = r'$S_2(T=0)$')
        mins = S2n - S2d                                       
        maxs = S2n + S2d                                       
        fill_between(Beta, mins,maxs, facecolor = colors[i%len(colors)], alpha =0.25, edgecolor='None')
        ax1.plot(Beta[0],S2n[0],
                color = colors[i%len(colors)],linewidth=4)#,label = r"$Data \, L=%2d$" %Lx) 
        (A,B,C) = (0.4,1,0.4)
        #flshift = np.max(np.where(Beta<1.14)) 
        
        
        minB = CutOff(j_z)
        if (Beta[0] > Lx) and (j_z>0): flshift = 0# np.max(np.where(Beta<Lx)) 
        else:             
            if Beta[0] < minB: flshift = np.max(np.where(Beta<minB))
            else:              flshift = 0
        fhshift = len(Beta)# np.min(np.where(Beta>245)) 
        
            
        ZBeta = np.linspace(Beta[flshift],Beta[fhshift-1],1000) 

        #coeff, var_matrix = curve_fit(RenyiCorrection_C,np.array(Beta)[flshift:fhshift],S2n[flshift:fhshift],p0=(C))
        #(C) = coeff
        #errs = np.sqrt(var_matrix.diagonal())
        #S2pred  = RenyiCorrection_C(ZBeta,C)
        #ax1.plot(ZBeta, S2pred, color=colors[i%len(colors)], linewidth = 1.5, label = r"$\mathrm{L=%2d \, J_z=%0.2f }$" %(Lx,j_z))

        #coeff, var_matrix = curve_fit(RenyiCorrection_DeltaC,np.array(Beta)[flshift:fhshift],S2n[flshift:fhshift],p0=(delta_fit,C))
        #(Delta,C) = coeff
        #errs = np.sqrt(var_matrix.diagonal())
        #print 'Length: %3d J_z: %4.3f Theory: %4.3f Fit: %4.3f +/- %4.3f' %(L, j_z, delta_fit,Delta,errs[0])
        #LDelta  += [Delta]
        #LDeltaE += [errs[0]]
        #S2pred  = RenyiCorrection_DeltaC(ZBeta,Delta,C)
        #ax1.plot(ZBeta, S2pred, color=colors[i%len(colors)], linewidth = 1.5, label = r"$\mathrm{L=%2d \, J_z=%0.2f }$" %(Lx,j_z))
        
        
        if  args.cutoff:
            #minB = 30
            flshift2 = flshift
            if Beta[0] < minB: flshift = np.max(np.where(Beta<minB))
            else:              flshift = 0
          
            leg = False
            ax2.plot([Beta[0],Beta[flshift]],[fDeltaJz(j_z),fDeltaJz(j_z)], ls = '-',color=colors[i%len(colors)], linewidth = 1.5)
            Bs   = []
            mins = []
            maxs = []
            chisqs = []
            cdfs   = []
            Bs2    = []
            for flshift in range(flshift):
                ZBeta = np.linspace(Beta[flshift],Beta[fhshift-1],1000) 
                coeff, var_matrix = curve_fit(RenyiCorrection_DeltaC,np.array(Beta)[flshift:fhshift],S2n[flshift:fhshift],p0=(delta_fit,C))
                (Delta,C) = coeff
                if type(var_matrix) is float: break
                errs = np.sqrt(var_matrix.diagonal())
                if  flshift!=flshift2: 
                    if  not(leg): 
                        ax2.errorbar(Beta[flshift],Delta,errs[0], color=colors[i%len(colors)], label = r"$\mathrm{L=%2d \, J_z=%0.2f }$" %(Lx,j_z))
                        leg = True   
                    else: 
                        mins += [Delta - errs[0]]
                        maxs += [Delta + errs[0]]
                        Bs   += [Beta[flshift]]
                        ax2.errorbar(Beta[flshift],Delta,errs[0], color=colors[i%len(colors)])
                else:       
                    ax2.errorbar(Beta[flshift],Delta,errs[0], color='black')
                
                Bs2    += [Beta[flshift]]
                dof     = len(np.array(Beta)[flshift:fhshift]) - len(coeff)
                chisqs += [sum(((S2n[flshift:fhshift]- RenyiCorrection_DeltaC(np.array(Beta)[flshift:fhshift],Delta,C))/np.array(S2d)[flshift:fhshift])**2)]
                cdfs   += [special.chdtrc(dof,chisqs[-1])]
            
            ax2.fill_between(Bs, mins,maxs, facecolor = colors[i%len(colors)], alpha =0.25, edgecolor='None')
            #ax4.plot(Bs2,chisqs, color=colors[i%len(colors)])
            ax4.plot(Bs2,cdfs, color=colors[i%len(colors)])
        i += 1

    ax1.set_title(r'$A_{\Delta} e^{-\beta*B_{\Delta}}-C$')
    #lg = ax1.legend(loc='best',ncol = 2,frameon=False)
    #lg.draggable(state=True)
    if args.cutoff:
        lg = ax2.legend(loc='best',ncol = 2,frameon=False)
        lg.draggable(state=True)
        ax2.set_ylim([0,1]) 
        ax4.set_ylim([0,1]) 
        #ax3.set_xlim([-1.05,1.05]) 

    ax3  = subplot2grid((2,2), (1,0),colspan=2)
    ax3.set_xlabel(r'$J_z$')
    ax3.set_ylabel(r'$\Delta$')
    ax3.set_ylim([0,0.5]) 
    ax3.set_xlim([-1.05,1.05]) 
    #LJz = np.arccos(np.array(LJz))
    #print LJz
    if  len(LJz)==len(LDelta):
        ax3.errorbar(np.array(LJz), LDelta, LDeltaE,
                     ls='',marker='.',color=colors[2],mec=colors[2],mfc='white', 
                     label=r"$\mathrm{QMC}$")

    Jzs    = np.linspace(-1.0,1.0,5000)
    deltas = fDeltaJz(Jzs)
    #Jzs = np.arccos(Jzs)
    ax3.plot(Jzs*1.0, deltas, color=colors[0], linewidth = 1.5, label = r"$\mathrm{Theory}$")
    lg = ax3.legend(loc='best',ncol = 1,frameon=False)
    lg.draggable(state=True)
    tight_layout()
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
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

        

