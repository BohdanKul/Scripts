from optparse import OptionParser
import pyutils as py
import numpy as np
from uncertainties   import unumpy
import scipy.optimize as optimize
import matplotlib.pylab as pl
import matplotlib as mpl
from scipy.optimize import curve_fit
import loadgmt,kevent
from scipy import special
import ssexyhelp
import argparse

def HalfLog(x,a,c):
    return a+0.5*np.log(x)/(1.0*x)+c/(1.0*x)

def OneLog(x,a,c):
    return a+1.0*np.log(x)/(1.0*x)+c/(1.0*x)

def ZeroLog(x,a,c):
    return a+c/(1.0*x)

def FreeLog1(x,a,b,c):
    return (a + b*np.log(x)/(1.0*x) + c/(1.0*x))#+ d/(1.0*x*x))*0.5

def FreeLog2(x,a,b):
    return a + b*np.log(x)/(1.0*x)

def main():

    parser = argparse.ArgumentParser(description='Scaling analysis of MI as function of Beta')
    parser.add_argument('fileNames',             help='MI files', nargs='+')
    args = parser.parse_args() 


    nfit  = 3
    cfits = ['a','b','c','d'] 
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
    Z = [[0,0],[0,0]]
    nhs = 6
    min, max = (0, 8)
    mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
    levels = np.linspace(min,max,nhs)
    CS3 = pl.contourf(Z, levels, cmap=mymap)
    pl.clf()
    Clb = pl.colorbar(CS3) # using the colorbar info I got from contourf
    Clb.ax.set_ylabel(r'$\mathrm{\beta \, [K^{-1}]}$')
    
    ax = pl.subplot(1,1,1)
    ax.set_xlabel(r'$\mathrm{L_x}$')
    ax.set_ylabel(r'$\mathrm{MI/(2L_x)}$')
    
    f,axs = pl.subplots(nfit)
    pl.connect('key_press_event',kevent.press)
    
    n = len(open(args.fileNames[0]).readlines())-1
    Lxs  = np.zeros(len(args.fileNames))
    MIs  = np.zeros((n,len(args.fileNames)))
    dMIs = np.zeros((n,len(args.fileNames)))
    for i,fileName in enumerate(args.fileNames):
        scalarhelp = ssexyhelp.ScalarReduce(fileName)
        parmap = scalarhelp.getParmap()
        Lx     = int(parmap['Lx'])
        Lxs[i] = Lx
        
        scalarhelp = ssexyhelp.ScalarReduce(fileName)
        MIs[:,i], dMIs[:,i] = scalarhelp.getAverages('MI')
    rb, Beta = scalarhelp.getrParams()
    
    sorted_indices = np.argsort(Lxs)
    Lxs = Lxs[sorted_indices]

    
    LogCoeff  = np.zeros((n,nfit))
    dLogCoeff = np.zeros((n,nfit))
    #chisqs    = np.zeros(n)
    #cdfs      = np.zeros(n)
    for i,b in enumerate(Beta):
        MI  =  MIs[i,:][sorted_indices]
        dMI = dMIs[i,:][sorted_indices]
 
        LogCoeff[i,:], var_matrix = curve_fit(FreeLog1,Lxs,MI,np.ones(nfit))
        dLogCoeff[i,:] = np.sqrt(np.diagonal(var_matrix))
        #dof       = len(Lcs) - len(LogCoeff[i,:])
        #chisqs[i] = sum(((Mis-FreeLog1(Lxs,*LogCoeff))/sigma)**2)
        #cdfs      = special.chdtrc(dof,chisq)

    for i in range(nfit):
        axs[i].errorbar(Beta,LogCoeff[:,i],dLogCoeff[:,i], label=r"$\mathrm{a+b*log(L)/L+c/L}$")
        #axs[i].set_xlabel(r'$\mathrm{\beta}[K^{-1}]$')
        axs[i].set_xlabel(r'$\mathrm{\beta}[K^{-1}]$')
        axs[i].set_ylabel(r'$\mathrm{%s}$' %cfits[i])
    axs[i].legend(loc=4)
    
    # Using contourf to provide my colorbar info, then clearing the figure
    lookup_betas = range(-len(Beta),-1,1)
    #lookup_betas = range(-20,-1,1)
    
    step = float(max-min)/len(lookup_betas)
    for j,i in enumerate(lookup_betas):
        r = (float(step*j)-min)/(max-min)
        g = 0
        b = 1-r
        ax.plot(Lxs,MIs[i,:][sorted_indices],marker = 'o',linestyle='',color = (r,g,b))#,label=r'$\mathrm{\beta=%0.3f}$' %Beta[i])
        ax.plot(Lxs,FreeLog1(Lxs,*LogCoeff[i,:]),color = (r,g,b))
        #ax.plot(Beta,MIs[:,0],color = colors[j%len(colors)])
    
    pl.tight_layout()
    ax.legend()
    pl.show()
if __name__ == "__main__":
   main()
