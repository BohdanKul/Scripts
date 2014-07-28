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
from numpy import array
from uncertainties       import unumpy
from utils import *

def main():

    parser = argparse.ArgumentParser() 
    parser.add_argument("--const", type = str, help='Reduce file with a constant increament')
    parser.add_argument("--var",   type = str, help='Reduce file with a variable increament')
    parser.add_argument('fileNames', help='Reduce scalar estimator', nargs='+')
    
    #Get paths to two reduce files
    args = parser.parse_args() 
    fconst = args.const
    fvar   = args.var

    #Set up the figure
    rcParams.update(mplrc.aps['params'])
    figure(1,(8,6))
    ax  = subplot(1,1,1)
    ax.set_xlabel(r'$\mathrm{A_{ext}}$')
    ax.set_ylabel(r'$\mathrm{Z[A_{i+1}]/Z[A_i]}$')
    #ax.set_ylabel(r'$\mathrm{dZ[A]/Z[A=0]}$')
    #ax.set_ylabel(r'$\mathrm{S_2^A}$')
    connect('key_press_event',kevent.press)

    #Create helper objects for both files
    #for k,fconst in enumerate(args.fileNames):
    #    consthelp  = ssexyhelp.ScalarReduce(fconst)
    #    fileparams = consthelp.getParmap()
    #    dA         = fileparams['dA']
    #    consthelp.loadData()
    #    rconst, drconst  = consthelp.getAverages('ALRatio')
    #    if len(rconst.shape)<1: rconst = unumpy.uarray(rconst,drconst)
    #    else:                   rconst = unumpy.uarray(rconst,drconst)[0]
    #    #rconst = -unumpy.log(rconst)
    #    #if k==0: ax.errorbar(int(dA), unumpy.nominal_values(rconst), unumpy.std_devs(rconst), ls = '', marker='s',color='r', label=r'$\mathrm{ALRT \, with  \, \Delta A=A}$') 
    #    #else   : ax.errorbar(int(dA), unumpy.nominal_values(rconst), unumpy.std_devs(rconst), ls = '', marker='s',color='r') 
    #    if k==0: ax.plot(int(dA), unumpy.std_devs(rconst),  ls = '', marker='s',color='r', label=r'$\mathrm{ALRT \, with  \, \Delta A=A$}') 
    #    else   : ax.plot(int(dA), unumpy.std_devs(rconst),  ls = '', marker='s',color='r') 
    #    #if k==0: ax.plot(int(dA), unumpy.std_devs(rconst)/unumpy.nominal_values(rconst),  ls = '', marker='s',color='r', label=r'$\mathrm{ALRT \, with  \, \Delta A=A$}') 
    #    #else   : ax.plot(int(dA), unumpy.std_devs(rconst)/unumpy.nominal_values(rconst),  ls = '', marker='s',color='r') 

        
    consthelp  = ssexyhelp.ScalarReduce(args.const)
    varhelp    = ssexyhelp.ScalarReduce(args.var)
    
    #Get their physical parameters from the filename
    #fileparams = varhelp.getParmap()
    #Lx         = fileparams['Lx']
    #Beta       = fileparams['b']
    
    #Load the files
    consthelp.loadData()
    varhelp.loadData()
    
    #Get the range of region-A sizes
    rb, As = consthelp.getrParams()
    rb, xconst = consthelp.getrParams()

    #Get the ratios from the advanced loop column
    rconst, drconst  = consthelp.getAverages('ALRatio')
    rconst = unumpy.uarray(rconst,drconst)
    #rconst = 1.0/rconst

    rvar, drvar = varhelp.getAverages('ALRatio')
    rvar        = unumpy.uarray(rvar,drvar)
    rb, xvar    = varhelp.getrParams()
    #rvar1, drvar1 = varhelp.getAverages('nAext')
    #rvar1         = unumpy.uarray(rvar1,drvar1)
    #rvar2, drvar2 = varhelp.getAverages('nAred')
    #rvar2         = unumpy.uarray(rvar2,drvar2)
    #rvar          = rvar1/rvar2
    #Get the products from the const increament data
    #for i in reversed(range(0,len(rconst)-1)):
    #    rconst[i] = rconst[i]*rconst[i+1]
    
    #for i in range(1,len(rconst)):
    #    rconst[i] = rconst[i]*rconst[i-1]

   
    
    #rvar   = -unumpy.log(rvar[:10])
    #rconst = -unumpy.log(rconst)
    #for i in reversed(range(0,len(rvar)-1)):
    #    rvar[i] = rvar[i]*rvar[i+1]
    #Plot the results
    #x = range(1,len(rconst)+1)
    #ax.plot(x,  unumpy.std_devs(unumpy.log(rconst))/unumpy.nominal_values(unumpy.log(rconst)),  ls = '', marker='s',color='r', label=r'$\mathrm{dA = const}$') 
    #ax.plot(x,  unumpy.std_devs(unumpy.log(rconst))/unumpy.nominal_values(unumpy.log(rconst)),  ls = '', marker='s',color='r', label=r'$\mathrm{dA = const}$') 
    #ax.plot(x,  unumpy.std_devs(rvar)/unumpy.nominal_values(rvar),    ls = '', marker='o',color='b', label=r'$\mathrm{dA = var}$') 
    x = As+(As[1]-As[0])
    #x  = As
    ax.errorbar(x,  unumpy.nominal_values(rconst),  unumpy.std_devs(rconst),  ls = ':', marker='s',color='b', label=r'$\mathrm{fixed A,\, moving \, \Delta A}$') 
    #ax.plot(x,    unumpy.std_devs(rconst)/unumpy.nominal_values(rconst),  ls = '', marker='s',color='b', label=r'$\mathrm{ALRT \, with \, \Delta A=8}$') 
    #ax.plot(x, unumpy.std_devs(rconst),  ls = '', marker='s',color='b', label=r'$\mathrm{\Delta A = const}$') 
    
    #x = 2*np.arange(1,len(rvar)+1)
    #rvar = unumpy.uarray(np.ones_like(rvar)*0.1319,np.ones_like(rvar)*0.0011)
    ax.errorbar(x,  unumpy.nominal_values(rvar),    unumpy.std_devs(rvar),    ls = '-', marker='o',color='r', label=r'$\mathrm{Expected \, value}$') 
    #ax.plot(x,  unumpy.nominal_values(rvar-rconst)/unumpy.std_devs(rvar-rconst),    ls = '', marker='o',color='b', label=r'$\mathrm{dA = var}$') 
    #ax.plot(x, unumpy.std_devs(rvar),    ls = '', marker='o',color='b', label=r'$\mathrm{SRT}$') 

    tight_layout()
    lgd = legend(loc='best',frameon=False)
    lgd.draggable(state=True)
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

