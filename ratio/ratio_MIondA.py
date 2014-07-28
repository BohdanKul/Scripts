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
import numpy as np

def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('fileNames', help='Reduce scalar estimator', nargs='+')
    args = parser.parse_args() 

   
    #Set up the figure
    rcParams.update(mplrc.aps['params'])
    figure(1,(8,6))
    ax  = subplot(1,1,1)
    ax.set_xlabel(r'$\mathrm{A_{ext}/L}$')
    #ax.set_xlabel(r'$\mathrm{\beta}$')
    ax.set_ylabel(r'$\mathrm{Z[A}]/Z[A_{i-1}]}$')
    #ax.set_ylabel(r'$\mathrm{Z[A}]/Z[A=0]}$')
    #ax.set_ylabel(r'$\mathrm{\frac{dZ[A]/Z[A=0]}{Z[A]/Z[A=0]}}$')
    #ax.set_ylabel(r'$\mathrm{abs(\frac{-dlog(Z[A]/Z[A=0]}{Z[A]/Z[A=0]}))}$')
    #ax.set_ylabel(r'$\mathrm{\frac{dS_2^A}{S_2^A}}$')
    #ax.set_ylabel(r'$\mathrm{dS_2^A}$')
    ax.set_ylabel(r'$\mathrm{S_2^A}$')
    #ax.set_ylabel(r'$\mathrm{MI/2L}$')
    #ax.set_ylabel(r'$\mathrm{(MI_{\beta=184}-MI_{\beta=92})/2L}$')
    #ax.set_ylabel(r'$\mathrm{(S_2^{\beta=184}-S_2^{\beta=92})/L}$')
    #ax.set_ylabel(r'$\mathrm{S_2/L - MI/2L}$')
    #ax.set_xlim((-1,65))
    connect('key_press_event',kevent.press)


    for j,fileName in enumerate(args.fileNames):
        scalarhelp = ssexyhelp.ScalarReduce(fileName)
    
        # Get parameters of a simulation
        fileparams = scalarhelp.getParmap()
        #dA = 3
        dA         = fileparams['dA']
        Lx         = fileparams['Lx']
        Ly         = fileparams['Ly']
        Beta       = fileparams['b']
    
        #Load the files
        scalarhelp.loadData()
    
        #Get the range of region-A sizes
        rp, As = scalarhelp.getrParams()

        #Get the ratios from the advanced loop column
        ratios, dratios  = scalarhelp.getAverages('ALRatio')
        #ratios  = np.insert(ratios,0,0)
        #dratios = np.insert(dratios,0,0)
        #As      = np.insert(As,0,0)


        #dratios *= 1.0/(2.0**j)
        ratios = unumpy.uarray(ratios,dratios)

        #Get the products from the const increament data
        if len(ratios.shape)==1:
           for i in range(1,len(ratios)):
               ratios[i] = ratios[i-1]*ratios[i]

        #Plot the results
        ratios = -unumpy.log(ratios)/(float(Lx)**2)
        ratios = np.insert(ratios,0,0)
        #ratios = ratios[:] + ratios[-1::-1]
        #ratios = ratios - ratios[-1]
        #ratios /=2.0
        #ratios = ratios[:-1] + ratios[-2::-1]
        colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
        if len(ratios.shape)==1: 
            x = (As+(As[1]-As[0]))
            x /=float(x[-1]) 
            x = np.insert(x,0,0)

            #x = As
            #ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
            #if j%2 == 0:
            #    tratios = ratios
            #else:
            #    ratios = ratios - tratios
            #    ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  color=colors[j],ls = '-', marker=['s','o','d','>'][j//2], label=r'$\mathrm{L = %2.0d}$' %(int(Lx))) 
            #ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  ls = '-', color=colors[j//2], marker=['s','o','d','>'][j%2], label=r'$\mathrm{L = %2.0d \, \beta=%0.0f}$' %(int(Lx),float(Beta))) 
            tratios = ratios[:] + ratios[-1::-1]
            tratios = tratios - ratios[-1]
            tratios /=2.0
            ratios = ratios - tratios
            ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  color=colors[j//4],ls = '-', marker=['s','o','d','>'][j%4], label=r'$\mathrm{L = %2.0d \, \beta=%0.0f}$' %(int(Lx),float(Beta))) 
            #ratios = ratios[4:-5]
            #x      = x[4:-5] 
            #for i,ratio in enumerate(ratios):
            #    if j==0: ax.errorbar(float(Beta),  unumpy.nominal_values(ratio),  unumpy.std_devs(ratio), color=colors[i],  ls = '-',ms=7, marker=['s','o'][i%2],label='A=%2d' %x[i]) 
            #    else:    ax.errorbar(float(Beta),  unumpy.nominal_values(ratio),  unumpy.std_devs(ratio), color=colors[i],  ls = '-',ms=7, marker=['s','o'][i%2])
            #ratios = 1.0/ratios
            #ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
        else:
            x = int(dA)
            ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  ls = '-', marker='s',color='orange') 
        #ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
        
        #if j==0: ref = ratios
        #mask = array(range(0,64+1,2**j))[1:] - 1
        #ax.plot(x, unumpy.std_devs(ratios)/unumpy.std_devs(ref[mask]),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
        
        #ax.plot(x, unumpy.std_devs(ratios)/unumpy.std_devs(ref[mask]),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
        
        #ax.plot(x, unumpy.std_devs(ratios),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
        #ax.plot(x, unumpy.std_devs(ratios)/unumpy.nominal_values(ratios),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
    
    ax.set_xlim(-0.02,1.02)
    tight_layout()
    legend(loc='best',frameon=False)
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

