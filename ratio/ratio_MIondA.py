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
from utils import *
import numpy as np

def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('fileNames', help='Reduce scalar estimator', nargs='+')
    args = parser.parse_args() 

   
    #Set up the figure
    rcParams.update(mplrc.aps['params'])
    figure(1,(8,5))
    #figure(1,(3,2.25))
    ax  = subplot(1,1,1)
    #ax.set_xlabel(r'$\mathrm{A/(L_xL_y)}$')
    #ax.set_ylabel(r'$\mathrm{I_2^A}$')
    #ax.set_xlabel(r'$\mathrm{\beta = J/T \, [K^{-1}]}$')
    #ax.set_ylabel(r'$\mathrm{I_2^{half}/L_x}$')
    #minorLocator   = AutoMinorLocator(5)                                        
    #ax.xaxis.set_minor_locator(minorLocator)                                    
    #minorLocator   = AutoMinorLocator(5)                                        
    #ax.yaxis.set_minor_locator(minorLocator)
    connect('key_press_event',kevent.press)
    #ax.yaxis.set_ticks([0,0.05,0.10,0.15])
    #ax.yaxis.set_ticks([0,0.15])
    #ax.xaxis.set_ticks([0.995,1.00])
    
    Xs     = []
    Values = []
    Errors = []
    used = []
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
        ratios = -unumpy.log(ratios)#/(float(Lx))
        if  len(As.shape) == 1:
            ratios = np.insert(ratios,0,0)
            ratios = ratios[:] + ratios[-1::-1]
            ratios = ratios - ratios[-1]
            ratios /=2.0
        colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
       
        if len(As.shape)==1: 
            x = (As+(As[1]-As[0]))
            x /=float(x[-1]) 
            x = np.insert(x,0,0)
            ratios = ratios/float(Lx)
            ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  
                        color = colors[j], ls = '-', marker='s', 
                        label=r'$\mathrm{dA=%02d}$' %int(dA)) 
                        #label=r'$\mathrm{L_x x L_y=%02d x %02d}$' %(int(Lx),int(Lx))) 
            #ax.plot(x,  unumpy.std_devs(ratios),  
            #            color = colors[j], ls = '-', marker='s', 
            #            label=r'$\mathrm{dA=%02d}$' %int(dA)) 
           
            #d = 8.0/3.0
            #d = 8.0/1.0
            #d = 8.0/2.0
            d = 8.0/4.0
            print x[len(x)//d]
            Xs.append(Beta)
            value = unumpy.nominal_values(ratios)[len(x)//d]
            error = unumpy.std_devs(ratios)[len(x)//d]
            Values.append(value)
            Errors.append(error)
            #ax.errorbar(1.0/d,value,error,color='b',marker='s')
            #if j%2 == 0:
            #    tratios = ratios
            #else:
            #    ratios = ratios - tratios
            #ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  color=colors[j],ls = '-', marker=['s','o','d','>'][j//2], label=r'$\mathrm{L = %2.0d}$' %(int(Lx))) 
            #ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  ls = '-', color=colors[j//2], marker=['s','o','d','>'][j%2], label=r'$\mathrm{L = %2.0d \, \beta=%0.0f}$' %(int(Lx),float(Beta))) 
            #tratios = ratios[:] + ratios[-1::-1]
            #tratios = tratios - ratios[-1]
            #tratios /=2.0
            #ratios = ratios - tratios
            #ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  color=colors[j//4],ls = '-', marker=['s','o','d','>'][j%4], label=r'$\mathrm{L = %2.0d \, \beta=%0.0f}$' %(int(Lx),float(Beta))) 
            #ratios = ratios[4:-5]
            #x      = x[4:-5] 
            #if not(int(Lx) in used):
            #    ax.errorbar(float(Beta),  value, error, 
            #                color=colors[int(Lx)//4 -2],  ls = '-', marker='s',
            #                label=r'$\mathrm{L_xxL_y=%02d x %02d}$' %(int(Lx),int(Lx))) 
            #    used.append(int(Lx))
            #else:
            #    ax.errorbar(float(Beta),  value, error, 
            #                color=colors[int(Lx)//4 -2],  ls = '-', marker='s')
            #ax.errorbar(float(Beta),  unumpy.nominal_values(ratio),  unumpy.std_devs(ratio), color=colors[i],  ls = '-',ms=7, marker=['s','o'][i%2])
            #ratios = 1.0/ratios
            #ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
        else:
            x = float(dA)/float(dA)
            ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  
                        ls = '-', marker='s',color='orange',
                        label=r'$\mathrm{dA=%02d}$' %int(dA)) 
        #ax.errorbar(x,  unumpy.nominal_values(ratios),  unumpy.std_devs(ratios),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
        
        #if j==0: ref = ratios
        #mask = array(range(0,64+1,2**j))[1:] - 1
        #ax.plot(x, unumpy.std_devs(ratios)/unumpy.std_devs(ref[mask]),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
        
        #ax.plot(x, unumpy.std_devs(ratios)/unumpy.std_devs(ref[mask]),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
        
        #ax.plot(x, unumpy.std_devs(ratios),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
        #ax.plot(x, unumpy.std_devs(ratios)/unumpy.nominal_values(ratios),  ls = '-', marker='s', label=r'$\mathrm{dA=%2d}$' %int(dA)) 
    
    print "Xs:     ", Xs
    print "Values: ", Values
    print "Errors: ", Errors
    #print len(Xs),' ',len(Values),' ',len(Errors)
    #ax.set_xlim(-.005,1.25)
    #ax.set_xlim(-.005,1.005)
    #ax.set_xlim(-.005,3700)
    #ax.set_xlim(0.995,1.001)
    #ax.set_ylim(0,.15)
    tight_layout()
    lgd = legend(loc='best',ncol=1,frameon=False)
    lgd.draggable(state=True)
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

