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
from scipy import fft
from scipy import fftpack

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

def GetZ_Spline(uspline,Beta,offset,nbetas):        
    
    if nbetas == 0: nbetas = 200
    minB   =  Beta[offset]
    maxB   = Beta[-1]

    Beta   = np.linspace(minB,maxB,nbetas)
    Z      = np.zeros_like(Beta)
    dZ     = np.zeros_like(Beta)

    ispline = uspline.antiderivative()
    a = ispline(Beta[0])
    for i,b in enumerate(Beta):
        Z[i] = ispline(b) - a
    
    return unumpy.uarray(-Z,dZ),Beta


def CheckUniq(l):
    l = list(l)
    size = len(l)
    lo = []
    for i in range(size):
        elem = l.pop()
        if  elem in l:
            #print "Not unique temperature ", elem
            lo.append(size - i -1)

    return lo

def CheckEqual(l1,l2):
    
    l1 = list(l1)
    l2 = list(l2)
    if  len(l1) > len(l2): 
        bigger  = Set(l1)
        smaller = Set(l2)                                    
    else:                 
        bigger  = Set(l2)
        smaller = Set(l1)                                    
    return bigger - smaller
    


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
    parser.add_argument('fileNames',   help='Scalar estimator files', nargs='+')
    parser.add_argument('--energy',    action='store_true', default=False, help = "Show energy plots")
    parser.add_argument('--renyi',     action='store_true', default=False, help = "Show energy plots")
    parser.add_argument('--partition', action='store_true', default=False, help = "Show Renyi entropy plots plots")
    parser.add_argument('--mutual',    action='store_true', default=False, help = "Show Renyi entropy plots plots")
    parser.add_argument('--smart',     action='store_true', default=False, help = "Combine results from r=1,2 full replicas")
    parser.add_argument('--save',      action='store_true', default=False, help = "Save MI to hard drive")
    parser.add_argument('--check',     action='store_true', default=False, help = "Turn on consistency beta check")
    parser.add_argument('--spline',    type=int,            default=0,     help = "Set to > 0 to perform a spline fit")
    args = parser.parse_args() 

    rcParams.update(mplrc.aps['params'])
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
    figure(1,(7,6))
    connect('key_press_event',kevent.press)
    nplots = args.renyi + args.mutual + args.energy
    iplot  = nplots
    if  args.mutual:
        ax  = subplot(nplots,1,iplot)
        iplot -= 1
        ax.set_ylabel(r'$\mathrm{I/(L_x)}$')
        ax.set_xlabel(r'$\mathrm{%s}[K^{-1}]$' %parMap['b'])
    if  args.renyi:    
        ax1  = subplot(nplots,1,iplot)
        iplot -= 1
        ax1.set_ylabel(r'$\mathrm{(S_2-S_2[\beta_{max}])/L}$')
        ax1.set_xlabel(r'$\mathrm{%s}[K^{-1}]$' %parMap['b'])
        #ax1.set_yscale('log')
    if  args.energy or args.partition:    
        ax2  = subplot(nplots,1,iplot)
        iplot -= 1
        #ax2.set_xlabel(r'$\mathrm{%s}$' %parMap['b'])
        #ax2.set_xlabel(r'$\mathrm{\beta [K^{-1}]}$')
        #ax2.set_ylabel(r'$\mathrm{(E_{fit}-E_{data})/dE_{data}}$')
        #ax2.set_ylabel(r'$\mathrm{(E_{fit}-E_{data})/dE_{data}}$')
       


    order = ['Lx','Ly','T','b']
    
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

    if  (args.renyi or args.mutual): 
        for Lx,dic in Lxs.items():
            if  not( ('full',1) in dic):
                print "One replica reduce file was not found for ", Lx
                sys.exit(0)
            else:
                t    = dic[('full',1)]
                del dic[('full',1)]
                Lxs[Lx] = OrderedDict(dic)
                Lxs[Lx][('full',1)] = t
    
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
            print " Geometry: ", geom, " Data length:  ", len(Beta)
           
            
            if (Beta[0]!=0):
                print " Inserting value at infinity: Beta[0] = %0.4f" %Beta[0]
                Beta = np.insert(Beta,0,0)
                E    = np.insert(E,0,0)
                dE   = np.insert(dE,0,0)
   
            if  args.check:
                print " Checking beta uniqueness"
                mask = CheckUniq(Beta)
                if len(mask)>0:
                    print "   Deleting repeated betas: ", Beta[mask]
                    Beta = np.delete(Beta,mask)
                    E    = np.delete(E,mask)
                    dE   = np.delete(dE,mask)
            
            Beta = Beta[shift:]
            #if j>0: CheckEqual(Beta2,Beta)
            #Beta2 = Beta
            
            Norm = int(Lx)**2 
            E    = E[shift:]*Norm
            dE   = dE[shift:]*Norm
           
            LxData[geom] = {}
            LxData[geom]['Beta'] = Beta
            LxData[geom]['E']    = E
            LxData[geom]['dE']   = dE
      
        if  args.check:
            print " Checking beta correspondence between full2 and full1"
            diff = CheckEqual(LxData[('full',2)]['Beta'],LxData[('full',1)]['Beta'])
            if  len(diff)>0:
                print "   full 2 <> full 1"
                if  (not(args.spline) and not(args.smart)):
                    exit(0)
            #print " Checking beta correspondence between full2 and half2"
            #diff = CheckEqual(LxData[('full',2)]['Beta'],LxData[('half',2)]['Beta'])
            #if  len(diff)>0:
            #    print "   full 2 <> half 2"
            #    if  not(args.spline):
            #        exit(0)


         
        if  args.smart:
            print " Smart mode"
            # Combine results for r = 1 and 2 full replicas based on the relation:
            # n_r2(Beta) = n_r1(2*Beta)
            # E_r2(Beta) = 2*E_r1(2*Beta)
            Beta2r = LxData[('full',2)]['Beta'] 
            Beta1r = LxData[('full',1)]['Beta']
            E2r    = LxData[('full',2)]['E']
            E1r    = LxData[('full',1)]['E']
            dE2r   = LxData[('full',2)]['dE']
            dE1r   = LxData[('full',1)]['dE']

            # if 2r has more data points
            dBeta2r = np.array([])
            indices = []
            print "   Length 2r: ", len(Beta2r), " 1r: ",len(Beta1r)
            if  len(Beta2r)>len(Beta1r):
                print "   Exporting betas from 2r to 1r" 
                dBeta2r = 2*Beta2r

                # Find the extra betas that can be used
                for index,b in enumerate(dBeta2r):
                    if  (not(b in Beta1r) and (b in Beta2r)):
                        indices.append(index)
            print "   ",len(dBeta2r[indices]), ' elements added'
            nBeta1r = np.concatenate((Beta1r,dBeta2r[indices]))
            nE1r    = np.concatenate((E1r,  0.5*E2r[indices]))
            ndE1r   = np.concatenate((dE1r,0.5*dE2r[indices]))
            
            # Sort 
            sort_indices = np.argsort(nBeta1r)
            nBeta1r = nBeta1r[sort_indices]
            nE1r    =    nE1r[sort_indices]
            ndE1r   =   ndE1r[sort_indices]
            
            LxData[('full',1)]['Beta'] = nBeta1r
            LxData[('full',1)]['E']    = nE1r
            LxData[('full',1)]['dE']   = ndE1r



# ====================== Process data =============================

        cutoff = 0
        for geom, data in reversed(dic.items()):
            # Take the processed data
            Beta  = LxData[geom]['Beta'][cutoff:]  
            E     = LxData[geom]['E'][cutoff:]       
            dE    = LxData[geom]['dE'][cutoff:]      
            if  cutoff==0:
                dE[0] = dE[1] 

            # Generate a new dataset if the spline interpolation is turned on
            if  (args.spline!=0):
                nBeta  = np.linspace(min(Beta),max(Beta),args.spline)
                nE,ndE = SplineGenerate(Beta,E,dE,nBeta)
            else:
                (nBeta,nE,ndE) = (Beta,E,dE)

            uspline = interpolate.UnivariateSpline(Beta[cutoff:],E[cutoff:], w=1./dE[cutoff:],k=3)

            # ------------Energy plots-------------
            if  args.energy:
                    ax2.errorbar(Beta,E,dE,
                                color = colors[i], ls='',
                                label = ((r'$\mathrm{Lx=%2.0d}$' %Lx) +'\n'+ (r'$\mathrm{S_{2%s}}$' %geom[0] )))
                    if  args.spline:
                        FancyErrorbar(ax2,nBeta, unumpy.uarray(nE,ndE), 
                                      color=colors[i],label="")

            # Stop here if no further analysis is required
            if  not(args.mutual or args.renyi or args.partition):
                continue

            # Otherwise, integrate the energy curve
            nZ, dZ  = Bootstrap(SplineIntegrate,nBeta,nE,ndE,250)
            (nZ,dZ) = (-nZ,dZ)
            Z       = unumpy.uarray(nZ,dZ)

            # ------------Partition plots----------
            if args.partition:
                        FancyErrorbar(ax2,nBeta,Z,\
                              color = colors[i],
                              label = r'$\mathrm{r=%0.2d;\, A=%s}$' %(geom[1],geom[0]))

            
            # ------------Entropy plots------------
            if  args.renyi or args.mutual:
                if  geom == ('full',1): Zf1 = Z 
                else:
                    S2  = -Z + 2*Zf1  
                    if   geom == ('half',2): S2A  = S2B = S2 = S2 + np.log(2)*(int(Lx)**2)*0.5  
                    elif geom == ('A',2):    S2A  = S2 = S2  + np.log(2)*9.0
                    elif geom == ('B',2):    S2B  = S2 = S2  + np.log(2)*7.0
                    elif geom == ('full',2): S2AB = S2 = S2  + np.log(2)*(int(Lx)**2)

                    if  args.renyi:
                        S2n = unumpy.nominal_values(S2)#-S2[-1])
                        S2d = unumpy.std_devs(S2)#-S2[-1])
                        FancyErrorbar(ax1,nBeta,S2,
                                 color = colors[i%len(colors)],
                                 label = ((r'$\mathrm{Lx=%2.0d}$' %Lx) +'\n'+ (r'$\mathrm{S_{2%s}}$' %geom[0] )))

            j += 1
        #T, MI = loadtxt("MI_ED.dat", unpack=True)
        #ax.plot(1.0/np.array(T), MI/float(Lx),\
        #        marker='',color=colors[-1],\
        #        label=r'$\mathrm{ED}$')
            
            
        # ------------MI plots----------
        if args.mutual:

            I  = (S2A+S2B-S2AB)/float(Lx)
            FancyErrorbar(ax,nBeta,I,\
                          color = colors[(i)%len(colors)],
                          label = (r'$\mathrm{r=%0.2d;\, A=%s}$' %(geom[1],geom[0])))

            if  args.save:
                filename = 'reduce-%s_Lx-%02d.dat' %('mutual',Lx)
                outFile  = open(filename,'w')
                outFile.write('#%15s%16s%16s\n' %('Beta','MI','+/-'))
                (nI,dI) = (unumpy.nominal_values(I),
                           unumpy.std_devs(I))
                for k in range(len(I)):
                    outFile.write('%16.8E%16.8E%16.8E\n' %(nBeta[k],nI[k],dI[k]))
                outFile.close()
        i += 1

    if  args.mutual:
        ax.legend(loc='best',frameon=False)
    if args.renyi:
       ax1.legend(loc='best',frameon=False)
    if args.energy or args.partition:
       ax2.legend(loc='best',frameon=False)
    tight_layout()
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

