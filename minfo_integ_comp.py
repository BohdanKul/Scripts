import os,sys,glob
import loadgmt,kevent
import ssexyhelp
import MCstat
from optparse import OptionParser
import argparse
from pylab import *
import mplrc
import numpy as np
from matplotlib.ticker import MaxNLocator
from scipy import interpolate
from sets import Set
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
    

def MultiplePlot(axs, toplot, x,y,dy, rx, ry,**kwargs):
    iplot = 0
    if dy[0] == 0: dy[0] = dy[1]
    #print dy
    if  toplot[0]: 
        axs[iplot].plot(rx,ry,ls='--',**kwargs)
        FancyErrorbar(axs[iplot],x, unumpy.uarray(y,dy), **kwargs)
        iplot += 1

    if  toplot[1]: 
        axs[iplot].plot(x,(y-ry)/dy,**kwargs)
        iplot += 1
    
    if  toplot[2]: 
        axs[iplot].plot(x,dy,**kwargs)

    return 0 



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
    parser.add_argument('--diff',     action='store_true', default=False, help = "Turn on consistency beta check")
    parser.add_argument('--error',     action='store_true', default=False, help = "Turn on consistency beta check")
    parser.add_argument('--value',     action='store_true', default=False, help = "Turn on consistency beta check")
    parser.add_argument('--spline',    type=int,            default=0,     help = "Set to > 0 to perform a spline fit")
    args = parser.parse_args() 

    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
    
    # Array of columns and raws to display
    acolumn = np.array([args.energy,args.partition,args.renyi,args.mutual])
    araw    = np.array([args.value,args.diff,args.error])
    
    # The number of columns and raws to display
    ncolumn = len(np.where(araw)[0])
    nraw    = len(np.where(acolumn)[0])
    
    # Determine indices of subplots to display
    mcolumn = np.tile(acolumn.reshape(4,1),(1,3))
    mraw    = np.tile(araw,(4,1))
    index  = mcolumn*mraw

    # Subplots labelling
    ylabels = np.array([[r'$\mathrm{E}$',      r'$\mathrm{(E-E_r)/dE}$',             r'$\mathrm{dE}$'],
                        [r'$\mathrm{ln(Z)}$',  r'$\mathrm{(ln(Z)-ln(Z)_r)/dln(Z)}$', r'$\mathrm{dln(Z)}$'],
                        [r'$\mathrm{S_2/L_x}$',  r'$\mathrm{(S_2-S_{2_r})/dS_{2_r}}$', r'$\mathrm{dS_{2_r}/L_x}$'],
                        [r'$\mathrm{I/(L_x)}$',r'$\mathrm{(I-I_r)/dI}$',             r'$\mathrm{dI/(L_x)}$' ]])   
    ylabels = ylabels[index].reshape(nraw,ncolumn)
    
    # Create subplots grid
    fig,axm  = subplots(nraw,ncolumn)
    connect('key_press_event',kevent.press)
    rcParams.update(mplrc.aps['params'])
    
    # Label subplots grid
    i = 0
    for xy, ax in np.ndenumerate(axm):
        if  i == 0: 
            legendax = ax
            i += 0
        ax.set_ylabel(ylabels[xy])
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        if  (xy[0] == nraw-1): 
            ax.set_xlabel(r'$\mathrm{%s}[K^{-1}]$' %parMap['b'])

    # Reshape the subplots grid matrix
    axmf = axm.flatten()
    axm  = np.zeros_like(index,dtype=type(ax))
    i = 0
    for xy,toplot in np.ndenumerate(index):
        if  toplot:
            axm[xy] = axmf[i]
            i += 1
    non0index = index

    # Get the data
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

        cutoff     = 0
        highcutoff = -1
        for geom, data in reversed(dic.items()):
            # Take the processed data
            Beta  = LxData[geom]['Beta'][cutoff:highcutoff]  
            E     = LxData[geom]['E'][cutoff:highcutoff]       
            dE    = LxData[geom]['dE'][cutoff:highcutoff]      
            if  cutoff==0:
                dE[0] = 1.0e-12#dE[1] 

            # "True" E and its integral splines
            rspline  = interpolate.UnivariateSpline(Beta,E, w=1./dE,k=3)
            rispline = rspline.antiderivative()
           
            tBeta  = np.linspace(min(Beta),max(Beta),args.spline)
            tE,tdE = SplineGenerate(Beta,E,dE,tBeta)

            
            
            # Generate a new dataset if the spline interpolation is turned on
            if  (args.spline!=0):
                #nBeta  = np.linspace(min(Beta),max(Beta),args.spline)
                nBeta  = tBeta
                nE,ndE = SplineGenerate(tBeta,tE,tdE, nBeta)
            else:
                (nBeta,nE,ndE) = (Beta,E,dE)

            # ------------Energy plots-------------
            if  args.energy:
                    #MultiplePlot(axm[0,:], np.array([False,False,True]), 
                    #             Beta,nE,dE, Beta, dE,
                    #             color=colors[j%len(colors)],label='yo')
                    MultiplePlot(axm[0,:], non0index[0,:], nBeta,nE,ndE, tBeta, tE,
                                 color=colors[j%len(colors)],
                                 label = ((r'$\mathrm{Lx=%2.0d}$' %Lx) +'\n'+ (r'$\mathrm{S_{2%s}}$' %geom[0] )))

            # Stop here if no further analysis is required
            if  not(args.mutual or args.renyi or args.partition):
                continue

            # "True" energy integral
            tZ = -(rispline(tBeta) - rispline(tBeta[0]))/float(Lx)

            # Otherwise, integrate the energy curve
            nZ, dZ  = Bootstrap(SplineIntegrate,nBeta,nE,ndE,250)
            (nZ,dZ) = (-nZ/float(Lx),dZ/float(Lx))
            Z       = unumpy.uarray(nZ,dZ)

            #nZ, dZ  = Bootstrap(SimpsonIntegrate,nBeta,nE,ndE,250)
            #(nZ,dZ) = (-nZ/float(Lx),dZ/float(Lx))
            #Z       = unumpy.uarray(nZ,dZ)
            
            #nE =  nE+ np.random.normal(0,1.0,len(nE))*ndE # Generate noise
            #nZ, dZ  = SimpsonIntegrateError(nBeta,nE,ndE)
            #(nZ,dZ) = (-nZ/float(Lx),dZ/float(Lx))
            #Z       = unumpy.uarray(nZ,dZ)
            
            # ------------Partition plots----------
            if args.partition:
                    MultiplePlot(axm[1,:], non0index[1,:], nBeta,nZ,dZ, tBeta, tZ,
                    color=colors[j%len(colors)],
                    label=r'$\mathrm{r=%0.2d;\, A=%s}$' %(geom[1],geom[0]))

            
            # ------------Entropy plots------------
            if  args.renyi or args.mutual:
                if  geom == ('full',1): 
                    Zf1  = Z 
                    tZf1 = tZ
                else:
                    S2   = -Z  + 2*Zf1  
                    tS2  = -tZ + 2*tZf1  
                    if   geom == ('half',2): 
                         S2A   = S2B  = S2
                         tS2A  = tS2B = tS2
                    elif geom == ('A',2):    
                         S2A   = S2   #+ np.log(2)*(int(Lx)**2)*0.5  
                         tS2A  = tS2  
                    elif geom == ('B',2):    
                         S2B   = S2   #+ np.log(2)*(int(Lx)**2)*0.5
                         tS2B  = tS2   
                    elif geom == ('full',2): 
                         S2AB  = S2   #+ np.log(2)*(int(Lx)**2)
                         tS2AB = tS2   

                    if  args.renyi:
                        S2n = -unumpy.nominal_values(S2)#-S2[-1])
                        S2d = unumpy.std_devs(S2)#-S2[-1])
                        MultiplePlot(axm[2,:], non0index[2,:], nBeta,S2n,S2d, tBeta, -tS2,
                                     color=colors[j%len(colors)],
                                     label = ((r'$\mathrm{Lx=%2.0d}$' %Lx) +'\n'+ (r'$\mathrm{S_{2%s}}$' %geom[0] )))

            j += 1
        T, MI = loadtxt("MI_ED.dat", unpack=True)
        axm[3,0].plot(1.0/np.array(T), MI/float(Lx),\
                marker='',color=colors[-1],\
                label=r'$\mathrm{ED}$')
            
            
        # ------------MI plots----------
        if args.mutual:
            I   = S2A  + S2B  - S2AB
            tI  = tS2A + tS2B - tS2AB
            In = unumpy.nominal_values(I)
            Id = unumpy.std_devs(I)
            MultiplePlot(axm[3,:], non0index[3,:], nBeta,In,Id, tBeta, tI,color=colors[j%len(colors)],label='yo')
            #FancyErrorbar(ax,nBeta,I,\
            #              colors[(i)%len(colors)],
            #              (r'$\mathrm{r=%0.2d;\, A=%s}$' %(geom[1],geom[0])))

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

    if  args.mutual or args.renyi or args.energy or args.partition:
        legendax.legend(loc='upper right',frameon=False)
    fig.tight_layout()
    subplots_adjust(wspace = 0.31,hspace = .25, bottom = 0.08)
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

