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

def GetZ(Beta,E,dE, offset):        
    
    #spline = interpolate.splrep(Beta,E, w=1./dE)
    #Beta  = np.arange(min(Beta),max(Beta),0.001)
    

    ##E = np.random.normal(E,dE)
    #spline = interpolate.splrep(Beta,E)
    #Beta  = np.arange(0,max(Beta),0.01)
    #E   = interpolate.splev(Beta,spline)
    
    Sleng = len(Beta) - offset 
    Z  = np.zeros(Sleng)
    dZ = np.zeros(Sleng)
    for i in range(offset,len(Beta)):
        #print i, E[i], Beta[i]
        
        Z[i-offset]  = integrate.simps((-1)*E[:i],Beta[:i], even='first')
        if  i<3:
            dZ[i-offset] = 0    
        else:
            dZ[i-offset] = simps_error(dE[:i], Beta[:i], even='first')
    Beta = Beta[offset:]
    return unumpy.uarray(Z,dZ),Beta

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
    parser.add_argument('fileNames', help='Scalar estimator files', nargs='+')
    parser.add_argument('--energy',    action='store_true', default=False, help = "Show energy plots")
    parser.add_argument('--renyi',     action='store_true', default=False, help = "Show energy plots")
    parser.add_argument('--partition', action='store_true', default=False, help = "Show Renyi entropy plots plots")
    parser.add_argument('--mutual', action='store_true', default=False, help = "Show Renyi entropy plots plots")
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
        ax.set_ylabel(r'$\mathrm{I/L}$')
        ax.set_xlabel(r'$\mathrm{%s}[K^{-1}]$' %parMap['b'])
    if  args.renyi:    
        ax1  = subplot(nplots,1,iplot)
        iplot -= 1
        ax1.set_ylabel(r'$\mathrm{(S_2-S_2[\beta_{max}])/L}$')
        ax1.set_xlabel(r'$\mathrm{%s}[K^{-1}]$' %parMap['b'])
        #ax1.set_yscale('log')
    if  args.energy:    
        ax2  = subplot(nplots,1,iplot)
        iplot -= 1
        #ax2.set_xlabel(r'$\mathrm{%s}$' %parMap['b'])
        #ax2.set_xlabel(r'$\mathrm{\beta [K^{-1}]}$')
        #ax2.set_ylabel(r'$\mathrm{(E_{fit}-E_{data})/dE_{data}}$')
        #ax2.set_ylabel(r'$\mathrm{(E_{fit}-E_{data})/dE_{data}}$')
       


    order = ['Lx','Ly','T','b']
    
    Lxs = {}
    offset = 1

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
        for geom, rscalar in reversed(dic.items()):
            #print geom
            
            rscalar.loadData()
            rb, Beta = rscalar.getrParams()
            E, dE    = rscalar.getAverages('ET')
            
            #Beta = Beta[::2]
            #E    = E[::2]
            #dE   = dE[::2]
            #if (1.0/Beta[0]) > 20:
            #    print "Inserting value at infinity"
            #    Beta = np.insert(Beta,0,0)
            #    E    = np.insert(E,0,0)
            #    dE   = np.insert(dE,0,0)
   
            mask = CheckUniq(Beta)
            if len(mask)>0:
                print "Deleting repeating betas: "#, Beta[mask]
                Beta = np.delete(Beta,mask)
                E    = np.delete(E,mask)
                dE   = np.delete(dE,mask)
            
            Beta = Beta[shift:]
            if j>0: CheckEqual(Beta2,Beta)
            Beta2 = Beta
            
            Norm = int(Lx)**2 
            E    = E[shift:]*Norm
            dE   = dE[shift:]*Norm
            
            #Z,BetaZ = GetZ(Beta,E,dE,offset)
        
            if args.energy:
                        #ax2.plot(Beta, E, color = colors[i])#, label = r'$\mathrm{fit}$')
                        #Beta  = np.arange(min(Beta),max(Beta),0.001)
                        #if i==0: ax2.plot(Beta[1:], E, color = colors[i], marker='o', label=r'$\mathrm{r=%0.2d;\, A=%s}$' %(geom[1],geom[0]))
                        #else:    ax2.plot(Beta[1:], E, color = colors[i], marker='o')
                        #if i==0: ax2.errorbar(Beta, E, dE, linestyle = 'None', color = colors[i], marker='o', label=r'$\mathrm{r=%0.2d;\, A=%s}$' %(geom[1],geom[0]))
                        #else:    ax2.errorbar(Beta, E, dE, linestyle = 'None', color = colors[i], marker='o')
                        spline = interpolate.splrep(Beta[1:],E[1:], w=1./dE[1:])
                        Ediff   = (interpolate.splev(Beta[1:],spline) - E[1:])/dE[1:]
                        #n, bins, patches = ax2.hist(Ediff, 20,normed=1,  facecolor='green', alpha=0.6,label=r'$\mathrm{(E_{fit}-E_{data})/dE_{data}}$')
                        #ax2.plot(Beta[1:],E, color = colors[i], label = r'$\mathrm{Original}$')
                        
                        Ehigh = Ediff[np.where(Beta[1:]> 2.529)]
                        BetaFFT = Beta[1:][np.where(Beta[1:]> 2.529)]
                        #ax2.plot(BetaFFT,E, color = colors[i], label = r'$\mathrm{FFT}$')
                    
                        main_sigh = Ehigh
                        nFs = 5
                        #for j in range(nFs):
                        #    sample_freq = fftpack.fftfreq(main_sigh.size, d=0.01)
                        #    pidxs = np.where(sample_freq > 0)
                        #    freqs = sample_freq[pidxs]
                        #    period = 1.0/freqs
                        #    EFFT = fftpack.fft(main_sigh)
                        #    power = np.abs(EFFT[pidxs])
                        #    maxfreq = freqs[power.argmax()]
                        #    #print maxfreq,1.0/maxfreq
                        #    #if j==0: ax2.plot(freqs,power, color = colors[i], label = r'$\mathrm{FFT}$')
                        #    #ax2.plot(maxfreq,power[freqs == maxfreq],marker='s', color = colors[i+1])
                        #    
                        #    EFFT[np.abs(sample_freq) == maxfreq] = 0
                        #    main_sigh = fftpack.ifft(EFFT)

                        #ax2.plot(BetaFFT,main_sigh, color = colors[i+1])

                        Elow = Ediff[np.where(Beta[1:]< 2.53)]
                        BetaFFT = Beta[1:][np.where(Beta[1:]< 2.53)]
                        #ax2.plot(BetaFFT,E, color = colors[i], label = r'$\mathrm{FFT}$')
                    
                        main_sig = Elow 
                        for j in range(nFs):
                            sample_freq = fftpack.fftfreq(main_sig.size, d=0.004)
                            pidxs = np.where(sample_freq > 0)
                            freqs = sample_freq[pidxs]
                            period = 1.0/freqs
                            EFFT = fftpack.fft(main_sig)
                            power = np.abs(EFFT[pidxs])
                            maxfreq = freqs[power.argmax()]
                            #print maxfreq,1.0/maxfreq
                            if j==0: ax2.plot(1.0/freqs,power, color = colors[i%len(colors)], label = r'$\mathrm{FFT}$')
                            ax2.plot(1.0/maxfreq,power[freqs == maxfreq],marker='s', color = colors[(i+1)%len(colors)])
                            
                            EFFT[np.abs(sample_freq) == maxfreq] = 0
                            main_sig = fftpack.ifft(EFFT)

                        #ax2.plot(BetaFFT,main_sig, color = colors[i+1], label = r'$\mathrm{5 \, Fs \, \  removed}$')



                        #ax2.plot(freqs,abs(EFFT), color = colors[i], label = r'$\mathrm{FFT}$')
                        #ax2.plot(E, color = colors[i+1])#, label = r'$\mathrm{fit}$')
                        #ax2.plot(Beta[1:], E, color = colors[i])#, label = r'$\mathrm{fit}$')
                        #RHT0, REedges = np.histogram(, xhbinN,normed=True)
                        #n, bins, patches = ax2.hist(np.concatenate((main_sig,main_sigh)), 20,normed=1,  facecolor='blue', alpha=0.15,label=r'$\mathrm{Removed \, 4 \, Fs}$')
                        #n, bins, patches = ax2.hist(E, 20,normed=1,  facecolor='green', alpha=0.15,label=r'$\mathrm{(E_{fit}-E_{data})/dE_{data}}$')
                        #E = np.insert(interpolate.splev(Beta[1:],spline) - np.concatenate((main_sig,main_sigh))*dE[1:],0,E[0])
                        #ax2.plot(Beta, E1-E, color = colors[i+1])#, label = r'$\mathrm{fit}$')

            Z,BetaZ = GetZ(Beta,E,dE,offset)
            if args.partition:
                        Zn = unumpy.nominal_values(Z)
                        Zd = unumpy.std_devs(Z)
                        #if i==0: ax2.errorbar(BetaZ, Zn, Zd, linestyle = '-', color = colors[i], label=r'$\mathrm{r=%0.2d;\, A=%s}$' %(geom[1],geom[0]))
                        #else:    ax2.errorbar(BetaZ, Zn, Zd, linestyle = '-', color = colors[i])
                        FancyErrorbar(ax2,BetaZ,Z,\
                              colors[i],r'$\mathrm{r=%0.2d;\, A=%s}$' %(geom[1],geom[0]))
                        print Z[-1] 
            if  args.renyi or args.mutual:
                if  geom == ('full',1): Zf1 = Z 
                else:
                    S2  = -Z + 2*Zf1  
                    if   geom == ('half',2): S2A  = S2B = S2
                    elif geom == ('A',2):    S2A  = S2   #+ np.log(2)*(int(Lx)**2)*0.5  
                    elif geom == ('B',2):    S2B  = S2   #+ np.log(2)*(int(Lx)**2)*0.5
                    elif geom == ('full',2): S2AB = S2   #+ np.log(2)*(int(Lx)**2)

                    #FancyErrorbar(ax1,np.array(BetaZ),((S2-S2[-1])*1.0)/float(Lx),\
                    #      colors[i],((r'$\mathrm{Lx=%2.0d}$' %Lx) +'\n'+ (r'$\mathrm{S_{2%s}}$' %geom[0] )))
                    if  args.renyi:
                        S2n = unumpy.nominal_values(((S2-S2[-1])*1.0)/float(Lx))
                        S2d = unumpy.std_devs(((S2-S2[-1])*1.0)/float(Lx))
                        ax1.errorbar(np.array(BetaZ),S2n,S2d,
                              color = colors[i],label = ((r'$\mathrm{Lx=%2.0d}$' %Lx) +'\n'+ (r'$\mathrm{S_{2%s}}$' %geom[0] )))

            j += 1
            i += 1
        #T, MI = loadtxt("MI_ED.dat", unpack=True)
        #ax.plot(1.0/np.array(T), MI/float(Lx),\
        #        marker='',color=colors[-1],\
        #        label=r'$\mathrm{ED}$')
            
            
       
        if args.renyi:
            dd = 2+2
        if args.mutual:
            I = S2A+S2B-S2AB
            FancyErrorbar(ax,np.array(BetaZ),(I*1.0)/float(Lx),\
                          colors[i],r'$\mathrm{L=%2.0d}$' %Lx)
            rstd = I*1.0/float(Lx)
            #ax.plot(BetaZ,unumpy.std_devs(rstd)/unumpy.nominal_values(rstd), marker='s', color = colors[i],label = r'$\mathrm{Lx=%2.0d}$' %Lx)
            print Lx, BetaZ[61], (unumpy.std_devs(rstd)[1:]/unumpy.nominal_values(rstd)[1:])[57], (unumpy.std_devs(rstd)[1:]/unumpy.nominal_values(rstd)[1:])[-1]
        i += 1

    if  args.mutual:
        ax.legend(loc='best',frameon=False)
    if args.renyi:
       ax1.legend(loc='best',frameon=False)
    if args.energy:
       ax2.legend(loc='best',frameon=False)
    tight_layout()
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

