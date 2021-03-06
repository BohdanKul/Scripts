from optparse import OptionParser
import pyutils as py
import numpy as np
import scipy.optimize as optimize
import matplotlib.pylab as pl
import matplotlib
from scipy.optimize import curve_fit
from matplotlib.ticker import AutoMinorLocator
import loadgmt,kevent
from scipy import special
from utils import Bootstrap
from mplrc import *
import mplrc

def HalfLog(x,a,c):
    rho_s      = 0.26974
    c_sw       = 1.1348
    #gamma_free = c
    #gamma = gamma_free + 0.5*np.log(2.0*np.pi) 
    #c = 0.5*np.log(2.0*rho_s/c_sw)+gamma
    return a+0.5*np.log(rho_s*x/c_sw)/(1.0*x)+c/(1.0*x)

def HalfLog2(x,a,c,d):
    rho_s      = 0.26974
    c_sw       = 1.1348
    #gamma_free = 0.014
    #gamma_free = c
    #gamma = gamma_free + 0.5*np.log(2.0*np.pi) 
    #c = 0.5*np.log(2.0*rho_s/c_sw)+gamma
    #return a+0.5*np.log(x)/(1.0*x)+c/(1.0*x) + d/(1.0*x*x)
    return a + 0.5*np.log(rho_s*x/c_sw)/(1.0*x) + c/(1.0*x) + d/(1.0*x*x)

def OneLog(x,a,c):
    return a+1.0*np.log(x)/(1.0*x)+c/(1.0*x)

def ZeroLog(x,a,c):
    return a+c/(1.0*x)

def FreeLog(x,a,b,c):
    b = b/2.0
    return a + b*np.log(x)/(1.0*x) + c/(1.0*x)

def FreeLog2(x,a,b,c,d):
    rho_s      = 0.26974
    c_sw       = 1.1348
    #gamma_free = 0.014
    #gamma_free = c
    #gamma = gamma_free + 0.5*np.log(2*np.pi) 
    #c = b*np.log(rho_s/c_sw)+gamma
    #return a + b*np.log(x)/(1.0*x) + c/(1.0*x) + d/(1.0*x*x)
    print 'yo'
    return a + b*np.log(rho_s*x/c_sw)/(1.0*x) + c/(1.0*x) + d/(1.0*x*x)

def MetaFreeLog(x,a,c):
    global bglob
    return a + bglob*np.log(x)/(1.0*x) + c/(1.0*x)

def main():
    parser = OptionParser(description='Fitting to a noisy data generated by a known function')
    parser.add_option("--npoints", type="int",   help="number of data points") 
    parser.add_option("--low",     type="float", help="smallest data point") 
    parser.add_option("--high",    type="float", help="highest data point") 
    parser.add_option("--sigma",   type="float", help="std of noise") 
    (options, args) = parser.parse_args() 

    pl.figure(1,(8,3))
    pl.rcParams.update(mplrc.aps['params'])
    ax = pl.subplot(1,1,1)
    ax.set_xlabel(r'$L_y$')
    ax.set_ylabel(r'$I_2^{A}/L_y$')
    minorLocator   = AutoMinorLocator(5)                                        
    ax.xaxis.set_minor_locator(minorLocator)                                    
    minorLocator   = AutoMinorLocator(5)                                        
    ax.yaxis.set_minor_locator(minorLocator)

    choice = 4 #4#0 
    off  = 0

    FitFuncs = [FreeLog,FreeLog2,OneLog,HalfLog,HalfLog2,ZeroLog]
    FitEqs   = [r'0.5N_G log(L)/L+',r'b \, log(L)/L+',r'log(L)/L+',r'0.5log(L)/L+',r'0.5log(\frac{\rho_s}{c}L)/L+','']
    for i,fit in enumerate(FitEqs):
        if   choice == 1: FitEqs[i] = r'a+'+fit+r'c(b,\gamma_{free})/L'
        elif choice == 3: FitEqs[i] = r'a+'+fit+r'c(\gamma_{free})/L'
        elif choice == 4: FitEqs[i] = r'a+'+fit+r'\gamma_{ord}/L'
        else:             FitEqs[i] = r'a+'+fit+r'd/L'
    FitEqs[1] = FitEqs[1] + r'+d/L^2'
    FitEqs[4] = FitEqs[4] + r'+d/L^2'
    
    FitFunc = FitFuncs[choice]
    FitEq  = FitEqs[choice]
    if  choice == 0: 
        clab  = ['a','N_G','d']
        guess = np.ones(3)
    elif choice == 1:
         clab  = ['a','b',r'\gamma_{free}','d']
         guess = np.ones(4)
    elif choice == 3:
         clab  = ['a',r'\gamma_{free}']
         guess = np.ones(2)
    elif choice == 4:
         clab  = ['a',r'\gamma_{ord}','d']
         guess = np.ones(3)
    else:
        clab = ['a','c']
        guess = np.ones(2)
    ax.set_title(r'$Fit\, to\, the\, form:\, %s$' %FitEq)
    pl.connect('key_press_event',kevent.press)
    
    sigma = options.sigma    

    off = 1
    offh = 8
    beta = 184

  
 
    gammae_exp = []
    gammav_exp = []
    #----------------------------------------------------------------------------------------------------------------------------------------------
    #1/8-------------------------------------------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------------------------
    off = 0
    Ls   = np.array([8,16,24,32])
    #MIs  = np.array([0.20170987004291954, 0.16195953475521396, 0.14235010902049608, 0.13113092651426075]) 
    #dMIs = np.array([0.00019319572927642786, 0.00023402606005673606, 0.00024543616376892271, 0.00019744156940868488])
    MIs  = np.array([0.20170987004291954, 0.16195953475521396, 0.14235010902049608, 0.13107297667495113]) 
    dMIs = np.array([0.00019319572927642786, 0.00023402606005673606, 0.00024543616376892271, 0.00016942024648351241])
    
    print Ls
    coeff, var_matrix = curve_fit(FitFunc,Ls,MIs,guess,sigma=1.0/(dMIs**2))
    err = np.sqrt(np.diagonal(var_matrix))
    dof     = len(Ls) - len(coeff)
    chisq   = sum(((MIs-FitFunc(Ls,*coeff))/dMIs)**2)
    cdf     = special.chdtrc(dof,chisq)
    lab = 'Fit: '
    if  choice == 4: 
        gammav_exp.append(coeff[-2])
        gammae_exp.append(err[-2])
    for i in range(len(guess)):
        lab += r'$%s = %0.3f(%0.3f);\,$' %(clab[i],coeff[i],err[i])
    lab += r'$\chi^2/DOF=%0.2f;\,$' %(chisq/float(dof))
    lab += r'$\, p-value = %0.2f$' %cdf
    #pl.errorbar(Ls,MIs,dMIs,ls='',label=r'$\beta =%3.0f$'%beta,color='y')
    pl.errorbar(Ls,MIs,dMIs,ls='',label=r'$Data\, for\, A/(L_xL_y)=1/8$',color=fcolors(4))
    nLs = np.linspace(Ls[0],Ls[-1],100)
    pl.plot(nLs,FitFunc(nLs,*coeff),label=lab,color=fcolors(4))

    #----------------------------------------------------------------------------------------------------------------------------------------------
    #2/8------------------------------------------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------------------------
    off = 1
    Ls   = np.array([4,8, 12, 16, 20, 24, 28,32])[off:]
    #MIs  = np.array([ 0.27140758767953033, 0.22113400442581177, 0.19095315323315157, 0.17249654628602423, 0.15946872869885381, 0.14956906776848633, 0.14271734664170935, 0.13648367829389965])[off:]
    #dMIs = np.array([ 0.00035063105703284086, 0.00028390662712597328, 0.00030945156199167442, 0.00032588203068154991, 0.00021749486902703514, 0.00034962902498649281, 0.00033862852457444855,0.00027553624634053158])[off:]
    MIs  = np.array([ 0.27140758767953033, 0.22113400442581177, 0.19095315323315157, 0.17249654628602423, 0.15946872869885381, 0.14956906776848633, 0.14271734664170935, 0.13639189179498032])[off:]
    dMIs = np.array([ 0.00035063105703284086, 0.00028390662712597328, 0.00030945156199167442, 0.00032588203068154991, 0.00021749486902703514, 0.00034962902498649281, 0.00033862852457444855,0.0002417771731823463])[off:]

    coeff, var_matrix = curve_fit(FitFunc,Ls,MIs,guess,sigma=1.0/(dMIs**2))
    err = np.sqrt(np.diagonal(var_matrix))
    dof     = len(Ls) - len(coeff)
    chisq   = sum(((MIs-FitFunc(Ls,*coeff))/dMIs)**2)
    cdf     = special.chdtrc(dof,chisq)
    lab = 'Fit: '
    if  choice == 4: 
        gammav_exp.append(coeff[-2])
        gammae_exp.append(err[-2])
    for i in range(len(guess)):
        lab += r'$%s = %0.3f(%0.3f);\,$' %(clab[i],coeff[i],err[i])
    lab += r'$\chi^2/DOF=%0.2f;\,$' %(chisq/float(dof))
    lab += r'$\, p-value = %0.2f$' %cdf
    #pl.errorbar(Ls,MIs,dMIs,ls='',label=r'$\beta =%3.0f$'%beta,color='b')
    pl.errorbar(Ls,MIs,dMIs,ls='',label=r'$Data\, for\, A/(L_xL_y)= 1/4$',color=fcolors(6))
    nLs = np.linspace(Ls[0],Ls[-1],100)
    pl.plot(nLs,FitFunc(nLs,*coeff),label=lab,color=fcolors(6))



    #----------------------------------------------------------------------------------------------------------------------------------------------
    #3/8-------------------------------------------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------------------------
    off = 0
    Ls   = np.array([8,16,24,32])
    #MIs  = np.array([0.22996345403546045, 0.17674766273907488, 0.15299083379548556, 0.13889845232701539]) 
    #dMIs = np.array([0.00035023175479636717, 0.00039822181187702241, 0.00037804621953715169, 0.00032578195392620591])
    MIs  = np.array([0.22996345403546045, 0.17674766273907488, 0.15299083379548556, 0.13883803379406248]) 
    dMIs = np.array([0.00035023175479636717, 0.00039822181187702241, 0.00037804621953715169, 0.00029383814706680241])
    
    coeff, var_matrix = curve_fit(FitFunc,Ls,MIs,guess,sigma=1.0/(dMIs**2))
    err = np.sqrt(np.diagonal(var_matrix))
    dof     = len(Ls) - len(coeff)
    chisq   = sum(((MIs-FitFunc(Ls,*coeff))/dMIs)**2)
    cdf     = special.chdtrc(dof,chisq)
    lab = 'Fit: '
    if  choice == 4: 
        gammav_exp.append(coeff[-2])
        gammae_exp.append(err[-2])
    for i in range(len(guess)):
        lab += r'$%s = %0.3f(%0.3f);\,$' %(clab[i],coeff[i],err[i])
    lab += r'$\chi^2/DOF=%0.2f;\,$' %(chisq/float(dof))
    lab += r'$\, p-value = %0.2f$' %cdf
    #pl.errorbar(Ls,MIs,dMIs,ls='',label=r'$\beta =%3.0f$'%beta,color='y')
    pl.errorbar(Ls,MIs,dMIs,ls='',label=r'$Data\, for\, A/(L_xL_y)=3/8$',color=fcolors(3))
    nLs = np.linspace(Ls[0],Ls[-1],100)
    pl.plot(nLs,FitFunc(nLs,*coeff),label=lab,color=fcolors(3))

    #----------------------------------------------------------------------------------------------------------------------------------------------
    # 4/8 ---------------------------------------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------------------------
    Ls   = np.array([8,10, 12, 16, 20, 24,28,32])
    #MIs  = np.array([0.23251280150358983, 0.21337256591279627,0.19797531796757287, 0.17775003649902879, 0.16411815945015565, 0.15374365325070799, 0.14583169547075336, 0.13937802905214949])
    #dMIs = np.array([0.00040395552486499571, 0.00062070973509060503,0.00043835040610511769, 0.00046104739832369675, 0.00030739428405506897, 0.00043982489905168917, 0.00047582622823638226, 0.00037807728176096899])
    MIs  = np.array([0.23251280150358983, 0.21337256591279627,0.19797531796757287, 0.17775003649902879, 0.16411815945015565, 0.15374365325070799, 0.14583169547075336, 0.13931301552362246])
    dMIs = np.array([0.00040395552486499571, 0.00062070973509060503,0.00043835040610511769, 0.00046104739832369675, 0.00030739428405506897, 0.00043982489905168917, 0.00047582622823638226, 0.00033882833155090007])
    coeff, var_matrix = curve_fit(FitFunc,Ls,MIs,guess,sigma=1.0/(dMIs**2))
    err = np.sqrt(np.diagonal(var_matrix))
    dof     = len(Ls) - len(coeff)
    chisq   = sum(((MIs-FitFunc(Ls,*coeff))/dMIs)**2)
    cdf     = special.chdtrc(dof,chisq)
    lab = 'Fit: '
    if  choice == 4: 
        gammav_exp.append(coeff[-2])
        gammae_exp.append(err[-2])
    for i in range(len(guess)):
        lab += r'$%s = %0.3f(%0.3f);\,$' %(clab[i],coeff[i],err[i])
    lab += r'$\chi^2/DOF=%0.2f;\,$' %(chisq/float(dof))
    lab += r'$\, p-value = %0.2f$' %cdf
    #pl.errorbar(Ls,MIs,dMIs,ls='',label=r'$\beta =%3.0f$'%beta,color='y')
    pl.errorbar(Ls,MIs,dMIs,ls='',label=r'$Data\, for\, A/(L_xL_y)=0.50$',color=fcolors(2))
    nLs = np.linspace(Ls[0],Ls[-1],100)
    pl.plot(nLs,FitFunc(nLs,*coeff),label=lab,color=fcolors(2))
    pl.tight_layout()

    #ax.set_ylim([0.13,0.24]) 
    #ax.set_xlim([7.5,32.5]) 
    lgd = pl.legend()
    lgd.draw_frame(False)    


    if  choice == 4:
        geom_exp  = [0.125,0.250,0.375,0.500]
        geom_thr  = [0.100,0.125,0.1500,0.2000,0.2500,0.3000,0.3500,0.3750,0.4000,0.4500,0.5000]
        gamma_thr = [0.598,0.672,0.7255,0.8008,0.8512,0.8866,0.9116,0.9209,0.9283,0.9379,0.9410]
        print np.array(gamma_thr) - 0.5*np.log(2.0*np.pi)
        pl.figure(2,(8,3))
        pl.connect('key_press_event',kevent.press)
        pl.rcParams.update(mplrc.aps['params'])
        ax = pl.subplot(1,1,1)
        ax.set_xlabel(r'$A/(L_xL_y)$')
        ax.set_ylabel(r'$\gamma_{ord}$')
        minorLocator   = AutoMinorLocator(5)                                        
        ax.xaxis.set_minor_locator(minorLocator)                                    
        minorLocator   = AutoMinorLocator(5)                                        
        ax.yaxis.set_minor_locator(minorLocator)
        ax.plot(geom_thr,gamma_thr,color=colors[1],label=r"$Theory$")
        ax.errorbar(geom_exp,gammav_exp,gammae_exp,
                    ls='', color=colors[2], label=r"$MC$")
        lgd = pl.legend()
        lgd.draw_frame(False)    
        lgd.draggable(state=True)
        ax.set_xlim([0.1,0.51]) 
        pl.tight_layout()
     
    pl.show()
if __name__ == "__main__":
   main()
