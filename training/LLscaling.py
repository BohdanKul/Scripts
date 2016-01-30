# Visualize best trained KL divergence between classical and quantum models. 
# KL-divergence is extracted from the corresponding training files.

import Hfile, bmachine, kevent, ssexyhelp
import argparse, collections
from uncertainties import unumpy
import matplotlib as mpl
import numpy as np
import pylab as pl

from matplotlib import pyplot as plt

# ----------------------------------------------------------------------
def main(): 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--quant',  '-Q',  help='Quantum training file',   nargs='+')
    parser.add_argument('--squant1', '-S', help='Fixed delta quantum training file',   nargs='+')
    parser.add_argument('--aquant1',  help='Approximate quantum training file',   nargs='+')
    parser.add_argument('--squant2', help='Fixed delta quantum training file',   nargs='+')
    parser.add_argument('--squant3', help='Fixed delta quantum training file',   nargs='+')
    parser.add_argument('--aquant2',  help='Approximate quantum training file',   nargs='+')
    parser.add_argument('--aquant3',  help='Approximate quantum training file',   nargs='+')
    parser.add_argument('--class','-C', help='Classical training file', nargs='+')

    args = vars(parser.parse_args())

    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]
    fig = pl.figure(1, figsize=(10,5))
    pl.connect('key_press_event',kevent.press)
    ax  = pl.subplot(111)

    
    # ----------------------------------------------------------------------
    if 'class' in args.keys():
        cdata = {}
        i = -1
        for filename in args['class']:
            data = Hfile.cleanData(np.loadtxt(filename))
            
            LL      = data[-1,0]
            fparams = ssexyhelp.getReduceParamMap(filename)
            N       = int(fparams['N'])
            
            if not(N in cdata.keys()): cdata[N] = []
            cdata[N] += [LL]
        
        (Ns, vLLs, eLLs) = (cdata.keys(), [], [])
        Ns.sort()
        for N in Ns:
            vLLs += [np.mean(np.array(cdata[N]))]
            eLLs += [np.std(np.array(cdata[N]))/float(np.sqrt(len(cdata[N])))]
        
        cLLs = unumpy.uarray((np.array(vLLs), np.array(eLLs)))
        #ax.errorbar(Ns, unumpy.nominal_values(cLLs), unumpy.std_devs(cLLs), color=colors[1], lw=1, marker='o', ms=4, ls='', label=r'$classical$')

    # ----------------------------------------------------------------------
    if 'quant' in args.keys():
        qdata = {}
        i = -1
        for filename in args['quant']:
            data = Hfile.cleanData(np.loadtxt(filename))
            
            LL      = data[-1,0]
            fparams = ssexyhelp.getReduceParamMap(filename)
            N       = int(fparams['N'])
            
            if not(N in qdata.keys()): qdata[N] = []
            qdata[N] += [LL]

        (Ns, vLLs, eLLs) = (qdata.keys(), [], [])
        Ns.sort()
        for N in Ns:
            vLLs += [np.mean(np.array(qdata[N]))]
            eLLs += [np.std(np.array(qdata[N]))/float(np.sqrt(len(qdata[N])))]

        eqLLs = cLLs - unumpy.uarray((np.array(vLLs), np.array(eLLs)))
        ax.errorbar(Ns, unumpy.nominal_values(eqLLs), unumpy.std_devs(eqLLs), color=colors[0], lw=1, ls=':', marker='s', ms=4, label=r'$exact \, quantum$')
        #ax.plot(Ns, unumpy.nominal_values(eqLLs), color=colors[2], lw=1, ls='-', marker='s', ms=4, label=r'$exact \, quantum$')


    # ----------------------------------------------------------------------
    if 'aquant1' in args.keys():
        aqdata = {}
        i = -1
        for filename in args['aquant1']:
            data = Hfile.cleanData(np.loadtxt(filename))
            
            LL      = data[-1,0]
            fparams = ssexyhelp.getReduceParamMap(filename)
            N       = int(fparams['N'])
            
            if not(N in aqdata.keys()): aqdata[N] = []
            aqdata[N] += [LL]

        (Ns, vLLs, eLLs) = (aqdata.keys(), [], [])
        Ns.sort()
        for N in Ns:
            vLLs += [np.mean(np.array(aqdata[N]))]
            eLLs += [np.std(np.array(aqdata[N]))/float(np.sqrt(len(aqdata[N])))]
        
        aqLLs = cLLs - unumpy.uarray((np.array(vLLs), np.array(eLLs)))
        ax.errorbar(Ns, unumpy.nominal_values(aqLLs), unumpy.std_devs(aqLLs), color=colors[1], lw=1, ls='-', marker='p', ms=4, label=r'$appr. \,quantum \, \Delta=1$')
        #ax.plot(Ns, unumpy.nominal_values(aqLLs), color=colors[3], lw=1, ls='-', marker='p', ms=4, label=r'$appr. \,quantum \, \Delta=2$')

    # ----------------------------------------------------------------------
    if 'squant1' in args.keys():
        qdata = {}
        i = -1
        for filename in args['squant1']:
            print filename
            data = Hfile.cleanData(np.loadtxt(filename))
            
            LL      = data[-1,0]
            fparams = ssexyhelp.getReduceParamMap(filename)
            N       = int(fparams['N'])
            
            if not(N in qdata.keys()): qdata[N] = []
            qdata[N] += [LL]

        (Ns, vLLs, eLLs) = (qdata.keys(), [], [])
        Ns.sort()
        for N in Ns:
            vLLs += [np.mean(np.array(qdata[N]))]

            eLLs += [np.std(np.array(qdata[N]))/float(np.sqrt(len(qdata[N])))]

        eqLLs = cLLs[:-1] - unumpy.uarray((np.array(vLLs), np.array(eLLs)))
        ax.errorbar(Ns, unumpy.nominal_values(eqLLs), unumpy.std_devs(eqLLs), color=colors[1], lw=1, ls='--', marker='s', ms=4, label=r'$fixed \, quantum \, \Delta=1$')
        #ax.plot(Ns, unumpy.nominal_values(eqLLs), color=colors[5], lw=1, ls='-', marker='s', ms=4, label=r'$fixed \, quantum \, \Delta =1$')

 
    # ----------------------------------------------------------------------
    if 'aquant2' in args.keys():
        aqdata = {}
        i = -1
        for filename in args['aquant2']:
            data = Hfile.cleanData(np.loadtxt(filename))
            
            LL      = data[-1,0]
            fparams = ssexyhelp.getReduceParamMap(filename)
            N       = int(fparams['N'])
            
            if not(N in aqdata.keys()): aqdata[N] = []
            aqdata[N] += [LL]

        (Ns, vLLs, eLLs) = (aqdata.keys(), [], [])
        Ns.sort()
        for N in Ns:
            vLLs += [np.mean(np.array(aqdata[N]))]
            eLLs += [np.std(np.array(aqdata[N]))/float(np.sqrt(len(aqdata[N])))]
        
        aqLLs = cLLs - unumpy.uarray((np.array(vLLs), np.array(eLLs)))
        ax.errorbar(Ns, unumpy.nominal_values(aqLLs), unumpy.std_devs(aqLLs), color=colors[2], lw=1, ls='-', marker='p', ms=4, label=r'$appr. \,quantum \, \Delta=2$')
        #ax.plot(Ns, unumpy.nominal_values(aqLLs), color=colors[3], lw=1, ls='-', marker='p', ms=4, label=r'$appr. \,quantum \, \Delta=2$')

    # ----------------------------------------------------------------------
    if 'squant2' in args.keys():
        qdata = {}
        i = -1
        for filename in args['squant2']:
            data = Hfile.cleanData(np.loadtxt(filename))
            
            LL      = data[-1,0]
            fparams = ssexyhelp.getReduceParamMap(filename)
            N       = int(fparams['N'])
            
            if not(N in qdata.keys()): qdata[N] = []
            qdata[N] += [LL]

        (Ns, vLLs, eLLs) = (qdata.keys(), [], [])
        Ns.sort()
        for N in Ns:
            vLLs += [np.mean(np.array(qdata[N]))]

            eLLs += [np.std(np.array(qdata[N]))/float(np.sqrt(len(qdata[N])))]

        eqLLs = cLLs[:-1] - unumpy.uarray((np.array(vLLs), np.array(eLLs)))
        ax.errorbar(Ns, unumpy.nominal_values(eqLLs), unumpy.std_devs(eqLLs), color=colors[2], lw=1, ls='--', marker='s', ms=4, label=r'$fixed \, quantum \, \Delta=3$')
        #ax.plot(Ns, unumpy.nominal_values(eqLLs), color=colors[5], lw=1, ls='-', marker='s', ms=4, label=r'$fixed \, quantum \, \Delta =1$')

 
    # ----------------------------------------------------------------------
    if 'aquant3' in args.keys():
        aqdata = {}
        i = -1
        for filename in args['aquant3']:
            data = Hfile.cleanData(np.loadtxt(filename))
            
            LL      = data[-1,0]
            fparams = ssexyhelp.getReduceParamMap(filename)
            N       = int(fparams['N'])
            
            if not(N in aqdata.keys()): aqdata[N] = []
            aqdata[N] += [LL]

        (Ns, vLLs, eLLs) = (aqdata.keys(), [], [])
        Ns.sort()
        for N in Ns:
            vLLs += [np.mean(np.array(aqdata[N]))]
            eLLs += [np.std(np.array(aqdata[N]))/float(np.sqrt(len(aqdata[N])))]
        
        aqLLs = cLLs - unumpy.uarray((np.array(vLLs), np.array(eLLs)))
        ax.errorbar(Ns, unumpy.nominal_values(aqLLs), unumpy.std_devs(aqLLs), color=colors[3], lw=1, ls='-', marker='p', ms=4, label=r'$appr. \,quantum \, \Delta=3$')
        #ax.plot(Ns, unumpy.nominal_values(aqLLs), color=colors[3], lw=1, ls='-', marker='p', ms=4, label=r'$appr. \,quantum \, \Delta=2$')

    # ----------------------------------------------------------------------
    if 'squant3' in args.keys():
        qdata = {}
        i = -1
        for filename in args['squant3']:
            data = Hfile.cleanData(np.loadtxt(filename))
            
            LL      = data[-1,0]
            fparams = ssexyhelp.getReduceParamMap(filename)
            N       = int(fparams['N'])
            
            if not(N in qdata.keys()): qdata[N] = []
            qdata[N] += [LL]

        (Ns, vLLs, eLLs) = (qdata.keys(), [], [])
        Ns.sort()
        for N in Ns:
            vLLs += [np.mean(np.array(qdata[N]))]

            eLLs += [np.std(np.array(qdata[N]))/float(np.sqrt(len(qdata[N])))]

        eqLLs = cLLs[:-1] - unumpy.uarray((np.array(vLLs), np.array(eLLs)))
        ax.errorbar(Ns, unumpy.nominal_values(eqLLs), unumpy.std_devs(eqLLs), color=colors[3], lw=1, ls='--', marker='s', ms=4, label=r'$fixed \, quantum \, \Delta=3$',capthick=2.0, capsize=3.0, markersize='5')
        #ax.plot(Ns, unumpy.nominal_values(eqLLs), color=colors[5], lw=1, ls='-', marker='s', ms=4, label=r'$fixed \, quantum \, \Delta =1$')

    
    pl.xlabel(r'$System \, size$')
    pl.ylabel(r'$\Delta LL$')
    ax.set_xlim([3,11])
    ax.set_yscale('log')
    lgd = pl.legend(loc = 'best')
    lgd.draggable(state=True)
    lgd.draw_frame(False)
    pl.tight_layout()
    pl.show()
   
    #diff = np.array(cLL) - np.array(qLL)
    #plt.hist(diff, 10, normed=1, alpha=0.75)
    #plt.show()
   
# ------------------------------------------------------------ ----------
# ------------------------------------------------------------ ----------
if __name__ == "__main__": 
    main()

