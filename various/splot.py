#!/usr/bin/python
# pimcplot.py
# Adrian Del Maestro
# 07.20.2009
# 
# Plot rough estimators vs. MC Bins for files supplied as input

import matplotlib
#matplotlib.use('TKAgg')

from uncertainties import ufloat 
from uncertainties import umath 
import MCstat
import zipfile
import os,sys
import pyutils
import loadgmt,kevent
from pylab import *
import argparse
import pimchelp
import mplrc 
import numpy as np
import MCstat

# ----------------------------------------------------------------------
def getStats(data,dim=0):
    ''' Get the average and error of all columns in the data matrix. '''

    if ndim(data) > dim:
        numBins  = size(data,dim) 
        dataAve  = average(data,dim) 
        dataAve2 = average(data*data,dim) 
        bins = MCstat.bin(data) 
        dataErr = amax(bins,axis=0)
        dataErr2 = sqrt( abs(dataAve2-dataAve**2)/(1.0*numBins-1.0) ) 

#        try:
#            bins = MCstat.bin(data) 
#            dataErr = amax(bins,axis=0)
#        except:
#            dataErr   = sqrt( abs(dataAve2-dataAve**2)/(1.0*numBins-1.0) ) 
    else:
        dataAve = data
        dataErr = 0.0*data

    return dataAve,dataErr




# -----------------------------------------------------------------------------
# Averaging functions
# -----------------------------------------------------------------------------
def cumulativeMovingAverage(data):
    '''Compute the cumulative mean as a function of bin index.'''
    CMA = zeros_like(data)
    CMA[0] = data[0]
    for n in range(len(data)-1):
        CMA[n+1] = (data[n+1] + n*CMA[n])/(1.0*(n+1))
    return CMA

# -----------------------------------------------------------------------------
def simpleMovingAverage(period,data):
    assert period == int(period) and period > 0, "Period must be an integer >0"
  
    #avoid averaging over a larger number of elements than there is 
    period = min(size(data,0),period)

    import numpy as np
    weightings = np.repeat(1.0, period) / period 
    return np.convolve(data, weightings, 'valid')
    #return array([sum(data[i:i+period])/(1.0*period) for i in xrange(len(data)-period+1)])

def GetFileParams(fileName):
    keys = ['Nx','Ny','T']
    [Nx,Ny,T] = fileName.split('-')[1:4] 
    return dict(zip(keys,[int(Nx),int(Ny),float(T)]))

def GetHeaderNumber(fileName,header):
    ffile = open(fileName,'r')
    headers = ffile.readline().lstrip('#').split()
    col = list([headers.index(header)])
    return col

def H(p):
    return -p*np.log(p)/np.log(2.0) - (1.0-p)*np.log(1.0-p)/np.log(2.0)
# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='Plot Raw MC Equilibration Data for Scalar Estimators.')
    parser.add_argument('fileNames', help='Scalar estimator files', nargs='+')
    parser.add_argument('--estimator','-e', help='A list of estimator names that \
                        are to be plotted.', type=str)
    parser.add_argument('--skip','-s', help='Number of measurements to be skipped \
                        in the average plot.', type=int, default=0)
    parser.add_argument('--period','-p', help='Period of the simple moving \
                        average. default=50', type=int, default=50)

    args = parser.parse_args()
    fileNames = args.fileNames
    if len(fileNames) < 1:
        parser.error("Need to specify at least one scalar estimator file")
    
    toSort = True
    for fileName in fileNames:
        if fileName.find('/') != -1:
           toSort = False
    if toSort:     
       fileNames.sort()
    
    ffile = open(fileNames[0],'r')
    headers = ffile.readline().lstrip('#').split()
    headers += ['srt']
    print headers

    # If we don't choose an estimator, provide a list of possible ones
    if (not args.estimator) or (args.estimator not in headers):
        errorString = "Need to specify one of:\n"      
        for head,index in headers.iteritems():
            errorString += "\"%s\"" % head + "   " 
        parser.error(errorString)



    numFiles = len(fileNames)

    #colors  = loadgmt.getColorList('cw/1','cw1-029',max(numFiles,2))
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
 
    fig = figure(1,figsize=(13,6))
    ax = subplot(111)
    connect('key_press_event',kevent.press)
   
    #rcParams.update(mplrc.aps['params'])
    n = 0
    E0 = 0
    Es  = []
    dEs = []
    srt = False
    if args.estimator=='srt': srt=True
    for i,fileName in enumerate(fileNames):
        dataFile  = open(fileName,'r');
        dataLines = dataFile.readlines();
        ffile = open(fileName,'r')
        headers = ffile.readline().lstrip('#').split()
        headers += ['srt']
        if len(dataLines) > 2:
            #params = GetFileParams(fileName)
            if  srt and not('ALRatio' in headers):
                col    = GetHeaderNumber(fileName,'nAred')
                dataR = loadtxt(fileName,usecols=col)
                col    = GetHeaderNumber(fileName,'nAext')
                dataE = loadtxt(fileName,usecols=col)
                    
                daveR = amax(MCstat.bin(dataR[args.skip:]),axis=0)
                aveR  = average(dataR[args.skip:])
                daveE = amax(MCstat.bin(dataE[args.skip:]),axis=0)
                aveE  = average(dataE[args.skip:])
                #S2 = -umath.log(ufloat(aveE,daveE)/ufloat(aveR,daveR))
                S2 = -1.0*umath.log(ufloat(aveE,daveE)/ufloat(aveR,daveR))
                print 'Entropy = ', S2
                Es  += [S2.n]
                dEs += [S2.s]
            else:
                if srt: col    = GetHeaderNumber(fileName,'ALRatio')
                else  : col    = GetHeaderNumber(fileName,args.estimator)
                data = loadtxt(fileName,usecols=col)

                ID = fileName[-14:-4]
                if  size(data) > 1:
                    sma = simpleMovingAverage(args.period,data[args.skip:])
                    ax.plot(sma,color=colors[i%len(colors)],linewidth=3,linestyle='-')
                    bins = MCstat.bin(data[args.skip:]) 
                    dataErr = amax(bins,axis=0)
                    dataAve = average(data[args.skip:])
                    if srt or (args.estimator=='ALRatio'):
                        S2 = -1.0*umath.log(ufloat(dataAve,dataErr)) 
                        print 'Entropy = ', S2
                        Es  += [S2.n]
                        dEs += [S2.s]
                    #ax.plot(bins,color=colors[i%len(colors)],linewidth=1,marker='None',linestyle='-')
                    #if E0 == 0: E0 = dataAve
                    print '%0.6f +/- %0.6f ' %(dataAve-E0,dataErr)
                    #print dataAve-E0
        else:
            print '%s contains no measurements' %fileName
    xlabel('MC bins (p=%s)' %args.period)
    ylabel(args.estimator)
    #tight_layout()
    print Es
    print dEs
    #legend() 
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()


