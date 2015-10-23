#!/usr/bin/python 
#
# rsubmit.py
# Adrian Del Maestro
# 09.04.2009
#
# Generate a torque submit script for the pimc code which creates
# a pbs file for various sets of parameters.  This restart script
# reads all log files in a current directory by default, or will
# take a list of id numbers from the command line

import sys,os,clusters
import pimchelp
from optparse  import OptionParser
import ssexyhelp

# -----------------------------------------------------------------------------
def getPIMCommand(fname):
    ''' Get the command string from the file name.'''

    params = ssexyhelp.getParamMap(fname)    
    command = './ssexy.e '
    for par,value in params.iteritems():
        if isinstance(value,int): command += "-%s %d "    %(par,value)
        else:                     command += "-%s %0.3f " %(par,value)
    command += "-%s %d " %('t',params['a']+16)
    return command

# -----------------------------------------------------------------------------
# Begin Main Program
# -----------------------------------------------------------------------------
def main():

    # setup the command line parser options 
    parser = OptionParser() 
    parser.add_option("-T", "--temperature", dest="T", type="float",
                      help="simulation temperature in Kelvin") 
    parser.add_option("-b", "--beta", dest="b", type="float",
                      help="number of particles") 
    parser.add_option("-d", "--delta", dest="d", type="float",
                      help="strength of SzSz interaction") 
    parser.add_option("-r", "--replica", dest = "r",type="int",
                      help="number of replica copies") 
    parser.add_option("-x", "--Lx", dest="x", type="int",
                      help="lattice width") 
    parser.add_option("-y", "--Ly", dest="y", type="int",
                      help="lattice height") 
    parser.add_option('-f', dest="folder", 
            help="path to the folder containing state files to be resubmitted. Default = ./OUTPUT")
    parser.add_option("--run",type=str,
            help="optional id that will be added to the submit script's name")
    parser.add_option("-e", "--exclude", action="append", dest="exID", type="int",\
            help="a list of PIMC ID numbers to exclude")
    parser.add_option("-c", "--cluster", dest="cluster", choices=['westgrid','sharcnet','scinet','clumeq','bluemoon'],\
            help="target cluster: [westgrid,sharcnet,scinet,clumeq,bluemoon]") 

    parser.set_defaults(load=False)
    parser.set_defaults(folder='OUTPUT')
    parser.set_defaults(run='')

    # parse the command line options and get the reduce flag
    (options, args) = parser.parse_args() 

    # remember the current path but switch temporary to the files location folder
    CurrentFolder  = os.getcwd() 
    options.folder = options.folder.rstrip('/')
    os.chdir(options.folder)

    if (not options.cluster):
        parser.error("need to specify a cluster")


    # Get the data string and create the pimc helper object
    ssexy = ssexyhelp.SSEXYHelp(options)
    ssexy.getSimulationParameters()
    if  (not ssexy.id): 
        print "No filenames detected satisfying those criteria"
        sys.exit()
    
    
    estFileNames  = ssexy.getFileList('estimator',idList=ssexy.id) 

    # If we have excluded any ID's we remove them from the list
    if options.exID:
        for id in options.exID:
            for n,fname in enumerate(logFileNames):
                if int(id) == ssexy.getID(fname):
                    logFileNames.pop(n)

    #Load sumbission command lines from the log files 
    SubmitData = []
    i = 0
    N = 0
    nN= 0

    #print As
    for estFile in estFileNames:   
        stateFile =  'state'+estFile[9:]
        fsize = len(open(estFile,'r').readlines())
        #if os.path.getsize(stateFile) != 0: 
        if fsize > 100:
            PIMCommand = getPIMCommand(estFile)
            PIMCommand += '-m %d ' %(1500)
            PIMCommand += '-s %s/%s' %(options.folder,stateFile)
            SubmitData.append(PIMCommand)

    #switch back to folder we used to be 
    os.chdir(CurrentFolder)
     
    # Now create the submission files
    if options.cluster == 'westgrid':
        clusters.westgrid(SubmitData,options.run)

    if options.cluster == 'sharcnet':
        clusters.sharcnet(SubmitData,options.run)

    if options.cluster == 'scinet':
        clusters.scinet(SubmitData,options.run)

    if options.cluster == 'clumeq':
        clusters.clumeq(SubmitData,options.run)

    if options.cluster == 'bluemoon':
        clusters.bluemoon(SubmitData,options.run)

# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
