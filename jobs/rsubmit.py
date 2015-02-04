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
from optparse import OptionParser

# -----------------------------------------------------------------------------
def getPIMCommand(fname):
    ''' Get the command string from the log file.'''
    # Open the log file and get the command string
    logFile = open(fname,'r')
    logLines = logFile.readlines()
    line = logLines[2]

    return line[2:]

# -----------------------------------------------------------------------------
# Begin Main Program
# -----------------------------------------------------------------------------
def main():

    # setup the command line parser options 
    parser = OptionParser() 
    parser.add_option("-T", "--temperature", dest="T", type="float", \
            help="simulation temperature in Kelvin") 
    parser.add_option("-N", "--number-particles", dest="N", type="int",\
            help="number of particles") 
    parser.add_option("-n", "--density", dest="n", type="float",\
            help="number density in Angstroms^{-d}")
    parser.add_option("-P", "--number-time-slices", dest="P", type="int",\
            help="number of time slices")
    parser.add_option("-u", "--chemical-potential", dest="mu", type="float",\
            help="chemical potential in Kelvin") 
    parser.add_option("-L", "--length", dest="L", type="float",\
            help="length in Angstroms") 
    parser.add_option("-t", "--imag-time-step", dest="tau", type="float",\
            help="imaginary time step") 
    parser.add_option("-V", "--volume", dest="V", type="float",
            help="volume in Angstroms^d") 
    parser.add_option("-R", "--radius", dest="R", type="float",
            help="radius in Angstroms") 
    parser.add_option('-f', dest="folder", 
            help="path to the folder containing state files to be resubmitted. Default = ./OUTPUT")
    parser.add_option("-r",type=str, dest='run',
            help="optional run number that will be added to the script's name")
    parser.add_option("--canonical", action="store_true", dest="canonical", 
            help="are we in the canonical ensemble?")
    parser.add_option("-i", "--id", action="append", dest="pimcID", type="int",\
            help="a list of PIMC ID numbers to include")
    parser.add_option("-e", "--exclude", action="append", dest="exID", type="int",\
            help="a list of PIMC ID numbers to exclude")
    parser.add_option("-c", "--cluster", dest="cluster", choices=['westgrid','sharcnet','scinet','clumeq','bluemoon'],\
            help="target cluster: [westgrid,sharcnet,scinet,clumeq,bluemoon]") 
    parser.add_option("--load", action="store_true", 
            help="do we want to load with -s flag? Default: false")

    parser.set_defaults(canonical=False)
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

    # Check that we are in the correct ensemble
    pimchelp.checkEnsemble(options.canonical)

    # Get the data string and create the pimc helper object
    dataName = pimchelp.getFileString(options,reduce=False)
    pimc = pimchelp.PimcHelp(dataName,options.canonical)

    # Fill up the simulation parameters maps
    pimc.getSimulationParameters()

    # Delete those pimcIDs that do not satysfy parameters 
    # that are not contained in the pimc output filanames' structure
    #"implicit" parameters
    pimc.ApplyImplicitParameters()

    # We get either all the log files in the requested folder, or just the
    # requested files by their ID number
    if not options.pimcID:
        logFileNames   = pimc.getFileList('log',pimc.id)
    else:
        logFileNames  = pimc.getFileList('log',idList=options.pimcID) 


    # If we have excluded any ID's we remove them from the list
    if options.exID:
        for id in options.exID:
            for n,fname in enumerate(logFileNames):
                if int(id) == pimc.getID(fname):
                    logFileNames.pop(n)

    #Load sumbission command lines from the log files 
    SubmitData = []
    i = 0
    for logFile in logFileNames:   
        if os.path.getsize(logFile) != 0: 
           PIMCommand = getPIMCommand(logFile).rstrip('\n')
           #in case we want to start a new job loading it with state files
           #we need to replace -R with -s flag
           if options.load:
              PIMCid = logFile[-13:-4]
              #get the name of the state file corresponding to logFile
              lsCommand = 'ls -1 *state*%s*' % (PIMCid)
              stateFile = os.popen(lsCommand).read().split('\n')[0]
              #print 'state %s\n' % stateFile
              #modify the submit command
              PIMCommand = PIMCommand.replace('-R '+PIMCid,'-s '+options.folder+'/'+stateFile)
              #print 'submit command %s\n' % PIMCommand
           SubmitData.append(PIMCommand)
        else:
            print 'Cannot resubmit PIMCid with an empty %s' %logFile  



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
