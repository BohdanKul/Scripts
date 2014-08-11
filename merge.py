# merge.py
# Adrian Del Maestro
# 10.19.2009
# 
# Merge the results of parallel SSEXY output files, and move the 
# originals to an archive

import os,sys,glob,shutil
import ssexyhelp
from optparse import OptionParser
from pylab import *
import numpy as np
import MCstat

def ReadSkipping(fileName,skip=0):
    '''Read measurements from fileName while skipping first ${skip} of them'''
    
    # if fileName doesnt exist, exit
    if not(len(glob.glob(fileName)) > 0):
       return 0

    # Take the content from the file
    inFile = open(fileName,'r');
    inLines = inFile.readlines();
    inFile.close()

    # Strip any comment lines
    if len(inLines) != 0:
       inLines.pop(0)

    #Number of measurements after skipping
    numLines = len(inLines)-skip

    # check in case we want to skip more lines than there are
    if numLines > 0:
       for i in range(skip):
           inLines.pop(0)   
    else:
         inLines =[]
   
    return inLines, numLines      

def CreateFileSkipping(FromfileName,oldSSEXYID,newSSEXYID,skip=0,average=False,waverage=False):
    '''Create a new file whose oldSSEXYID is replaced with newSSEXYID and
      skip ${skip} first measurements
      Set average to True in order to write out just one average'''

    #Read and adapt the header
    inFile = open(FromfileName,'r');
    #firstLine  = inFile.readline().replace(str(oldSSEXYID),str(newSSEXYID));
    secondLine = inFile.readline()
    inFile.close();

    #Read the measurements we want to keep
    inList,numLines = ReadSkipping(FromfileName,skip)
    if len(inList) == 0:
       return '',0
    # get the output file name and open the file for writing
    outName = FromfileName.replace(str(oldSSEXYID),str(newSSEXYID))
    outFile = open('MERGED/' + outName,'w');
    print('To: MERGED/%-80s' % outName)
 
    #write the header
    #outFile.write(firstLine)
    secondLine = secondLine.rstrip('\n')
    if waverage:
       headers = secondLine.split()
       temp = '#%15s%16s' %(headers[1],'+/-')
       for header in headers[2:]:
           temp += '%16s%16s' %(header,'+/-')
       secondLine = temp 
    if average:
       secondLine += '%16s' %('Bins')
    
    
    outFile.write(secondLine+'\n')
    #if there is more measurements than we want to skip
    if numLines > 0:
       #write the content
       if average or waverage:
          data = np.loadtxt(inList)
          aves = np.average(data,0)
          print "From: %s" %FromfileName
          if  average:
              for ave in aves:
                  outFile.write('%16.8E' %ave)
              outFile.write('%16d' %len(data[0:]))    
          if  waverage:
              errs = amax(MCstat.bin(data),axis=0)
              for (ave,err) in zip(aves,errs):
                  outFile.write('%16.8E%16.8E' %(ave,err))
             
          outFile.write('\n')
       else:
           outFile.writelines(inList) 

    return outFile,numLines


def addFromToSkipping(FromfileName,outFile,skip=0,average=False,waverage=False):
    '''add SSEXY measurements from FromfileName to toFile, while skipping skip first measurements'''
     
    #Read the measurements we want to keep
    inList,numLines = ReadSkipping(FromfileName,skip)

    if len(inList) == 0:
       return 0

    #Write them
    if len(inList) != 0:
       if average or waverage:
          data = np.loadtxt(inList)
          aves = np.average(data,0)
          print "From: %s" %FromfileName
          if  average:
              for ave in aves:
                  outFile.write('%16.8E' %ave)
              outFile.write('%16d' %len(data[0:]))    
          if  waverage:
              errs = amax(MCstat.bin(data),axis=0)
              for (ave,err) in zip(aves,errs):
                  outFile.write('%16.8E%16.8E' %(ave,err))
             
          outFile.write('\n')
       else:
           outFile.writelines(inList) 

   
    return numLines



# -----------------------------------------------------------------------------
def mergeData(ssexy,type,newID,skip=0,isRestarted=False,ssexyIds=None,average=False,waverage=False):
    ''' Merge the results of the SSEXY data files to a single file. '''

    #file-names to merge
    fileNames = ssexy.getFileList(type,ssexyIds)
    print 'Merging files #: ',len(fileNames)
    print fileNames
    #Total number of measurements in the merged file
    numLines = 0
    #Have we created an output file?
    fileExist = False
  
    for i,fname in enumerate(fileNames):
        # Does the file exist?
        if len(glob.glob(fname)) > 0:
           #if we havent yet an output file, create one
           if not fileExist:
               print fname
               outFile,numLines = CreateFileSkipping(fname,ssexyIds[0],newID,skip,average,waverage)
               fileExist = (outFile != '') 
               #taken it is a restarted job, if there is less measurements than we want to skip in 
               #the first file we will need to skip more in the files that follow
               if (isRestarted) and (numLines < 0):
                  skip = -numLines   
               #if it is not a restarted job, but the file contains less than enough measurements
               if (not isRestarted) and (numLines < 0):
                  numLines = 0 
           else:
               #if we are merging files that got started from a parent state file
               #we skip measurements only from the parent  
               if isRestarted:
                   #if the parent has more measurements that we want to skip
                   if not(numLines < 0):
                       numLines += addFromToSkipping(fname,outFile,average,waverage)  
                   #if not, we will have to skip some measurements in children
                   else:
                       numLines = addFromToSkipping(fname,outFile,skip,average,waverage)   
                       skip = -numLines     
               else:
                   numLines += addFromToSkipping(fname,outFile,skip,average,waverage)   
                   if numLines < 0:
                      numLines = 0         
    outFile.close() 
    print('%10d' %numLines)


    return

#--------------------------------------------------------------------------------------------------------------------

def getNewSSEXYID(ssexy):
    ''' Find a new SSEXYID which is the average of the ones to merge, and make sure it doesn't already exist'''
    newID = 0
    for id in ssexy.id:
        newID += int(id)
    newID = int(newID/(1.0*len(ssexy.id)))
    # Now we keep incrementing the ID number until we are sure it is unique
    while ( (len(glob.glob('*est*%09d*' % newID)) > 0) or
            (len(glob.glob('MERGED/*est*%09d*' % newID)) > 0) or 
            (len(glob.glob('../*est*%09d*' % newID)) > 0) ):
        newID += 1
    return newID

#--------------------------------------------------------------------------------------------------------------------

def getMergeSets(ssexy,VarP):
    '''Create a dictionnary of sets of ssexyids to be merged'''
    
    #dictionnary key   = a value of the varying paramater
    #dicitonnary value = a list of ssexyids with this value
    MergeSets= {}
    
    for SSEXYid,param in ssexy.params.items():
        if MergeSets.has_key(param[VarP]):
            MergeSets[param[VarP]].append(SSEXYid)  
        else:
            MergeSets[param[VarP]] = [SSEXYid] 
    
    #sort each subset
    for sets in MergeSets.iterkeys():
        MergeSets[sets].sort()     

    return MergeSets  

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = OptionParser() 
    parser.add_option("-T", "--temperature", dest="T", type="float",
                      help="simulation temperature in Kelvin") 
    parser.add_option("-b", "--beta", dest="B", type="float",
                      help="simulation temperature in inverse Kelvin") 
    parser.add_option("-x",  dest="x", type="float",
                      help="x length") 
    parser.add_option("-y",  dest="y", type="float",
                      help="y length") 
    parser.add_option("-r",  dest="r", type="float",
                      help="replicas number") 
    parser.add_option("-v", "--varp", dest="varp",
                      choices=['T','B','x','y','r','p','a'], 
                      help="varying parameter, one of [T,B,x,y,r,p,a]") 
    parser.add_option("--restarted", action="store_true", dest="restarted",
                      help="are we merging an ancestry chain of ssexys that got restarted from a each other?")
    parser.add_option("-s", "--skip", dest="skip", type="int",
                      help="how many input lines should we skip?")
    parser.add_option("--average", action="store_true",
                      help="Do we want to average over each file's data?")
    parser.add_option("--waverage", action="store_true",
                      help="Do we want to perform a weighted average?")
    parser.set_defaults(skip=0)

    parser.set_defaults(restarted=False)
    parser.set_defaults(average  =False)
    parser.set_defaults(waverage  =False)

    # parse the command line options and get the reduce flag
    (options, args) = parser.parse_args() 
    if len(args) > 0: 
        parser.error("incorrect number of arguments")

    # We check if we have a MERGED directory, if not create it
    if os.path.exists('MERGED') == False: 
       os.makedirs('MERGED')
    

    # Create the SSEXY analysis helper
    ssexy = ssexyhelp.SSEXYHelp(options)
 
    # Fill up the simulation parameters maps
    ssexy.getSimulationParameters()

    #if there is not need to merge with a varying parameter
    if (not options.varp):
        #Create new ssexyID
        newID = getNewSSEXYID(ssexy)
    
        # Merge all the output files
        print('Merged data files:')
        for type in ssexy.dataType:
            mergeData(ssexy,type,newID,options.skip,options.restarted,ssexy.id,average=options.average,waverage=options.waverage)
    

     #with a varying parameter, one needs to group corresponding ssexyIds
    else:
         #group ssexyIds with the same varying parameter  
         MergeSets = getMergeSets(ssexy,options.varp)
         for varp in sorted(MergeSets.iterkeys()):
             mergeSet = MergeSets[varp]
             print('\nMerged data files for %s=%s:\n' %(options.varp,varp))
             print('SSEXYids to merge: %s' %mergeSet)


             #if there is only one ssexyId with a varp, the just copy the files
             if (len(mergeSet) == 1):
                lsCommand = 'ls *log*%s*' %mergeSet[0]
                LogName = os.popen(lsCommand).read().split('\n')[0]
                shutil.copyfile(LogName,'MERGED/'+LogName) 
                
                for type in ssexy.dataType: 
                    lsCommand = "ls *%s*%s*" %(type,mergeSet[0]) 
                    fileName = os.popen(lsCommand).read().split('\n')[0]
                    outFile,numLines = CreateFileSkipping(fileName,mergeSet[0],mergeSet[0],options.skip,average = options.average,waverage=options.waverage)
                    print('%10d' %numLines)    
                    outFile.close
    
                lsCommand = "ls CYLINDER/*%s*" %mergeSet[0]
                fileNames = os.popen(lsCommand).read().split('\n')
                fileNames.pop()
                for files in fileNames:
                    outFile,numLines = CreateFileSkipping(files,mergeSet[0],mergeSet[0],options.skip)
                    print('%10d' %numLines)    
                    outFile.close
             #otherwise we need to be careful what files do we merge together
             else:
                 #Create new ssexyID
                 newID = getNewSSEXYID(ssexy)      
                 
                 lsCommand = 'ls *state*%s*' %mergeSet[0]
                 #oldLogName = os.popen(lsCommand).read().split('\n')[0]
                 #newLogName = oldLogName.replace(str(mergeSet[0]),str(newID))
                 #shutil.copyfile(oldLogName,'MERGED/'+newLogName)
                 
                 for type in ['estimator']: #ssexy.dataType:
                     mergeData(ssexy,type,newID,options.skip,options.restarted,mergeSet,average=options.average,waverage=options.waverage)
                 #os.system('cp %s %s' % (oldLogName,'MERGED/'+newLogName))    




# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

