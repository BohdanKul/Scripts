import sys,os,shutil
from optparse import OptionParser

def GetPIMCs(fileName):
    '''Creates a list of PIMC ids from a job output file generated on a cluster'''
    inLines = []
    try:
        inFile = open(fileName,'r')
        inLines = inFile.readlines();
    except:
          return 0
    PIMCs = []
    Mstate = False
    Mstate = True
    for line in inLines:
        #if  line.find('Legs') != -1:
        #    Mstate = True
        #    continue
        if  Mstate:
            if  line.find('Measurement taken') != -1:
                pimcid = line[:9]
                PIMCs.append(pimcid)
                break
    inFile.close
    return PIMCs

def GetUnique(PIMCs):
    '''Returns unique entries in PIMCs list'''
    res = [PIMCs[0]] 
    for PIMC in PIMCs[1:]:
        newflag = bool(1)
        for pimcid in res:
            if pimcid == PIMC:
               newflag = bool(0)
               break
        if newflag == True:  
           res.append(PIMC)
    return res 

def DetectDiff(inList, lookFor):
    '''Detects which elements of lookFor are not contained in inList'''
          
    absElements = []
    for n,pimcid in enumerate(lookFor):
        existsflag = False
        for elements in inList:
            if elements.find(pimcid) != -1:
               existsflag = True
               break
        if existsflag == False:
           absElements.append(n)    
    return absElements

def CopyDataFiles(BaseDir, newfolderName,PIMCs,fileList):
   
    for pimcid in PIMCs: 
        for fileName in fileList:
            if fileName.find(pimcid) != -1:    
               shutil.copyfile(BaseDir+fileName, newfolderName+fileName)

def MoveDataFiles(BaseDir, newfolderName,PIMCs,fileList):

    for pimcid in PIMCs: 
        for fileName in fileList:
            if fileName.find(pimcid) != -1:    
               shutil.move(BaseDir+fileName, newfolderName+fileName)

def DeleteDataFiles(BaseDir,PIMCs,fileList):

    for pimcid in PIMCs: 
        for fileName in fileList:
            if fileName.find(pimcid) != -1:    
               os.remove(BaseDir+fileName)

def main():

    parser = OptionParser() 
    parser.add_option("-r", "--run", help="run number") 
    parser.add_option("-m", action="store_true", dest="move",    default=False,\
                      help="Check to move pimc data files to the 'run-${run}' folder")
    parser.add_option("-d", action="store_true", dest="delete",  default=False,\
                      help="Check to delete pimc data files")
    parser.add_option("-c", action="store_true", dest="copy",    default=False,\
                      help="Check to copy pimc data files to the 'run-${run}' folder")
    parser.add_option("-a", action="store_true", dest="archive", default=False,\
                      help="Check to archive pimc data files to the archive folder")
    parser.add_option("--out",    dest = 'outDir',\
                      help="Folder containing pimc redirected output files. Default = ''", default = './')
    parser.add_option("--OUTPUT", dest = 'OUTPUTDir',\
                      help="Folder containing pimc data files. Default = '../OUTPUT/'", default = '../OUTPUT/')     
    
    (options, args) = parser.parse_args()

    options.outDir = options.outDir.rstrip('/')+'/'
    options.OUTPUTDir = options.OUTPUTDir.rstrip('/')+'/'

    if options.run == None:
       parser.error("Must specify flag -r to determine relevant output files to proceed")

    if options.delete and options.copy:
       parser.error("Cant delete and copy at the same time")
    if options.move and options.copy:
       parser.error("Cand move and copy at the same time")
#---------------------------------------------------------------------------------------------#
#---------------Get unique pimcids from the redirected output files---------------------------#
#---------------------------------------------------------------------------------------------#

    fileList = os.listdir(options.outDir)
    PIMCs = []
    outputFiles = [] 

    for files in fileList:
        if files.find(options.run) != -1:
           outputFiles.append(files)
           temp = GetPIMCs(options.outDir+files)
           if temp == 0:
              print '%s contains no measurements' %options.outDir+files
           else:
               PIMCs.extend(temp)
    print PIMCs
    if len(outputFiles) == 0:
       print 'No redirected output files for the run %s detected' %options.run
       sys.exit() 
    
    if len(PIMCs) == 0:
       print 'No measurements for the run %s detected' %options.run
       sys.exit() 

#Store unique pimcids
    PIMCs = GetUnique(PIMCs)
#Keep an extra copy for the Cylinder folder 
    cPIMCs = PIMCs[:] 

#---------------------------------------------------------------------------------------------#
#-----------------------------------------Normal output---------------------------------------#
#---------------------------------------------------------------------------------------------#

    fileList = os.listdir(options.OUTPUTDir)
    delpimcid = DetectDiff(fileList, PIMCs)
    for n in reversed(delpimcid):
        print('%s doesnt longer exist') %PIMCs[n]
        PIMCs.pop(n)
    if len(PIMCs) == 0:
       print 'None of the files cant be found'
       sys.exit()
#-------------------------- archive data-files-----------------------
    if (options.archive):
        #Create an archive folder
        if not(os.path.exists(options.OUTPUTDir+'archive/')):
           os.makedirs(options.OUTPUTDir+'archive/')
           for pimcid in PIMCs: 
               for fileName in fileList:
                   if fileName.find(pimcid) != -1:    
                      shutil.copyfile(options.OUTPUTDir+fileName, options.OUTPUTDir+'archive/'+fileName)

#----------------------------------------------------------------------
    #Create an output folder 
    if options.copy or options.move:
       folderName = options.OUTPUTDir+'run-'+options.run+'/'
       if os.path.exists(folderName): 
           print("Folder %s already exists" %(folderName))
           sys.exit()
       else:
           os.makedirs(folderName)
           for outFile in outputFiles:
               shutil.copyfile(options.outDir+outFile,folderName+outFile)
#-------------------------- Copy data-files----------------------------
    if (options.copy):
       CopyDataFiles(options.OUTPUTDir,folderName,PIMCs,fileList)

#-------------------------- Move data-files----------------------------
    if (options.move):
       MoveDataFiles(options.OUTPUTDir,folderName,PIMCs,fileList) 

#-------------------------- Delete data-files--------------------------
    if (options.delete):
       DeleteDataFiles(options.OUTPUTDir,PIMCs,fileList) 

#---------------------------------------------------------------------------------------------#
#--------------------------------------Cylindrical output-------------------------------------#
#---------------------------------------------------------------------------------------------#
    if os.path.exists(options.OUTPUTDir+'CYLINDER'): 
       cfileList = os.listdir(options.OUTPUTDir+'CYLINDER') 
    else:
         print 'No cylindrical ouput'
         sys.exit()


    for n,fileName in enumerate(cfileList):
        cfileList[n] = 'CYLINDER/' + fileName
  
#-------------------------- archive data-files-----------------------
    if (options.archive):
        #Create an archive folder
        if not(os.path.exists(options.OUTPUTDir+'archive/CYLINDER')):
           os.makedirs(options.OUTPUTDir+'archive/CYLINDER')
           for pimcid in PIMCs: 
               for fileName in cfileList:
                   if fileName.find(pimcid) != -1:    
                      shutil.copyfile(options.OUTPUTDir+fileName, options.OUTPUTDir+'archive/'+fileName) 

    if options.copy or options.move: 
       os.makedirs(folderName+'CYLINDER/')

#-------------------------- Copy data-files----------------------------
    if (options.copy):
       CopyDataFiles(options.OUTPUTDir,folderName,PIMCs,cfileList)

#-------------------------- Move data-files----------------------------
    if (options.move):
       MoveDataFiles(options.OUTPUTDir, folderName,PIMCs,cfileList) 

#-------------------------- Delete data-files--------------------------
    if (options.delete):
       DeleteDataFiles(options.OUTPUTDir, PIMCs,cfileList) 


if __name__ == "__main__": 
    main()

