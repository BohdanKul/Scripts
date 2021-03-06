''' PimcHelp - helper methods for analzing pimc output data.
'''
import os
from operator import itemgetter, attrgetter
from numpy import *

# -----------------------------------------------------------------------------
def getHeadersFromFile(fileName): 
    ''' Get the data column headers from a PIMC output file. '''

    inFile = open(fileName,'r');
    inLines = inFile.readlines();
    headers = inLines[0].split()
    headers.pop(0)
    inFile.close()
    return headers

## ----------------------------------------------------------------------------
def getParamMap(fname):
    '''Get the parameters from the output filename.  
    
    We need to be careful due to the problems associated with splitting at the
    '-' character when there are potential minus signs.
    '''

    fileParts  = fname.rstrip('.dat').split('-')
    pIndex = []
    for n,part in enumerate(fileParts):
        if part == '':
            fileParts[n+1] = '-' + fileParts[n+1]
            pIndex.append(n)
        elif ((len(part)==1) and (part.isalpha())):
            fileParts[n+1] = part + '-' + fileParts[n+1]
            fileParts.pop(n)
    for n in pIndex:
        fileParts.pop(n)

    
    params = {'r': int(fileParts[1]),
              'x': int(fileParts[2]),       
              'y': int(fileParts[3])}
    if 't' in fileParts[4]: params['T']    = float(fileParts[4][1:])
    else:                   params['b']    = float(fileParts[4][1:])
    k = 5
    if 'd' in fileParts[5]:
        params['d']    = float(fileParts[5][1:])
        k = 6 
    if fileParts[k].isdigit(): params['a'] = int(fileParts[k])
    else: 
        if 'p' in fileParts[k]: params['p'] = int(fileParts[k][1:])
    return params 

#---------------------------------------------------------------------------
def getReduceParamMap(fname):
    '''Get the parameters from the output filename.  
    '''

    fileParts  = os.path.split(fname)[1].rstrip('.dat').split('_')
    fileParts.pop(0)

    paramMap = {}
    for part in fileParts:
        if  '-' in part:
            if  '--' in part: 
                (item,temp,value)=part.split('-')
                value = '-' + value
            else:
                (item,value)=part.split('-')
            #paramMap.update(dict([part.split('-')]))
            if '.' in value: paramMap[item] = float(value)
            else:            paramMap[item] = value
        else:
            paramMap[part] = ''

    return paramMap 

# -------------------------------------------------------------------------------
def getWildCardString(options):
    ''' Using the command line flags, form the input file string that will
        be used to open all data files. '''

    out = ''
    if options.b is not None: flagB = "%06.3f" % options.b; out += '_b-'+flagB
    else:                     flagB = "*"

    if options.T is not None: flagT = "%06.3f" % options.T; out += '_T-'+flagT
    else:                     flagT = "*"

    if options.x is not None: flagx = "%03d" % options.x;   out += '_Lx-'+flagx
    else:                     flagx = "*"

    if options.y is not None: flagy = "%03d" % options.y;   out += '_Ly-'+flagy
    else:                     flagy = "*"
    
    if options.r is not None: flagr = "%02d" % options.r;   out += '_r-'+flagr
    else:                     flagr = "*"
    
    if options.d is not None: flagd = "%06.3f" % options.d;   out += '_delta-'+flagd
    else:                     flagd = "*"

    if    options.T: dataName = '%s-%s-%s-t%s' % (flagr,flagx,flagy,flagT)
    elif  options.b: dataName = '%s-%s-%s-b%s' % (flagr,flagx,flagy,flagB)
    else:            dataName = '%s-%s-%s-%s'  % (flagr,flagx,flagy,"*")

    if    options.d: dataName += '-d%s' % (flagd)

    dataName += '*.dat'
    return dataName, out



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
class ScalarReduce:

     def __init__(self, fileName):

        self.fname    = fileName
        self.paramMap = getReduceParamMap(fileName)
        self.headers  = getHeadersFromFile(fileName)
        self.rvar     = self.headers[0]
        self.loadData()
## -------------------------------------------------------------------------------
     def getTupleIdstr(self,order):
            
         tupleIdstr = ''  
         for item in order:
             if item in self.paramMap.keys():
                if self.paramMap[item]:
                   tupleIdstr += r'%s=%s\,' %(item,self.paramMap[item])
                else:    
                   tupleIdstr += r'%s\,' %(item)
         return tupleIdstr

#-------------------------------------------------------------------------------
     def getTupleId(self,order):
            
         tupleId = ()  
         for item in order:
             if self.paramMap[item]:
                TupleId += self.paramMap[item]
         return TupleId

# -------------------------------------------------------------------------------
     def getHeaders(self):
         return self.headers

# -------------------------------------------------------------------------------
     def getParmap(self):
         return self.paramMap

# -------------------------------------------------------------------------------
     def loadData(self):
         self.data     = loadtxt(self.fname)
         
# -------------------------------------------------------------------------------
     def getData(self):
         return self.data
# -------------------------------------------------------------------------------
     def getrParams(self):
         if len(self.data.shape)>1: return self.rvar, self.data[:,0]
         else:                      return self.rvar, self.data[0]
# -------------------------------------------------------------------------------
     def getAverages(self,estimator):
         if not(estimator in self.headers):
             return [],[]
         else:   
             index = self.headers.index(estimator)
             if len(self.data.shape)>1: return self.data[:,index], self.data[:,index+1]
             else:                       return self.data[index],   self.data[index+1]



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
class SSEXYHelp:
    ''' Helper methods for analzing SSEXY output data. '''

    # ----------------------------------------------------------------------
    def __init__(self,options,baseDir=''):

        self.baseDir  = baseDir
        self.dataName, out = getWildCardString(options) 
        self.dataType = ['estimator']

    # -----------------------------------------------------------------------------
    def getID(self,fileName): 
        ''' Return the ID number corresponding to a given filename. '''
        ID = int(fileName.rstrip('.dat').split('-')[-1])
        return ID
          
    # -----------------------------------------------------------------------------
    def getSimulationParameters(self): 
        '''Get the full list of parameter maps for each input file and a list of
           ID numbers. '''

        # Get the list files
        fileNames = self.getFileList("estimator")

        self.params = {}
        self.id = []
        for fname in fileNames:
            #ID = self.getID(fname)
            ID = fname
            self.id.append(ID)
            self.params[ID] = getParamMap(fname)

    # -----------------------------------------------------------------------------
    def getFileList(self,type,idList=None):
        ''' Get a list of input files based on their type, or possibly a number
            of unique ID's'''

        fileNames = []

        # We want all the file names here
        if not idList:
            lsCommand = 'ls -1 %s%s-%s' % (self.baseDir,type,self.dataName)
            fileNames = os.popen(lsCommand).read().split('\n')
            fileNames.pop()

        # Otherwise we just go through and get the ID's we need
        else:
            for id in idList: 
                #lsCommand = 'ls -1 %s%s-*%s.dat' % (self.baseDir,type,id)
                lsCommand = 'ls -1 %s%s' % (self.baseDir, id)
                fileNames.extend(os.popen(lsCommand).read().split('\n'))
                fileNames.pop() 

        return fileNames


