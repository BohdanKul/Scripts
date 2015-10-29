from numpy import *


colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", 
          "#E17597", "#C1B546",(0,1,0),(1,0,0),(1,1,0),(1,0,1),
          (0,0,0),(0,0,1), (1,0.5,0.5),(0.5,0.5,0.5),
          (0.5,0,0),(1,0.5,0)]

# Generate filename based on system parameters
def GetFileName(estimator,OBC,x):
    filename = ''
    filename += estimator+'-'
    filename += 'OBC'*OBC+'PBC'*(not OBC)+'-'
    filename += '%02d.dat' %x
    #filename += '%1.2f' %H

    return filename

# Retrieve system parameters from a filename
def GetParams(filename):
    ParamsDic = {}
    filename = filename.strip('.dat').split('-')
    ParamsDic['estimator'] = filename[0]
    ParamsDic['OBC']       = filename[1]=='OBC'
    ParamsDic['PBC']       = filename[1]=='PBC'
    ParamsDic['x']         = int(filename[2])
    ParamsDic['y']         = int(filename[3])
    #ParamsDic['H']         = float(filename[4])

    return ParamsDic

# Get scalar data from a file
def GetFileData(filename):
    fdata = loadtxt(filename)
    
    return fdata

# Retrieve data for many clusters geometries 
def GetPBCData(filenames):
    Data = {}
    dim = 1
    for filename in filenames:
        fParams = GetParams(filename)
        Data[(fParams['x'],fParams['y'])] = GetFileData(filename)
        if fParams['y']==2: dim = 2

    return Data,dim


