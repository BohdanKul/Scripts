import argparse 
import ssexyhelp
import numpy as np
# ----------------------------------------------------------------------
def getInfoSpecifiers(fileparams):
    allspecs = []
    for elem in fileparams:
        infospecs = {}
        for key, value in elem.iteritems():
            if value!='': 
                infospecs[key] = value
        allspecs += [infospecs]

    return allspecs 

# ----------------------------------------------------------------------
def getAllSpecifiers(fileparams):
    allkeys = []
    for elem in fileparams:
        for key, value in elem.iteritems():
            if not(key in allkeys):
                allkeys += [key]

    return allkeys

# ----------------------------------------------------------------------
def augmentMissingSpecifiers(fileparams):
    allSpec = getAllSpecifiers(fileparams)
    for elem in fileparams:
        for spec in allSpec:
            if not(spec in elem.keys()):
                elem[spec] = 0.0

    return fileparams

# ----------------------------------------------------------------------
def cleanData(view):
    if np.isnan(view[-1, 0]) or (view[-2, 0] < view[-1, 0]):
        view = view[:-1,:]

    return view

# ----------------------------------------------------------------------
def computeEnergies(header, data):
    split1 = header.index('<Z0>')
    split2 = header.index('<X0>')      - split1
    split3 = header.index('<ZZ(0,1)>') - split1
    
    lEs = data[:, 1:split1]*data[:, split1:]
    
    QE = np.sum(lEs[:, :split2], axis=1) + np.sum(lEs[:, split3:], axis=1)
    CE = np.sum(lEs[:,  split2 : split3], axis=1) 
    QE = QE[np.newaxis].T
    CE = CE[np.newaxis].T

    return QE, CE

def getHeaders(filename):
    inFile = open(filename,'r');
    inLine = inFile.readlines();
    inLine[0] = inLine[0].rstrip()
    oheaders = inLine[0].split('   ')
    oheaders.pop(0)
    nheaders = []
    for head in oheaders:
        if head!='': nheaders += [head.replace(' ','')]
    inFile.close()

    return nheaders



# ----------------------------------------------------------------------
def main(): 

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('filenames', help='Training output files', nargs='+')

    args = parser.parse_args()
    filenames = args.filenames
    specs = []
    for i,filen in enumerate(filenames):
        fparams = ssexyhelp.getReduceParamMap(filen)
        fparams['filename'] = filen
        fparams['index']    = i
        specs += [fparams]

    specs = getInfoSpecifiers(specs)
    specs = augmentMissingSpecifiers(specs)
    
    for spec in specs:
        filename = spec.pop('filename')
        
        # Read off the list of headers
        # Load and check the data
        fdata = np.loadtxt(filename) 
        fdata = cleanData(fdata.view())
       
        headers = getHeaders(filename)

        # Augment the data set with new quantities
        QE, CE = computeEnergies(headers, fdata)
        nheaders += ['E_Q', 'E_C']
        fdata = np.hstack([fdata, QE, CE])
        spec['data'] = fdata.view(dtype= zip(nheaders, [str(fdata.dtype)]*len(nheaders))).copy()

# ----------------------------------------------------------------------
def getCategory(L, key, value):
    cat = []
    for dic in L:
        if dic[key] == value:
            cat += [dic]

    return cat

# ----------------------------------------------------------------------
def convert2Dic(L):
    D = {}
    for key in L[0].keys():
        if key!= 'data':
            D[key] = []
        else:
            for header in L[0][key].dtype.names:
                D[header] = []
    for l in L:
        for key, value in l.iteritems():
            if key!= 'data':
                D[key] += [value]
        else:
            for header in value.dtype.names:
                D[header] = [value]

    return D

# ----------------------------------------------------------------------
def getData(L, *args):
    data = []
    for 

# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()


