#!/usr/bin/python
# rename.py
# Adrian Del Maestro
# 09.28.2009
# 
# Rename a series of PIMC output files based on an input file

import os,sys,glob
import argparse

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

	# setup the command line parser options 

        parser = argparse.ArgumentParser(description='Detect id collisions')
        parser.add_argument('fileNames', help='Estimator files', nargs='+')
        parser.add_argument('-r',dest='rename',action='store_true',help='Set to rename problematic files')
        parser.set_defaults(feature=False)
        # parse the command line options and get the file name
        args = parser.parse_args()
	
        fileType = ['state','estimator']
        dic = {}
        for fileName in args.fileNames:

            fileParts = fileName.split('estimator')

            # Now break up the data name into descriptive parts
            dataName = fileParts[-1]
            dataName = dataName.rstrip('.dat')
            dataParts = dataName.split('-')
            oldID = int(dataParts[-1])
            if  not(oldID in dic):
                dic[oldID] = fileName
            else:
                print 'Collision:\n', dic[oldID], '\n', fileName

                if args.rename:
                    newID = oldID + 1
                    while len(glob.glob('*estimator*-%09d*' % newID)) > 0:
                            newID += 1

                    dataName = ''
                    for lab in dataParts[:-1]:
                            dataName += lab + '-'
                    dataName += '%09d.dat' % newID
                    for type in fileType:
                            oldName = fileParts[0] + type + fileParts[-1]
                            newName = fileParts[0] + type + dataName
                            os.popen('mv %s %s' % (oldName,newName))
                    print 'Rename to: ', newName
                             

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
	main()
