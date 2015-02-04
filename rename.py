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

            print fileName
            # Now break up the data name into descriptive parts
            dataFile  = open(fileName,'r');
            dataLines = dataFile.readlines()[:500];
            dataFile.close()
            dataFile  = open(fileName,'w');
            dataFile.writelines(dataLines)
                        

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
	main()
