
import os,sys
import argparse
# ----------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='Plot Raw MC Equilibration Data for Scalar Estimators.')
    parser.add_argument('fileNames', help='Scalar estimator files', nargs='+')
    args = parser.parse_args()

    fileNames = args.fileNames
   
    ffile = open(fileNames[0],'r')
    for line in ffile.readlines():
        if  line.find('command')!=-1:
            print 'sqsub -q NRAP_893 -o out/run-b-184_saw_01-%J --mpp=2G -r 7d ', line[line.find('./ssexy'):],

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()



