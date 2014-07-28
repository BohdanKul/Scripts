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
    abad = [0] + range(32,1600,32)
    for line in ffile.readlines():
        #if  line.find('jobid')!=-1:
        #    jobid = line[line.find('jobid:')+6:].lstrip()
        #if  line.find('command')!=-1:
        #    if  line.find('-x 40 -y 40'):
        #        a = line[line.find('-a')+2:line.find('-t')-1]
        #        a = int(a)
        #        if  a in abad:
        #            print jobid.rstrip('\n'),' ',
        if (line.find('serial')!=-1) or (line.find('NRAP')!=-1) or (line.find('bkulchyt')!=-1):
           if  line.find('-x 40 -y 40')!=-1:
               a = line[line.find('-a')+2:line.find('-t')-1]
               a = int(a)
               if  a in abad:
                   print line[0:8],' ',
            
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()



