import os,sys
# -----------------------------------------------------------------------------
def graham(commandLines, run, memory):
    ''' Write a pbs submit script for CLUMEQ. '''

   # Open the pbs file and write its header
    numOptions = len(commandLines)
    fileName = 'run_%s.sh' % run

    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#!/bin/bash
#SBATCH --time 168:00:00
#SBATCH --output=run_%s.out
#SBATCH --mem-per-cpu=%s
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1\n''' %(memory, run))
    pbsFile.write('''
case $SLURM_ARRAY_TASK_ID in\n''')
    # Create the command string and make the case structure
    for n in range(numOptions):
        pbsFile.write('%d)\n' % (n))
        pbsFile.write('%s &\n' %  (commandLines[n]))
        pbsFile.write(';;\n')
    pbsFile.write('esac\nwait\n')


    pbsFile.close();

    print '\nSubmit jobs with: sbatch --array 1-%d %s\n' % (numOptions, fileName)


# -----------------------------------------------------------------------------
def clumeq(commandLines,run):
    ''' Write a pbs submit script for CLUMEQ. '''

   # Open the pbs file and write its header
    numOptions = len(commandLines)
    fileName = 'submit-pimc_%s_clumeq.pbs' % run
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash
#PBS -A zwu-174-aa
#PBS -l walltime=48:00:00
#PBS -V MOAB_JOBARRAYINDEX
#PBS -l nodes=1:ppn=8\n''')
    pbsFile.write('#PBS -N run-%s\n#PBS -e out/run-%s-%%J-%%I\n#PBS -o out/run-%s-%%J-%%I\n' % (run,run,run))
    if (numOptions%8==0):
        numNodes = numOptions//8
    else:
        numNodes = numOptions//8+1
    pbsFile.write('#PBS -t [1-%d]\n ' % (numNodes))
    pbsFile.write('''
cd ${PBS_O_WORKDIR}
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bkulchytskyy/gcc/lib:/opt/torque/lib:/opt/moab/lib:
case $MOAB_JOBARRAYINDEX in\n''')
    # Create the command string and make the case structure
    if (numOptions%8 == 0):
        for n in range(1,numNodes+1):
            pbsFile.write('%d)\n' % (n))
            pbsFile.write('sleep %d\n' % (20*(n-1)))
            for i in range(0,8):
                pbsFile.write('%s &\nsleep 2\n' %  (commandLines[(n-1)*8+i]))
            pbsFile.write(';;\n')
    else:
        for n in range(1,numNodes):
            pbsFile.write('%d)\n' % (n))
            pbsFile.write('sleep %d\n' % (20*(n-1)))
            for i in range(0,8):
                pbsFile.write('%s &\nsleep 2\n' % (commandLines[(n-1)*8+i]))
            pbsFile.write(';;\n')
        pbsFile.write('%d)\n' % (numNodes))
        pbsFile.write('sleep %d\n' % (20*(numNodes-1)))
        for i in range(0,numOptions%8):
            pbsFile.write('%s &\nsleep 2\n' % (commandLines[(numNodes-1)*8+i]))
        pbsFile.write(';;\n')
    pbsFile.write('esac\nwait\n')


    pbsFile.close();

    print '\nSubmit jobs with: msub %s\n' % (fileName)

# -----------------------------------------------------------------------------
def westgrid(commandLines,run):
    ''' Write a pbs submit script for westgrid. '''

    numOptions = len(commandLines)
    # Open the pbs file and write its header
    fileName = 'submit-pimc_%s_westgrid.pbs' % run
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash\n
#PBS -l walltime=120:00:00
#PBS -N PIMC
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_JOBID}\n
# Do not send email
#PBS -M agdelma@phas.ubc.ca
#PBS -m n\n
# Start job script
cd $PBS_O_WORKDIR
echo \"Starting run at: `date`\"

case ${PBS_ARRAYID} in\n''')

    # Create the command string and make the case structure
    for n in range(numOptions):
        pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,commandLines[n]))

    pbsFile.write('esac\nwait\n')
    pbsFile.close();

    print '\nSubmit jobs with: qsub -t 0-%d %s\n' % (numOptions-1,fileName)

# -----------------------------------------------------------------------------
def sharcnet(commandLines,run):
    ''' Write a submit script for sharcnet. '''

    numOptions = len(commandLines)
    # Open the script file and write its header
    fileName = 'submit-pimc_%s_sharcnet' % run
    scriptFile = open(fileName,'w')
    scriptFile.write('''#!/bin/bash
# Sharcnet pimc submit script\n\n''')

    # Create the command string and output to submit file
    for n in range(numOptions):
        name = '''out/run-%s-%%J''' % run
        #scriptFile.write('sleep %d\nsqsub -q NRAP_893 -o %s --mpp=2G -r 1h %s\n' % (2,name,commandLines[n]))
        scriptFile.write('sqsub -q NRAP_893 -o %s --mpp=1G -r 7d %s\n' % (name,commandLines[n]))

    scriptFile.close();
    os.system('chmod u+x %s'%fileName)
    print '\nSubmit jobs with: ./%s\n' % (fileName)

# -----------------------------------------------------------------------------
def scinet(commandLines,run):
    ''' Write a pbs submit script for scinet. '''

    numOptions = len(commandLines)
    if numOptions != 8:
        print 'For scinet, must submit in multiples of 8'
        return

    # Open the pbs file and write its header
    fileName = 'submit-pimc_%s_scinet' % run
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash
# MOAB/Torque submission script for multiple serial jobs on SciNet GPC
#
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N serialx8_pimc
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_JOBID}
# Do not send email
#PBS -M agdelma@phas.ubc.ca
#PBS -m n

# Start job script
cd $PBS_O_WORKDIR
echo \"Starting run at: `date`\"\n\n''')

    # Create the command string and output to submit file
    for n in range(numOptions):
        pbsFile.write('(sleep %02d; %s) &\n' % (10,commandLines[n]))
    pbsFile.write('wait')
    pbsFile.close();
    print '\nSubmit job with: qsub %s\n'%fileName

# -----------------------------------------------------------------------------
def bluemoon(commandLines,run):
    ''' Write a pbs submit script for bluemoon '''

    numOptions = len(commandLines)
    # Open the pbs file and write its header
    fileName = 'submit-pimc%s_bluemoon.pbs' % run
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash\n
#PBS -l nodes=1:ppn=1
#PBS -l walltime=30:00:00
#PBS -N PIMC
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_JOBID}\n
#PBS -M Max.Graves@uvm.edu
#PBS -m n\n
# Start job script
cd $HOME/PIMC
echo \"Starting run at: `date`\"

case ${PBS_ARRAYID} in\n''')

    # Create the command string and make the case structure
    for n in range(numOptions):
        pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,commandLines[n]))

    pbsFile.write('esac\nwait\necho \"Finished run at: `date`\"')
    pbsFile.close();

    print '\nSubmit jobs with: qsub -t 0-%d %s\n' % (numOptions-1,fileName)

