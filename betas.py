from numpy import *
import sys,os,shutil

#lowL1 = 0
#lowL2 = 2.52
#lowStep = 0.004
#lb  = arange(lowL1,lowL2,lowStep)
#lb2 = lb*2
#
#highL1 = 2.52
#highL2 = 8.02
#highStep = 0.01
#hb  = arange(highL1,highL2,highStep)
#hb2 = hb*2

#usedb = [str('%0.3f ') %e for e in sorted(concatenate((lb2,hb2)))]
#oldb  = [str('%0.3f ') %e for e in sorted(concatenate((lb,hb)))]

lowL1 = 0
lowL2 = 8
lowStep = 0.04
lb  = arange(lowL1,lowL2,lowStep)

highL1 = 2.53
highL2 = 8
highStep = 0.01
hb1  = arange(highL1,highL2,highStep)

highL1 = 0
highL2 = 2.52
highStep = 0.004
hb2  = arange(highL1,highL2,highStep)

fileList = os.listdir('.')
usedb = []
for files in fileList:
    if  files.find('estimator-01')==0:
        usedb.append(str('%0.3f') %float(files.split('-')[-3][1:])) 

#print ' '.join(e for e in sorted(usedb))


#usedb = [str('%0.3f ') %e for e in sorted(lb)]
oldb  = [str('%0.3f') %e for e in sorted(concatenate((hb1,hb2)))]

#print ' '.join(e for e in sorted(oldb))
print type(oldb)
newb = []
for b in oldb:
    if  not(b in usedb):
        newb.append(b)


print ' '.join(e for e in newb)

