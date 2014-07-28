from numpy import *

lowL1 = 0
lowL2 = 2.52
lowStep = 0.004
lb  = arange(lowL1,lowL2,lowStep)
lb2 = lb*2

highL1 = 2.52
highL2 = 8.02
highStep = 0.01
hb  = arange(highL1,highL2,highStep)
hb2 = hb*2


usedb = [str('%0.3f ') %e for e in sorted(concatenate((lb2,hb2)))]
oldb  = [str('%0.3f ') %e for e in sorted(concatenate((lb,hb)))]

newb = []
for b in oldb:
    if  not(b in usedb):
        newb.append(b)


print ' '.join(e for e in newb)

