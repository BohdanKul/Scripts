from numpy import *

def RWstats(coins,Nwalks,Lbound):
    # Random walk
    cprobs = array([0.3,0.5])
    Ubound  = len(coins)
    Nwins   = 0
    ipos    = Lbound
    WinsStats  = zeros(len(coins),dtype=float)
    for i in range(Nwalks):
        pos = ipos
        bWinsStats = zeros(len(coins),dtype=bool)
        while (pos in range(Lbound,Ubound)):
            if not(bWinsStats[pos]): 
                bWinsStats[pos] = True
                WinsStats[pos] += 1.0
            if  random.rand(1)[0]<cprobs[coins[pos]]:
                pos += 1.0
            else:
                pos -= 1.0
        #if (pos == Ubound): Nwins += 1
    
    return WinsStats/float(Nwalks)
    #return float(Nwalks-Nwins)/float(Nwalks)
    #print "# wins:   ", Nwins
    #print "# looses: ", Nwalks-Nwins
    #print "# total:  ", Nwalks
    #print "Probability of ruin: ", float(Nwalks-Nwins)/float(Nwalks)

# Generate a distribution of coins
def GenerateCoins(N,ndiff):
    coins = zeros(N,dtype=int)
    while ndiff>0:
          cindex = random.randint(Lbound,len(coins))
          if coins[cindex]==0:
             if random.rand(1)<0.5:
                coins[cindex] = 1
                ndiff -= 1

    return coins

N     = 5
ndiff = 3
Nwalks  = 1000
Lbound  = 0
Nstats  = 10


print "# coins: ", N, " # diff coins: ", ndiff, "# walks: ", Nwalks," # stats: ", Nstats


print "Determined coins distribution: "
stats  = zeros((Nstats,N),dtype=float)
coins = GenerateCoins(N,ndiff) 
for i in range(Nstats):
    stats[i,:] = RWstats(coins,Nwalks,Lbound)
print "Probability of a win: ", mean(stats,axis=0), " +/- ", std(stats,axis=0)

print "Random coins distribution: "
stats  = zeros((Nstats,N),dtype=float)
for i in range(Nstats):
    coins = GenerateCoins(N,ndiff)
    stats[i,:] = RWstats(coins,Nwalks,Lbound)
print "Probability of a win: ", mean(stats,axis=0), " +/- ", std(stats,axis=0)


