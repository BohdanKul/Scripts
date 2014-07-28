import matplotlib as mpl
import matplotlib.pyplot as plt
from helper import *
from numpy  import *
from sys    import exit
import pylab as pl
import argparse
import mplrc
import loadgmt,kevent

#----------------------------------------------------------
#------------------Cluster expansion 2-d-------------------
#----------------------------------------------------------
def ClusterExpansion2d(Pn):
    # Subclusters multiplicity matrix
    Multiplicity2d = { 1: array([ 0, 0, 0, 0, 0, 0,0,0,0,0]),
                       2: array([ 2, 0, 0, 0, 0, 0,0,0,0,0]),
                       3: array([ 3, 2, 0, 0, 0, 0,0,0,0,0]),
                       4: array([ 4, 3, 2, 0, 0, 0,0,0,0,0]),
                       5: array([ 4, 4, 0, 0, 0, 0,0,0,0,0]),
                       6: array([ 6, 7, 2, 0, 2, 0,0,0,0,0]),
                       7: array([ 8,10, 4, 2, 3, 2,0,0,0,0]),
                       8: array([ 9,12, 6, 0, 4, 4,0,0,0,0]),
                       9: array([12,17,10, 3, 6, 7,2,2,0,0]),
                      10: array([16,24,16, 8, 9,12,6,4,4,0])
                   }
    L = array([1,2,2,2,1,2,2,1,2,1])

    # Weights
    Wn = zeros(len(Pn))
    for n in range(1,len(Pn)+1):
        Wn[n-1] = Pn[n-1] - sum(Wn*Multiplicity2d[n][:len(Pn)])
    
    # Observables
    for n in range(1,len(Pn)+1):
        Pn[n-1] = sum(Wn[:n-1]*L[:n-1])

    return Wn,Pn
#----------------------------------------------------------
#----------------------------------------------------------


#----------------------------------------------------------
#------------------Cluster expansion 1-d-------------------
#----------------------------------------------------------
def ClusterExpansion1d(Pn):
    #Subclusters multiplicity matrix
    Multiplicity1d = { 1: array([ 0, 0,0, 0,0,0,0,0,0,0,0]),
                       2: array([ 2, 0,0, 0,0,0,0,0,0,0,0]),
                       3: array([ 3, 2,0, 0,0,0,0,0,0,0,0]),
                       4: array([ 4, 3,2, 0,0,0,0,0,0,0,0]),
                       5: array([ 5, 4,3, 2,0,0,0,0,0,0,0]),
                       6: array([ 6, 5,4, 3,2,0,0,0,0,0,0]),
                       7: array([ 7, 6,5, 4,3,2,0,0,0,0,0]),
                       8: array([ 8, 7,6, 5,4,3,2,0,0,0,0]),
                       9: array([ 9, 8,7, 6,5,4,3,2,0,0,0]),
                      10: array([10, 9,8, 7,6,5,4,3,2,0,0]),
                      11: array([11,10,9, 8,7,6,5,4,3,2,0]),
                      12: array([12,11,10,9,8,7,6,5,4,3,2])
                    }
# Weights
    Wn = zeros(len(Pn))
    Wn[1-1] = Pn[1-1]
    Wn[2-1] = Pn[2-1] - 2*Pn[1-1]
    for n in range(3,len(Pn)+1):
        Wn[n-1] = Pn[n-1] - 2*Pn[(n-1)-1] + Pn[(n-2)-1]

    #Observables
    for n in range(1,len(Pn)+1):
        Pn[n-1] = sum(Wn[0:n-1]) #+ sum(Wn*Multiplicity1d[n][:len(Pn)])

    return Wn, Pn
#----------------------------------------------------------
#----------------------------------------------------------



#----------------------------------------------------------
#---------------------------Main---------------------------
#----------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='Cluster expansion for TFIM')
    parser.add_argument('fileNames', help='Scalar estimator files', nargs='+')
    args = parser.parse_args() 

    #----------------------------------------------------------
    #---------------------------Plotting-----------------------
    #----------------------------------------------------------

    nhs = 41
    min, max = (0, 4)
    
    # Setting up a colormap that's a simple transtion
    mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])

    # Using contourf to provide my colorbar info, then clearing the figure
    Z = [[0,0],[0,0]]
    levels = linspace(min,max,nhs)
    CS3 = plt.contourf(Z, levels, cmap=mymap)
    plt.clf()


    
    pl.figure(1,(6,6))
    pl.connect('key_press_event',kevent.press)
    pl.rcParams.update(mplrc.aps['params'])
    axEW = pl.subplot(121)
    axEW.set_ylabel(r'$\mathrm{W_E}$')
    axEW.set_xlabel(r'$\mathrm{Cluster \, size}$')

    axMW = pl.subplot(122)
    axMW.set_ylabel(r'$\mathrm{W_M}$')
    axMW.set_xlabel(r'$\mathrm{Cluster \, size}$')
    Clb = plt.colorbar(CS3) # using the colorbar info I got from contourf
    Clb.ax.set_ylabel(r'$\mathrm{h/J}$')
    
    #axE = pl.subplot(121)
    #axE.set_xlabel(r'$\mathrm{h/J}$')
    #axE.set_ylabel(r'$\mathrm{E/N}$')
    #
    #axM = pl.subplot(122)
    #axM.set_xlabel(r'$\mathrm{h/J}$')
    #axM.set_ylabel(r'$\mathrm{M/N}$')

    #----------------------------------------------------------
    #----------------------------------------------------------


    N2Cluster2d = array([(1,1),(2,1),(3,1),(4,1),
                         (2,2),(3,2),(4,2),
                         (3,3),(4,3),
                         (4,4)
                       ])

    # Get all subclusters data from files
    ClustersData, dim = GetPBCData(args.fileNames)

    # Number of subclusters
    size = len(args.fileNames)

    step = float(max-min)/len(ClustersData[(1,1)][::,0])
    for i,H in enumerate(ClustersData[(1,1)][::,0]):
        En  = zeros(size) 
        Mn  = zeros(size)
        En2 = zeros(size) 
        Mn2 = zeros(size)
        print H
        # Put E and M data in the right order
        for (x,y) in ClustersData.keys():
            if dim == 2: index = N2Cluster2d.tolist().index([x,y])
            else:        index = x-1

            #print (x,y),index
            En[index]  = ClustersData[(x,y)][i,1]*x*y
            Mn[index]  = ClustersData[(x,y)][i,2]*x*y
            En2[index] = ClustersData[(x,y)][i,3]*x*y
            Mn2[index] = ClustersData[(x,y)][i,4]*x*y
        # Perform cluster expansion 
        if  dim == 2: 
            M_Wn, M_Pn   = ClusterExpansion2d(Mn)
            E_Wn, E_Pn   = ClusterExpansion2d(En)
            M2_Wn, M2_Pn = ClusterExpansion2d(Mn2)
            E2_Wn, E2_Pn = ClusterExpansion2d(En2)
        else:
            M_Wn, M_Pn = ClusterExpansion1d(Mn)
            E_Wn, E_Pn = ClusterExpansion1d(En)

        r = (float(step*i)-min)/(max-min)
        g = 0
        b = 1-r
        axEW.plot(range(1,len(E_Wn)+1),E_Wn, ls = '-', marker='s', color=(r,g,b))
        axMW.plot(range(1,len(M_Wn)+1),M_Wn, ls = '-', marker='s', color=(r,g,b),label=r"$\mathrm{H = %0.2f}$" %H)
        #axM.plot(-H,M_Pn[-1], ls = '--', marker='s', color=colors[0])
        #axE.plot(-H,E_Pn[-1], ls = '-', marker='s', color=colors[0])
        #axM.plot(-H,M2_Pn[-1], ls = '-', marker='s', color=colors[1])
        #axE.plot(-H,E2_Pn[-1], ls = '-', marker='s', color=colors[1])
        if H == -1: 
            print M_Wn
            #axM.plot(-H,M_Pn[-1], ls = '--', marker='s', color=colors[0],label=r'$\mathrm{M_0}$')
            #axE.plot(-H,E_Pn[-1], ls = '-', marker='s', color=colors[0],label=r'$\mathrm{E_0}$')
            #axM.plot(-H,M2_Pn[-1], ls = '-', marker='s', color=colors[1],label=r'$\mathrm{M_1}$')
            #axE.plot(-H,E2_Pn[-1], ls = '-', marker='s', color=colors[1],label=r'$\mathrm{E_1}$')
    #lgd = axE.legend(loc='best')
    #lgd.draggable(state=True)
    #lgd = pl.legend(loc='best')
    #lgd.draggable(state=True)
    pl.tight_layout()
    pl.show()
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()


