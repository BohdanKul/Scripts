import numpy as np
import pylab as pl
from helper import *
from scipy.sparse.linalg import eigsh
from scipy.sparse import lil_matrix 
import argparse
import mplrc
import loadgmt,kevent

#A simple diagonalization for the 1D TFIM with OBCs

def main():

    parser = argparse.ArgumentParser(description='Exact diagonilization for TSFI')
    parser.add_argument('-x',   type=int, default=8, help = "Lattice x-dimension")
    parser.add_argument('-a',   type=int, default=4, help = "Number of sites in region A")
    parser.add_argument('-z',   type=float, default=1, help = "Jz coupling")
    parser.add_argument('--OBC', action='store_true', default=False, help = "Open boundary conditions?")
    parser.add_argument('--lanczos',action='store_true', default=False, help = "Use Lanczos method")
    args = parser.parse_args() 

#------------------------------------------------------------Plotting-----------------------
#----------------------------------------------------------

    #pl.rcParams.update(mplrc.aps['params'])
    pl.connect('key_press_event',kevent.press)
    
    pl.figure(1,(8,6))
    nplots = 1
    axE = pl.subplot(nplots,1,1)
    axE.set_xlabel(r'$\mathrm{\Delta}$')
    axE.set_ylabel(r'$\mathrm{(E-E_0)/N}$')
    pl.title(r'$L=%d$' %(args.x)) 

    #axG = pl.subplot(nplots,1,2)
    #axG.set_xlabel('x')
    #axG.set_ylabel('G')

#----------------------------------------------------------
#-----------------------Constants--------------------------
#----------------------------------------------------------
    nlevels = 10                   # Number of lowest energy eigenstates
    #Deltas = np.array([3.044])

    #Deltas = np.array([1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0])
    #Deltas = np.array([-0.4])
    N = args.x
    A = args.a
    Jz  = args.z
    #Jzs = np.arange(-1.05,-0.90,0.05)
    Jzs = np.arange(-1.05,1.05,0.05)
    npoints = len(Jzs)
    HDim = 2**N  #Hilbert Space Dimension
    rDim = 2**A  #Density matrix dimension
    tDim = 2**(N-A)
    OBC = args.OBC  

#----------------------------------------------------------
#-----------------------Variables--------------------------
#----------------------------------------------------------
    # Energy of eigen states
    Es  = np.zeros((nlevels,npoints))

    ## Energy utput file`
    #filename = GetFileName('estimator',args.OBC,args.x)    
    #outefile = open(filename,'w')
    #headers  = '#%15s' %('Jz')
    #for i in range(nlevels):
    #    headers += '%16s' %('e'+str(i))
    #headers += '\n'
    #outefile.write(headers)

    ## Ground state state
    #filename = GetFileName('ground',args.OBC,args.x)    
    #outgfile = open(filename,'w')
    #headers  = '#%15s' %('Jz')
    #for i in range(HDim):
    #    headers += '%16s' %(str(i))
    #headers += '\n'
    #outgfile.write(headers)



    for l,Jz in enumerate(Jzs):
        #print 'Jz = ', Jz
    #----------------------------------------------------------
    #--------------------Diagonal term-------------------------
    #----------------------------------------------------------
        #print "Computing diagonal term"

        #Compute diagonal Hamiltonian term first
        Ham_diag = np.zeros(HDim) 
        #Loop through the hilbert space
        for bra in range (0,HDim):
            #Next, loop through all the Hamiltonian elements
            Ham_diag[bra] = 0
            tbra = bra

            # PBC, go through all bra spins 
            # This operation depends on dimensionnality

            # For 1-site, diagonal contribution is zero
            if N > 1:
                for s in range (0,N-OBC):    
                    S0 = 2*((tbra>>s)&1) - 1      
                    S1 = 2*((tbra>>((s+1)%N))&1) - 1  
                    Ham_diag[bra] += Jz*S0*S1/4.0    # Energy of a neighbouring pair


    #----------------------------------------------------------
    #----------------Off-diagonal term-------------------------
    #----------------------------------------------------------
        #print "Computing off-diagonal term"

        #Create Hamiltonian matrix with filled diagonal terms
        #For lanczos method, use a sparse matrix
        if args.lanczos: Ham = lil_matrix((HDim,HDim))
        else:            Ham = np.zeros((HDim,HDim))


        #Set the diagonal
        if args.lanczos: Ham.setdiag(Ham_diag)

        #Loop through the hilbert space
        for bra in range(0,HDim):
        
            #Diagonal matrix elemnt 
            if not(args.lanczos): Ham[bra,bra] = Ham_diag[bra]  
            
            #Now let's do the off-diagonal terms
            tbra = bra
            for s in range (0,N-OBC):
                S0 = 2*((tbra>>s)&1) - 1      
                S1 = 2*((tbra>>((s+1)%N))&1) - 1  
                if  (S0*S1 == -1):
                    sy = (2**s + 2**((s+1)%N))  #Bits to be flipped
                    ket = tbra^sy               #Flip those bits in the bra
                    Ham[bra,ket] = 0.5
        
        #print Ham
        #this returns the eigenvalues w and eigenvectors v
        if  args.lanczos: 
            w, v = eigsh(Ham,k=nlevels,which='SA',maxiter=25000)  
        else:            w, v = np.linalg.eig(Ham)    
       
        idx = w.argsort()  #this sorts the eigenvalues and eigenvectors 
        w = w[idx]
        v=  v[:,idx]
        #gstate  = v[:,0]
        #gstate /= np.vdot(gstate,gstate)
        #
        #
        #density = np.zeros((rDim,rDim))
        #for t in range(0,tDim):
        #    tbra = t*rDim
        #    for a in range(0,rDim):
        #        ket   = tbra+a
        #        kCoef = gstate[ket]
        #        if  (np.absolute(kCoef) != 0):
        #            for b in range(0,rDim):
        #                bra   = tbra+b
        #                bCoef = gstate[bra]
        #                if  (np.absolute(bCoef) != 0):
        #                    density[b,a] += np.vdot(bCoef,kCoef)
        #                    #density[a,b] += np.vdot(kCoef,bCoef)
       
        #density = density/np.trace(density)
        ##print density
        #s2 = -np.log(np.trace(np.linalg.matrix_power(density,2)))
        ##w2, v2 = np.linalg.eig(np.dot(density,density))    
        ##s2 = -np.log(np.sum(w2))
        #print 'Entropy: ', s2
        # Compute energy of the lowest lying eigenstates
        #outefile.write('%16.8E' %H)
        #print 'Energies: '
        for j in range(min(nlevels,len(w))):
            Es[j,l] = w[j]/float(N)
            #print 'E ',j, ' = ',Es[j,l]
    #    outefile.write('%16.8E' %(Es[j,i])
        #print ('Jz = %0.3f' %Jz) + (' Energy = %0.7f' %Es[0,l])+ (' Entropy = %0.7f' %s2)
    #outefile.write('\n')

    for j in range(nlevels):
        axE.plot(Jzs, Es[nlevels-j-1,:]-Es[0,:], ls = '-', marker='.', color=colors[(nlevels-j-1)%len(colors)], label=r"$\mathrm{E_{%d}}$" %(nlevels-j-1))
    # Save the ground state
    #outgfile.write('%16.8E' %Jz)
    #for k in range(HDim):
    #    outgfile.write('%16.8E' %np.real(v[k,0]))
    #outgfile.write('\n')
    #outefile.close()
    #outgfile.close()
#----------------------------------------------------------
#----------------------------------------------------------




#----------------------------------------------------------
#---------------------------Plotting-----------------------
#----------------------------------------------------------

    #for j in range(npoints):
    #    axG.plot(gstate, ls = '-', marker='s', color=colors[j%len(colors)], label=r"$\mathrm{J_z = %0.2f}$" %Jzs[j])

    pl.legend(loc='best',ncol=2)
    pl.tight_layout()
    pl.show()



#----------------------------------------------------------
#----------------------------------------------------------
if __name__ == "__main__": 
    main()

