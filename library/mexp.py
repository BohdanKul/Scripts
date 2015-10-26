#import scipy.sparse.eye as eye
from numpy        import eye
from numpy        import zeros
from numpy        import frexp 
from numpy.linalg import norm 
from numpy.linalg import inv
from math         import ceil

import numpy as np

def expmchk():
# EXPMCHK Check the class of input A and
#    initialize M_VALS and THETA accordingly.
    m_vals = np.array([3, 5, 7, 9, 13])
    theta  = np.array([ 
                        1.495585217958292e-002,  # m_vals = 3
                        2.539398330063230e-001,  # m_vals = 5
                        9.504178996162932e-001,  # m_vals = 7
                        2.097847961257068e+000,  # m_vals = 9
                        5.371920351148152e+000
                     ])# m_vals = 13
    
    return m_vals, theta

def getPadeCoefficients(m):
# GETPADECOEFFICIENTS Coefficients of numerator P of Pade approximant
#    C = GETPADECOEFFICIENTS returns coefficients of numerator
#    of [M/M] Pade approximant, where M = 3,5,7,9,13.
    if   (m==3):  c = np.array([120, 60, 12, 1])
    elif (m==5):  c = np.array([30240, 15120, 3360, 420, 30, 1])
    elif (m==7):  c = np.array([17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1])
    elif (m==9):  c = np.array([17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1])
    elif (m==13): c = np.array([
                                 64764752532480000, 32382376266240000, 7771770303897600, 
                                 1187353796428800,  129060195264000,   10559470521600, 
                                 670442572800,      33522128640,       1323241920,
                                 40840800,          960960,            16380,  182,  1
                               ])
        
    return c      

def PadeApproximantOfDegree(A,m):
#PADEAPPROXIMANTOFDEGREE  Pade approximant to exponential.
#   F = PADEAPPROXIMANTOFDEGREE(M) is the degree M diagonal
#   Pade approximant to EXP(A), where M = 3, 5, 7, 9 or 13.
#   Series are evaluated in decreasing order of powers, which is
#   in approx. increasing order of maximum norms of the terms.

    n = max(A.shape)
    c = getPadeCoefficients(m)

    # Evaluate Pade approximant.
    if (m == 13):
        # For optimal evaluation need different formula for m >= 12.
        A2 = A*A
        A4 = A2*A2 
        A6 = A2*A4
        U = A * (A6*(c[13]*A6 + c[11]*A4 + c[9]*A2) + c[7]*A6 + c[5]*A4 + c[3]*A2 + c[1]*eye(n,n) )
        V =      A6*(c[12]*A6 + c[10]*A4 + c[8]*A2) + c[6]*A6 + c[4]*A4 + c[2]*A2 + c[0]*eye(n,n)
                
        F = inv(V-U)*(V+U)

    else: # m == 3, 5, 7, 9
        #Apowers = []*int(ceil((m+1)/2.0))
        Apowers = []

        Apowers.append(eye(n,n))
        Apowers.append(A*A)
        
        for j in range(2, int(ceil((m+1)/2.0))):
            Apowers.append(Apowers[j-1]*Apowers[1])
    
        U = zeros((n,n)) 
        V = zeros((n,n))
        
        # Not sure about the right point of the range
        for j in  range(m+1, 1, -2):     
            U = U + c[j-1]*Apowers[(j-1)//2]
    
        U = A*U
        
        for j in  range(m, 0, -2):     
            V = V + c[j-1]*Apowers[(j+1-1)//2]
    
        F = inv(V-U)*(V+U)
    
    return F

def expm(A):
# EXPM   Matrix exponential.
#   EXPM(X) is the matrix exponential of X.  EXPM is computed using
#   a scaling and squaring algorithm with a Pade approximation.
#
# Julia implementation closely based on MATLAB code by Nicholas Higham
#

# Initialization
    m_vals, theta = expmchk()

    normA = norm(A,1)

    if normA <= theta[-1]:
        # no scaling and squaring is required.
        for i in range(len(m_vals)):
            if normA <= theta[i]:
                F = PadeApproximantOfDegree(A,m_vals[i])
                break
    else:
        t,s = frexp(normA/float(theta[-1]))
        s = s - (t == 0.5) # adjust s if normA/theta(end) is a power of 2.
        A = (A/2.0)**s          # Scaling
        F = PadeApproximantOfDegree(A,m_vals[-1])
        
        for i in range(s):
            F = F*F   # Squaring

    return F
# End of expm
