#import scipy.sparse.eye as eye
from numpy import frexp 
from math  import ceil
from scipy.sparse import issparse

import numpy               as np
import scipy.sparse        as sp
import numpy.linalg        as nlin
import scipy.sparse.linalg as slin

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
        if issparse(A):
            A2 = A*A
            A4 = A2*A2 
            A6 = A2*A2
            U =  A*(A6*(c[13]*A6 + c[11]*A4 + c[9]*A2) + c[7]*A6 + c[5]*A4 + c[3]*A2 + c[1]*sp.eye(n))
            V =     A6*(c[12]*A6 + c[10]*A4 + c[8]*A2) + c[6]*A6 + c[4]*A4 + c[2]*A2 + c[0]*sp.eye(n)
                    
            F =  (slin.inv(V-U))*(V+U)

             
        else:
            A2 = np.dot(A,A)
            A4 = np.dot(A2, A2) 
            A6 = np.dot(A2, A4)
            U =  np.dot(A, np.dot(A6, c[13]*A6 + c[11]*A4 + c[9]*A2) + c[7]*A6 + c[5]*A4 + c[3]*A2 + c[1]*np.eye(n))
            V =            np.dot(A6, c[12]*A6 + c[10]*A4 + c[8]*A2) + c[6]*A6 + c[4]*A4 + c[2]*A2 + c[0]*np.eye(n)
                    
            F =  np.dot(nlin.inv(V-U), (V+U))

    else: # m == 3, 5, 7, 9
        #Apowers = []*int(ceil((m+1)/2.0))
        Apowers = []

        if issparse(A):
            Apowers.append(sp.eye(n,n))
            Apowers.append(A*A)
        else:
            Apowers.append(np.eye(n,n))
            Apowers.append(np.dot(A,A))
        
        for j in range(2, int(ceil((m+1)/2.0))):
            if issparse(A): Apowers.append(Apowers[j-1] * Apowers[1])
            else:           Apowers.append(np.dot(Apowers[j-1], Apowers[1]))
    
        if issparse(A):
            U = sp.csc_matrix((n,n)) 
            V = sp.csc_matrix((n,n))
        else:
            U = np.zeros((n,n)) 
            V = np.zeros((n,n))
        
        # Not sure about the right point of the range
        for j in  range(m+1, 1, -2):     
            U = U + c[j-1]*Apowers[(j-1)//2]
    
        if issparse(A):  U = A*U
        else: U = np.dot(A, U)
        
        for j in  range(m, 0, -2):     
            V = V + c[j-1]*Apowers[(j+1-1)//2]
    
        if issparse(A): F =       (slin.inv(V-U)*(V+U))
        else:           F = np.dot(nlin.inv(V-U),(V+U))
   
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

    normA = 0
    if issparse(A): normA = np.amax((A.multiply(A.sign())).sum(0)) 
    else:           normA = nlin.norm(A,1) 
    
    if normA <= theta[-1]:
        # no scaling and squaring is required.
        for i in range(len(m_vals)):
            if normA <= theta[i]:
                F = PadeApproximantOfDegree(A, m_vals[i])
                break
    else:
        t,s = frexp(normA/float(theta[-1]))
        s = s - (t == 0.5) # adjust s if normA/theta(end) is a power of 2.
        A = A/(2.0**s)     # Scaling
        F = PadeApproximantOfDegree(A, m_vals[-1])
        
        for i in range(s):
            if issparse(A): F = F*F
            else:           F = np.dot(F,F)   

    return F
# End of expm
