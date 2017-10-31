# A polished version of grad ascent implementation with the help LBFGS. It is capable of training
# a semi-restricted Boltzmann machines in classical, quantum and semi-quantum modes. Not fully tested.

import Hfile, Hbuilder, semiRBM
import argparse, os
import numpy.random as rm
import numpy as np
from scipy.optimize import fmin_l_bfgs_b as LBFGS

record_prog = True

#------------------------------------------------------------------------------
# Call back function in the LBFGS routine called every time a new step is performed.
# It is used to track the progress of the optimization.
#------------------------------------------------------------------------------
def prog_tracker(inter):

    global record_prog
    record_prog = True

    return 0


#------------------------------------------------------------------------------
# The function to optimize with the LBFGS routine. It evaluates the
# Log-Likelihood and its gradient.
#------------------------------------------------------------------------------
def comp_LL_grad(params, *args):

    # Unpack all the arguments
    (N_v, N_h, beta, bonds, data, weights, train_mode, file_name) = args

    # Initialize the number of spins, the data set length and the number of bonds
    N_s = N_h + N_v
    N_d = len(data)
    N_b = len(bonds)



    # Initialize Boltzmann machine
    couplings = unpack_params(params, N_s, train_mode)
    rbm = semiRBM.Semi_Restricted_Boltzmann_Machine(N_v, N_h, beta, couplings)


    # Compute the log-likelihood
    LL = 0
    for index, data_vector in enumerate(data[:-1]):

        # Clamp the visible units
        rbm.set_projector(data_vector)

        # Add up the contribution of this vector
        P   = np.real(rbm.eval_projector())
        LL -= np.log(P) * weights[index]



    # Compute local averages
    aves  = np.zeros((N_d, N_s*2  + N_b))

    for index, data_vector in enumerate(data):
        # Convert the data vector from the basis position to its z-basis representation
        # Note: 0 -> 1, 1 -> -1 in the conversion
        data_vector_z = -1.0*(np.array(data_vector)*2-1)

        # Duplicate the vector. The new vector will also contain hidden units
        # averages when their number is non-zero
        vector_z = data_vector_z

        # If we can compute the expectation values easily, we do it
        if (train_mode in ['class', 'quant_approx']) and (index!=N_d-1):

           # Compute 1-body Z operators for the visible units
           aves[index, :N_v] =  data_vector_z* beta



           # Similar computation for the hidden units is more involved
           if (N_h > 0):

               # Get J couplings between the hidden and visible layers
               J = rbm.get_rest_layer_couplings()

               # Select longitudinal fields over the hidden units
               h = couplings[N_v : N_s]

               # Select transverse field over the hidden units
               d = couplings[N_s + N_v : N_s + N_s]

               # Compute the effective longitudinal field over the hidden units
               h_eff = h + np.dot(J, data_vector_z.reshape(-1, 1)).T[0]

               # Compute D ref. Eq 36
               D = np.sqrt(d*d + h_eff*h_eff)

               # Finally, compute Z operator expectations over the hidden units

               hidden_vector_z  =  h_eff/D * np.tanh(D)
               aves[index, N_v : N_s] = hidden_vector_z * beta

               # Treat the expectation values over the hidden units as observed values
               vector_z = np.hstack([data_vector_z, hidden_vector_z])

               #print vector_z
           # 2-body ZZ operators

           # Compute ZZ expectation for all the bonds
           for bindex, bond in enumerate(bonds):

               (s_0, s_1) = bond
               # Note - sign cancels out during the conversion
               aves[index, 2*N_s + bindex] = vector_z[s_0] * vector_z[s_1] * beta


        # otherwise we resort to the computationally intensive method
        else:

            # Clamp the visible units
            rbm.set_projector(data_vector)

            aves[index,:] = rbm.comp_local_averages()

    gradient = np.sum((aves*weights),  axis=0)

    # Record progress to file once in a while
    global record_prog
    if  record_prog:

        # Build up the string to be recorded
        new_line = '    %16.8E' %LL
        for val in couplings:  new_line+= '    %16.8E' %val
        for val in aves[-1,:]: new_line+= '    %16.8E' %val
        new_line += '\n'

        # Write it to the file
        f = open(file_name, 'a')
        f.write(new_line)
        f.close()

    # Prevent writing to the file at each call to the function
    record_prog = False

    grad_LL = pack_params(gradient, N_s, train_mode)

    return LL, grad_LL


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def unpack_params(params, N_s, train_mode):

    # classical
    if  train_mode in ["class"]:
        (H, J)  = np.split(params, [N_s])
        D       = np.zeros(N_s)

    # quantum
    elif train_mode in ["quant"]:
        (H, delta, J)  = np.split(params, [N_s, N_s+1])
        D = np.ones(N_s)*delta

    # quantum without the need to optimize over the transverse field
    elif train_mode in ["quant_fixed", "quant_approx"]:
        (H, J) = np.split(params, [N_s])
        global delta
        D = np.ones(N_s)*delta

    return np.hstack([H, D, J])


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def pack_params(gradient, N_s, train_mode):

    (Z, X, ZZ) = np.split(gradient, [N_s, N_s*2])
    if   train_mode in ["class", "quant_fixed", 'quant_approx']:
         return np.hstack([Z, ZZ])

    elif train_mode == "quant":
         return np.hstack([Z, np.sum(X), ZZ])






#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def bitfield(n):
    return np.array([1 if digit=='1' else 0 for digit in bin(n)[2:]])



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def main():

    modes  =  ['class', 'quant', 'quant_fixed', 'quant_approx']

    # setup the command line parser options
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--seed',  help='Seed used to generate data', type=int, default=0)
    parser.add_argument('--data',  help='Data file',                  type=str)
    parser.add_argument('--delta', help='Transverse field',           type=float)
    parser.add_argument('--beta',  help='Inverse temperature ',       type=float)
    parser.add_argument('--Nv',    help='Number of visible units',    type=int)
    parser.add_argument('--Nh',    help='Number of hidden units',    type=int)
    parser.add_argument('--mode',  help='Training mode', choices=modes)

    args = vars(parser.parse_args())
    rm.seed(args['seed'])

    train_mode = args['mode']
    beta       = args['beta']
    N_v        = args['Nv']
    N_h        = args['Nh']

    global delta
    if args['delta']!= None:
       delta=args['delta']


    N_s = N_v + N_h


    # Load the training data
    data          = np.loadtxt(args['data'], dtype='int32')
    data, weights = np.unique(data, return_counts=True, axis=0)
    weights       = weights/float(np.sum(weights))

    # Compute the data set entropy
    entropy = -1.0*np.sum(weights*np.log(weights))
    print 'Data set entropy: %0.2f' %entropy

    # Add an empty list to the data set. It is equivalent for an unclamped machine.
    data  = data.tolist()
    data += [[]]                        # empty vector represents unclamped visible units
    weights = np.append(weights, -beta) # and its weight

    # Record the size of the data set and transpose the weights vector for future manipulations
    N_d   = len(data)
    weights = weights.reshape((N_d,1))



    # We need to know the connectivity of the machine for some computations outside of the RBM
    bonds, offset = Hbuilder.semiRBM_Connected(N_v, N_h)
    N_b   = len(bonds)


    # Create the log file
    file_name = 'train_'

    if args['mode'] == 'class': file_name += 'mode-%s_' %train_mode
    else:                       file_name += 'mode-%s_delta-%04.2f_' %(train_mode, delta)

    file_name += os.path.split(args['data'])[1][5:]
    f = open(file_name, 'w')

    # Record the header
    header = '#%19s' %('LL')
    for i in range(N_s):   header += '%20s' %('H'+str(i))
    for i in range(N_s):   header += '%20s' %('D'+str(i))
    for bond in bonds:     header += '%20s' %('J(%d, %d)' %(bond[0], bond[1]))

    for i in range(N_s):   header += '%20s' %('<Z'+str(i)+'>')
    for i in range(N_s):   header += '%20s' %('<X'+str(i)+'>')
    for bond in bonds:     header += '%20s' %('<ZZ(%d, %d)>' %(bond[0], bond[1]))

    header +='\n'
    f.write(header)



    # Generate an initial guess for Hamiltonian parameters
    H  = (rm.rand(N_s)-0.5)*0.1
    J  = (rm.rand(N_b)-0.5)*0.1

    if train_mode in ['class', 'quant_fixed', 'quant_approx']:
         params = np.hstack([H, J])

    elif train_mode == 'quant':
         params = np.hstack([H, delta, J])



    # Perform the gradient descent with the help of LBFGS procedure
    print '---------Gradient descent----------'
    fx, fLL, info = LBFGS(comp_LL_grad, x0=params, args = (N_v, N_h, beta, bonds, data, weights, train_mode, file_name), iprint = 1, pgtol=1e-05, factr=1e13/2.2/10000.0, maxiter=50, callback=prog_tracker)
    #fx, fLL, info = LBFGS(LL, x0=params, fprime=grad, args = (Ns, beta, bonds, data, weights, args['mode']), iprint = 1, pgtol=1e-08, factr=1e13/2.2/10000.0, maxiter=500, callback=progtracker)

    print fx
    print info

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()



