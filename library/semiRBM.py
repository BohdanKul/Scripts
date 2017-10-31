import numpy         as np
import scipy.sparse  as sparse
import scipy.linalg  as slin
import Hbuilder
import Hfile
import mexp
from   scipy.integrate import simps, trapz


#--------------------------------------------------------------------------------------------------
# [0,0,1] --> 1
#--------------------------------------------------------------------------------------------------
def Bits2Int(bitlist):
    out = 0
    for bit in bitlist:
        out = (out << 1) | bit

    return out

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
class Semi_Restricted_Boltzmann_Machine:
    def __init__(self, N_v, N_h, beta, couplings):

        # The size of hilbert space associated with the hidden units
        self.hidden_dim  = 2**N_h

        # Spind numbers
        self.N_spins = N_v + N_h
        self.N_h     = N_h

        # Fields are stored as bonds acting over a single site
        self.bonds = [[e] for e in (range(self.N_spins) + range(self.N_spins))]

        # Add a set of semi-restricted couplings
        semi_rest, self.fully_connected_offset = Hbuilder.semiRBM_Connected(N_v, N_h)
        self.bonds = self.bonds + semi_rest

        self.couplings   = couplings
        self.N_couplings = len(self.bonds)

        # Introduce the matrix form of local elements in the Hamiltonian
        I = np.eye(2)
        X = np.array([[0, 1], [1, 0]])
        Z = np.array([[1, 0], [0, -1]])


        # The Hamiltonian
        self.H = np.zeros((2**self.N_spins, 2**self.N_spins))

        # The list of local operators we are interested in.
        # The elements of the list are tuples (operator_representation, is_diagonal_flag)
        self.oper_list = []

        # Build and store those objects at the same time
        for index, (coupling, bond) in enumerate(zip(couplings, self.bonds[:])):

            # The operator type can be determined by its index
            if  index in range(self.N_spins, self.N_spins*2):
                oper = Hbuilder.EmbedOper(X, self.N_spins, bond+[1.0])

                # for the off-diagonal operators, store both row and column indices of non-zero elements
                self.oper_list  += [(oper.nonzero(), False)]

            else:
                oper = Hbuilder.EmbedOper(Z, self.N_spins, bond+[1.0])

                # keep track only of the diagonal (non-zero) elements for the diagonal operators
                self.oper_list += [(oper.diagonal(), True)]

            self.H = self.H + oper*coupling


        # Compute the full unnormalized density matrix
        self.beta = beta
        self.rho  = np.matrix(mexp.expm(-self.H*beta), copy=False)


        # Store separately its diagonal elements and normalize them
        self.norm_rho_diag  = np.diagonal(self.rho).copy()
        self.normalization  = np.sum(self.norm_rho_diag)
        self.norm_rho_diag /= self.normalization

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
    def get_bonds(self):

        return self.bonds


#--------------------------------------------------------------------------------------------------
# Return couplings between hidden and visible units as a matrix. Each row represents
# a set of couplings between a hidden unit and all of the visible ones.
#--------------------------------------------------------------------------------------------------
    def get_rest_layer_couplings(self):

        interlayer_couplings = self.couplings[2 * self.N_spins + self.fully_connected_offset:]
        interlayer_couplings = interlayer_couplings.reshape(self.N_h, -1)

        return interlayer_couplings


#--------------------------------------------------------------------------------------------------
#   Clamp visible units to a particular state
#--------------------------------------------------------------------------------------------------
    def set_projector(self, cbits):

        if  len(cbits)>0:
            self.is_clamped = True
            self.proj_index = Bits2Int(cbits + [0]*self.N_h)
            #print cbits, self.proj_index
        else:
            self.is_clamped = False

#--------------------------------------------------------------------------------------------------
#   Evaluate the probability of generating a particular state of visible units from the QBM
#--------------------------------------------------------------------------------------------------
    def eval_projector(self):

        if self.is_clamped: return np.sum( self.norm_rho_diag[self.proj_index : self.proj_index + self.hidden_dim] )
        else:               return 0


#--------------------------------------------------------------------------------------------------
#   Compute the expectation values of local operators. It works for clamped and unclamped QBMs
#--------------------------------------------------------------------------------------------------
    def comp_local_averages(self):
        if   self.is_clamped: return self.comp_local_averages_clamped()
        else:                 return self.comp_local_averages_unclamped()


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
    def comp_local_averages_unclamped(self):

        # Normalize the imaginary time evolution operator
        U = self.rho/self.normalization

        # Initialize to 0 the storage for the local expectation values
        aves = np.zeros(self.N_couplings)

        # Compute all the local expectation values
        for index, oper_repr in enumerate(self.oper_list):
            aves[index] = self.comp_local_operator(oper_repr, U)

        return aves

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
    def comp_local_averages_clamped(self):

        # Define the imaginary-time propagation operator
        U  = np.matrix(np.zeros((2**self.N_spins, 2**self.N_spins)), copy=False)

        # Apply the projector operator from the right
        U[:, self.proj_index : self.proj_index + self.hidden_dim]  = self.rho[:,  self.proj_index : self.proj_index + self.hidden_dim ]

        # Normalize U by Tr_{h} <projector, h| e^{-H \beta} |projector, h>
        clamped_norm = np.sum( np.diagonal(self.rho)[self.proj_index : self.proj_index + self.hidden_dim ])
        U /= clamped_norm

        # Define discrete imaginary-time step propagators
        N_grid    = 12
        tau       = 1.0*self.beta/(2.0*N_grid)

        Ub0 = np.matrix( mexp.expm( self.H*tau), copy=False)
        Uf0 = np.matrix( mexp.expm(-self.H*tau), copy=False)

        # Evaluate time-evolved operators on a discrete grid
	grid_aves = np.zeros((N_grid+1, self.N_couplings))

        for step in range(N_grid+1):

            # At tau = 0, both propagators are just identities
            if  step>0:
                U = Ub0 * U * Uf0

            # compute all the local expectation values
            for index, oper_repr in enumerate(self.oper_list):
                grid_aves[step][index] = self.comp_local_operator(oper_repr, U)



        # Numerically integrate
        aves = np.zeros(self.N_couplings)
        xs   = np.linspace(0.0, N_grid+1-1, N_grid+1)*tau

        for index in range(self.N_couplings):
            aves[index] = trapz(grid_aves[:, index], xs)*2.0

        return aves


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
    def comp_local_operator(self, oper_repr, U):
        '''
            The second element in oper_repr tuple indicates whether the operator is diagonal.
            If it is diagonal:     the first element are the operator's diagonal elements
            If it is off-diagonal: the first element are indices of the operator's non-zero elements
                                   expressed as a tuple (rows, cols)
        '''

        is_diagonal = oper_repr[1]
        if is_diagonal:
            # Only diagonal elements of U are needed to compute <ZZ>
            # (all others yield zero when multiplied with ZZ)

            return np.sum(oper_repr[0] * U.diagonal().T )
        else:
            # The sum of diagonal elements resulting from the multiplication
            # of X with U can be computed by selecting the elements in U.T
            # with indices such that the same elements in X are non-zero
            # and summing them over.
            (X_rows, X_cols) = oper_repr[0]

            # notice the column and row indices are reversed (U.T)
            return np.sum(U[X_cols, X_rows])


