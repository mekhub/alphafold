from partition_helpers import *
from output_helpers import *

##################################################################################################
class AlphaFoldParams:
    '''
    Paramters that Define the statistical mechanical model for RNA folding
    '''
    def __init__( self ):
        # Four parameter model
        self.Kd_BP  = 0.001;
        self.C_init = 1          # a bit like exp(a) in multiloop
        self.l      = 0.5        # a bit like exp(b) in multiloop
        self.l_BP   = 0.1        # a bit like exp(c) in multiloop
        self.C_std = 1; # 1 M. drops out in end (up to overall scale factor).
        self.C_init_BP = self.C_init * (self.l_BP/self.l) # 0.2
        self.min_loop_length = 1

    def get_variables( self ):
        return ( self.Kd_BP, self.C_init, self.l, self.l_BP, self.C_std, self.C_init_BP, self.min_loop_length )

##################################################################################################
def partition( sequences, params = AlphaFoldParams(), circle = False, verbose = False ):
    '''
    Wrapper function into Partition() class
    '''
    p = Partition( sequences, params )
    p.circle  = circle
    p.run()
    if verbose: p.show_matrices()
    p.show_results()
    return ( p.Z_final[0], p.bpp, p.dZ_final[0] )

##################################################################################################
class Partition:
    '''
    Statistical mechanical model for RNA folding, testing a bunch of extensions and with lots of cross-checks
    (C) R. Das, Stanford University, 2018
    '''
    def __init__( self, sequences, params ):
        '''
        Required user input.
        sequences = string with sequence, or array of strings (sequences of interacting strands)
        params    = AlphaFoldParams object
        '''
        self.sequences = sequences
        self.params = params
        self.verbose = False # user can update later --> full matrix output
        self.circle = False  # user can update later --> circularize sequence
        return

    ##############################################################################################
    def run( self ):
        '''
        Do the dynamic programming to fill partition function matrices
        '''
        initialize_sequence_information( self ) # N, sequence, is_cutpoint, any_cutpoint
        initialize_dynamic_programming_matrices( self ) # ( Z_BP, C_eff, Z_linear, dZ_BP, dC_eff, dZ_linear )

        (Kd_BP, C_init, l, l_BP, C_std, C_init_BP, min_loop_length, N, \
         sequence, is_cutpoint, any_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear ) = unpack_variables( self )

        # do the dynamic programming
        for offset in range( 1, N ): #length of subfragment
            for i in range( N ): #index of subfragment
                j = (i + offset) % N;  # N cyclizes
                update_Z_BP( self, i, j )
                update_C_eff( self, i, j )
                update_Z_linear( self, i, j )

        get_Z_final( self )    # (Z, dZ_final)
        get_bpp_matrix( self ) # fill base pair probability matrix
        run_cross_checks( self ) # cross-check
        return

    ##############################################################################################
    def show_results( self ):
        print 'sequence =', self.sequence
        cutpoint = ''
        for i in range( self.N ):
            if self.is_cutpoint[ i ]: cutpoint += 'X'
            else: cutpoint += '-'
        print 'cutpoint =', cutpoint
        print 'circle   = ', self.circle
        print 'Z =',self.Z_final[0]
        return

    ##################################################################################################
    def show_matrices( self ):
        output_DP( "Z_BP", self.Z_BP )
        output_DP( "C_eff", self.C_eff, self.Z_final )
        output_DP( "Z_linear", self.Z_linear )
        output_square( "BPP", self.bpp );


##################################################################################################
# Following three functions hold ALL THE GOOD STUFF.
##################################################################################################
def update_Z_BP( self, i, j ):
    '''
    Z_BP is the partition function for all structures that base pair i and j.
    Relies on previous Z_BP, C_eff, Z_linear available for subfragments.
    TODO:  Will soon generalize to arbitrary number of base pairs beyond C-G
    '''
    (Kd_BP, C_init, l, l_BP, C_std, C_init_BP, min_loop_length, N, \
     sequence, is_cutpoint, any_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear ) = unpack_variables( self )
    offset = ( j - i ) % N

    # Residues that are base paired must *not* bring together contiguous loop with length less than min_loop length
    if (( sequence[i] == 'C' and sequence[j] == 'G' ) or ( sequence[i] == 'G' and sequence[j] == 'C' )) and \
          ( any_cutpoint[i][j] or ( ((j-i-1) % N)) >= min_loop_length  ) and \
          ( any_cutpoint[j][i] or ( ((i-j-1) % N)) >= min_loop_length  ):

        # base pair closes a loop
        if (not is_cutpoint[ i ]) and (not is_cutpoint[ (j-1) % N]):
            Z_BP[i][j]  += (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N] * l * l)
            dZ_BP[i][j] += (1.0/Kd_BP ) * ( dC_eff[(i+1) % N][(j-1) % N] * l * l)

        # base pair brings together two strands that were previously disconnected
        for c in range( i, i+offset ):
            if is_cutpoint[c % N]:
                # strand 1  (i --> c), strand 2  (c+1 -- > j)
                Z_comp1  = Z_comp2  = 1
                dZ_comp1 = dZ_comp2 = 0
                if c != i :
                    Z_comp1  = Z_linear [(i+1) % N][c % N]
                    dZ_comp1 = dZ_linear[(i+1) % N][c % N]
                if (c+1)%N != j:
                    Z_comp2  = Z_linear [(c+1) % N][(j-1) % N]
                    dZ_comp2 = dZ_linear[(c+1) % N][(j-1) % N]
                Z_product  = Z_comp1 * Z_comp2
                dZ_product = dZ_comp1 * Z_comp2 + Z_comp1 * dZ_comp2

                Z_BP [i][j] += (C_std/Kd_BP) * (l/l_BP) * Z_product
                dZ_BP[i][j] += (C_std/Kd_BP) * (l/l_BP) * dZ_product

    # key 'special sauce' for derivative w.r.t. Kd_BP
    dZ_BP[i][j] += -(1.0/Kd_BP) * Z_BP[i][j]

    return

##################################################################################################
def update_C_eff( self, i, j ):
    '''
    C_eff tracks the effective molarity of a loop starting at i and ending at j
    Assumes a model where each additional element multiplicatively reduces the effective molarity, by
      the variables l, l_BP, etc.
    Relies on previous Z_BP, C_eff, Z_linear available for subfragments.
    Relies on Z_BP being already filled out for i,j
    TODO: In near future, will include possibility of multiple C_eff terms, which combined together will
      allow for free energy costs of loop closure to scale apprpoximately log-linearly rather than
      linearly with loop size.
    '''
    offset = ( j - i ) % self.N

    (Kd_BP, C_init, l, l_BP, C_std, C_init_BP, min_loop_length, N, \
     sequence, is_cutpoint, any_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear ) = unpack_variables( self )

    # j is not base paired: Extension by one residue from j-1 to j.
    if not is_cutpoint[(j-1) % N]:
        C_eff[i][j]  += C_eff[i][(j-1) % N] * l
        dC_eff[i][j] += dC_eff[i][(j-1) % N] * l

    # j is base paired, and its partner is i
    C_eff[i][j]  += C_init_BP * Z_BP[i][j]
    dC_eff[i][j] += C_init_BP * dZ_BP[i][j]

    # j is base paired, and its partner is k > i
    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            C_eff[i][j]  += C_eff[i][(k-1) % N] * Z_BP[k % N][j] * l_BP
            dC_eff[i][j] += ( dC_eff[i][(k-1) % N] * Z_BP[k % N][j] + C_eff[i][(k-1) % N] * dZ_BP[k % N][j] ) * l_BP


##################################################################################################
def update_Z_linear( self, i, j ):
    '''
    Z_linear tracks the total partition function from i to j, assuming all intervening residues are covalently connected (or base-paired).
    Relies on previous Z_BP, C_eff, Z_linear available for subfragments.
    Relies on Z_BP being already filled out for i,j
    '''
    offset = ( j - i ) % self.N

    (Kd_BP, C_init, l, l_BP, C_std, C_init_BP, min_loop_length, N, \
     sequence, is_cutpoint, any_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear ) = unpack_variables( self )

    # j is not base paired: Extension by one residue from j-1 to j.
    if not is_cutpoint[(j-1) % N]:
        Z_linear[i][j]  += Z_linear[i][(j - 1) % N]
        dZ_linear[i][j] += dZ_linear[i][(j - 1) % N]

    # j is base paired, and its partner is i
    Z_linear[i][j]  += Z_BP[i][j]
    dZ_linear[i][j] += dZ_BP[i][j]

    # j is base paired, and its partner is k > i
    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            Z_linear[i][j]  += Z_linear[i][(k-1) % N] * Z_BP[k % N][j]
            dZ_linear[i][j] += ( dZ_linear[i][(k-1) % N] * Z_BP[k % N][j] + Z_linear[i][(k-1) % N] * dZ_BP[k % N][j] )

##################################################################################################
def get_Z_final( self ):
    # Z_final is total partition function, and is computed at end of filling dynamic programming arrays
    # Get the answer (in N ways!) --> so final output is actually Z_final(i), an array.
    # Equality of the array is tested in run_cross_checks()
    Z_final  = []
    dZ_final = []

    (Kd_BP, C_init, l, l_BP, C_std, C_init_BP, min_loop_length, N, \
     sequence, is_cutpoint, any_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear ) = unpack_variables( self )
    for i in range( N ):
        Z_final.append( 0 )
        dZ_final.append( 0 )

        if self.is_cutpoint[(i + N - 1) % N]:
            Z_final[i]  += Z_linear[i][(i-1) % N]
            dZ_final[i] += dZ_linear[i][(i-1) % N]
        else:
            # Scaling Z_final by Kd_lig/C_std to match previous literature conventions
            Z_final[i]  += C_eff[i][(i - 1) % N] * l / C_std
            dZ_final[i] += dC_eff[i][(i - 1) % N] * l / C_std

            for c in range( i, i + N - 1):
                if self.is_cutpoint[c % N]:
                    #any split segments, combined independently
                    Z_final[i]  += Z_linear[i][c % N] * Z_linear[(c+1) % N][ i - 1 ]
                    dZ_final[i] += ( dZ_linear[i][c % N] * Z_linear[(c+1) % N][ i - 1 ] + Z_linear[i][c % N] * dZ_linear[(c+1) % N][ i - 1 ] )

    self.Z_final = Z_final
    self.dZ_final = dZ_final

##################################################################################################
def get_bpp_matrix( self ):
    '''
    Getting base pair probability matrix.
    Gets carried out pretty fast since we've already computed the sum over structures in i..j encapsulated by a pair (i,j), as well
      as structures in j..i encapsulated by those pairs.
    So: it becomes easy to calculate partition function over all structures with base pair (i,j), and then divide by total Z.
    '''

    # base pair probability matrix
    self.bpp = initialize_zero_matrix( self.N );
    for i in range( self.N ):
        for j in range( self.N ):
            self.bpp[i][j] = self.Z_BP[i][j] * self.Z_BP[j][i] * self.params.Kd_BP * (self.params.l_BP / self.params.l) / self.Z_final[0]


##################################################################################################
def run_cross_checks( self ):
    # stringent test that partition function is correct -- all the Z(i,i) agree.
    for i in range( self.N ):
        assert( abs( ( self.Z_final[i] - self.Z_final[0] ) / self.Z_final[0] ) < 1.0e-5 )
        assert( abs( ( self.dZ_final[i] - self.dZ_final[0] ) / self.dZ_final[0] ) < 1.0e-5 )

    # calculate bpp_tot = -dlog Z_final /dlog Kd_BP in two ways! wow cool test
    bpp_tot = 0.0
    for i in range( self.N ):
        for j in range( self.N ):
            bpp_tot += self.bpp[i][j]/2.0 # to avoid double counting (i,j) and (j,i)
    bpp_tot_based_on_deriv = -self.dZ_final[0] * self.params.Kd_BP / self.Z_final[0]
    assert( abs( ( bpp_tot - bpp_tot_based_on_deriv )/bpp_tot ) < 1.0e-5 )

##################################################################################################
def initialize_sequence_information( self ):
    '''
    Create sequence information from sequences of strands:

    INPUT:
    sequences = sequences of interacting strands (array of strings)
    circle    = user asks for nucleotides N and 1 to be ligated ('circularized') (bool)

    OUTPUT:
    sequence     = concatenated sequence (string, length N)
    is_cutpoint  = is cut ('nick','chainbreak') or not (Array of bool, length N)
    any_cutpoint = any cutpoint exists between i and j (N X N)
    '''
    # initialize sequence
    if isinstance( self.sequences, str ): self.sequence = self.sequences
    else:
        self.sequence = ''
        for i in range( len( self.sequences ) ): self.sequence += self.sequences[i]
    self.N = len( self.sequence )

    # initialize cutpoint information
    self.is_cutpoint = [False] * self.N
    if isinstance( self.sequences, list ):
        L = 0
        for i in range( len(self.sequences)-1 ):
            L = L + len( self.sequences[i] )
            self.is_cutpoint[ L-1 ] = True
    if not self.circle: self.is_cutpoint[ self.N-1 ] = True

    self.any_cutpoint = initialize_any_cutpoint( self.is_cutpoint )

##################################################################################################
def initialize_dynamic_programming_matrices( self ):
    '''
    A bunch of zero matrices. Only non-trivial thing is
    initialization of (i,i) [diagonal]:
         Z_BP(i,i)     = 0
         C_eff(i,i)    = C_init (units of M)
         Z_linear(i,i) = 1
    '''
    N = self.N
    # initialize dynamic programming matrices
    self.C_eff    = initialize_zero_matrix( N );
    self.Z_BP     = initialize_zero_matrix( N );
    self.Z_linear = initialize_zero_matrix( N );
    for i in range( N ): #length of fragment
        self.C_eff[ i ][ i ] = self.params.C_init
        self.Z_linear[ i ][ i ] = 1

    # first calculate derivatives with respect to Kd_BP
    self.dC_eff    = initialize_zero_matrix( N );
    self.dZ_BP     = initialize_zero_matrix( N );
    self.dZ_linear = initialize_zero_matrix( N );


##################################################################################################
def unpack_variables( self ):
    '''
    This helper function just lets me write out equations without
    using "self" which obscures connection to my handwritten equations
    '''
    return (self.params.Kd_BP, self.params.C_init, self.params.l, self.params.l_BP, self.params.C_std, self.params.C_init_BP, \
            self.params.min_loop_length, \
            self.N, self.sequence, self.is_cutpoint, self.any_cutpoint,  \
            self.Z_BP, self.dZ_BP, self.C_eff, self.dC_eff, self.Z_linear, self.dZ_linear )

