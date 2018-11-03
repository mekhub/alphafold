from output_helpers import _show_results, _show_matrices
from copy import deepcopy
from alphafold.secstruct import *
from alphafold.explicit_recursions import *
from alphafold.dynamic_programming import *
from alphafold.backtrack import mfe

##################################################################################################
class AlphaFoldParams:
    '''
    Parameters that Define the statistical mechanical model for RNA folding
    '''
    def __init__( self ):
        # Seven parameter model
        self.C_init = 1.0     # Effective molarity for starting each loop (units of M)
        self.l      = 0.5     # Effective molarity penalty for each linkages in loop (dimensionless)
        self.Kd_BP  = 0.0002  # Kd for forming base pair (units of M )
        self.l_BP   = 0.2     # Effective molarity penalty for each base pair in loop (dimensionless)
        self.C_eff_stacked_pair = 1e4 # Effective molarity for forming stacked pair (units of M)
        self.K_coax = 100     # coax bonus for contiguous helices (dimensionless). Set to 0 to turn off coax (dimensionless)
        self.l_coax = 200     # Effective molarity bonus for each coaxial stack in loop. Initial guess: C_eff_stacked_pair / (C_init*l*K_coax)
        self.C_std = 1.0      # 1 M. drops out in end (up to overall scale factor).
        self.min_loop_length = 1 # Disallow apical loops smaller than this size (integer)
        self.allow_strained_3WJ = False # Prevent strained three-way-junctions with two helices coaxially stacked and no spacer nucleotides to other helix.

    def get_variables( self ):
        return ( self.C_init, self.l, self.Kd_BP, self.l_BP, self.C_eff_stacked_pair, self.K_coax, self.l_coax, self.C_std, self.min_loop_length, self.allow_strained_3WJ )

##################################################################################################
def partition( sequences, params = AlphaFoldParams(), circle = False, verbose = False, calc_deriv = False, use_explicit_recursions = True, backtrack = False ):
    '''
    Wrapper function into Partition() class
    '''
    #if use_explicit_recursions: from explicit_recursions import *  # overrides simpler objected-oriented formulae in recursions.py
    p = Partition( sequences, params, calc_deriv )
    p.circle  = circle
    p.run()
    if backtrack: p.calc_mfe()
    if verbose: p.show_matrices()
    p.show_results()
    p.run_cross_checks()

    return ( p.Z_final[0].Q, p.bpp, p.bps_MFE, p.Z_final[0].dQ )

##################################################################################################
class Partition:
    '''
    Statistical mechanical model for RNA folding, testing a bunch of extensions and with lots of cross-checks.
    TODO: complete expressions for derivatives (only doing derivatives w.r.t. Kd_BP right now)
    TODO: replace dynamic programming matrices with a class that auto-updates derivatives, caches each contribution for backtracking, and automatically does the modulo N wrapping
    (C) R. Das, Stanford University, 2018
    '''
    def __init__( self, sequences, params, calc_deriv = False ):
        '''
        Required user input.
        sequences = string with sequence, or array of strings (sequences of interacting strands)
        params    = AlphaFoldParams object
        '''
        self.sequences = sequences
        self.params = params
        self.circle = False  # user can update later --> circularize sequence
        self.calc_deriv = calc_deriv
        self.calc_contrib = False
        self.bps_MFE = []
        return

    ##############################################################################################
    def run( self ):
        '''
        Do the dynamic programming to fill partition function matrices
        '''
        initialize_sequence_information( self ) # N, sequence, is_cutpoint, any_intervening_cutpoint
        initialize_dynamic_programming_matrices( self ) # ( Z_BP, C_eff, Z_linear, Z_cut, Z_coax; dZ_BP, dC_eff, dZ_linear, dZ_cut, dZ_coax )
        initialize_base_pair_types( self )

        # do the dynamic programming
        for offset in range( 1, self.N ): #length of subfragment
            for i in range( self.N ):     #index of subfragment
                j = (i + offset) % self.N;  # N cyclizes
                # some preliminary helpers
                update_Z_cut( self, i, j )
                # base pairs and co-axial stacks
                for base_pair_type in self.base_pair_types: update_Z_BPq( self, base_pair_type, i, j )
                update_Z_BP( self, i, j )
                update_Z_coax( self, i, j )
                # C_eff makes use of information on Z_BP, so compute last
                update_C_eff_basic( self, i, j )
                update_C_eff_no_BP_singlet( self, i, j )
                update_C_eff_no_coax_singlet( self, i, j )
                update_C_eff( self, i, j )
                update_Z_linear( self, i, j )

        get_Z_final( self )    # (Z, dZ_final)
        get_bpp_matrix( self ) # fill base pair probability matrix
        return

    def calc_mfe( self ): _calc_mfe( self )
    # boring member functions -- defined later.
    def show_results( self ): _show_results( self )
    def show_matrices( self ): _show_matrices( self )
    def run_cross_checks( self ): _run_cross_checks( self )

##################################################################################################
def get_bpp_matrix( self ):
    '''
    Getting base pair probability matrix.
    Gets carried out pretty fast since we've already computed the sum over structures in i..j encapsulated by a pair (i,j), as well
      as structures in j..i encapsulated by those pairs.
    So: it becomes easy to calculate partition function over all structures with base pair (i,j), and then divide by total Z.
    '''

    # base pair probability matrix
    self.bpp = initialize_zero_matrix( self.N )
    for i in range( self.N ):
        for j in range( self.N ):
            self.bpp[i][j] = self.Z_BP[i][j].Q * self.Z_BP[j][i].Q * self.params.Kd_BP / self.Z_final[0].Q

##################################################################################################
def _run_cross_checks( self ):
    # stringent test that partition function is correct -- all the Z(i,i) agree.
    for i in range( self.N ):
        assert( abs( ( self.Z_final[i].Q - self.Z_final[0].Q ) / self.Z_final[0].Q ) < 1.0e-5 )
        if self.calc_deriv and self.Z_final[0].dQ > 0:
            assert( self.Z_final[0].dQ == 0 or  abs( ( self.Z_final[i].dQ - self.Z_final[0].dQ ) / self.Z_final[0].dQ ) < 1.0e-5 )

    # calculate bpp_tot = -dlog Z_final /dlog Kd_BP in two ways! wow cool test
    if self.calc_deriv:
        bpp_tot = 0.0
        for i in range( self.N ):
            for j in range( self.N ):
                bpp_tot += self.bpp[i][j]/2.0 # to avoid double counting (i,j) and (j,i)
        bpp_tot_based_on_deriv = -self.Z_final[0].dQ * self.params.Kd_BP / self.Z_final[0].Q
        if bpp_tot > 0: assert( abs( ( bpp_tot - bpp_tot_based_on_deriv )/bpp_tot ) < 1.0e-5 )

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
    any_intervening_cutpoint = any cutpoint exists between i and j (N X N)
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

    self.any_intervening_cutpoint = initialize_any_intervening_cutpoint( self.is_cutpoint )

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
    self.Z_BP     = DynamicProgrammingMatrix( N );
    self.Z_linear = DynamicProgrammingMatrix( N, diag_val = 1.0 );
    self.Z_cut    = DynamicProgrammingMatrix( N );
    self.Z_coax   = DynamicProgrammingMatrix( N );
    self.C_eff_basic           = DynamicProgrammingMatrix( N, diag_val = self.params.C_init );
    self.C_eff_no_coax_singlet = DynamicProgrammingMatrix( N, diag_val = self.params.C_init );
    self.C_eff_no_BP_singlet   = DynamicProgrammingMatrix( N, diag_val = self.params.C_init );
    self.C_eff                 = DynamicProgrammingMatrix( N, diag_val = self.params.C_init );

##################################################################################################
def initialize_zero_matrix( N ):
    X = WrappedArray( N )
    for i in range( N ): X[ i ] = WrappedArray( N )
    return X

##################################################################################################
class BasePairType:
    def __init__( self, nt1, nt2, Kd_BP, N ):
        '''
        Uh, a little weird to have Z in here.
        '''
        self.nt1 = nt1
        self.nt2 = nt2
        self.Kd_BP = Kd_BP
        self.Z_BP  = DynamicProgrammingMatrix( N );
        self.match_lowercase = ( nt1 == '' and nt2 == '' )

##################################################################################################
def initialize_base_pair_types( self ):
    N = self.N
    self.base_pair_types = []
    self.base_pair_types.append( BasePairType( 'C', 'G', self.params.Kd_BP, N ) )
    self.base_pair_types.append( BasePairType( 'G', 'C', self.params.Kd_BP, N ) )
    self.base_pair_types.append( BasePairType( 'A', 'U', self.params.Kd_BP, N ) )
    self.base_pair_types.append( BasePairType( 'U', 'A', self.params.Kd_BP, N ) )
    self.base_pair_types.append( BasePairType( '', '', self.params.Kd_BP, N ) ) # generic match

##################################################################################################
def initialize_any_intervening_cutpoint( is_cutpoint ):
    N = len( is_cutpoint )
    any_intervening_cutpoint = [[]]*N
    for i in range( N ): any_intervening_cutpoint[i] = [False]*N
    for i in range( N ): #index of subfragment
        found_cutpoint = False
        any_intervening_cutpoint[ i ][ i ] = False
        for offset in range( N ): #length of subfragment
            j = (i + offset) % N;  # N cyclizes
            any_intervening_cutpoint[ i ][ j ] = found_cutpoint
            if is_cutpoint[ j ]: found_cutpoint = True
    return any_intervening_cutpoint

##################################################################################################
def _calc_mfe( self ):
    N = self.N
    p_MFE = [0.0]*N
    bps_MFE = [[]]*N
    self.calc_contrib = True
    get_Z_final( self )
    for i in range( 1 ):
        (bps_MFE[i], p_MFE[i] ) = mfe( self, self.Z_final[i].contribs )
        assert( abs( ( p_MFE[i] - p_MFE[0] ) / p_MFE[0] ) < 1.0e-5 )
    print
    print 'Doing backtrack to get minimum free energy structure:'
    print  secstruct(bps_MFE[0],N), "   ", p_MFE[0], "[MFE]"
    print
    self.bps_MFE = bps_MFE



