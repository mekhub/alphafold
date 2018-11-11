from output_helpers import _show_results, _show_matrices
from copy import deepcopy
from alphafold.secstruct import *
from alphafold.wrapped_array import WrappedArray
from alphafold.backtrack import mfe, boltzmann_sample, enumerative_backtrack
from alphafold.parameters import AlphaFoldParams

##################################################################################################
def partition( sequences, params = AlphaFoldParams(), circle = False, verbose = False,  mfe = False, calc_bpp = False,
               n_stochastic = 0, do_enumeration = False, calc_deriv = False, use_simple_recursions = False  ):
    '''
    Wrapper function into Partition() class
    Returns Partition object p which holds results like:

      p.Z   = final partition function (where unfolded state has unity weight)
      p.bpp = matrix of base pair probabilities (if requested by user with calc_bpp = True)
      p.struct_MFE = minimum free energy secondary structure in dot-parens notation
      p.bps_MFE  = minimum free energy secondary structure as sorted list of base pairs
      p.dZ  = derivative of Z w.r.t. Kd_BP (will later generalize) (if requested by user with calc_deriv = True)

    '''
    p = Partition( sequences, params, calc_deriv, calc_all_elements = calc_bpp )
    p.use_simple_recursions = use_simple_recursions
    p.circle  = circle
    p.run()
    if calc_bpp:         p.get_bpp_matrix()
    if mfe:              p.calc_mfe()
    if n_stochastic > 0: p.stochastic_backtrack( n_stochastic )
    if do_enumeration:   p.enumerative_backtrack()
    if verbose:          p.show_matrices()
    p.show_results()
    p.run_cross_checks()

    return p

##################################################################################################
class Partition:
    '''
    Statistical mechanical model for RNA folding, testing a bunch of extensions and with lots of cross-checks.
    TODO: complete expressions for derivatives (only doing derivatives w.r.t. Kd_BP right now)
    TODO: replace dynamic programming matrices with a class that auto-updates derivatives, caches each contribution for backtracking, and automatically does the modulo N wrapping
    (C) R. Das, Stanford University, 2018
    '''
    def __init__( self, sequences, params, calc_deriv = False, calc_all_elements = False ):
        '''
        Required user input.
        sequences = string with sequence, or array of strings (sequences of interacting strands)
        params    = AlphaFoldParams object
        '''
        self.sequences = sequences
        self.params = params
        self.circle = False  # user can update later --> circularize sequence
        self.options = PartitionOptions( calc_deriv = calc_deriv )
        self.use_simple_recursions = False
        self.calc_all_elements     = calc_all_elements
        self.calc_bpp = False

        # for output:
        self.Z       = 0
        self.bpp     = []
        self.bps_MFE = []
        self.struct_MFE = ''
        self.struct_stochastic = []
        self.struct_enumerate  = []
        return

    ##############################################################################################
    def run( self ):
        '''
        Do the dynamic programming to fill partition function matrices
        '''
        initialize_sequence_information( self ) # N, sequence, ligated, all_ligated
        initialize_base_pair_types( self ) # C-G base pair, etc.
        initialize_dynamic_programming_matrices( self ) # ( Z_BP, C_eff, Z_linear, Z_cut, Z_coax, etc. )

        # do the dynamic programming
        for offset in range( 1, self.N ): #length of subfragment
            for i in range( self.N ):     #index of subfragment
                if (not self.calc_all_elements) and ( i + offset ) >= self.N: continue
                j = (i + offset) % self.N;  # N cyclizes
                for Z in self.Z_all: Z.update( self, i, j )

        for i in range( self.N): self.Z_final.update( self, i )

        self.Z  = self.Z_final.val(0)
        self.dZ = self.Z_final.deriv(0)

    # boring member functions -- defined later.
    def get_bpp_matrix( self ): _get_bpp_matrix( self ) # fill base pair probability matrix
    def calc_mfe( self ): _calc_mfe( self )
    def stochastic_backtrack( self, N ): _stochastic_backtrack( self, N )
    def enumerative_backtrack( self ): _enumerative_backtrack( self )
    def show_results( self ): _show_results( self )
    def show_matrices( self ): _show_matrices( self )
    def run_cross_checks( self ): _run_cross_checks( self )

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
    sequences = self.sequences
    if isinstance( sequences, str ): self.sequence = sequences
    else:
        self.sequence = ''
        for i in range( len( sequences ) ): self.sequence += sequences[i]
    self.N = len( self.sequence )
    N = self.N

    ligated = [True] * N
    if self.use_simple_recursions: ligated = WrappedArray( N, True )

    if isinstance( sequences, list ):
        L = 0
        for i in range( len(sequences)-1 ):
            L = L + len( sequences[i] )
            ligated[ L-1 ] = False
    if not self.circle: ligated[ N-1 ] = False

    self.ligated = ligated
    self.all_ligated = initialize_all_ligated( ligated )

##################################################################################################
class PartitionOptions:
    def __init__( self, calc_deriv = False ):
        self.calc_deriv   = calc_deriv
        self.calc_contrib = False

##################################################################################################
class BasePairType:
    def __init__( self, nt1, nt2, Kd_BP, match_lowercase = False ):
        '''
        Two sequence characters that get matched to base pair, e.g., 'C' and 'G';
           and the dissociation constant K_d (in M) associated with initial pair formation.
        match_lowercase means match 'x' to 'x', 'y' to 'y', etc.
        TODO: also store chemical modification info.
        '''
        self.nt1 = nt1
        self.nt2 = nt2
        self.Kd_BP = Kd_BP
        self.match_lowercase = ( nt1 == '*' and nt2 == '*' and match_lowercase )

##################################################################################################
def initialize_base_pair_types( self ):
    N = self.N
    self.base_pair_types = []
    self.base_pair_types.append( BasePairType( 'C', 'G', self.params.Kd_BP ) )
    self.base_pair_types.append( BasePairType( 'G', 'C', self.params.Kd_BP ) )
    self.base_pair_types.append( BasePairType( 'A', 'U', self.params.Kd_BP ) )
    self.base_pair_types.append( BasePairType( 'U', 'A', self.params.Kd_BP ) )
    self.base_pair_types.append( BasePairType( '*', '*', self.params.Kd_BP, match_lowercase = True ) )

##################################################################################################
def initialize_dynamic_programming_matrices( self ):
    '''
    A bunch of zero matrices. Only non-trivial thing is
    initialization of (i,i) [diagonal]:
         Z_BP(i,i)     = 0
         C_eff(i,i)    = C_init (units of M)
         Z_linear(i,i) = 1
    And: the order of initialization in the list Z_all will
      determine the actual order of updates during dynamic programmming at each (i,j).
    '''

    from explicit_recursions import update_Z_BPq, update_Z_BP, update_Z_cut, update_Z_coax, update_C_eff_basic, update_C_eff_no_BP_singlet, update_C_eff_no_coax_singlet, update_C_eff, update_Z_final, update_Z_linear
    from explicit_dynamic_programming import DynamicProgrammingMatrix, DynamicProgrammingList
    if self.use_simple_recursions: # over-ride with simpler recursions that are easier for user to input.
        from recursions import update_Z_BPq, update_Z_BP, update_Z_cut, update_Z_coax, update_C_eff_basic, update_C_eff_no_BP_singlet, update_C_eff_no_coax_singlet, update_C_eff, update_Z_final, update_Z_linear
        from dynamic_programming import DynamicProgrammingMatrix, DynamicProgrammingList

    N = self.N

    # Collection of all N X N dynamic programming matrices -- order in this list will
    #  determine order of updates.
    self.Z_all = Z_all = []

    # some preliminary helpers
    self.Z_cut    = DynamicProgrammingMatrix( N, DPlist = Z_all, update_func = update_Z_cut, options = self.options );

    # base pairs and co-axial stacks
    self.Z_BPq = {}
    for base_pair_type in self.base_pair_types:
        # the bpt = base_pair_type holds the base_pair_type info in the lambda (Python FAQ)
        update_func = lambda partition,i,j,bpt=base_pair_type: update_Z_BPq(partition,i,j,bpt)
        self.Z_BPq[ base_pair_type ] = DynamicProgrammingMatrix( N, DPlist = Z_all,
                                                                 update_func = update_func, options = self.options )
    self.Z_BP     = DynamicProgrammingMatrix( N, DPlist = Z_all, update_func = update_Z_BP, options = self.options );
    self.Z_coax   = DynamicProgrammingMatrix( N, DPlist = Z_all, update_func = update_Z_coax, options = self.options );

    # C_eff makes use of information on Z_BP, so compute last
    C_init = self.params.C_init
    self.C_eff_basic           = DynamicProgrammingMatrix( N, diag_val = C_init, DPlist = Z_all, update_func = update_C_eff_basic, options = self.options );
    self.C_eff_no_BP_singlet   = DynamicProgrammingMatrix( N, diag_val = C_init, DPlist = Z_all, update_func = update_C_eff_no_BP_singlet, options = self.options );
    self.C_eff_no_coax_singlet = DynamicProgrammingMatrix( N, diag_val = C_init, DPlist = Z_all, update_func = update_C_eff_no_coax_singlet, options = self.options );
    self.C_eff                 = DynamicProgrammingMatrix( N, diag_val = C_init, DPlist = Z_all, update_func = update_C_eff, options = self.options );

    self.Z_linear = DynamicProgrammingMatrix( N, diag_val = 1.0, DPlist = Z_all, update_func = update_Z_linear, options = self.options );

    # Last DP 1-D list (not a 2-D N x N matrix)
    self.Z_final = DynamicProgrammingList( N, update_func = update_Z_final, options = self.options  )

##################################################################################################
def initialize_zero_matrix( N ):
    from alphafold.dynamic_programming import WrappedArray
    X = WrappedArray( N )
    for i in range( N ): X[ i ] = WrappedArray( N )
    return X

##################################################################################################
def initialize_all_ligated( ligated ):
    '''
    all_ligated is needed to keep track of whether an apical loop is long enough
    to be 'closed' into a hairpin by base pair formation.
    TODO: alternatively could create a 'hairpin_OK' matrix -- that
          would be more analogous to pre-scanning for protein/ligand binding sites too.
    '''
    N = len( ligated )
    all_ligated = initialize_zero_matrix( N )
    for i in range( N ): #index of subfragment
        found_cutpoint = False
        all_ligated[ i ][ i ] = True
        for offset in range( N ): #length of subfragment
            j = (i + offset) % N;  # N cyclizes
            all_ligated[ i ][ j ] = ( not found_cutpoint )
            if not ligated[ j ]: found_cutpoint = True
    return all_ligated

##################################################################################################
def _get_bpp_matrix( self ):
    '''
    Getting base pair probability matrix.
    Gets carried out pretty fast since we've already computed the sum over structures in i..j encapsulated by a pair (i,j), as well
      as structures in j..i encapsulated by those pairs.
    So: it becomes easy to calculate partition function over all structures with base pair (i,j), and then divide by total Z.
    '''
    assert( self.calc_all_elements )
    self.bpp = initialize_zero_matrix( self.N )
    for i in range( self.N ):
        for j in range( self.N ):
            self.bpp[i][j] = self.Z_BP.val(i,j) * self.Z_BP.val(j,i) * self.params.Kd_BP / self.Z_final.val(0)

##################################################################################################
def _calc_mfe( self ):
    #
    # Wrapper into mfe(), written out in backtrack.py
    #
    N = self.N
    p_MFE   = [0.0]*N
    bps_MFE = [[]]*N

    # there are actually numerous ways to calculate MFE if we did all N^2 elements -- let's check.
    n_test = N if self.calc_all_elements else 1

    for i in range( n_test ):
        (bps_MFE[i], p_MFE[i] ) = mfe( self, self.Z_final.get_contribs(self,i) )
        assert( abs( ( p_MFE[i] - p_MFE[0] ) / p_MFE[0] ) < 1.0e-5 )
        assert( bps_MFE[i] == bps_MFE[0] )

    print
    print 'Doing backtrack to get minimum free energy structure:'
    print  secstruct(bps_MFE[0],N), "   ", p_MFE[0], "[MFE]"
    print
    self.bps_MFE = bps_MFE[0]
    self.struct_MFE = secstruct( bps_MFE[0], N)

##################################################################################################
def _stochastic_backtrack( self, N_backtrack ):
    #
    # Get stochastic, Boltzmann-weighted structural samples from partition function
    #
    print 'Doing',N_backtrack,'stochastic backtracks to get Boltzmann-weighted ensemble'
    for i in range( N_backtrack ):
        bps, p = boltzmann_sample( self, self.Z_final.get_contribs(self,0) )
        print secstruct(bps,self.N), "   ", p, "[stochastic]"
        self.struct_stochastic.append( secstruct(bps,self.N) )
    print

    return

##################################################################################################
def _enumerative_backtrack( self ):
    #
    # Enumerate all structures, and track their probabilities
    #
    p_bps = enumerative_backtrack( self )
    for (p,bps) in p_bps:
        print secstruct(bps,self.N), "   ", p, "[enumerative]"
        self.struct_enumerate.append( secstruct(bps,self.N) )
    p_tot = sum( p_bp[0] for p_bp in p_bps )
    print 'p_tot = ',p_tot
    assert( abs(p_tot - 1.0) < 1.0e-5 )
    return

##################################################################################################
def _run_cross_checks( self ):
    # stringent test that partition function is correct -- all the Z(i,i) agree.
    if self.calc_all_elements:
        for i in range( self.N ):
            assert( abs( ( self.Z_final.val(i) - self.Z_final.val(0) ) / self.Z_final.val(0) ) < 1.0e-5 )

    if self.calc_all_elements and self.options.calc_deriv and self.Z_final.deriv(0) > 0:
        for i in range( self.N ):
            assert( self.Z_final.deriv(0) == 0 or  abs( ( self.Z_final.deriv(i) - self.Z_final.deriv(0) ) / self.Z_final.deriv(0) ) < 1.0e-5 )

    # calculate bpp_tot = -dlog Z_final /dlog Kd_BP in two ways! wow cool test
    if self.options.calc_deriv and len(self.bpp)>0:
        bpp_tot = 0.0
        for i in range( self.N ):
            for j in range( self.N ):
                bpp_tot += self.bpp[i][j]/2.0 # to avoid double counting (i,j) and (j,i)
        bpp_tot_based_on_deriv = -self.Z_final.deriv(0) * self.params.Kd_BP / self.Z_final.val(0)
        print 'bpp_tot',bpp_tot,'bpp_tot_based_on_deriv',bpp_tot_based_on_deriv
        if bpp_tot > 0: assert( abs( ( bpp_tot - bpp_tot_based_on_deriv )/bpp_tot ) < 1.0e-5 )
