from output_helpers import _show_results, _show_matrices
from copy import deepcopy

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
def partition( sequences, params = AlphaFoldParams(), circle = False, verbose = False, calc_deriv = False ):
    '''
    Wrapper function into Partition() class
    '''
    p = Partition( sequences, params, calc_deriv )
    p.circle  = circle
    p.run()
    if verbose: p.show_matrices()
    p.show_results()
    p.run_cross_checks()

    return ( p.Z_final[0], p.bpp, p.dZ_final[0] )

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
        return

    ##############################################################################################
    def run( self ):
        '''
        Do the dynamic programming to fill partition function matrices
        '''
        initialize_sequence_information( self ) # N, sequence, is_cutpoint, any_intervening_cutpoint
        initialize_dynamic_programming_matrices( self ) # ( Z_BP, C_eff, Z_linear, Z_cut, Z_coax; dZ_BP, dC_eff, dZ_linear, dZ_cut, dZ_coax )

        # do the dynamic programming
        for offset in range( 1, self.N ): #length of subfragment
            for i in range( self.N ):     #index of subfragment
                j = (i + offset) % self.N;  # N cyclizes
                # some preliminary helpers
                update_Z_cut( self, i, j )
                # base pairs and co-axial stacks
                update_Z_BP( self, i, j )
                update_Z_coax( self, i, j )
                # C_eff makes use of information on Z_BP, so compute last
                update_C_eff( self, i, j )
                update_Z_linear( self, i, j )

        get_Z_final( self )    # (Z, dZ_final)
        get_bpp_matrix( self ) # fill base pair probability matrix
        return

    # boring member functions -- defined later.
    def show_results( self ): _show_results( self )
    def show_matrices( self ): _show_matrices( self )
    def run_cross_checks( self ): _run_cross_checks( self )

##################################################################################################
# Following four functions hold ALL THE GOOD STUFF.
##################################################################################################
def update_Z_cut( self, i, j ):
    '''
    Z_cut is the partition function for independently combining one contiguous/bonded segment emerging out of i to a cutpoint c, and another segment that goes from c+1 to j.
    Useful for Z_BP and Z_final calcs below.
    Analogous to 'exterior' Z in Mathews calc & Dirks multistrand calc.
    '''
    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, is_cutpoint, any_intervening_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear, Z_cut, dZ_cut, Z_coax, dZ_coax ) = unpack_variables( self )
    offset = ( j - i ) % N
    for c in range( i, i+offset ):
        if is_cutpoint[c % N]:
            # strand 1  (i --> c), strand 2  (c+1 -- > j)
            Z_seg1  = Z_seg2  = 1
            dZ_seg1 = dZ_seg2 = 0
            if c != i :
                Z_seg1  = Z_linear [(i+1) % N][c % N]
                dZ_seg1 = dZ_linear[(i+1) % N][c % N]
            if (c+1)%N != j:
                Z_seg2  = Z_linear [(c+1) % N][(j-1) % N]
                dZ_seg2 = dZ_linear[(c+1) % N][(j-1) % N]
            Z_cut [i][j] += Z_seg1 * Z_seg2
            dZ_cut[i][j] += dZ_seg1 * Z_seg2 + Z_seg1 * dZ_seg2
            #Z_cut_contrib[i][j].append( Z_linear_contrib

##################################################################################################
def update_Z_BP( self, i, j ):
    '''
    Z_BP is the partition function for all structures that base pair i and j.
    Relies on previous Z_BP, C_eff, Z_linear available for subfragments.
    '''
    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, is_cutpoint, any_intervening_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear, Z_cut, dZ_cut, Z_coax, dZ_coax ) = unpack_variables( self )
    offset = ( j - i ) % N

    ( C_eff_for_coax, dC_eff_for_coax, C_eff_for_BP, dC_eff_for_BP ) = \
    (C_eff, dC_eff, C_eff, dC_eff) if allow_strained_3WJ else (self.C_eff_no_BP_singlet.X, self.C_eff_no_BP_singlet.dX, self.C_eff_no_coax_singlet.X, self.C_eff_no_coax_singlet.dX)

    # minimum loop length -- no other way to penalize short segments.
    if ( not any_intervening_cutpoint[i][j] and ( ((j-i-1) % N)) < min_loop_length ): return
    if ( not any_intervening_cutpoint[j][i] and ( ((i-j-1) % N)) < min_loop_length ): return

    for base_pair_type in self.base_pair_types:
        if base_pair_type.match_lowercase:
            if not (sequence[i].islower() and sequence[j].islower() and sequence[i]==sequence[j] ): continue
        else:
            if not ( sequence[i] == base_pair_type.nt1 and sequence[ j ] == base_pair_type.nt2 ): continue

        Z_BPq  = base_pair_type.Z_BP
        dZ_BPq = base_pair_type.dZ_BP
        Kd_BPq = base_pair_type.Kd_BP

        if (not is_cutpoint[ i ]) and (not is_cutpoint[ (j-1) % N]):
            # base pair closes a loop
            #
            #    ~~~~~~
            #   ~      ~
            # i+1      j-1
            #   \       /
            #    i ... j
            #
            Z_BPq[i][j]  += (1.0/Kd_BPq ) * ( C_eff_for_BP [(i+1) % N][(j-1) % N] * l * l * l_BP)
            dZ_BPq[i][j] += (1.0/Kd_BPq ) * ( dC_eff_for_BP[(i+1) % N][(j-1) % N] * l * l * l_BP)
            #Z_BP_contrib[i][j].append( [[C_eff_for_BP, i+1, j-1, (1.0/Kd_BPq)*l*l*l_BP]] )

            # base pair forms a stacked pair with previous pair
            #      ___
            #     /   \
            #  i+1 ... j-1
            #    |     |
            #    i ... j
            #
            # TODO: generalize C_eff_stacked_pair to be function of base pairs q (at i,j) and r (at i+1,j-1)
            Z_BPq[i][j]  += (1.0/Kd_BPq ) * C_eff_stacked_pair * Z_BP[(i+1) % N][(j-1) % N]
            dZ_BPq[i][j] += (1.0/Kd_BPq ) * C_eff_stacked_pair * dZ_BP[(i+1) % N][(j-1) % N]

        # base pair brings together two strands that were previously disconnected
        #
        #   \       /
        #    i ... j
        #
        Z_BPq [i][j] += (C_std/Kd_BPq) * Z_cut[i][j]
        dZ_BPq[i][j] += (C_std/Kd_BPq) * dZ_cut[i][j]

        if (not is_cutpoint[i]) and (not is_cutpoint[j-1]):

            # coaxial stack of bp (i,j) and (i+1,k)...  "left stack",  and closes loop on right.
            #      ___
            #     /   \
            #  i+1 ... k - k+1 ~
            #    |              ~
            #    i ... j - j-1 ~
            #
            for k in range( i+2, i+offset-1 ):
                if not is_cutpoint[k % N]:
                    Z_BPq [i][j] += Z_BP[(i+1) % N][k % N] * C_eff_for_coax[(k+1) % N][(j-1) % N] * l**2 * l_coax * K_coax / Kd_BPq
                    dZ_BPq[i][j] += (dZ_BP[(i+1) % N][k % N] * C_eff_for_coax[(k+1) % N][(j-1) % N] +
                                    Z_BP[(i+1) % N][k % N] * dC_eff_for_coax[(k+1) % N][(j-1) % N] ) * l**2 * l_coax * K_coax / Kd_BPq

            # coaxial stack of bp (i,j) and (k,j-1)...  close loop on left, and "right stack"
            #            ___
            #           /   \
            #  ~ k-1 - k ... j-1
            # ~              |
            #  ~ i+1 - i ... j
            #
            for k in range( i+2, i+offset-1 ):
                if not is_cutpoint[(k-1) % N]:
                    Z_BPq [i][j] += C_eff_for_coax[(i+1) % N][(k-1) % N] * Z_BP[k % N][(j-1) % N] * l**2 * l_coax * K_coax / Kd_BPq
                    dZ_BPq[i][j] += (dC_eff_for_coax[(i+1) % N][(k-1) % N] * Z_BP[k % N][(j-1) % N] +
                                    C_eff_for_coax[(i+1) % N][(k-1) % N] * dZ_BP[k % N][(j-1) % N] ) * l**2 * l_coax * K_coax / Kd_BPq

        # "left stack" but no loop closed on right (free strands hanging off j end)
        #      ___
        #     /   \
        #  i+1 ... k -
        #    |
        #    i ... j -
        #
        if not is_cutpoint[ i ]:
            for k in range( i+2, i+offset ):
                Z_BPq[i][j] += Z_BP[(i+1) % N][k % N] * Z_cut[k % N][j] * C_std * K_coax / Kd_BPq
                dZ_BPq[i][j] += (dZ_BP[(i+1) % N][k % N] * Z_cut[k % N][j] + Z_BP[(i+1) % N][k % N] * dZ_cut[k % N][j] ) * C_std * K_coax / Kd_BPq

        # "right stack" but no loop closed on left (free strands hanging off i end)
        #       ___
        #      /   \
        #   - k ... j-1
        #           |
        #   - i ... j
        #
        if not is_cutpoint[(j-1) % N]:
            for k in range( i, i+offset-1 ):
                Z_BPq[i][j] += Z_cut[i][k % N] * Z_BP[k % N][(j-1) % N] * C_std * K_coax / Kd_BPq
                dZ_BPq[i][j]+= ( dZ_cut[i][k % N] * Z_BP[k % N][(j-1) % N]  + Z_cut[i][k % N] * dZ_BP[k % N][(j-1) % N] ) * C_std * K_coax / Kd_BPq

        # key 'special sauce' for derivative w.r.t. Kd_BP
        dZ_BPq[i][j] += -(1.0/Kd_BPq) * Z_BPq[i][j]

        Z_BP[i][j]  = Z_BPq[i][j]
        dZ_BP[i][j] = dZ_BPq[i][j]
        #Z_BP_contrib[i][j] += Z_BPq_contrib[i][j]

##################################################################################################
def update_Z_coax( self, i, j ):
    '''
    Z_coax(i,j) is the partition function for all structures that form coaxial stacks between (i,k) and (k+1,j) for some k
    '''
    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, is_cutpoint, any_intervening_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear, Z_cut, dZ_cut, Z_coax, dZ_coax ) = unpack_variables( self )
    offset = ( j - i ) % N

    #  all structures that form coaxial stacks between (i,k) and (k+1,j) for some k
    #
    #       -- k - k+1 -
    #      /   :    :   \
    #      \   :    :   /
    #       -- i    j --
    #
    for k in range( i+1, i+offset-1 ):
        if not is_cutpoint[ k % N ]:
            Z_coax[i][j]  += Z_BP[i][k % N] * Z_BP[(k+1) % N][j] * K_coax
            dZ_coax[i][j] += (dZ_BP[i][k % N] * Z_BP[(k+1) % N][j] + Z_BP[i][k % N] * dZ_BP[(k+1) % N][j]) * K_coax

##################################################################################################
def update_C_eff( self, i, j ):
    '''
    C_eff tracks the effective molarity of a loop starting at i and ending at j
    Assumes a model where each additional element multiplicatively reduces the effective molarity, by
      the variables l, l_BP, C_eff_stacked_pair, K_coax, etc.
    Relies on previous Z_BP, C_eff, Z_linear available for subfragments.
    Relies on Z_BP being already filled out for i,j
    TODO: In near future, will include possibility of multiple C_eff terms, which combined together will
      allow for free energy costs of loop closure to scale apprpoximately log-linearly rather than
      linearly with loop size.
    '''
    offset = ( j - i ) % self.N

    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, is_cutpoint, any_intervening_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear, Z_cut, dZ_cut, Z_coax, dZ_coax ) = unpack_variables( self )

    exclude_strained_3WJ = (not allow_strained_3WJ) and (offset == N-1) and (not is_cutpoint[j] )

    # j is not base paired or coaxially stacked: Extension by one residue from j-1 to j.
    #
    #    i ~~~~~~ j-1 - j
    #
    if not is_cutpoint[(j-1) % N]:
        C_eff[i][j]  += C_eff[i][(j-1) % N] * l
        dC_eff[i][j] += dC_eff[i][(j-1) % N] * l
        #C_eff_contrib[i][j].append( [[C_eff, i, j-1, l]] )

    # j is base paired, and its partner is k > i. (look below for case with i and j base paired)
    #                 ___
    #                /   \
    #    i ~~~~k-1 - k...j
    #
    (C_eff_for_BP, dC_eff_for_BP) = (self.C_eff_no_coax_singlet.X, self.C_eff_no_coax_singlet.dX ) if exclude_strained_3WJ else (C_eff, dC_eff)
    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            C_eff[i][j]  += C_eff_for_BP[i][(k-1) % N] * l * Z_BP[k % N][j] * l_BP
            dC_eff[i][j] += ( dC_eff_for_BP[i][(k-1) % N] * Z_BP[k % N][j] + C_eff_for_BP[i][(k-1) % N] * dZ_BP[k % N][j] ) * l * l_BP

    # j is coax-stacked, and its partner is k > i.  (look below for case with i and j coaxially stacked)
    #               _______
    #              / :   : \
    #              \ :   : /
    #    i ~~~~k-1 - k   j
    #
    (C_eff_for_coax, dC_eff_for_coax) = (self.C_eff_no_BP_singlet.X, self.C_eff_no_BP_singlet.dX) if exclude_strained_3WJ else (C_eff, dC_eff)
    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            C_eff[i][j]  += C_eff_for_coax[i][(k-1) % N] * Z_coax[k % N][j] * l * l_coax
            dC_eff[i][j] += (dC_eff_for_coax[i][(k-1) % N] * Z_coax[k % N][j] + C_eff_for_coax[i][(k-1) % N] * dZ_coax[k % N][j]) * l * l_coax

    # some helper arrays that prevent closure of any 3WJ with a single coaxial stack and single helix with not intervening loop nucleotides
    self.C_eff_no_coax_singlet.X[i][j] =  C_eff[i][j] + C_init *  Z_BP[i][j] * l_BP
    self.C_eff_no_coax_singlet.dX[i][j] = dC_eff[i][j] + C_init * dZ_BP[i][j] * l_BP

    self.C_eff_no_BP_singlet.X[i][j] =  C_eff[i][j] + C_init *  Z_coax[i][j] * l_coax
    self.C_eff_no_BP_singlet.dX[i][j] = dC_eff[i][j] + C_init * dZ_coax[i][j] * l_coax

    # j is base paired, and its partner is i
    #      ___
    #     /   \
    #  i+1 ... j-1
    #    |     |
    #    i ... j
    #
    C_eff[i][j]  += C_init * Z_BP[i][j] * l_BP
    dC_eff[i][j] += C_init * dZ_BP[i][j] * l_BP
    #C_eff_contrib[i][j].append( [[Z_BP, i, j, C_init * l_BP]] )

    # j is coax-stacked, and its partner is i.
    #       ------------
    #      /   :    :   \
    #      \   :    :   /
    #       -- i    j --
    #
    C_eff[i][j]  += C_init * Z_coax[i][j] * l_coax
    dC_eff[i][j] += C_init * dZ_coax[i][j] * l_coax

##################################################################################################
def update_Z_linear( self, i, j ):
    '''
    Z_linear tracks the total partition function from i to j, assuming all intervening residues are covalently connected (or base-paired).
    Relies on previous Z_BP, C_eff, Z_linear available for subfragments.
    Relies on Z_BP being already filled out for i,j
    '''
    offset = ( j - i ) % self.N

    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, is_cutpoint, any_intervening_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear, Z_cut, dZ_cut, Z_coax, dZ_coax ) = unpack_variables( self )

    # j is not base paired: Extension by one residue from j-1 to j.
    #
    #    i ~~~~~~ j-1 - j
    #
    if not is_cutpoint[(j-1) % N]:
        Z_linear[i][j]  += Z_linear[i][(j - 1) % N]
        dZ_linear[i][j] += dZ_linear[i][(j - 1) % N]

    # j is base paired, and its partner is i
    #                 ___
    #                /   \
    #    i ~~~~k-1 - k...j
    #
    Z_linear[i][j]  += Z_BP[i][j]
    dZ_linear[i][j] += dZ_BP[i][j]

    # j is base paired, and its partner is k > i
    #                 ___
    #                /   \
    #    i ~~~~k-1 - k...j
    #
    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            Z_linear[i][j]  += Z_linear[i][(k-1) % N] * Z_BP[k % N][j]
            dZ_linear[i][j] += ( dZ_linear[i][(k-1) % N] * Z_BP[k % N][j] + Z_linear[i][(k-1) % N] * dZ_BP[k % N][j] )

    # j is coax-stacked, and its partner is i.
    #       ------------
    #      /   :    :   \
    #      \   :    :   /
    #       -- i    j --
    #
    Z_linear[i][j]  += Z_coax[i][j]
    dZ_linear[i][j] += dZ_coax[i][j]

    # j is coax-stacked, and its partner is k > i.
    #
    #               _______
    #              / :   : \
    #              \ :   : /
    #    i ~~~~k-1 - k   j
    #
    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            Z_linear[i][j]  += Z_linear[i][(k-1) % N] * Z_coax[k % N][j]
            dZ_linear[i][j] += dZ_linear[i][(k-1) % N] * Z_coax[k % N][j] + Z_linear[i][(k-1) % N] * dZ_coax[k % N][j]

##################################################################################################
def get_Z_final( self ):
    # Z_final is total partition function, and is computed at end of filling dynamic programming arrays
    # Get the answer (in N ways!) --> so final output is actually Z_final(i), an array.
    # Equality of the array is tested in run_cross_checks()
    Z_final  = []
    dZ_final = []

    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, is_cutpoint, any_intervening_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear, Z_cut, dZ_cut, Z_coax, dZ_coax ) = unpack_variables( self )

    for i in range( N ):
        Z_final.append( 0.0 )
        dZ_final.append( 0.0 )

        if self.is_cutpoint[(i + N - 1) % N]:
            #
            #      i ------- i-1
            #
            #     or equivalently
            #        ________
            #       /        \
            #       \        /
            #        i-1    i
            #
            Z_final[i]  += Z_linear[i][(i-1) % N]
            dZ_final[i] += dZ_linear[i][(i-1) % N]
            #Z_final_contrib[i].append( Z_linear, i, i-1, 1.0 )
        else:
            # Need to 'ligate' across i-1 to i
            # Scaling Z_final by Kd_lig/C_std to match previous literature conventions

            # Need to remove Z_coax contribution from C_eff, since its covered by C_eff_stacked_pair below.
            Z_final[i]  += self.C_eff_no_coax_singlet.X[i][(i - 1) % N] * l / C_std
            dZ_final[i] += self.C_eff_no_coax_singlet.dX[i][(i - 1) % N] * l / C_std

            for c in range( i, i + N - 1):
                if self.is_cutpoint[c % N]:
                    #any split segments, combined independently
                    #
                    #   c+1 --- i-1 - i --- c
                    #               *
                    Z_final[i]  += Z_linear[i][c % N] * Z_linear[(c+1) % N][(i-1) % N ]
                    dZ_final[i] += ( dZ_linear[i][c % N] * Z_linear[(c+1) % N][(i-1) % N ] + Z_linear[i][c % N] * dZ_linear[(c+1) % N][(i-1) % N ] )

            # base pair forms a stacked pair with previous pair
            #
            #   - j+1 - j -
            #      :    :
            #      :    :
            #   - i-1 - i -
            #         *
            for j in range( i+1, (i + N - 1) ):
                if not is_cutpoint[ j % N ]:
                    Z_final[i]  += C_eff_stacked_pair * Z_BP[i % N][j % N] * Z_BP[(j+1) % N][(i - 1) % N]
                    dZ_final[i] += C_eff_stacked_pair * (dZ_BP[i % N][j % N] * Z_BP[(j+1) % N][(i - 1) % N] + Z_BP[i % N][j % N] * dZ_BP[(j+1) % N][(i - 1) % N] )

            (C_eff_for_coax, dC_eff_for_coax ) = (C_eff, dC_eff ) if allow_strained_3WJ else (self.C_eff_no_BP_singlet.X,self.C_eff_no_BP_singlet.dX)

            # New co-axial stack might form across ligation junction
            for j in range( i + 1, i + N - 2):
                # If the two coaxially stacked base pairs are connected by a loop.
                #
                #       ~~~~
                #   -- k    j --
                #  /   :    :   \
                #  \   :    :   /
                #   - i-1 - i --
                #         *
                for k in range( j + 2, i + N - 1):
                    if is_cutpoint[j % N]: continue
                    if is_cutpoint[(k-1) % N]: continue
                    Z_final[i]  += Z_BP[i][j % N] * C_eff_for_coax[(j+1) % N][(k-1) % N] * Z_BP[k % N][(i-1) % N] * l * l * l_coax * K_coax
                    dZ_final[i] += (dZ_BP[i][j % N] *  C_eff_for_coax[(j+1) % N][(k-1) % N] *  Z_BP[k % N][(i-1) % N] +
                                     Z_BP[i][j % N] * dC_eff_for_coax[(j+1) % N][(k-1) % N] *  Z_BP[k % N][(i-1) % N] +
                                     Z_BP[i][j % N] *  C_eff_for_coax[(j+1) % N][(k-1) % N] * dZ_BP[k % N][(i-1) % N]) * l * l * l_coax * K_coax

                # If the two stacked base pairs are in split segments
                #
                #      \    /
                #   -- k    j --
                #  /   :    :   \
                #  \   :    :   /
                #   - i-1 - i --
                #         *
                for k in range( j + 1, i + N - 1):
                    Z_final[i]  += Z_BP[i][j % N] * Z_cut[j % N][k % N] * Z_BP[k % N][(i-1) % N] * K_coax
                    dZ_final[i] += (dZ_BP[i][j % N] * Z_cut[j % N][k % N] * Z_BP[k % N][(i-1) % N] +
                                    Z_BP[i][j % N] * dZ_cut[j % N][k % N] * Z_BP[k % N][(i-1) % N] +
                                    Z_BP[i][j % N] * Z_cut[j % N][k % N] * dZ_BP[k % N][(i-1) % N]) * K_coax

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
    self.bpp = DynamicProgrammingData( self.N );
    for i in range( self.N ):
        for j in range( self.N ):
            self.bpp.X[i][j] = self.Z_BP[i][j] * self.Z_BP[j][i] * self.params.Kd_BP / self.Z_final[0]

##################################################################################################
def _run_cross_checks( self ):
    # stringent test that partition function is correct -- all the Z(i,i) agree.
    for i in range( self.N ):
        assert( abs( ( self.Z_final[i] - self.Z_final[0] ) / self.Z_final[0] ) < 1.0e-5 )
        if self.calc_deriv and self.dZ_final[0] > 0:
            assert( self.dZ_final[0] == 0 or  abs( ( self.dZ_final[i] - self.dZ_final[0] ) / self.dZ_final[0] ) < 1.0e-5 )

    # calculate bpp_tot = -dlog Z_final /dlog Kd_BP in two ways! wow cool test
    if self.calc_deriv:
        bpp_tot = 0.0
        for i in range( self.N ):
            for j in range( self.N ):
                bpp_tot += self.bpp[i][j]/2.0 # to avoid double counting (i,j) and (j,i)
        bpp_tot_based_on_deriv = -self.dZ_final[0] * self.params.Kd_BP / self.Z_final[0]
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

###################################################################################################################33
class DynamicProgrammingData:
    '''
    Dynamic programming object, with derivs and contribution accumulation.
     X   = values (N x N)
     dX  = derivatives (N X N)
     X_contrib = contributions (coming soon)
    '''
    def __init__( self, N ):
        self.X = []
        for i in range( N ): self.X.append( [0.0]*N )

        self.dX = deepcopy( self.X ) # another zero matrix.

        self.X_contrib = []
        for i in range( N ): self.X_contrib.append( [[]]*N )

    def __getitem__( self, idx ):
        # overloaded []. warning: overhead! directly access object.X[ idx ] in inner loops.
        return self.X[ idx ]

    def __len__( self ): return len( self.X )

    def add( self, i, j, b ):
        #  trying out a function that might make code more readable,
        #  (could hide all contribution accumulation for backtracking -- and
        #   perhaps even derivatives -- inside class!)
        # but this kind of thing appears to take up too much overhead.
        self.X[i][j]  += b
        self.dX[i][j] += 0
        self.X_contrib[i][j].append( [i,j,b] )

##################################################################################################
class BasePairType:
    def __init__( self, nt1, nt2, Kd_BP, N ):
        '''
        Uh, a little weird to have Z in here.
        '''
        self.nt1 = nt1
        self.nt2 = nt2
        self.Kd_BP = Kd_BP
        self.Z_BP_DP  = DynamicProgrammingData( N );
        self.Z_BP  = self.Z_BP_DP.X
        self.dZ_BP = self.Z_BP_DP.dX;
        self.match_lowercase = ( nt1 == '' and nt2 == '' )

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
    self.Z_BP     = DynamicProgrammingData( N );
    self.Z_linear = DynamicProgrammingData( N );
    self.Z_cut    = DynamicProgrammingData( N );
    self.Z_coax   = DynamicProgrammingData( N );
    self.C_eff    = DynamicProgrammingData( N );
    for i in range( N ): #length of fragment
        self.Z_linear.X[i][i] = 1
        self.C_eff.X[i][i]                 = self.params.C_init
    self.C_eff_no_coax_singlet = deepcopy( self.C_eff )
    self.C_eff_no_BP_singlet   = deepcopy( self.C_eff )

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
def unpack_variables( self ):
    '''
    This helper function just lets me write out equations without
    using "self" which obscures connection to my handwritten equations
    In C++, will just use convention of object variables like N_, sequence_.
    '''
    return self.params.get_variables() + \
           ( self.N, self.sequence, self.is_cutpoint, self.any_intervening_cutpoint,  \
             self.Z_BP.X, self.Z_BP.dX, self.C_eff.X, self.C_eff.dX, self.Z_linear.X, self.Z_linear.dX, self.Z_cut.X, self.Z_cut.dX, self.Z_coax.X, self.Z_coax.dX )

