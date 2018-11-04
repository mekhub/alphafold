from dynamic_programming import DynamicProgrammingData

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
     sequence, ligated, all_ligated, Z_BP, C_eff_basic, C_eff_no_BP_singlet, C_eff_no_coax_singlet, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv, calc_contrib ) = unpack_variables( self )
    offset = ( j - i ) % N
    for c in range( i, i+offset ):
        if not ligated[c % N]:
            # strand 1  (i --> c), strand 2  (c+1 -- > j)
            Z_seg1  = DynamicProgrammingData( 1.0 )
            Z_seg2  = DynamicProgrammingData( 1.0 )
            if c != i :
                Z_seg1  = Z_linear[(i+1) % N][c % N]
            if (c+1)%N != j:
                Z_seg2  = Z_linear[(c+1) % N][(j-1) % N]
            Z_cut[i][j] += Z_seg1 * Z_seg2
            #Z_cut[i][j].contribs.append( Z_linear

##################################################################################################l
def update_Z_BPq( self, i, j, base_pair_type ):
    '''
    Z_BPq is the partition function for all structures that base pair i and j with base_pair_type
    Relies on previous Z contributions available for subfragments, and Z_cut for this fragment i,j
    '''

    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, ligated, all_ligated, Z_BP, C_eff_basic, C_eff_no_BP_singlet, C_eff_no_coax_singlet, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv, calc_contrib ) = unpack_variables( self )
    offset = ( j - i ) % N

    ( C_eff_for_coax, C_eff_for_BP ) = (C_eff, C_eff ) if allow_strained_3WJ else (self.C_eff_no_BP_singlet, self.C_eff_no_coax_singlet )

    # minimum loop length -- no other way to penalize short segments.
    if ( all_ligated[i][j] and ( ((j-i-1) % N)) < min_loop_length ): return
    if ( all_ligated[j][i] and ( ((i-j-1) % N)) < min_loop_length ): return

    if base_pair_type.match_lowercase:
        if not (sequence[i].islower() and sequence[j].islower() and sequence[i]==sequence[j] ): return
    else:
        if not ( sequence[i] == base_pair_type.nt1 and sequence[ j ] == base_pair_type.nt2 ): return

    (Z_BPq, Kd_BPq)  = ( self.Z_BPq[ base_pair_type ], base_pair_type.Kd_BP )

    if (ligated[ i ]) and (ligated[ (j-1) % N]):
        # base pair closes a loop
        #
        #    ~~~~~~
        #   ~      ~
        # i+1      j-1
        #   \       /
        #    i ... j
        #
        Z_BPq[i][j]  += (1.0/Kd_BPq ) * ( C_eff_for_BP [(i+1) % N][(j-1) % N] * l * l * l_BP)

        # base pair forms a stacked pair with previous pair
        #      ___
        #     /   \
        #  i+1 ... j-1
        #    |     |
        #    i ... j
        #
        # TODO: generalize C_eff_stacked_pair to be function of base pairs q (at i,j) and r (at i+1,j-1)
        Z_BPq[i][j]  += (1.0/Kd_BPq ) * C_eff_stacked_pair * Z_BP[(i+1) % N][(j-1) % N]

    # base pair brings together two strands that were previously disconnected
    #
    #   \       /
    #    i ... j
    #
    Z_BPq [i][j] += (C_std/Kd_BPq) * Z_cut[i][j]

    if (ligated[i]) and (ligated[j-1]):

        # coaxial stack of bp (i,j) and (i+1,k)...  "left stack",  and closes loop on right.
        #      ___
        #     /   \
        #  i+1 ... k - k+1 ~
        #    |              ~
        #    i ... j - j-1 ~
        #
        for k in range( i+2, i+offset-1 ):
            if ligated[k % N]: Z_BPq [i][j] += Z_BP[(i+1) % N][k % N] * C_eff_for_coax[(k+1) % N][(j-1) % N] * l**2 * l_coax * K_coax / Kd_BPq

        # coaxial stack of bp (i,j) and (k,j-1)...  close loop on left, and "right stack"
        #            ___
        #           /   \
        #  ~ k-1 - k ... j-1
        # ~              |
        #  ~ i+1 - i ... j
        #
        for k in range( i+2, i+offset-1 ):
            if ligated[(k-1) % N]: Z_BPq [i][j] += C_eff_for_coax[(i+1) % N][(k-1) % N] * Z_BP[k % N][(j-1) % N] * l**2 * l_coax * K_coax / Kd_BPq

    # "left stack" but no loop closed on right (free strands hanging off j end)
    #      ___
    #     /   \
    #  i+1 ... k -
    #    |
    #    i ... j -
    #
    if ligated[ i ]:
        for k in range( i+2, i+offset ): Z_BPq[i][j] += Z_BP[(i+1) % N][k % N] * Z_cut[k % N][j] * C_std * K_coax / Kd_BPq

    # "right stack" but no loop closed on left (free strands hanging off i end)
    #       ___
    #      /   \
    #   - k ... j-1
    #           |
    #   - i ... j
    #
    if ligated[(j-1) % N]:
        for k in range( i, i+offset-1 ): Z_BPq[i][j] += Z_cut[i][k % N] * Z_BP[k % N][(j-1) % N] * C_std * K_coax / Kd_BPq

    # key 'special sauce' for derivative w.r.t. Kd_BP
    if calc_deriv: Z_BPq[i][j].dQ += -(1.0/Kd_BPq) * Z_BPq[i][j].Q


##################################################################################################
def update_Z_BP( self, i, j ):
    '''
    Z_BP is the partition function for all structures that base pair i and j.
    All the Z_BPq (partition functions for each base pair type) must have been
    filled in already for i,j.
    '''
    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, ligated, all_ligated, Z_BP, C_eff_basic, C_eff_no_BP_singlet, C_eff_no_coax_singlet, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv, calc_contrib ) = unpack_variables( self )

    for base_pair_type in self.base_pair_types:
        Z_BP[i][j]  += self.Z_BPq[ base_pair_type ][i][j]

##################################################################################################
def update_Z_coax( self, i, j ):
    '''
    Z_coax(i,j) is the partition function for all structures that form coaxial stacks between (i,k) and (k+1,j) for some k
    '''
    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, ligated, all_ligated, Z_BP, C_eff_basic, C_eff_no_BP_singlet, C_eff_no_coax_singlet, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv, calc_contrib ) = unpack_variables( self )
    offset = ( j - i ) % N

    #  all structures that form coaxial stacks between (i,k) and (k+1,j) for some k
    #
    #       -- k - k+1 -
    #      /   :    :   \
    #      \   :    :   /
    #       -- i    j --
    #
    for k in range( i+1, i+offset-1 ):
        if ligated[ k % N ]: Z_coax[i][j]  += Z_BP[i][k % N] * Z_BP[(k+1) % N][j] * K_coax

##################################################################################################
def update_C_eff_basic( self, i, j ):
    '''
    C_eff tracks the effective molarity of a loop starting at i and ending at j
    Assumes a model where each additional element multiplicatively reduces the effective molarity, by
      the variables l, l_BP, C_eff_stacked_pair, K_coax, etc.
    Relies on previous Z_BP, C_eff_basic, C_eff_no_BP_singlet, C_eff_no_coax_singlet, C_eff, Z_linear available for subfragments.
    Relies on Z_BP being already filled out for i,j
    TODO: In near future, will include possibility of multiple C_eff terms, which combined together will
      allow for free energy costs of loop closure to scale approximately log-linearly rather than
      linearly with loop size.
    '''
    offset = ( j - i ) % self.N

    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, ligated, all_ligated, Z_BP, C_eff_basic, C_eff_no_BP_singlet, C_eff_no_coax_singlet, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv, calc_contrib ) = unpack_variables( self )

    exclude_strained_3WJ = (not allow_strained_3WJ) and (offset == N-1) and (ligated[j] )

    # j is not base paired or coaxially stacked: Extension by one residue from j-1 to j.
    #
    #    i ~~~~~~ j-1 - j
    #
    if ligated[(j-1) % N]: C_eff_basic[i][j]  += C_eff[i][(j-1) % N] * l

    # j is base paired, and its partner is k > i. (look below for case with i and j base paired)
    #                 ___
    #                /   \
    #    i ~~~~k-1 - k...j
    #
    C_eff_for_BP = self.C_eff_no_coax_singlet if exclude_strained_3WJ else C_eff
    for k in range( i+1, i+offset):
        if ligated[ (k-1) % N]: C_eff_basic[i][j]  += C_eff_for_BP[i][(k-1) % N] * l * Z_BP[k % N][j] * l_BP

    # j is coax-stacked, and its partner is k > i.  (look below for case with i and j coaxially stacked)
    #               _______
    #              / :   : \
    #              \ :   : /
    #    i ~~~~k-1 - k   j
    #
    C_eff_for_coax = self.C_eff_no_BP_singlet if exclude_strained_3WJ else C_eff
    for k in range( i+1, i+offset):
        if ligated[ (k-1) % N]: C_eff_basic[i][j]  += C_eff_for_coax[i][(k-1) % N] * Z_coax[k % N][j] * l * l_coax

##################################################################################################
def update_C_eff_no_coax_singlet( self, i, j ):
    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, ligated, all_ligated, Z_BP, C_eff_basic, C_eff_no_BP_singlet, C_eff_no_coax_singlet, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv, calc_contrib ) = unpack_variables( self )

    # some helper arrays that prevent closure of any 3WJ with a single coaxial stack and single helix with not intervening loop nucleotides
    C_eff_no_coax_singlet[i][j] += C_eff_basic[i][j]
    C_eff_no_coax_singlet[i][j] += C_init *  Z_BP[i][j] * l_BP

##################################################################################################
def update_C_eff_no_BP_singlet( self, i, j ):
    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, ligated, all_ligated, Z_BP, C_eff_basic, C_eff_no_BP_singlet, C_eff_no_coax_singlet, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv, calc_contrib ) = unpack_variables( self )

    C_eff_no_BP_singlet[i][j] += C_eff_basic[i][j]
    C_eff_no_BP_singlet[i][j] += C_init *  Z_coax[i][j] * l_coax

##################################################################################################
def update_C_eff( self, i, j ):
    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, ligated, all_ligated, Z_BP, C_eff_basic, C_eff_no_BP_singlet, C_eff_no_coax_singlet, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv, calc_contrib ) = unpack_variables( self )

    C_eff[i][j] += C_eff_basic[i][j]

    # j is base paired, and its partner is i
    #      ___
    #     /   \
    #  i+1 ... j-1
    #    |     |
    #    i ... j
    #
    C_eff[i][j] += C_init * Z_BP[i][j] * l_BP

    # j is coax-stacked, and its partner is i.
    #       ------------
    #      /   :    :   \
    #      \   :    :   /
    #       -- i    j --
    #
    C_eff[i][j] += C_init * Z_coax[i][j] * l_coax

##################################################################################################
def update_Z_linear( self, i, j ):
    '''
    Z_linear tracks the total partition function from i to j, assuming all intervening residues are covalently connected (or base-paired).
    Relies on previous Z_BP, C_eff_basic, C_eff_no_BP_singlet, C_eff_no_coax_singlet, C_eff, Z_linear available for subfragments.
    Relies on Z_BP being already filled out for i,j
    '''
    offset = ( j - i ) % self.N

    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, ligated, all_ligated, Z_BP, C_eff_basic, C_eff_no_BP_singlet, C_eff_no_coax_singlet, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv, calc_contrib ) = unpack_variables( self )

    # j is not base paired: Extension by one residue from j-1 to j.
    #
    #    i ~~~~~~ j-1 - j
    #
    if ligated[(j-1) % N]: Z_linear[i][j] += Z_linear[i][(j - 1) % N]

    # j is base paired, and its partner is i
    #     ___
    #    /   \
    #    i...j
    #
    Z_linear[i][j] += Z_BP[i][j]

    # j is base paired, and its partner is k > i
    #                 ___
    #                /   \
    #    i ~~~~k-1 - k...j
    #
    for k in range( i+1, i+offset):
        if ligated[ (k-1) % N]: Z_linear[i][j] += Z_linear[i][(k-1) % N] * Z_BP[k % N][j]

    # j is coax-stacked, and its partner is i.
    #       ------------
    #      /   :    :   \
    #      \   :    :   /
    #       -- i    j --
    #
    Z_linear[i][j] += Z_coax[i][j]

    # j is coax-stacked, and its partner is k > i.
    #
    #               _______
    #              / :   : \
    #              \ :   : /
    #    i ~~~~k-1 - k   j
    #
    for k in range( i+1, i+offset):
        if ligated[ (k-1) % N]: Z_linear[i][j] += Z_linear[i][(k-1) % N] * Z_coax[k % N][j]

##################################################################################################
def update_Z_final( self, i ):
    # Z_final is total partition function, and is computed at end of filling dynamic programming arrays
    # Get the answer (in N ways!) --> so final output is actually Z_final(i), an array.
    # Equality of the array is tested in run_cross_checks()
    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, ligated, all_ligated, Z_BP, C_eff_basic, C_eff_no_BP_singlet, C_eff_no_coax_singlet, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv, calc_contrib ) = unpack_variables( self )

    Z_final = self.Z_final
    if not ligated[(i - 1) % N]:
        #
        #      i ------- i-1
        #
        #     or equivalently
        #        ________
        #       /        \
        #       \        /
        #        i-1    i
        #
        Z_final[i] += Z_linear[i][(i-1) % N]
    else:
        # Need to 'ligate' across i-1 to i
        # Scaling Z_final by Kd_lig/C_std to match previous literature conventions

        # Need to remove Z_coax contribution from C_eff, since its covered by C_eff_stacked_pair below.
        Z_final[i] += self.C_eff_no_coax_singlet[i][(i - 1) % N] * l / C_std

        #any split segments, combined independently
        #
        #   c+1 --- i-1 - i --- c
        #               *
        for c in range( i, i + N - 1):
            if not ligated[c % N]: Z_final[i] += Z_linear[i][c % N] * Z_linear[(c+1) % N][(i-1) % N ]

        # base pair forms a stacked pair with previous pair
        #
        #   - j+1 - j -
        #      :    :
        #      :    :
        #   - i-1 - i -
        #         *
        for j in range( i+1, (i + N - 1) ):
            if ligated[ j % N ]: Z_final[i] += C_eff_stacked_pair * Z_BP[i % N][j % N] * Z_BP[(j+1) % N][(i - 1) % N]

        C_eff_for_coax = C_eff if allow_strained_3WJ else self.C_eff_no_BP_singlet

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
                if not ligated[j % N]: continue
                if not ligated[(k-1) % N]: continue
                Z_final[i] += Z_BP[i][j % N] * C_eff_for_coax[(j+1) % N][(k-1) % N] * Z_BP[k % N][(i-1) % N] * l * l * l_coax * K_coax

            # If the two stacked base pairs are in split segments
            #
            #      \    /
            #   -- k    j --
            #  /   :    :   \
            #  \   :    :   /
            #   - i-1 - i --
            #         *
            for k in range( j + 1, i + N - 1):
                Z_final[i] += Z_BP[i][j % N] * Z_cut[j % N][k % N] * Z_BP[k % N][(i-1) % N] * K_coax


##################################################################################################
def unpack_variables( self ):
    '''
    This helper function just lets me write out equations without
    using "self" which obscures connection to my handwritten equations
    In C++, will just use convention of object variables like N_, sequence_.
    '''
    return self.params.get_variables() + \
           ( self.N, self.sequence, self.ligated, self.all_ligated,  \
             self.Z_BP,self.C_eff_basic,self.C_eff_no_BP_singlet,self.C_eff_no_coax_singlet,self.C_eff,\
             self.Z_linear,self.Z_cut,self.Z_coax,\
             self.options.calc_deriv, self.options.calc_contrib )

