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
     sequence, is_cutpoint, any_intervening_cutpoint, Z_BP, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv ) = unpack_variables( self )
    offset = ( j - i ) % N
    for c in range( i, i+offset ):
        if is_cutpoint[c % N]:
            # strand 1  (i --> c), strand 2  (c+1 -- > j)
            Z_seg1  = Z_seg2  = 1
            if calc_deriv: dZ_seg1 = dZ_seg2 = 0
            if c != i :
                Z_seg1  = Z_linear [(i+1) % N][c % N]
                if calc_deriv: dZ_seg1 = Z_linear.dQ[(i+1) % N][c % N]
            if (c+1)%N != j:
                Z_seg2  = Z_linear [(c+1) % N][(j-1) % N]
                if calc_deriv: dZ_seg2 = Z_linear.dQ[(c+1) % N][(j-1) % N]
            Z_cut [i][j] += Z_seg1 * Z_seg2
            if calc_deriv: Z_cut.dQ[i][j] += dZ_seg1 * Z_seg2 + Z_seg1 * dZ_seg2
            #Z_cut_contrib[i][j].append( Z_linear_contrib

##################################################################################################
def update_Z_BP( self, i, j, calc_contrib = False ):
    '''
    Z_BP is the partition function for all structures that base pair i and j.
    Relies on previous Z_BP, C_eff, Z_linear available for subfragments.
    '''
    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, is_cutpoint, any_intervening_cutpoint, Z_BP, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv ) = unpack_variables( self )
    offset = ( j - i ) % N

    ( C_eff_for_coax, C_eff_for_BP ) = (C_eff, C_eff ) if allow_strained_3WJ else (self.C_eff_no_BP_singlet, self.C_eff_no_coax_singlet )
    if calc_contrib: Z_BP.contrib[i][j] = []

    # minimum loop length -- no other way to penalize short segments.
    if ( not any_intervening_cutpoint[i][j] and ( ((j-i-1) % N)) < min_loop_length ): return
    if ( not any_intervening_cutpoint[j][i] and ( ((i-j-1) % N)) < min_loop_length ): return

    for base_pair_type in self.base_pair_types:
        if base_pair_type.match_lowercase:
            if not (sequence[i].islower() and sequence[j].islower() and sequence[i]==sequence[j] ): continue
        else:
            if not ( sequence[i] == base_pair_type.nt1 and sequence[ j ] == base_pair_type.nt2 ): continue

        (Z_BPq, Kd_BPq)  = ( base_pair_type.Z_BP_DP, base_pair_type.Kd_BP )

        if (not is_cutpoint[ i ]) and (not is_cutpoint[ (j-1) % N]):
            # base pair closes a loop
            #
            #    ~~~~~~
            #   ~      ~
            # i+1      j-1
            #   \       /
            #    i ... j
            #
            Z_BPq.Q[i][j]  += (1.0/Kd_BPq ) * ( C_eff_for_BP.Q [(i+1) % N][(j-1) % N] * l * l * l_BP)
            if calc_deriv: Z_BPq.dQ[i][j] += (1.0/Kd_BPq ) * ( C_eff_for_BP.dQ[(i+1) % N][(j-1) % N] * l * l * l_BP)
            if calc_contrib: self.Z_BP.contrib[i][j].append( ( (1.0/Kd_BPq ) * ( C_eff_for_BP.Q [(i+1) % N][(j-1) % N] * l * l * l_BP), [(C_eff_for_BP, i+1, j-1)] ) )

            # base pair forms a stacked pair with previous pair
            #      ___
            #     /   \
            #  i+1 ... j-1
            #    |     |
            #    i ... j
            #
            # TODO: generalize C_eff_stacked_pair to be function of base pairs q (at i,j) and r (at i+1,j-1)
            Z_BPq.Q[i][j]  += (1.0/Kd_BPq ) * C_eff_stacked_pair * Z_BP.Q[(i+1) % N][(j-1) % N]
            if calc_deriv: Z_BPq.dQ[i][j] += (1.0/Kd_BPq ) * C_eff_stacked_pair * Z_BP.dQ[(i+1) % N][(j-1) % N]

        # base pair brings together two strands that were previously disconnected
        #
        #   \       /
        #    i ... j
        #
        Z_BPq .Q[i][j] += (C_std/Kd_BPq) * Z_cut.Q[i][j]
        if calc_deriv: Z_BPq.dQ[i][j] += (C_std/Kd_BPq) * Z_cut.dQ[i][j]

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
                    Z_BPq .Q[i][j] += Z_BP.Q[(i+1) % N][k % N] * C_eff_for_coax.Q[(k+1) % N][(j-1) % N] * l**2 * l_coax * K_coax / Kd_BPq
                    if calc_deriv: Z_BPq.dQ[i][j] += (Z_BP.dQ[(i+1) % N][k % N] * C_eff_for_coax.Q[(k+1) % N][(j-1) % N] +
                                    Z_BP.Q[(i+1) % N][k % N] * C_eff_for_coax.dQ[(k+1) % N][(j-1) % N] ) * l**2 * l_coax * K_coax / Kd_BPq

            # coaxial stack of bp (i,j) and (k,j-1)...  close loop on left, and "right stack"
            #            ___
            #           /   \
            #  ~ k-1 - k ... j-1
            # ~              |
            #  ~ i+1 - i ... j
            #
            for k in range( i+2, i+offset-1 ):
                if not is_cutpoint[(k-1) % N]:
                    Z_BPq .Q[i][j] += C_eff_for_coax.Q[(i+1) % N][(k-1) % N] * Z_BP.Q[k % N][(j-1) % N] * l**2 * l_coax * K_coax / Kd_BPq
                    if calc_deriv: Z_BPq.dQ[i][j] += (C_eff_for_coax.dQ[(i+1) % N][(k-1) % N] * Z_BP.Q[k % N][(j-1) % N] +
                                    C_eff_for_coax.Q[(i+1) % N][(k-1) % N] * Z_BP.dQ[k % N][(j-1) % N] ) * l**2 * l_coax * K_coax / Kd_BPq

        # "left stack" but no loop closed on right (free strands hanging off j end)
        #      ___
        #     /   \
        #  i+1 ... k -
        #    |
        #    i ... j -
        #
        if not is_cutpoint[ i ]:
            for k in range( i+2, i+offset ):
                Z_BPq.Q[i][j] += Z_BP.Q[(i+1) % N][k % N] * Z_cut.Q[k % N][j] * C_std * K_coax / Kd_BPq
                if calc_deriv: Z_BPq.dQ[i][j] += (Z_BP.dQ[(i+1) % N][k % N] * Z_cut.Q[k % N][j] + Z_BP.Q[(i+1) % N][k % N] * Z_cut.dQ[k % N][j] ) * C_std * K_coax / Kd_BPq

        # "right stack" but no loop closed on left (free strands hanging off i end)
        #       ___
        #      /   \
        #   - k ... j-1
        #           |
        #   - i ... j
        #
        if not is_cutpoint[(j-1) % N]:
            for k in range( i, i+offset-1 ):
                Z_BPq.Q[i][j] += Z_cut.Q[i][k % N] * Z_BP.Q[k % N][(j-1) % N] * C_std * K_coax / Kd_BPq
                if calc_deriv: Z_BPq.dQ[i][j]+= ( Z_cut.dQ[i][k % N] * Z_BP.Q[k % N][(j-1) % N]  + Z_cut.Q[i][k % N] * Z_BP.dQ[k % N][(j-1) % N] ) * C_std * K_coax / Kd_BPq

        # key 'special sauce' for derivative w.r.t. Kd_BP
        if calc_deriv: Z_BPq.dQ[i][j] += -(1.0/Kd_BPq) * Z_BPq.Q[i][j]

        Z_BP.Q[i][j]  += Z_BPq.Q[i][j]
        if calc_deriv: Z_BP.dQ[i][j] += Z_BPq.dQ[i][j]
        #Z_BP_contrib[i][j] += Z_BPq_contrib[i][j]

##################################################################################################
def update_Z_coax( self, i, j ):
    '''
    Z_coax(i,j) is the partition function for all structures that form coaxial stacks between (i,k) and (k+1,j) for some k
    '''
    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, is_cutpoint, any_intervening_cutpoint, Z_BP, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv ) = unpack_variables( self )
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
            Z_coax.Q[i][j]  += Z_BP.Q[i][k % N] * Z_BP.Q[(k+1) % N][j] * K_coax
            Z_coax.dQ[i][j] += (Z_BP.dQ[i][k % N] * Z_BP.Q[(k+1) % N][j] + Z_BP.Q[i][k % N] * Z_BP.dQ[(k+1) % N][j]) * K_coax

##################################################################################################
def update_C_eff( self, i, j, calc_contrib = False ):
    '''
    C_eff tracks the effective molarity of a loop starting at i and ending at j
    Assumes a model where each additional element multiplicatively reduces the effective molarity, by
      the variables l, l_BP, C_eff_stacked_pair, K_coax, etc.
    Relies on previous Z_BP, C_eff, Z_linear available for subfragments.
    Relies on Z_BP being already filled out for i,j
    TODO: In near future, will include possibility of multiple C_eff terms, which combined together will
      allow for free energy costs of loop closure to scale approximately log-linearly rather than
      linearly with loop size.
    '''
    offset = ( j - i ) % self.N

    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, is_cutpoint, any_intervening_cutpoint, Z_BP, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv ) = unpack_variables( self )

    if ( calc_contrib ): self.C_eff.contrib[i][j] = []

    exclude_strained_3WJ = (not allow_strained_3WJ) and (offset == N-1) and (not is_cutpoint[j] )

    # j is not base paired or coaxially stacked: Extension by one residue from j-1 to j.
    #
    #    i ~~~~~~ j-1 - j
    #
    if not is_cutpoint[(j-1) % N]:
        C_eff.Q[i][j]  += C_eff.Q[i][(j-1) % N] * l
        if calc_deriv: C_eff.dQ[i][j] += C_eff.dQ[i][(j-1) % N] * l
        if calc_contrib: self.C_eff.contrib[i][j].append( (C_eff.Q[i][(j-1) % N] * l, [(self.C_eff,i,j-1)] ) )

    # j is base paired, and its partner is k > i. (look below for case with i and j base paired)
    #                 ___
    #                /   \
    #    i ~~~~k-1 - k...j
    #
    C_eff_for_BP = self.C_eff_no_coax_singlet if exclude_strained_3WJ else C_eff
    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            C_eff.Q[i][j]  += C_eff_for_BP.Q[i][(k-1) % N] * l * Z_BP.Q[k % N][j] * l_BP
            if calc_deriv: C_eff.dQ[i][j] += ( C_eff_for_BP.dQ[i][(k-1) % N] * Z_BP.Q[k % N][j] + C_eff_for_BP.Q[i][(k-1) % N] * Z_BP.dQ[k % N][j] ) * l * l_BP

    # j is coax-stacked, and its partner is k > i.  (look below for case with i and j coaxially stacked)
    #               _______
    #              / :   : \
    #              \ :   : /
    #    i ~~~~k-1 - k   j
    #
    C_eff_for_coax = self.C_eff_no_BP_singlet if exclude_strained_3WJ else C_eff
    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            C_eff.Q[i][j]  += C_eff_for_coax.Q[i][(k-1) % N] * Z_coax.Q[k % N][j] * l * l_coax
            if calc_deriv: C_eff.dQ[i][j] += (C_eff_for_coax.dQ[i][(k-1) % N] * Z_coax.Q[k % N][j] + C_eff_for_coax.Q[i][(k-1) % N] * Z_coax.dQ[k % N][j]) * l * l_coax

    # some helper arrays that prevent closure of any 3WJ with a single coaxial stack and single helix with not intervening loop nucleotides
    self.C_eff_no_coax_singlet.Q[i][j] =  C_eff.Q[i][j]  + C_init *  Z_BP.Q[i][j] * l_BP
    self.C_eff_no_coax_singlet.dQ[i][j] = C_eff.dQ[i][j] + C_init * Z_BP.dQ[i][j] * l_BP
    if calc_contrib: self.C_eff_no_coax_singlet.contrib[i][j] = self.C_eff.contrib[i][j] + [ (C_init * Z_BP.Q[i][j] * l_BP, [(self.Z_BP,i,j)] )]

    self.C_eff_no_BP_singlet.Q[i][j] =  C_eff.Q[i][j] + C_init *  Z_coax.Q[i][j] * l_coax
    self.C_eff_no_BP_singlet.dQ[i][j] = C_eff.dQ[i][j] + C_init * Z_coax.dQ[i][j] * l_coax


    # j is base paired, and its partner is i
    #      ___
    #     /   \
    #  i+1 ... j-1
    #    |     |
    #    i ... j
    #
    C_eff.Q[i][j]  += C_init * Z_BP.Q[i][j] * l_BP
    if calc_deriv: C_eff.dQ[i][j] += C_init * Z_BP.dQ[i][j] * l_BP
    #C_eff_contrib[i][j].append( [[Z_BP, i, j, C_init * l_BP]] )

    # j is coax-stacked, and its partner is i.
    #       ------------
    #      /   :    :   \
    #      \   :    :   /
    #       -- i    j --
    #
    C_eff.Q[i][j]  += C_init * Z_coax.Q[i][j] * l_coax
    if calc_deriv: C_eff.dQ[i][j] += C_init * Z_coax.dQ[i][j] * l_coax

##################################################################################################
def update_Z_linear( self, i, j, calc_contrib = False ):
    '''
    Z_linear tracks the total partition function from i to j, assuming all intervening residues are covalently connected (or base-paired).
    Relies on previous Z_BP, C_eff, Z_linear available for subfragments.
    Relies on Z_BP being already filled out for i,j
    '''
    offset = ( j - i ) % self.N

    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, is_cutpoint, any_intervening_cutpoint, Z_BP, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv ) = unpack_variables( self )

    if calc_contrib: self.Z_linear.contrib[i][j] = []

    # j is not base paired: Extension by one residue from j-1 to j.
    #
    #    i ~~~~~~ j-1 - j
    #
    if not is_cutpoint[(j-1) % N]:
        Z_linear.Q[i][j]  += Z_linear.Q[i][(j - 1) % N]
        if calc_deriv: Z_linear.dQ[i][j] += Z_linear.dQ[i][(j - 1) % N]

    # j is base paired, and its partner is i
    #     ___
    #    /   \
    #    i...j
    #
    Z_linear.Q[i][j]  += Z_BP.Q[i][j]
    if calc_deriv: Z_linear.dQ[i][j] += Z_BP.dQ[i][j]
    if calc_contrib: self.Z_linear.contrib[i][j].append( (Z_BP.Q[i][j], [(self.Z_BP,i,j)]) )

    # j is base paired, and its partner is k > i
    #                 ___
    #                /   \
    #    i ~~~~k-1 - k...j
    #
    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            Z_linear.Q[i][j]  += Z_linear.Q[i][(k-1) % N] * Z_BP.Q[k % N][j]
            if calc_deriv: Z_linear.dQ[i][j] += ( Z_linear.dQ[i][(k-1) % N] * Z_BP.Q[k % N][j] + Z_linear.Q[i][(k-1) % N] * Z_BP.dQ[k % N][j] )

    # j is coax-stacked, and its partner is i.
    #       ------------
    #      /   :    :   \
    #      \   :    :   /
    #       -- i    j --
    #
    Z_linear.Q[i][j]  += Z_coax.Q[i][j]
    if calc_deriv: Z_linear.dQ[i][j] += Z_coax.dQ[i][j]

    # j is coax-stacked, and its partner is k > i.
    #
    #               _______
    #              / :   : \
    #              \ :   : /
    #    i ~~~~k-1 - k   j
    #
    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            Z_linear.Q[i][j]  += Z_linear.Q[i][(k-1) % N] * Z_coax.Q[k % N][j]
            if calc_deriv: Z_linear.dQ[i][j] += Z_linear.dQ[i][(k-1) % N] * Z_coax.Q[k % N][j] + Z_linear.Q[i][(k-1) % N] * Z_coax.dQ[k % N][j]

##################################################################################################
def get_Z_final( self, calc_contrib = False ):
    # Z_final is total partition function, and is computed at end of filling dynamic programming arrays
    # Get the answer (in N ways!) --> so final output is actually Z_final(i), an array.
    # Equality of the array is tested in run_cross_checks()
    (C_init, l, Kd_BP, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ, N, \
     sequence, is_cutpoint, any_intervening_cutpoint, Z_BP, C_eff, Z_linear, Z_cut, Z_coax, calc_deriv ) = unpack_variables( self )

    Z_final =[0.0]*N
    dZ_final =[0.0]*N
    Z_final_contrib = []

    for i in range( N ):
        Z_final.append( 0.0 )
        dZ_final.append( 0.0 )
        Z_final_contrib.append( [] )

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
            Z_final[i]  += Z_linear.Q[i][(i-1) % N]
            if calc_deriv: dZ_final[i] += Z_linear.dQ[i][(i-1) % N]
            if calc_contrib: Z_final_contrib[i].append( (Z_linear.Q[i][(i-1) % N], [(self.Z_linear, i, i-1)]) )
        else:
            # Need to 'ligate' across i-1 to i
            # Scaling Z_final by Kd_lig/C_std to match previous literature conventions

            # Need to remove Z_coax contribution from C_eff, since its covered by C_eff_stacked_pair below.
            Z_final[i]  += self.C_eff_no_coax_singlet.Q[i][(i - 1) % N] * l / C_std
            if calc_deriv: dZ_final[i] += self.C_eff_no_coax_singlet.dQ[i][(i - 1) % N] * l / C_std

            for c in range( i, i + N - 1):
                if self.is_cutpoint[c % N]:
                    #any split segments, combined independently
                    #
                    #   c+1 --- i-1 - i --- c
                    #               *
                    Z_final[i]  += Z_linear.Q[i][c % N] * Z_linear.Q[(c+1) % N][(i-1) % N ]
                    if calc_deriv: dZ_final[i] += ( Z_linear.dQ[i][c % N] * Z_linear.Q[(c+1) % N][(i-1) % N ] + Z_linear.Q[i][c % N] * Z_linear.dQ[(c+1) % N][(i-1) % N ] )

            # base pair forms a stacked pair with previous pair
            #
            #   - j+1 - j -
            #      :    :
            #      :    :
            #   - i-1 - i -
            #         *
            for j in range( i+1, (i + N - 1) ):
                if not is_cutpoint[ j % N ]:
                    Z_final[i]  += C_eff_stacked_pair * Z_BP.Q[i % N][j % N] * Z_BP.Q[(j+1) % N][(i - 1) % N]
                    if calc_deriv: dZ_final[i] += C_eff_stacked_pair * (Z_BP.dQ[i % N][j % N] * Z_BP.Q[(j+1) % N][(i - 1) % N] + Z_BP.Q[i % N][j % N] * Z_BP.dQ[(j+1) % N][(i - 1) % N] )

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
                    if is_cutpoint[j % N]: continue
                    if is_cutpoint[(k-1) % N]: continue
                    Z_final[i]  += Z_BP.Q[i][j % N] * C_eff_for_coax.Q[(j+1) % N][(k-1) % N] * Z_BP.Q[k % N][(i-1) % N] * l * l * l_coax * K_coax
                    if calc_deriv: dZ_final[i] += (Z_BP.dQ[i][j % N] *  C_eff_for_coax.Q[(j+1) % N][(k-1) % N] *  Z_BP.Q[k % N][(i-1) % N] +
                                                   Z_BP.Q[i][j % N] * C_eff_for_coax.dQ[(j+1) % N][(k-1) % N] *  Z_BP.Q[k % N][(i-1) % N] +
                                                   Z_BP.Q[i][j % N] *  C_eff_for_coax.Q[(j+1) % N][(k-1) % N] * Z_BP.dQ[k % N][(i-1) % N]) * l * l * l_coax * K_coax

                # If the two stacked base pairs are in split segments
                #
                #      \    /
                #   -- k    j --
                #  /   :    :   \
                #  \   :    :   /
                #   - i-1 - i --
                #         *
                for k in range( j + 1, i + N - 1):
                    Z_final[i]  += Z_BP.Q[i][j % N] * Z_cut.Q[j % N][k % N] * Z_BP.Q[k % N][(i-1) % N] * K_coax
                    if calc_deriv: dZ_final[i] += (Z_BP.dQ[i][j % N] * Z_cut.Q[j % N][k % N] * Z_BP.Q[k % N][(i-1) % N] +
                                                   Z_BP.Q[i][j % N] * Z_cut.dQ[j % N][k % N] * Z_BP.Q[k % N][(i-1) % N] +
                                                   Z_BP.Q[i][j % N] * Z_cut.Q[j % N][k % N] * Z_BP.dQ[k % N][(i-1) % N]) * K_coax

    self.Z_final = Z_final
    self.dZ_final = dZ_final
    self.Z_final_contrib = Z_final_contrib

##################################################################################################
def unpack_variables( self ):
    '''
    This helper function just lets me write out equations without
    using "self" which obscures connection to my handwritten equations
    In C++, will just use convention of object variables like N_, sequence_.
    '''
    return self.params.get_variables() + \
           ( self.N, self.sequence, self.is_cutpoint, self.any_intervening_cutpoint,  \
             self.Z_BP,self.C_eff,self.Z_linear,self.Z_cut,self.Z_coax,\
             self.calc_deriv)

