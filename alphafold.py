#!/usr/bin/python
import argparse
import random
from util import *

C_init = 1
l    = 0.5
l_BP = 0.1
Kd_BP = 0.001;
C_std = 1; # 1 M
min_loop_length = 1

def partition( sequences, circle = False ):
    if isinstance( sequences, str ): sequence = sequences
    else:
        sequence = ''
        for i in range( len( sequences ) ): sequence += sequences[i]
    N = len( sequence )
    # Could also use numpy arrays, but
    # eventually I'd like to use linked lists to
    # simplify backtracking.
    C_eff = initialize_zero_matrix( N );
    Z_BP  = initialize_zero_matrix( N );
    Z_linear = initialize_zero_matrix( N );

    C_eff_contrib = initialize_contrib_matrix( N )
    Z_BP_contrib = initialize_contrib_matrix( N )
    Z_linear_contrib = initialize_contrib_matrix( N )

    # initialize
    for i in range( N ): #length of fragment
        C_eff[ i ][ i ] = C_init
        Z_linear[ i ][ i ] = 1

    is_cutpoint = [False]*N
    if isinstance( sequences, list ):
        L = 0
        for i in range( len(sequences)-1 ):
            L = L + len( sequences[i] )
            is_cutpoint[ L-1 ] = True
    if not circle: is_cutpoint[ N-1 ] = True

    any_cutpoint = initialize_zero_matrix( N )
    for i in range( N ): #index of subfragment
        found_cutpoint = False
        any_cutpoint[ i ][ i ] = False
        for offset in range( N ): #length of subfragment
            j = (i + offset) % N;  # N cyclizes
            any_cutpoint[ i ][ j ] = found_cutpoint
            if is_cutpoint[ j ]: found_cutpoint = True

    # do the dynamic programming
    for offset in range( 1, N ): #length of subfragment
        for i in range( N ): #index of subfragment
            j = (i + offset) % N;  # N cyclizes

            if (( sequence[i] == 'C' and sequence[j] == 'G' ) or ( sequence[i] == 'G' and sequence[j] == 'C' )) and \
                  ( any_cutpoint[i][j] or ( ((j-i-1) % N)) >= min_loop_length  ) and \
                  ( any_cutpoint[j][i] or ( ((i-j-1) % N)) >= min_loop_length  ):
                if (not is_cutpoint[ i ]) and (not is_cutpoint[ (j-1) % N]):
                    Z_BP[i][j] += (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N] * l * l * l_BP)
                    Z_BP_contrib[i][j].append( [ (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N] * l * l * l_BP), [[id(C_eff),i+1,j-1]] ] )

                for c in range( i, i+offset ):
                    if is_cutpoint[c % N]:
                        Z_product = 1
                        Z_product_contrib = []
                        if c != i :
                            Z_product *= Z_linear[i+1][c % N]
                            Z_product_contrib.append( [id(Z_linear),i+1,c] )
                        if (c+1)%N != j:
                            Z_product *= Z_linear[(c+1) % N][j-1]
                            Z_product_contrib.append( [id(Z_linear),c+1,j-1] )
                        Z_BP[i][j] += (C_std/Kd_BP) * Z_product
                        Z_BP_contrib[i][j].append( [(C_std/Kd_BP) * Z_product, Z_product_contrib] )

            if not is_cutpoint[(j-1) % N]:
                C_eff[i][j] += C_eff[i][(j-1) % N] * l
                C_eff_contrib[i][j].append( [C_eff[i][(j-1) % N] * l, [[id(C_eff), i, j-1]]] )

            C_eff[i][j] += C_init * Z_BP[i][j] * l_BP
            C_eff_contrib[i][j].append( [ C_init * Z_BP[i][j] * l_BP, [[id(Z_BP), i, j]] ] )

            for k in range( i+1, i+offset):
                if not is_cutpoint[ (k-1) % N]:
                    C_eff[i][j] += C_eff[i][(k-1) % N] * l * Z_BP[k % N][j] * l_BP
                    C_eff_contrib[i][j].append( [ C_eff[i][(k-1) % N] * l * Z_BP[k % N][j] * l_BP, \
                                                  [[ id(C_eff), i, k-1 ], [id(Z_BP), k, j ]] ] )

            if not is_cutpoint[(j-1) % N]:
                Z_linear[i][j] += Z_linear[i][(j - 1) % N]
                Z_linear_contrib[i][j].append( [ Z_linear[i][(j - 1) % N], [[id(Z_linear), i, j-1]]]  )

            Z_linear[i][j] += Z_BP[i][j]
            Z_linear_contrib[i][j].append( [ Z_BP[i][j], [ [id(Z_BP), i, j] ] ] )

            for k in range( i+1, i+offset):
                if not is_cutpoint[ (k-1) % N]:
                    Z_linear[i][j] += Z_linear[i][(k-1) % N] * Z_BP[k % N][j]
                    Z_linear_contrib[i][j].append( [ Z_linear[i][(k-1) % N] * Z_BP[k % N][j], \
                                                     [ [id(Z_linear), i, k-1], [id(Z_BP), k, j ] ] ] )

    # get the answer (in N ways!)
    Z_final = []
    Z_final_contrib = []
    for i in range( N ):
        Z_final.append( 0 )
        Z_final_contrib.append( [] )
        if not is_cutpoint[(i + N - 1) % N]:
            for c in range( i, i + N - 1):
                if is_cutpoint[c % N]:
                    Z_final[i] += Z_linear[i][c % N] * Z_linear[(c+1) % N][ i - 1 ] #any split segments, combined independently
                    Z_final_contrib[i].append( [ Z_linear[i][c % N] * Z_linear[(c+1) % N][ i - 1 ],\
                                                 [ [id(Z_linear), i, c ], [id(Z_linear), c+1, i-1] ] ] )

        if is_cutpoint[(i + N - 1) % N]:
            Z_final[i] += Z_linear[i][(i-1) % N]
            Z_final_contrib[i].append( [ Z_linear[i][(i-1) % N], [[id(Z_linear), i, i-1 ]] ] )
        else:
            # Scaling Z_final by Kd_lig/C_std to match previous literature conventions
            Z_final[i] += C_eff[i][(i - 1) % N] * l / C_std
            Z_final_contrib[i].append( [ C_eff[i][(i - 1) % N] * l / C_std, [[id(C_eff), i, i-1]] ] )

    # base pair probability matrix
    bpp = initialize_zero_matrix( N );
    for i in range( N ):
        for j in range( N ):
            bpp[i][j] = Z_BP[i][j] * Z_BP[j][i] * Kd_BP / Z_final[0]

    output_DP( "Z_BP", Z_BP )
    output_DP( "C_eff", C_eff, Z_final )
    output_DP( "Z_linear", Z_linear )
    output_square( "BPP", bpp );

    ######################################
    # Let's do MFE through backtracking
    def backtrack( contribs, bps, p = 1.0, mfe_mode = True ):
        if len( contribs ) == 0: return p
        contrib_sum = sum( contrib[0] for contrib in contribs )
        if mfe_mode:
            contrib     = max( contribs )
        else: # stochastic backtracking mode
            contrib_cumsum = [ contribs[0][0]/contrib_sum ]
            for contrib in contribs[1:]: contrib_cumsum.append( contrib_cumsum[-1] + contrib[0]/contrib_sum )
            assert( len( contrib_cumsum ) == len( contribs ) )
            r = random.random()
            for (idx,psum) in enumerate( contrib_cumsum ):
                if r < psum: break
            contrib = contribs[idx]

        p *= contrib[0]/contrib_sum
        for backtrack_info in contrib[1]:
            #print 'backtrack_info',backtrack_info, id( Z_BP), id( C_eff ), id( Z_linear)
            Z_backtrack_id = backtrack_info[0]
            i = backtrack_info[1]
            j = backtrack_info[2]
            backtrack_contrib = []
            if Z_backtrack_id == id(Z_BP):
                backtrack_contrib = Z_BP_contrib
                bps.append( (i%N,j%N) )
            elif Z_backtrack_id == id(C_eff):     backtrack_contrib = C_eff_contrib
            elif Z_backtrack_id == id(Z_linear):  backtrack_contrib = Z_linear_contrib
            p = backtrack( backtrack_contrib[i%N][j%N], bps, p, mfe_mode )
        return p

    ######################################
    # MFE tests
    p_MFE = [0.0]*N
    bps_MFE = [[]]*N
    for i in range( N ):
        bps_MFE[i] = []
        p_MFE[i] = backtrack( Z_final_contrib[i], bps_MFE[i] )
    for i in range( N ): assert( abs( ( p_MFE[i] - p_MFE[0] ) / p_MFE[0] ) < 1.0e-5 )
    print
    print  bps_MFE[0], "   ", p_MFE[0], "[MFE]"
    print

    ######################################
    # Stochastic backtrack tests
    N_backtrack = 10
    for i in range( N_backtrack ):
        bps = []
        p = backtrack( Z_final_contrib[0], bps, mfe_mode = False )
        print bps, "   ", p

    # stringent test that partition function is correct:
    for i in range( N ): assert( abs( ( Z_final[i] - Z_final[0] ) / Z_final[0] ) < 1.0e-5 )

    print 'sequence =', sequence
    cutpoint = ''
    for i in range( N ):
        if is_cutpoint[ i ]: cutpoint += 'X'
        else: cutpoint += ' '
    print 'cutpoint =', cutpoint
    print 'circle   = ', circle
    print 'Z =',Z_final[0]

    return ( Z_final[0], bpp )


if __name__=='__main__':
    parser = argparse.ArgumentParser( description = "Compute nearest neighbor model partitition function for RNA sequence" )
    parser.add_argument( "-s","-seq","--sequences",help="RNA sequences (separate by space)",nargs='*')
    parser.add_argument("--circle", action='store_true', default=False, help='Sequence is a circle')
    args     = parser.parse_args()
    sequences = args.sequences;
    circle   = args.circle;

    def output_test( Z, Z_ref = 0, bpp = [], bpp_idx= [], bpp_expected = 0):
        print 'Z =',Z_ref,' [expected]'
        assert( abs( (Z - Z_ref)/Z_ref )  < 1e-5 )
        print
        print 'bpp[0,4] = ',bpp[ bpp_idx[0] ][ bpp_idx[1] ]
        print 'bpp[0,4] = ',bpp_expected,' [expected]'
        assert( abs( (bpp[ bpp_idx[0] ][ bpp_idx[1] ] - bpp_expected)/bpp[ bpp_idx[0] ][ bpp_idx[1] ] )  < 1e-5 )
        print


    if sequences == None: # run tests
        # test of sequences where we know the final partition function.
        sequence = 'CAAAGAA'
        (Z, bpp) = partition( sequence, circle = True )
        output_test( Z, C_init  * (l**7) * (1 + (C_init * l_BP**2) / Kd_BP ) / C_std, \
                     bpp, [0,4], (C_init * l_BP**2/ Kd_BP) / ( 1 + C_init * l_BP**2/ Kd_BP) )

        sequence = 'CAG'
        (Z, bpp) = partition( sequence )
        output_test( Z, 1 + C_init * l**2 * l_BP / Kd_BP, \
                     bpp, [0,2], (C_init * l**2 * l_BP /Kd_BP)/( 1 + C_init * l**2 * l_BP/Kd_BP ) )

        sequences = ['C','G']
        (Z, bpp) = partition( sequences ) # note that Z sums over only base pair (not dissociated strands!)
        output_test( Z, C_std / Kd_BP, \
                     bpp, [0,1], 1.0 )

        sequences = ['GC','GC']
        (Z, bpp) = partition( sequences )
        output_test( Z, (C_std/Kd_BP)*(2 + l**2 * l_BP**2 *C_init/Kd_BP ),
                     bpp, [0,3], (1 + l**2 * l_BP**2 * C_init/Kd_BP )/(2 + l**2 * l_BP**2 *C_init/Kd_BP ) )

        sequence = 'CAGGC'
        (Z, bpp) = partition( sequence ) # note that Z sums over only base pair (not dissociated strands!)
        output_test( Z, 1 + C_init * l**2 *l_BP/Kd_BP * ( 2 + l ), \
                 bpp, [0,2], C_init*l**2*l_BP/Kd_BP /(  1+C_init*l**2*l_BP/Kd_BP * ( 2 + l )) )

        sequence = 'CGACG'
        (Z, bpp) = partition( sequence ) # note that Z sums over only base pair (not dissociated strands!)
        output_test( Z, 1 + C_init*l**2*l_BP/Kd_BP +
                     C_init*l**4*l_BP/Kd_BP  +
                     C_init**2 * (l_BP**3) * l**4 /Kd_BP /Kd_BP, \
                 bpp, [0,4], ( C_init*l**4*l_BP/Kd_BP  + C_init**2 * (l_BP**3) * l**4 /Kd_BP /Kd_BP ) / ( 1 + C_init*l**2*l_BP/Kd_BP + C_init*l**4*l_BP/Kd_BP  + C_init**2 * (l_BP**3) * l**4 /Kd_BP /Kd_BP )  )



    else:
        (Z, bpp ) = partition( sequences, circle )


