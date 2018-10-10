#!/usr/bin/python
import argparse
from util import *

C_init = 1
l    = 0.5
l_BP = 0.1
C_init_BP = C_init * (l_BP/l) # 0.2
Kd_BP = 0.001;
C_std = 1; # 1 M
Kd_lig = 1.0e-5 # drops out in final answer if connections/cutpoints are predefined
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
    Z_cut    = initialize_zero_matrix( N );

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
                    Z_BP[i][j] += (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N] * l * l)
                for c in range( i, i+offset ):
                    if is_cutpoint[c % N]:
                        Z_product = 1
                        if c != i : Z_product *= Z_linear[i+1][c % N]
                        if (c+1)%N != j:  Z_product *= Z_linear[(c+1) % N][j-1]
                        Z_BP[i][j] += (C_std/Kd_BP) * (l/l_BP) * Z_product

            if not is_cutpoint[(j-1) % N]: C_eff[i][j] += C_eff[i][(j-1) % N] * l
            C_eff[i][j] += C_init_BP * Z_BP[i][j]
            for k in range( i+1, i+offset):
                if not is_cutpoint[ (k-1) % N]: C_eff[i][j] += C_eff[i][(k-1) % N] * Z_BP[k % N][j] * l_BP

            if not is_cutpoint[(j-1) % N]: Z_linear[i][j] += Z_linear[i][(j - 1) % N]
            Z_linear[i][j] += Z_BP[i][j]
            for k in range( i+1, i+offset):
                if not is_cutpoint[ (k-1) % N]: Z_linear[i][j] += Z_linear[i][(k-1) % N] * Z_BP[k % N][j]

    # get the answer (in N ways!)
    Z_final = []
    for i in range( N ):
        Z_final.append( 0 )
        if not is_cutpoint[(i + N - 1) % N]:
            for c in range( i, i + N - 1):
                if is_cutpoint[c % N]: Z_final[i] += Z_linear[i][c % N] * Z_linear[(c+1) % N][ i - 1 ] #any split segments, combined independently

        if is_cutpoint[(i + N - 1) % N]:
            Z_final[i] += Z_linear[i][(i-1) % N]
        else:
            Z_close = C_eff[i][(i - 1) % N] * l / Kd_lig
            # Scaling Z_final by Kd_lig/C_std to match previous literature conventions
            Z_close *= Kd_lig / C_std
            Z_final[i] +=  Z_close

    # base pair probability matrix
    bpp = initialize_zero_matrix( N );
    for i in range( N ):
        for j in range( N ):
            bpp[i][j] = Z_BP[i][j] * Z_BP[j][i] * Kd_BP * (l_BP / l) / Z_final[0]

    output_DP( "Z_BP", Z_BP )
    output_DP( "C_eff", C_eff, Z_final )
    output_DP( "Z_linear", Z_linear )
    output_square( "BPP", bpp );

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
        output_test( Z, C_init  * (l**7) * (1 + C_init_BP / Kd_BP ) / C_std, \
                     bpp, [0,4], (C_init**2 * (l**3) * l_BP/ Kd_BP) / ( C_init * (l**4) + C_init**2 * (l**3) * l_BP/ Kd_BP) )

        sequence = 'CAG'
        (Z, bpp) = partition( sequence )
        output_test( Z, 1 + C_init * l**2 / Kd_BP, \
                     bpp, [0,2], (C_init * l**2/Kd_BP)/( 1 + C_init * l**2/Kd_BP ) )

        sequences = ['C','G']
        (Z, bpp) = partition( sequences ) # note that Z sums over only base pair (not dissociated strands!)
        output_test( Z, C_std * l / Kd_BP/l_BP, \
                     bpp, [0,1], 1.0 )

        sequences = ['GC','GC']
        (Z, bpp) = partition( sequences )
        output_test( Z, (C_std/Kd_BP)*(l/l_BP)*(2 + l*l_BP*C_init/Kd_BP ), \
                     bpp, [0,3], (1 + l*l_BP*C_init/Kd_BP )/(2 + l*l_BP*C_init/Kd_BP ) )

        sequence = 'CAGGC'
        (Z, bpp) = partition( sequence ) # note that Z sums over only base pair (not dissociated strands!)
        output_test( Z, 1+C_init*l**2/Kd_BP * ( 2 + l ), \
                     bpp, [0,2], C_init*l**2/Kd_BP /(  1+C_init*l**2/Kd_BP * ( 2 + l )) )



    else:
        (Z, bpp ) = partition( sequences, circle )


