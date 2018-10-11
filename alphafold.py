#!/usr/bin/python
import argparse
from util import *

# Four parameter model
Kd_BP  = 0.001;
C_init = 1          # a bit like exp(a) in multiloop
l      = 0.5        # a bit like exp(b) in multiloop
l_BP   = 0.1        # a bit like exp(c) in multiloop
params_default = [ Kd_BP, C_init, l, l_BP ]
C_std = 1; # 1 M. drops out in end (up to overall scale factor).

def partition( sequences, params = params_default, verbose = False, circle = False ):

    # unwrap the parameters of the model
    Kd_BP  = params[0]
    C_init = params[1]
    l      = params[2]
    l_BP   = params[3]
    C_init_BP = C_init * (l_BP/l) # 0.2
    min_loop_length = 1

    # initialize sequence and cutpoint info
    if isinstance( sequences, str ): sequence = sequences
    else:
        sequence = ''
        for i in range( len( sequences ) ): sequence += sequences[i]
    N = len( sequence )

    is_cutpoint = [False]*N
    if isinstance( sequences, list ):
        L = 0
        for i in range( len(sequences)-1 ):
            L = L + len( sequences[i] )
            is_cutpoint[ L-1 ] = True
    if not circle: is_cutpoint[ N-1 ] = True

    # initialize dynamic programming matrices
    # Could also use numpy arrays, but
    # eventually I'd like to use linked lists to
    # simplify backtracking.
    C_eff    = initialize_zero_matrix( N );
    Z_BP     = initialize_zero_matrix( N );
    Z_linear = initialize_zero_matrix( N );

    # first calculate derivatives with respect to Kd_BP
    dC_eff    = initialize_zero_matrix( N );
    dZ_BP     = initialize_zero_matrix( N );
    dZ_linear = initialize_zero_matrix( N );

    for i in range( N ): #length of fragment
        C_eff[ i ][ i ] = C_init
        Z_linear[ i ][ i ] = 1

    any_cutpoint = initialize_any_cutpoint( is_cutpoint )

    # do the dynamic programming
    for offset in range( 1, N ): #length of subfragment
        for i in range( N ): #index of subfragment
            j = (i + offset) % N;  # N cyclizes

            if (( sequence[i] == 'C' and sequence[j] == 'G' ) or ( sequence[i] == 'G' and sequence[j] == 'C' )) and \
                  ( any_cutpoint[i][j] or ( ((j-i-1) % N)) >= min_loop_length  ) and \
                  ( any_cutpoint[j][i] or ( ((i-j-1) % N)) >= min_loop_length  ):
                if (not is_cutpoint[ i ]) and (not is_cutpoint[ (j-1) % N]):
                    Z_BP[i][j]  += (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N] * l * l)
                    dZ_BP[i][j] += (1.0/Kd_BP ) * ( dC_eff[(i+1) % N][(j-1) % N] * l * l)
                for c in range( i, i+offset ):
                    if is_cutpoint[c % N]:
                        Z_comp1 = 1
                        Z_comp2 = 1
                        dZ_comp1 = 0
                        dZ_comp2 = 0
                        if c != i :
                            Z_comp1  = Z_linear[i+1][c % N]
                            dZ_comp1 = dZ_linear[i+1][c % N]
                        if (c+1)%N != j:
                            Z_comp2 = Z_linear[(c+1) % N][j-1]
                            dZ_comp2 = dZ_linear[(c+1) % N][j-1]
                        Z_product  = Z_comp1 * Z_comp2
                        dZ_product = dZ_comp1 * Z_comp2 + Z_comp1 * dZ_comp2

                        Z_BP[i][j] += (C_std/Kd_BP) * (l/l_BP) * Z_product
                        dZ_BP[i][j] += (C_std/Kd_BP) * (l/l_BP) * dZ_product
                dZ_BP[i][j] += -(1.0/Kd_BP) * Z_BP[i][j]

            if not is_cutpoint[(j-1) % N]:
                C_eff[i][j] += C_eff[i][(j-1) % N] * l
                dC_eff[i][j] += dC_eff[i][(j-1) % N] * l

            C_eff[i][j] += C_init_BP * Z_BP[i][j]
            dC_eff[i][j] += C_init_BP * dZ_BP[i][j]

            for k in range( i+1, i+offset):
                if not is_cutpoint[ (k-1) % N]:
                    C_eff[i][j] += C_eff[i][(k-1) % N] * Z_BP[k % N][j] * l_BP
                    dC_eff[i][j] += ( dC_eff[i][(k-1) % N] * Z_BP[k % N][j] + C_eff[i][(k-1) % N] * dZ_BP[k % N][j] ) * l_BP

            if not is_cutpoint[(j-1) % N]:
                Z_linear[i][j] += Z_linear[i][(j - 1) % N]
                dZ_linear[i][j] += dZ_linear[i][(j - 1) % N]

            Z_linear[i][j]  += Z_BP[i][j]
            dZ_linear[i][j] += dZ_BP[i][j]

            for k in range( i+1, i+offset):
                if not is_cutpoint[ (k-1) % N]:
                    Z_linear[i][j]  += Z_linear[i][(k-1) % N] * Z_BP[k % N][j]
                    dZ_linear[i][j] += ( dZ_linear[i][(k-1) % N] * Z_BP[k % N][j]+ Z_linear[i][(k-1) % N] * dZ_BP[k % N][j] )

    # get the answer (in N ways!)
    Z_final  = []
    dZ_final = []
    for i in range( N ):
        Z_final.append( 0 )
        dZ_final.append( 0 )

        if not is_cutpoint[(i + N - 1) % N]:
            for c in range( i, i + N - 1):
                if is_cutpoint[c % N]:
                    #any split segments, combined independently
                    Z_final[i]  += Z_linear[i][c % N] * Z_linear[(c+1) % N][ i - 1 ]
                    dZ_final[i] += ( dZ_linear[i][c % N] * Z_linear[(c+1) % N][ i - 1 ] + Z_linear[i][c % N] * dZ_linear[(c+1) % N][ i - 1 ] )

        if is_cutpoint[(i + N - 1) % N]:
            Z_final[i]  += Z_linear[i][(i-1) % N]
            dZ_final[i] += dZ_linear[i][(i-1) % N]
        else:
            # Scaling Z_final by Kd_lig/C_std to match previous literature conventions
            Z_final[i]  += C_eff[i][(i - 1) % N] * l / C_std
            dZ_final[i] += dC_eff[i][(i - 1) % N] * l / C_std

    # base pair probability matrix
    bpp = initialize_zero_matrix( N );
    bpp_tot = 0.0
    for i in range( N ):
        for j in range( N ):
            bpp[i][j] = Z_BP[i][j] * Z_BP[j][i] * Kd_BP * (l_BP / l) / Z_final[0]
            bpp_tot += bpp[i][j]/2.0 # to avoid double counting (i,j) and (j,i)

    if verbose:
        output_DP( "Z_BP", Z_BP )
        output_DP( "C_eff", C_eff, Z_final )
        output_DP( "Z_linear", Z_linear )
        output_square( "BPP", bpp );


    # stringent test that partition function is correct:
    for i in range( N ):
        assert( abs( ( Z_final[i] - Z_final[0] ) / Z_final[0] ) < 1.0e-5 )
        assert( abs( ( dZ_final[i] - dZ_final[0] ) / dZ_final[0] ) < 1.0e-5 )

    # calculate bpp_tot = -dlog Z_final /dlog Kd_BP in two ways! wow cool
    bpp_tot_based_on_deriv = -dZ_final[0] * Kd_BP / Z_final[0]
    assert( abs( ( bpp_tot - bpp_tot_based_on_deriv )/bpp_tot ) < 1.0e-5 )


    print 'sequence =', sequence
    cutpoint = ''
    for i in range( N ):
        if is_cutpoint[ i ]: cutpoint += 'X'
        else: cutpoint += ' '
    print 'cutpoint =', cutpoint
    print 'circle   = ', circle
    print 'Z =',Z_final[0]

    return ( Z_final[0], bpp, dZ_final[0] )

if __name__=='__main__':
    parser = argparse.ArgumentParser( description = "Compute nearest neighbor model partitition function for RNA sequence" )
    parser.add_argument( "-s","-seq","--sequences",help="RNA sequences (separate by space)",nargs='*')
    parser.add_argument("--circle", action='store_true', default=False, help='Sequence is a circle')
    args     = parser.parse_args()
    sequences = args.sequences;
    circle   = args.circle;

    if sequences == None: # run tests
        # test of sequences where we know the final partition function.
        sequence = 'CAAAGAA'
        (Z, bpp, dZ) = partition( sequence, circle = True )
        output_test( Z, C_init  * (l**7) * (1 + (C_init * l_BP/l) / Kd_BP ) / C_std, \
                     bpp, [0,4], (C_init**2 * (l**3) * l_BP/ Kd_BP) / ( C_init * (l**4) + C_init**2 * (l**3) * l_BP/ Kd_BP) )

        sequence = 'CAG'
        (Z, bpp, dZ) = partition( sequence )
        output_test( Z, 1 + C_init * l**2 / Kd_BP, \
                     bpp, [0,2], (C_init * l**2/Kd_BP)/( 1 + C_init * l**2/Kd_BP ) )

        sequences = ['C','G']
        (Z, bpp, dZ) = partition( sequences ) # note that Z sums over only base pair (not dissociated strands!)
        output_test( Z, C_std * l / Kd_BP/l_BP, \
                     bpp, [0,1], 1.0 )

        sequences = ['GC','GC']
        (Z, bpp, dZ) = partition( sequences )
        output_test( Z, (C_std/Kd_BP)*(l/l_BP)*(2 + l*l_BP*C_init/Kd_BP ), \
                     bpp, [0,3], (1 + l*l_BP*C_init/Kd_BP )/(2 + l*l_BP*C_init/Kd_BP ) )

        sequence = 'CAGGC'
        (Z, bpp, dZ) = partition( sequence ) # note that Z sums over only base pair (not dissociated strands!)
        output_test( Z, 1+C_init*l**2/Kd_BP * ( 2 + l ), \
                     bpp, [0,2], C_init*l**2/Kd_BP /(  1+C_init*l**2/Kd_BP * ( 2 + l )) )

        sequence = 'CGACG'
        (Z, bpp, dZ) = partition( sequence ) # note that Z sums over only base pair (not dissociated strands!)
        output_test( Z, 1 + C_init*l**2/Kd_BP + C_init*l**4/Kd_BP  + C_init * (l_BP/l) * l**4 /Kd_BP /Kd_BP , \
                     bpp, [0,4], ( C_init*l**4/Kd_BP  + C_init * (l_BP/l) * l**4 /Kd_BP /Kd_BP ) / ( 1 + C_init*l**2/Kd_BP + C_init*l**4/Kd_BP  + C_init * (l_BP/l) * l**4 /Kd_BP /Kd_BP )  )

        #################################################
        # let's do a numerical vs. analytic deriv test
        #################################################
        params = params_default
        delta = 1.0e-12
        params[0] += delta
        (Z_perturb, bpp_perturb, dZ_perturb) = partition( sequence, params ) # note that Z sums over only base pair (not dissociated strands!)
        dZ_numerical = (Z_perturb-Z)/delta
        print "dZ_dKd (numerical) =",dZ_numerical, ";  dZ_dKd (analytic) =",dZ

    else:
        (Z, bpp ) = partition( sequences, circle )


