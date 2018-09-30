#!/usr/bin/python
import argparse

C_init = 1
l    = 0.5
l_BP = 0.1
C_init_BP = C_init * (l_BP/l) # 0.2
Kd_BP = 0.001;
C_std = 1; # 1 M
Kd_lig = 1.0e-5 # drops out in final answer if connections/chainbreaks are predefined

def initialize_zero_matrix( N ):
    X = []
    for i in range( N ):
        X.append( [] )
        for j in range( N ): X[i].append( 0.0 )
    return X

def partition( sequence, circle ):
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

    is_chainbreak = [False]*N
    if not circle: is_chainbreak[ N-1 ] = True

    # do the dynamic programming
    for offset in range( 1, N ): #length of subfragment
        for i in range( N ): #index of subfragment
            j = (i + offset) % N;  # N cyclizes

            if ( sequence[i] == 'C' and sequence[j] == 'G' ) or ( sequence[i] == 'G' and sequence[j] == 'C' ):
                if not is_chainbreak[ (j-1) % N ]: Z_BP[i][j]     += ( 1.0 / Kd_BP ) * ( C_eff[i][ (j-1) % N] * l )
                for c in range( i, i+offset ):
                    if is_chainbreak[ c % N ]: Z_BP[i][j] += (C_std/Kd_BP) * (l/l_BP) * Z_linear[i][c % N] * Z_linear[(c+1) % N][j]

            if not is_chainbreak[ (j-1) % N ]: C_eff[ i ][ j ] += C_eff[i][ (j-1) % N] * l
            C_eff[ i ][ j ] += C_init_BP * Z_BP[i][j]
            for k in range( i+1, i+offset):
                C_eff[ i ][ j ] += C_eff[i][(k - 1 ) % N] * Z_BP[k % N][j] * l_BP

            if not is_chainbreak[ (j-1) % N ]: Z_linear[ i ][ j ] += Z_linear[ i ][ (j - 1) % N]
            Z_linear[ i ][ j ] += Z_BP[ i ][ j ]
            for k in range( i+1, i+offset):
                Z_linear[i][j] += Z_linear[i][(k - 1 ) % N] * Z_BP[k % N][j]

    # get the answer (in N ways!)
    Z_final = []
    for i in range( N ):
        Z_final.append( 0 )
        for c in range( i, i + N - 1):
            #any split segments, combined independently. connection does not affect boltzman weight.
            if is_chainbreak[ c % N ]: Z_final[i] += Z_linear[i][c % N] * Z_linear[(c+1) % N][j] #any split segments, combined independently

        if is_chainbreak[ (i + N - 1) % N ]:
            Z_final[ i ] += Z_linear[ i ][ (i-1) % N]
        else:
            Z_close = C_eff[ i ][ (i - 1) % N] * l / Kd_lig
            # Scaling Z_final by Kd_lig/C_std to match previous literature conventions
            Z_close *= Kd_lig / C_std
            Z_final[ i ] +=  Z_close

    # base pair probability matrix
    bpp = initialize_zero_matrix( N );
    for i in range( N ):
        for j in range( N ):
            bpp[ i ][ j ] = Z_BP[i][j] * Z_BP[j][i] * Kd_BP * (l_BP / l) / Z_final[0]

    # output of Z_linear
    print "Z_linear"
    for i in range( N ):
        for q in range( i ): print '         ', # padding to line up
        for j in range( N ):
            print ' %8.3f' % Z_linear[ i ][ (i + j) % N ],
        print

    # output of Z_linear
    print "Z_BP"
    for i in range( N ):
        for q in range( i ): print '         ', # padding to line up
        for j in range( N ):
            print ' %8.3f' % Z_BP[ i ][ (i + j) % N ],
        print

    # output of dynamic programming matrix
    print
    print "C_eff"
    for i in range( N ):
        for q in range( i ): print '         ', # padding to line up
        for j in range( N ):
            print ' %8.3f' % C_eff[ i ][ (i + j) % N ],
        print '==> %8.3f' % Z_final[ i ]

    print
    print "BPP"
    for i in range( N ):
        for j in range( N ):
            print ' %8.3f' % bpp[ i ][ j ],
        print

    # stringent test that partition function is correct:
    for i in range( 1, N ):
        assert( ( Z_final[i] - Z_final[0] ) / abs( Z_final[0] ) < 1.0e-5 )

    return ( Z_final[0], bpp )



parser = argparse.ArgumentParser( description = "Compute nearest neighbor model partitition function for RNA sequence" )
parser.add_argument( "-s","-seq","--sequence",default='',help="RNA sequence")
parser.add_argument("--circle", action='store_true', default=False, help='Sequence is a circle')
args     = parser.parse_args()
sequence = args.sequence;
circle   = args.circle;

if sequence == '': # run tests
    # test of sequences where we know the final partition function.
    sequence = 'CAAAGAA'
    (Z, bpp) = partition( sequence, circle = True )
    Z_ref = C_init  * (l**7) * (1 + C_init_BP / Kd_BP ) / C_std
    print 'sequence =', sequence
    print 'Z =',Z
    print 'Z =',Z_ref,' [expected]'
    assert( abs( (Z - Z_ref)/Z_ref )  < 1e-5 )
    print
    bpp_expected = (C_init**2 * (l**3) * l_BP/ Kd_BP) / ( C_init * (l**4) + C_init**2 * (l**3) * l_BP/ Kd_BP)
    print 'bpp[0,4] = ',bpp[0][4]
    print 'bpp[0,4] = ',bpp_expected,' [expected]'

    # test of sequences where we know the final partition function.
    sequence = 'CAG'
    Z_ref = 1 + C_init * l**2 / Kd_BP
    print 'Z =',Z_ref,' [expected]'
    (Z, bpp) = partition( sequence, circle = False )
    print 'sequence =', sequence
    print 'Z =',Z
    print 'Z =',Z_ref,' [expected]'
    assert( abs( (Z - Z_ref)/Z_ref )  < 1e-5 )
    print
    bpp_expected = 0 ### Need to supply!
    print 'bpp[0,2] = ',bpp[0][2]
    print 'bpp[0,2] = ',bpp_expected,' [expected]'


