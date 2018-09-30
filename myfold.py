#!/usr/bin/python
import argparse

C_init = 1
l    = 0.5
l_BP = 0.1
C_init_BP = C_init * (l_BP/l) # 0.2
Kd_BP = 0.001;
C_std = 1; # 1 M
Kd_lig = 1.0e-5 # drops out in final answer if connections/chainbreaks are predefined

def partition( sequence ):
    N = len( sequence )
    # Could also use numpy arrays, but
    # eventually I'd like to use linked lists to
    # simplify backtracking.
    C_eff = [];
    Z = [];

    # initialize
    for i in range( N ):
        C_eff.append( [] )
        Z.append( [] )
        for j in range( N ):
            C_eff[ i ].append( 0.0 )
            Z[     i ].append(  0.0 )
    for i in range( N ): #length of fragment
        C_eff[ i ][ i ] = C_init

    for offset in range( 1, N ): #length of subfragment
        for i in range( N ): #index of subfragment
            j = (i + offset) % N;  # N cyclizes

            if ( sequence[i] == 'C' and sequence[j] == 'G' ) or ( sequence[i] == 'G' and sequence[j] == 'C' ):
                Z[i][j]     = ( 1.0 / Kd_BP ) * ( C_eff[i][ (j-1) % N] * l )

            C_eff[ i ][ j ] = C_eff[i][ (j-1) % N] * l
            C_eff[ i ][ j ] += C_init_BP * Z[i][j]
            for k in range( i+1, i+offset):
                C_eff[i][j] += C_eff[i][(k - 1 ) % N] * Z[k % N][j] * l_BP

    Z_final = []
    for i in range( N ):
        Z_close = C_eff[ i ][ (i - 1) % N] * l / Kd_lig
        # Scaling Z_final by Kd_lig/C_std to match previous literature conventions
        Z_close *= Kd_lig / C_std
        Z_final.append( Z_close )


    for i in range( N ):
        for q in range( i ): print '         ',
        for j in range( N ):
            print ' %8.3f' % Z[ i ][ (i + j) % N ],
        print ' %8.3f' % Z_final[ i ],
        print

    # stringent test that partition function is correct:
    for i in range( 1, N ):
        assert( ( Z_final[i] - Z_final[0] ) / abs( Z_final[0] ) < 1.0e-5 )

    return Z_final[0]



sequence = 'CAGA'
sequence = 'CAAAGAA'
parser = argparse.ArgumentParser( description = "Compute nearest neighbor model partitition function for RNA sequence" )
parser.add_argument( "-s","-seq","--sequence", default='CAAAGAA',help="RNA sequence")
args = parser.parse_args()
sequence = args.sequence;
Z = partition( sequence )
print 'sequence =', sequence
print 'Z =',Z

if sequence == 'CAAAGAA':  # for testing
    Z_ref = C_init  * (l**7) * (1 + C_init_BP / Kd_BP ) / C_std
    print 'Z =',Z_ref,' [expected]'
    assert( abs( (Z - Z_ref)/Z_ref )  < 1e-5 )
