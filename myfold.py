#!/usr/bin/python


sequence = 'CAAAGAA'

C_init = 1;
l    = 0.5
l_BP = 0.1
N = len( sequence )
Kd = 0.001;

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

for offset in range( 1, N+1 ): #length of subfragment
    for i in range( N ): #index of subfragment
        j = (i + offset) % N;  # N cyclizes

        if ( sequence[i] == 'C' and sequence[j] == 'G' ) or ( sequence[i] == 'G' and sequence[j] == 'C' ):
            Z[i][j]     = ( 1.0 / Kd ) * ( C_eff[i][ (j-1) % N] * l )
            print "FOUND BP"

        C_eff[ i ][ j ] = C_eff[i][ (j-1) % N] * l
        C_eff[ i ][ j ] += C_init * Z[i][j] * l_BP/l
        for k in range( i+1, j):
            C_eff[i][j] += C_eff[i][(k - 1 ) % N] * Z[k % N][j] * l_BP

for i in range( N ):
    for j in range( N ):
        print ' %8.3f' % C_eff[ i ][ j ],
    print

print
for i in range( N ):
    for j in range( N ):
        print ' %8.3f' % Z[ i ][ j ],
    print
