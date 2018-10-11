def initialize_zero_matrix( N ):
    X = []
    for i in range( N ):
        X.append( [] )
        for j in range( N ): X[i].append( 0.0 )
    return X

def output_DP( tag, X, X_final = []):
    N = len( X )
    print
    print "-----", tag, "-----"
    for i in range( N ):
        for q in range( i ): print '          ', # padding to line up
        for j in range( N ):
            print ' %9.3f' % X[i][(i+j) % N],
        if len( X_final ) > 0: print '==> %9.3f' % X_final[i],
        print

def output_square( tag, X ):
    N = len( X )
    print
    print "-----", tag, "-----"
    for i in range( N ):
        for j in range( N ):
            print ' %9.3f' % X[i][j],
        print


def initialize_any_cutpoint( is_cutpoint ):
    N = len( is_cutpoint )
    any_cutpoint = initialize_zero_matrix( N )
    for i in range( N ): #index of subfragment
        found_cutpoint = False
        any_cutpoint[ i ][ i ] = False
        for offset in range( N ): #length of subfragment
            j = (i + offset) % N;  # N cyclizes
            any_cutpoint[ i ][ j ] = found_cutpoint
            if is_cutpoint[ j ]: found_cutpoint = True
    return any_cutpoint

def output_test( Z, Z_ref = 0, bpp = [], bpp_idx= [], bpp_expected = 0):
    print 'Z =',Z_ref,' [expected]'
    assert( abs( (Z - Z_ref)/Z_ref )  < 1e-5 )
    print
    print 'bpp[%d,%d] = ' % (bpp_idx[0],bpp_idx[1]),bpp[ bpp_idx[0] ][ bpp_idx[1] ]
    print 'bpp[%d,%d] = ' % (bpp_idx[0],bpp_idx[1]),bpp_expected,' [expected]'
    assert( abs( (bpp[ bpp_idx[0] ][ bpp_idx[1] ] - bpp_expected)/bpp[ bpp_idx[0] ][ bpp_idx[1] ] )  < 1e-5 )
    print

