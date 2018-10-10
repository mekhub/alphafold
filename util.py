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
