def assert_equal( x, y, tolerance = 1.0e-5 ):
    if abs(y) > 0:
        assert( abs( x - y )/ y  < tolerance )
    else:
        assert( abs( x ) == 0 )
    return
