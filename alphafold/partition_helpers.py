def initialize_zero_matrix( N ):
    X = []
    for i in range( N ):
        X.append( [] )
        for j in range( N ): X[i].append( 0.0 )
    return X

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

def initialize_sequence_information( sequences, circle ):
    # initialize sequence
    if isinstance( sequences, str ): sequence = sequences
    else:
        sequence = ''
        for i in range( len( sequences ) ): sequence += sequences[i]
    N = len( sequence )

    # initialize cutpoint information
    is_cutpoint = [False]*N
    if isinstance( sequences, list ):
        L = 0
        for i in range( len(sequences)-1 ):
            L = L + len( sequences[i] )
            is_cutpoint[ L-1 ] = True
    if not circle: is_cutpoint[ N-1 ] = True

    any_cutpoint = initialize_any_cutpoint( is_cutpoint )

    return ( sequence, is_cutpoint, any_cutpoint )


def initialize_dynamic_programming_matrices( N, C_init ):
    # initialize dynamic programming matrices
    C_eff    = initialize_zero_matrix( N );
    Z_BP     = initialize_zero_matrix( N );
    Z_linear = initialize_zero_matrix( N );
    for i in range( N ): #length of fragment
        C_eff[ i ][ i ] = C_init
        Z_linear[ i ][ i ] = 1

    # first calculate derivatives with respect to Kd_BP
    dC_eff    = initialize_zero_matrix( N );
    dZ_BP     = initialize_zero_matrix( N );
    dZ_linear = initialize_zero_matrix( N );

    return   ( Z_BP, C_eff, Z_linear, dZ_BP, dC_eff, dZ_linear )
