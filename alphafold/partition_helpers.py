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


def  update_Z_BP( idx, sequence, C_std, l, l_BP, Kd_BP, is_cutpoint, any_cutpoint, min_loop_length, Z_BP, C_eff, Z_linear, dZ_BP, dC_eff, dZ_linear ):
    i = idx[0]
    j = idx[1]
    N = len( Z_BP )
    offset = ( j - i ) % N

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

                Z_BP[i][j]  += (C_std/Kd_BP) * (l/l_BP) * Z_product
                dZ_BP[i][j] += (C_std/Kd_BP) * (l/l_BP) * dZ_product

        # key 'special sauce' for derivative w.r.t. Kd_BP
        dZ_BP[i][j] += -(1.0/Kd_BP) * Z_BP[i][j]
