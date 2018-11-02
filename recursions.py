from dynamic_programming import *

def update_Z( i, j, N, \
              sequence, is_cutpoint, any_cutpoint, \
              C_init, l, l_BP, Kd_BP, C_std, min_loop_length, \
              Z_BP, C_eff, Z_linear ):

    offset = ( j - i ) % N
    if (( sequence[i] == 'C' and sequence[j] == 'G' ) or ( sequence[i] == 'G' and sequence[j] == 'C' )) and \
       ( any_cutpoint[i][j] or ( ((j-i-1) % N)) >= min_loop_length  ) and \
       ( any_cutpoint[j][i] or ( ((i-j-1) % N)) >= min_loop_length  ):
        if (not is_cutpoint[ i ]) and (not is_cutpoint[ (j-1) % N]): Z_BP[i][j] += (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N] * l * l * l_BP)

        for c in range( i, i+offset ):
            if is_cutpoint[c % N]:
                Z_product = DynamicProgrammingData( 1.0 )
                if c != i :      Z_product *= Z_linear[i+1][c % N]
                if (c+1)%N != j: Z_product *= Z_linear[(c+1) % N][j-1]
                Z_BP[i][j] += (C_std/Kd_BP) * Z_product

    if not is_cutpoint[(j-1) % N]: C_eff[i][j] += C_eff[i][(j-1) % N] * l

    C_eff[i][j] += C_init * Z_BP[i][j] * l_BP

    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]: C_eff[i][j] += C_eff[i][(k-1) % N] * l * Z_BP[k % N][j] * l_BP

    if not is_cutpoint[(j-1) % N]: Z_linear[i][j] += Z_linear[i][(j - 1) % N]

    Z_linear[i][j] += Z_BP[i][j]

    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]: Z_linear[i][j] += Z_linear[i][(k-1) % N] * Z_BP[k % N][j]

def get_Z_final(N,
                sequence, is_cutpoint, any_cutpoint, \
                C_init, l, l_BP, Kd_BP, C_std, min_loop_length, \
                Z_BP, C_eff, Z_linear ):
    Z_final = []
    for i in range( N ):
        Z_final.append( DynamicProgrammingData() )
        if not is_cutpoint[(i + N - 1) % N]:
            for c in range( i, i + N - 1):
                if is_cutpoint[c % N]: Z_final[i] += Z_linear[i][c % N] * Z_linear[(c+1) % N][ i - 1 ] #any split segments, combined independently

        if is_cutpoint[(i + N - 1) % N]: Z_final[i] += Z_linear[i][(i-1) % N]
        else:                            Z_final[i] += C_eff[i][(i - 1) % N] * l / C_std

    return Z_final
