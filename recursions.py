from dynamic_programming import *

def update_Z( i, j, N, \
              sequence, is_cutpoint, any_cutpoint, \
              C_init, l, l_BP, Kd_BP, C_std, min_loop_length, \
              Z_BP, C_eff, Z_linear ):
    offset = ( j - i ) % N
    if (( sequence[i] == 'C' and sequence[j] == 'G' ) or ( sequence[i] == 'G' and sequence[j] == 'C' )) and \
                  ( any_cutpoint[i][j] or ( ((j-i-1) % N)) >= min_loop_length  ) and \
          ( any_cutpoint[j][i] or ( ((i-j-1) % N)) >= min_loop_length  ):
        if (not is_cutpoint[ i ]) and (not is_cutpoint[ (j-1) % N]):
            Z_BP[i][j].Q += (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N].Q * l * l * l_BP)
            Z_BP[i][j].contrib.append( [ (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N].Q * l * l * l_BP), [[id(C_eff),i+1,j-1]]] )

        for c in range( i, i+offset ):
            if is_cutpoint[c % N]:
                Z_product = DynamicProgrammingData()
                Z_product.Q = 1.0
                if c != i :
                    Z_product.Q *= Z_linear[i+1][c % N].Q
                    Z_product.contrib.append( [id(Z_linear),i+1,c] )
                if (c+1)%N != j:
                    Z_product.Q *= Z_linear[(c+1) % N][j-1].Q
                    Z_product.contrib.append( [id(Z_linear),c+1,j-1] )
                Z_BP[i][j].Q += (C_std/Kd_BP) * Z_product.Q
                if len( Z_product.contrib ) > 0: Z_BP[i][j].contrib.append( [(C_std/Kd_BP) * Z_product.Q, Z_product.contrib] )

    if not is_cutpoint[(j-1) % N]:
        C_eff[i][j] += C_eff[i][(j-1) % N] * l

    C_eff[i][j] += C_init * Z_BP[i][j] * l_BP

    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            C_eff[i][j].Q += C_eff[i][(k-1) % N].Q * l * Z_BP[k % N][j].Q * l_BP
            C_eff[i][j].contrib.append( [ C_eff[i][(k-1) % N].Q * l * Z_BP[k % N][j].Q * l_BP, \
                                          [[ id(C_eff), i, k-1 ], [id(Z_BP), k, j ]]] )

    if not is_cutpoint[(j-1) % N]:
        Z_linear[i][j].Q += Z_linear[i][(j - 1) % N].Q
        Z_linear[i][j].contrib.append( [ Z_linear[i][(j - 1) % N].Q, [[id(Z_linear), i, j-1]]]  )

    Z_linear[i][j].Q += Z_BP[i][j].Q
    Z_linear[i][j].contrib.append( [ Z_BP[i][j].Q, [ [id(Z_BP), i, j]]] )

    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            Z_linear[i][j].Q += Z_linear[i][(k-1) % N].Q * Z_BP[k % N][j].Q
            Z_linear[i][j].contrib.append( [ Z_linear[i][(k-1) % N].Q * Z_BP[k % N][j].Q, \
                                             [ [id(Z_linear), i, k-1], [id(Z_BP), k, j ]]] )

def get_Z_final(N,
                sequence, is_cutpoint, any_cutpoint, \
                C_init, l, l_BP, Kd_BP, C_std, min_loop_length, \
                Z_BP, C_eff, Z_linear ):
    Z_final = []
    for i in range( N ):
        Z_final.append( DynamicProgrammingData() )
        if not is_cutpoint[(i + N - 1) % N]:
            for c in range( i, i + N - 1):
                if is_cutpoint[c % N]:
                    Z_final[i].Q += Z_linear[i][c % N].Q * Z_linear[(c+1) % N][ i - 1 ].Q #any split segments, combined independently
                    Z_final[i].contrib.append( [ Z_linear[i][c % N].Q * Z_linear[(c+1) % N][ i - 1 ].Q,\
                                                [ [id(Z_linear), i, c ], [id(Z_linear), c+1, i-1] ] ] )

        if is_cutpoint[(i + N - 1) % N]:
            Z_final[i].Q += Z_linear[i][(i-1) % N].Q
            Z_final[i].contrib.append( [ Z_linear[i][(i-1) % N].Q, [[id(Z_linear), i, i-1 ]] ] )
        else:
            # Scaling Z_final by Kd_lig/C_std to match previous literature conventions
            Z_final[i].Q += C_eff[i][(i - 1) % N].Q * l / C_std
            Z_final[i].contrib.append( [ C_eff[i][(i - 1) % N].Q * l / C_std, [[id(C_eff), i, i-1]] ] )


    return Z_final
