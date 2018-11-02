from dynamic_programming import *

def update_Z( i, j, N, sequence, ligated, all_ligated, \
              C_init, l, l_BP, Kd_BP, C_std, min_loop_length, \
              Z_BP, C_eff, Z_linear ):
    offset = ( j - i ) % N
    if (( sequence[i] == 'C' and sequence[j] == 'G' ) or ( sequence[i] == 'G' and sequence[j] == 'C' )) and \
        ( not all_ligated[i][j] or ( ((j-i-1) % N)) >= min_loop_length  ) and \
        ( not all_ligated[j][i] or ( ((i-j-1) % N)) >= min_loop_length  ):
        if (ligated.array[ i ]) and (ligated.array[ (j-1) % N]):
            Z_BP[i][j].Q += (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N].Q * l * l * l_BP)
            Z_BP[i][j].contribs.append( [ (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N].Q * l * l * l_BP), [[C_eff,i+1,j-1]]] )

        for c in range( i, i+offset ):
            if not ligated.array[c % N]:
                Z_product = DynamicProgrammingData()
                Z_product.Q = 1.0
                if c != i :
                    Z_product.Q *= Z_linear[i+1][c % N].Q
                    Z_product.contribs.append( [Z_linear,i+1,c] )
                if (c+1)%N != j:
                    Z_product.Q *= Z_linear[(c+1) % N][j-1].Q
                    Z_product.contribs.append( [Z_linear,c+1,j-1] )
                Z_BP[i][j].Q += (C_std/Kd_BP) * Z_product.Q
                if len( Z_product.contribs ) > 0: Z_BP[i][j].contribs.append( [(C_std/Kd_BP) * Z_product.Q, Z_product.contribs] )

    if ligated.array[(j-1) % N]:
        C_eff[i][j].Q += C_eff[i][(j-1) % N].Q * l
        C_eff[i][j].contribs.append( [C_eff[i][(j-1) % N].Q * l, [[C_eff, i, j-1]]] )

    C_eff[i][j].Q += C_init * Z_BP[i][j].Q * l_BP
    C_eff[i][j].contribs.append( [ C_init * Z_BP[i][j].Q * l_BP, [[Z_BP, i, j]]] )

    for k in range( i+1, i+offset):
        if ligated.array[ (k-1) % N]:
            C_eff[i][j].Q += C_eff[i][(k-1) % N].Q * l * Z_BP[k % N][j].Q * l_BP
            C_eff[i][j].contribs.append( [ C_eff[i][(k-1) % N].Q * l * Z_BP[k % N][j].Q * l_BP, \
                                          [[ C_eff, i, k-1 ], [Z_BP, k, j ]]] )

    if ligated.array[(j-1) % N]:
        Z_linear[i][j].Q += Z_linear[i][(j - 1) % N].Q
        Z_linear[i][j].contribs.append( [ Z_linear[i][(j - 1) % N].Q, [[Z_linear, i, j-1]]]  )

    Z_linear[i][j].Q += Z_BP[i][j].Q
    Z_linear[i][j].contribs.append( [ Z_BP[i][j].Q, [ [Z_BP, i, j]]] )

    for k in range( i+1, i+offset):
        if ligated.array[ (k-1) % N]:
            Z_linear[i][j].Q += Z_linear[i][(k-1) % N].Q * Z_BP[k % N][j].Q
            Z_linear[i][j].contribs.append( [ Z_linear[i][(k-1) % N].Q * Z_BP[k % N][j].Q, \
                                             [ [Z_linear, i, k-1], [Z_BP, k, j ]]] )

def get_Z_final(N, sequence, ligated, all_ligated, \
                C_init, l, l_BP, Kd_BP, C_std, min_loop_length, \
                Z_BP, C_eff, Z_linear ):
    Z_final = []
    for i in range( N ):
        Z_final.append( DynamicProgrammingData() )
        if ligated.array[(i + N - 1) % N]:
            for c in range( i, i + N - 1):
                if not ligated.array[c % N]:
                    Z_final[i].Q += Z_linear[i][c % N].Q * Z_linear[(c+1) % N][ i - 1 ].Q #any split segments, combined independently
                    Z_final[i].contribs.append( [ Z_linear[i][c % N].Q * Z_linear[(c+1) % N][ i - 1 ].Q,\
                                                [ [Z_linear, i, c ], [Z_linear, c+1, i-1] ] ] )

        if not ligated.array[(i + N - 1) % N]:
            Z_final[i].Q += Z_linear[i][(i-1) % N].Q
            Z_final[i].contribs.append( [ Z_linear[i][(i-1) % N].Q, [[Z_linear, i, i-1 ]] ] )
        else:
            # Scaling Z_final by Kd_lig/C_std to match previous literature conventions
            Z_final[i].Q += C_eff[i][(i - 1) % N].Q * l / C_std
            Z_final[i].contribs.append( [ C_eff[i][(i - 1) % N].Q * l / C_std, [[C_eff, i, i-1]] ] )

    return Z_final
