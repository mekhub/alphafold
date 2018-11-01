def update_Z( i, j, N, \
              sequence, is_cutpoint, any_cutpoint, \
              C_init, l, l_BP, Kd_BP, C_std, min_loop_length, \
              Z_BP, Z_BP_contrib, C_eff, C_eff_contrib, Z_linear, Z_linear_contrib ):
    offset = ( j - i ) % N
    if (( sequence[i] == 'C' and sequence[j] == 'G' ) or ( sequence[i] == 'G' and sequence[j] == 'C' )) and \
                  ( any_cutpoint[i][j] or ( ((j-i-1) % N)) >= min_loop_length  ) and \
          ( any_cutpoint[j][i] or ( ((i-j-1) % N)) >= min_loop_length  ):
        if (not is_cutpoint[ i ]) and (not is_cutpoint[ (j-1) % N]):
            Z_BP[i][j] += (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N] * l * l * l_BP)
            Z_BP_contrib[i][j].append( [ (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N] * l * l * l_BP), [[id(C_eff),i+1,j-1]] ] )

        for c in range( i, i+offset ):
            if is_cutpoint[c % N]:
                Z_product = 1
                Z_product_contrib = []
                if c != i :
                    Z_product *= Z_linear[i+1][c % N]
                    Z_product_contrib.append( [id(Z_linear),i+1,c] )
                if (c+1)%N != j:
                    Z_product *= Z_linear[(c+1) % N][j-1]
                    Z_product_contrib.append( [id(Z_linear),c+1,j-1] )
                Z_BP[i][j] += (C_std/Kd_BP) * Z_product
                if len( Z_product_contrib ) > 0: Z_BP_contrib[i][j].append( [(C_std/Kd_BP) * Z_product, Z_product_contrib] )

    if not is_cutpoint[(j-1) % N]:
        C_eff[i][j] += C_eff[i][(j-1) % N] * l
        C_eff_contrib[i][j].append( [C_eff[i][(j-1) % N] * l, [[id(C_eff), i, j-1]]] )

    C_eff[i][j] += C_init * Z_BP[i][j] * l_BP
    C_eff_contrib[i][j].append( [ C_init * Z_BP[i][j] * l_BP, [[id(Z_BP), i, j]] ] )

    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            C_eff[i][j] += C_eff[i][(k-1) % N] * l * Z_BP[k % N][j] * l_BP
            C_eff_contrib[i][j].append( [ C_eff[i][(k-1) % N] * l * Z_BP[k % N][j] * l_BP, \
                                          [[ id(C_eff), i, k-1 ], [id(Z_BP), k, j ]] ] )

    if not is_cutpoint[(j-1) % N]:
        Z_linear[i][j] += Z_linear[i][(j - 1) % N]
        Z_linear_contrib[i][j].append( [ Z_linear[i][(j - 1) % N], [[id(Z_linear), i, j-1]]]  )

    Z_linear[i][j] += Z_BP[i][j]
    Z_linear_contrib[i][j].append( [ Z_BP[i][j], [ [id(Z_BP), i, j] ] ] )

    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            Z_linear[i][j] += Z_linear[i][(k-1) % N] * Z_BP[k % N][j]
            Z_linear_contrib[i][j].append( [ Z_linear[i][(k-1) % N] * Z_BP[k % N][j], \
                                             [ [id(Z_linear), i, k-1], [id(Z_BP), k, j ] ] ] )
