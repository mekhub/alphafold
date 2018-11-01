#!/usr/bin/python
import argparse
import random
from util import *
from copy import deepcopy
from secstruct import *
from recursions import *
#from explicit_recursions import *
from dynamic_programming import DynamicProgrammingData

C_init = 1.0
l    = 0.5
l_BP = 0.1
Kd_BP = 0.001;
C_std = 1.0; # 1 M
min_loop_length = 1

def partition( sequences, circle = False ):
    if isinstance( sequences, str ): sequence = sequences
    else:
        sequence = ''
        for i in range( len( sequences ) ): sequence += sequences[i]
    N = len( sequence )
    # Could also use numpy arrays, but
    # eventually I'd like to use linked lists to
    # simplify backtracking.
    C_eff = DynamicProgrammingMatrix( N );
    Z_BP  = DynamicProgrammingMatrix( N );
    Z_linear = DynamicProgrammingMatrix( N );
    set_ids( C_eff )
    set_ids( Z_BP )
    set_ids( Z_linear )


    # initialize
    for i in range( N ): #length of fragment
        C_eff[ i ][ i ].Q = C_init
        Z_linear[ i ][ i ].Q = 1.0

    is_cutpoint = [False]*N
    if isinstance( sequences, list ):
        L = 0
        for i in range( len(sequences)-1 ):
            L = L + len( sequences[i] )
            is_cutpoint[ L-1 ] = True
    if not circle: is_cutpoint[ N-1 ] = True

    any_cutpoint = initialize_zero_matrix( N )
    for i in range( N ): #index of subfragment
        found_cutpoint = False
        any_cutpoint[ i ][ i ] = False
        for offset in range( N ): #length of subfragment
            j = (i + offset) % N;  # N cyclizes
            any_cutpoint[ i ][ j ] = found_cutpoint
            if is_cutpoint[ j ]: found_cutpoint = True

    # do the dynamic programming
    for offset in range( 1, N ): #length of subfragment
        for i in range( N ): #index of subfragment
            j = (i + offset) % N;  # N cyclizes
            update_Z( i, j, N,
                      sequence, is_cutpoint, any_cutpoint, \
                      C_init, l, l_BP, Kd_BP, C_std, min_loop_length, \
                      Z_BP, C_eff, Z_linear )

    # get the answer (in N ways!)
    Z_final = get_Z_final(N,
                          sequence, is_cutpoint, any_cutpoint, \
                          C_init, l, l_BP, Kd_BP, C_std, min_loop_length, \
                          Z_BP, C_eff, Z_linear )

    # base pair probability matrix
    bpp = initialize_zero_matrix( N );
    for i in range( N ):
        for j in range( N ):
            bpp[i][j] =  Z_BP[i][j].Q * Z_BP[j][i].Q * Kd_BP  / Z_final[0].Q

    output_DP( "Z_BP", Z_BP )
    output_DP( "C_eff", C_eff, Z_final )
    output_DP( "Z_linear", Z_linear )
    output_square( "BPP", bpp );

    def get_random_contrib( contribs ):
        # Random sample weighted by probability. Must be a simple function for this.
        contrib_cumsum = [ contribs[0][0] ]
        for contrib in contribs[1:]: contrib_cumsum.append( contrib_cumsum[-1] + contrib[0] )
        r = random.random() * contrib_cumsum[ -1 ]
        for (idx,psum) in enumerate( contrib_cumsum ):
            if r < psum: return contribs[idx]

    def backtrack( contribs_input, mode = 'mfe' ):
        if len( contribs_input ) == 0: return []
        contrib_sum = sum( contrib[0] for contrib in contribs_input )
        if   mode == 'enumerative': contribs = deepcopy( contribs_input )
        elif mode == 'mfe':         contribs = [ max( contribs_input ) ]
        elif mode == 'stochastic' : contribs = [ get_random_contrib( contribs_input ) ]

        p_bps = [] # list of tuples of (p_structure, bps_structure) for each structure
        for contrib in contribs: # each option ('contribution' to this partition function of this sub-region)
            if ( contrib[0] == 0.0 ): continue
            p_contrib = contrib[0]/contrib_sum
            p_bps_contrib = [ [p_contrib,[]] ]

            for backtrack_info in contrib[1]: # each 'branch'
                ( Z_backtrack_id, i, j )  = backtrack_info
                if Z_backtrack_id == id(Z_BP):
                    backtrack_contrib = Z_BP[i%N][j%N].contrib
                    p_bps_contrib = [ [p_bp[0], p_bp[1]+[(i,j)] ] for p_bp in p_bps_contrib ]
                elif Z_backtrack_id == id(C_eff):     backtrack_contrib = C_eff[i%N][j%N].contrib
                elif Z_backtrack_id == id(Z_linear):  backtrack_contrib = Z_linear[i%N][j%N].contrib
                p_bps_component = backtrack( backtrack_contrib, mode )
                if len( p_bps_component ) == 0: continue
                # put together all branches
                p_bps_contrib_new = []
                for p_bps1 in p_bps_contrib:
                    for p_bps2 in p_bps_component:
                        p_bps_contrib_new.append( [p_bps1[0]*p_bps2[0], p_bps1[1]+p_bps2[1]] )
                p_bps_contrib = p_bps_contrib_new

            p_bps += p_bps_contrib
        return p_bps

    def mfe( Z_final_contrib ):
        p_bps = backtrack( Z_final_contrib, mode = 'mfe' )
        assert( len(p_bps) == 1 )
        return (p_bps[0][1],p_bps[0][0])

    def boltzmann_sample( Z_final_contrib ):
        p_bps = backtrack( Z_final_contrib, mode = 'stochastic' )
        assert( len(p_bps) == 1 )
        return (p_bps[0][1],p_bps[0][0])

    ######################################
    # MFE tests
    p_MFE = [0.0]*N
    bps_MFE = [[]]*N

    for i in range( N ): (bps_MFE[i], p_MFE[i] ) = mfe( Z_final[i].contrib )
    for i in range( N ): assert( abs( ( p_MFE[i] - p_MFE[0] ) / p_MFE[0] ) < 1.0e-5 )
    print
    print 'Doing backtrack to get minimum free energy structure:'
    print  secstruct(bps_MFE[0],N), "   ", p_MFE[0], "[MFE]"
    print

    ######################################
    # Stochastic backtrack tests
    N_backtrack = 10
    print 'Doing',N_backtrack,'stochastic backtracks to get Boltzmann-weighted ensemble'
    for i in range( N_backtrack ):
        (bps,p)= boltzmann_sample( Z_final[0].contrib )
        print secstruct(bps,N), "   ", p, "[stochastic]"
    print

    ######################################
    # Enumerative backtrack tests
    p_bps = backtrack( Z_final[0].contrib , 'enumerative' )
    for (p,bps) in p_bps:  print secstruct(bps,N), "   ", p, "[enumerative]"
    p_tot = sum( p_bp[0] for p_bp in p_bps )
    print 'p_tot = ',p_tot

    # stringent test that partition function is correct:
    for i in range( N ): assert( abs( ( Z_final[i].Q - Z_final[0].Q ) / Z_final[0].Q ) < 1.0e-5 )

    print 'sequence =', sequence
    cutpoint = ''
    for i in range( N ):
        if is_cutpoint[ i ]: cutpoint += 'X'
        else: cutpoint += ' '
    print 'cutpoint =', cutpoint
    print 'circle   = ', circle
    print 'Z =',Z_final[0].Q

    return ( Z_final[0].Q, bpp )

if __name__=='__main__':
    parser = argparse.ArgumentParser( description = "Compute nearest neighbor model partitition function for RNA sequence" )
    parser.add_argument( "-s","-seq","--sequences",help="RNA sequences (separate by space)",nargs='*')
    parser.add_argument("--circle", action='store_true', default=False, help='Sequence is a circle')
    args     = parser.parse_args()
    sequences = args.sequences;
    circle   = args.circle;

    def output_test( Z, Z_ref = 0, bpp = [], bpp_idx= [], bpp_expected = 0):
        print 'Z =',Z_ref,' [expected]'
        assert( abs( (Z - Z_ref)/Z_ref )  < 1e-5 )
        print
        print 'bpp[0,4] = ',bpp[ bpp_idx[0] ][ bpp_idx[1] ]
        print 'bpp[0,4] = ',bpp_expected,' [expected]'
        assert( abs( (bpp[ bpp_idx[0] ][ bpp_idx[1] ] - bpp_expected)/bpp[ bpp_idx[0] ][ bpp_idx[1] ] )  < 1e-5 )
        print

    if sequences == None: # run tests
        # test of sequences where we know the final partition function.
        sequence = 'CAAAGAA'
        (Z, bpp) = partition( sequence, circle = True )
        output_test( Z, C_init  * (l**7) * (1 + (C_init * l_BP**2) / Kd_BP ) / C_std, \
                     bpp, [0,4], (C_init * l_BP**2/ Kd_BP) / ( 1 + C_init * l_BP**2/ Kd_BP) )

        sequence = 'CAG'
        (Z, bpp) = partition( sequence )
        output_test( Z, 1 + C_init * l**2 * l_BP / Kd_BP, \
                     bpp, [0,2], (C_init * l**2 * l_BP /Kd_BP)/( 1 + C_init * l**2 * l_BP/Kd_BP ) )

        sequences = ['C','G']
        (Z, bpp) = partition( sequences ) # note that Z sums over only base pair (not dissociated strands!)
        output_test( Z, C_std / Kd_BP, \
                     bpp, [0,1], 1.0 )

        sequences = ['GC','GC']
        (Z, bpp) = partition( sequences )
        output_test( Z, (C_std/Kd_BP)*(2 + l**2 * l_BP**2 *C_init/Kd_BP ),
                     bpp, [0,3], (1 + l**2 * l_BP**2 * C_init/Kd_BP )/(2 + l**2 * l_BP**2 *C_init/Kd_BP ) )

        sequence = 'CAGGC'
        (Z, bpp) = partition( sequence ) # note that Z sums over only base pair (not dissociated strands!)
        output_test( Z, 1 + C_init * l**2 *l_BP/Kd_BP * ( 2 + l ), \
                 bpp, [0,2], C_init*l**2*l_BP/Kd_BP /(  1+C_init*l**2*l_BP/Kd_BP * ( 2 + l )) )

        sequence = 'CGACG'
        (Z, bpp) = partition( sequence ) # note that Z sums over only base pair (not dissociated strands!)
        output_test( Z, 1 + C_init*l**2*l_BP/Kd_BP +
                     C_init*l**4*l_BP/Kd_BP  +
                     C_init**2 * (l_BP**3) * l**4 /Kd_BP /Kd_BP, \
                 bpp, [0,4], ( C_init*l**4*l_BP/Kd_BP  + C_init**2 * (l_BP**3) * l**4 /Kd_BP /Kd_BP ) / ( 1 + C_init*l**2*l_BP/Kd_BP + C_init*l**4*l_BP/Kd_BP  + C_init**2 * (l_BP**3) * l**4 /Kd_BP /Kd_BP )  )



    else:
        (Z, bpp ) = partition( sequences, circle )

