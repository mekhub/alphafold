#!/usr/bin/python
import argparse
from alphafold.output_helpers import *
from alphafold.partition import *
from alphafold.parameters import get_minimal_params

def test_alphafold( verbose = False, use_simple_recursions = False ):
    test_params = get_minimal_params()
    (C_init, l, l_BP, C_eff_stacked_pair, K_coax, l_coax, C_std, min_loop_length, allow_strained_3WJ ) = test_params.get_variables()
    Kd_BP = test_params.base_pair_types[0].Kd_BP

    # test of sequences where we know the final partition function.
    sequence = 'CNNNGNN'
    p = partition( sequence, circle = True, params = test_params, calc_deriv = True, calc_bpp = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    output_test( p.Z, C_init  * (l**7) * (1 + (C_init * l_BP**2) / Kd_BP ) / C_std, \
                 p.bpp, [0,4], (C_init * l_BP**2/ Kd_BP) / ( 1 + C_init * l_BP**2/ Kd_BP) )

    structure= '(...)..'
    p = partition( sequence, circle = True, params = test_params, calc_deriv = True, calc_bpp = True, verbose = verbose, use_simple_recursions = use_simple_recursions, structure = structure )
    output_test( p.Z, C_init  * (l**7) * (C_init * l_BP**2) / Kd_BP / C_std, \
                 p.bpp, [0,4], 1.0 )

    sequence = 'CNG'
    p = partition( sequence, params = test_params, calc_deriv = True, mfe = True, calc_bpp = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    assert( p.bps_MFE == [(0,2)] )
    output_test( p.Z, 1 + C_init * l**2 * l_BP/ Kd_BP, \
                 p.bpp, [0,2], (C_init * l**2 * l_BP/Kd_BP)/( 1 + C_init * l**2 * l_BP/Kd_BP ) )

    sequences = ['C','G']
    p = partition( sequences, params = test_params, calc_deriv = True, calc_bpp = True, verbose = verbose, use_simple_recursions = use_simple_recursions ) # note that Z sums over only base pair (not dissociated strands!)
    output_test( p.Z, C_std/ Kd_BP, \
                 p.bpp, [0,1], 1.0 )

    sequences = ['GC','GC']
    p = partition( sequences, params = test_params, calc_deriv = True, calc_bpp = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    output_test( p.Z, (C_std/Kd_BP)*(2 + l**2 * l_BP**2 *C_init/Kd_BP + C_eff_stacked_pair/Kd_BP ), \
                 p.bpp, [0,3], (1 + l**2 * l_BP**2 * C_init/Kd_BP + C_eff_stacked_pair/Kd_BP )/(2 + l**2 * l_BP**2 *C_init/Kd_BP + C_eff_stacked_pair/Kd_BP ) )

    sequence = 'CNGGC'
    p = partition( sequence, params = test_params, calc_deriv = True, calc_bpp = True, verbose = verbose,  use_simple_recursions = use_simple_recursions )
    output_test( p.Z, 1 + C_init * l**2 *l_BP/Kd_BP * ( 2 + l ), \
                 p.bpp, [0,2], C_init*l**2*l_BP/Kd_BP /(  1+C_init*l**2*l_BP/Kd_BP * ( 2 + l )) )

    structure= '(..).'
    p = partition( sequence, params = test_params, structure = structure, calc_deriv = True, calc_bpp = True, verbose = verbose,  use_simple_recursions = use_simple_recursions )
    output_test( p.Z,  C_init * l**2 *l_BP/Kd_BP * l, \
                 p.bpp, [0,2], 0.0 )

    sequence = 'CGNCG'
    p = partition( sequence, params = test_params, calc_deriv = True, calc_bpp = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    output_test( p.Z, 1 + C_init*l**2*l_BP/Kd_BP +
                 C_init*l**4*l_BP/Kd_BP  +
                 C_init**2 * (l_BP**3) * l**4 /Kd_BP /Kd_BP +
                 C_init * l_BP * l**2 * C_eff_stacked_pair/Kd_BP /Kd_BP , \
                 p.bpp, [0,4], ( C_init*l**4*l_BP/Kd_BP  + C_init**2 * (l_BP**3) * l**4 /Kd_BP /Kd_BP  + C_init * l_BP * l**2 * C_eff_stacked_pair/Kd_BP /Kd_BP) / ( 1 + C_init*l**2*l_BP/Kd_BP + C_init*l**4*l_BP/Kd_BP  + C_init**2 * (l_BP**3) * l**4 /Kd_BP /Kd_BP + C_init * l_BP * l**2 * C_eff_stacked_pair/Kd_BP /Kd_BP )  )

    #################################################
    # let's do a numerical vs. analytic deriv test
    #################################################
    params_perturb = get_minimal_params()
    delta = 1.0e-10
    for base_pair_type in params_perturb.base_pair_types: base_pair_type.Kd_BP += delta
    p_perturb = partition( sequence, params = params_perturb ) # note that Z sums over only base pair (not dissociated strands!)
    dZ_numerical = (p_perturb.Z - p.Z)/delta
    print "dZ_dKd (numerical) =",dZ_numerical, ";  dZ_dKd (analytic) =",p.dZ
    assert( abs( dZ_numerical - p.dZ )/ abs( p.dZ ) < 1.0e-5 )
    print

    sequence = 'CNGCNG'
    p = partition( sequence, params = test_params, calc_deriv = True, calc_bpp = True, do_enumeration = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    output_test( p.Z, (1 + C_init * l**2 *l_BP/Kd_BP)**2  + C_init * l**5 * l_BP/Kd_BP + (C_init * l**2 *l_BP/Kd_BP)**2 * K_coax, \
                 p.bpp, [0,2], (C_init * l**2 *l_BP/Kd_BP*(1 + C_init * l**2 *l_BP/Kd_BP) + (C_init * l**2 *l_BP/Kd_BP)**2 * K_coax)/((1 + C_init * l**2 *l_BP/Kd_BP)**2  + C_init * l**5 * l_BP/Kd_BP + (C_init * l**2 *l_BP/Kd_BP)**2 * K_coax) )
    assert( set(p.struct_enumerate) == set(['......', '(.)...', '(....)', '...(.)', '(.)(.)']) )

    # stringent test of structure-constrained scores.
    Z_tot_ref = p.Z
    Z_enumerate = []
    structures = ['......', '(.)...', '(....)', '...(.)', '(.)(.)']
    Z_refs     = [1, C_init * l**2 *l_BP/Kd_BP, C_init * l**5 * l_BP/Kd_BP,  C_init * l**2 *l_BP/Kd_BP, (C_init * l**2 *l_BP/Kd_BP)**2 * (1+K_coax)]
    bpp_refs_0_2=[0,1,0,0,1]
    for n,structure in enumerate( structures ):
        p = partition( sequence, structure = structure, params = test_params, calc_deriv = True, calc_bpp = True, do_enumeration = False, verbose = verbose, use_simple_recursions = use_simple_recursions )
        output_test( p.Z, Z_refs[n], p.bpp, [0,2], bpp_refs_0_2[n] )
        Z_enumerate.append( p.Z )
    assert( abs( sum(Z_enumerate) - Z_tot_ref )/Z_tot_ref < 1.0e-6 )

    # testing extended alphabet & coaxial stacks
    sequence = ['xy','yz','zx']
    params_allow_strained_3WJ = get_minimal_params()
    params_allow_strained_3WJ.allow_strained_3WJ = True
    p = partition( sequence, params = params_allow_strained_3WJ, calc_deriv = True, calc_bpp = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref = 3*(C_std/Kd_BP)**2 * (1 + K_coax)  + \
            (C_std/Kd_BP)**2 * (C_init/Kd_BP) * l**3 * l_BP**3  + \
            3*(C_std/Kd_BP)**2 * (C_init/Kd_BP) * K_coax * l_coax*l**2 * l_BP
    bpp_ref = ( 2 * (C_std/Kd_BP)**2 * (1 + K_coax) + \
                (C_std/Kd_BP)**2 * (C_init/Kd_BP) * l**3 * l_BP**3 + \
                3*(C_std/Kd_BP)**2 * (C_init/Kd_BP) * K_coax * l_coax*l**2 * l_BP ) / Z_ref
    output_test( p.Z, Z_ref, p.bpp, [1,2], bpp_ref  )

    # testing extended alphabet & coaxial stacks
    sequence = ['xy','yz','zx']
    p = partition( sequence, params = test_params, calc_deriv = True, calc_bpp = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref = 3*(C_std/Kd_BP)**2 * (1 + K_coax)  + \
            (C_std/Kd_BP)**2 * (C_init/Kd_BP) * l**3 * l_BP**3
    bpp_ref = ( 2 * (C_std/Kd_BP)**2 * (1 + K_coax) + \
                (C_std/Kd_BP)**2 * (C_init/Kd_BP) * l**3 * l_BP**3 ) / Z_ref
    output_test( p.Z, Z_ref, p.bpp, [1,2], bpp_ref  )

    # test that caught a bug in Z_final
    sequence = 'NyNyxNx'
    p = partition( sequence, params = test_params, calc_deriv = True, calc_bpp = True, verbose = verbose, use_simple_recursions = use_simple_recursions )
    Z_ref = (1 + C_init * l**2 *l_BP/Kd_BP)**2  +(C_init * l**2 *l_BP/Kd_BP)**2 * K_coax
    bpp_ref = ( C_init * l**2 *l_BP/Kd_BP * (1 + C_init * l**2 *l_BP/Kd_BP)  + (C_init * l**2 *l_BP/Kd_BP)**2 * K_coax ) / Z_ref
    output_test( p.Z, Z_ref, p.bpp, [1,3], bpp_ref  )

    # test secstruct
    assert( secstruct( [(0,5),(1,4)],7 ) == '((..)).' )
    assert( bps(  '((..)).' ) == [(0,5),(1,4)] )
    assert( motifs( '.(((.)(.))).' ) == [[[1, 2], [9, 10]], [[2, 3], [5, 6], [8, 9]], [[3, 4, 5]], [[6, 7, 8]], [[10, 11, 0, 1]]] )
    assert( motifs( '(((.)(.))).'  ) == [[[0, 1], [8, 9]], [[1, 2], [4, 5], [7, 8]], [[2, 3, 4]], [[5, 6, 7]], [[9, 10, 0]]] )
    assert( motifs( '.(((.)(.)))'  ) == [[[1, 2], [9, 10]], [[2, 3], [5, 6], [8, 9]], [[3, 4, 5]], [[6, 7, 8]], [[10, 0, 1]]] )
    assert( motifs( '(((.)(.)))'   ) == [[[0, 1], [8, 9]], [[1, 2], [4, 5], [7, 8]], [[2, 3, 4]], [[5, 6, 7]], [[9, 0]]] )

if __name__=='__main__':
    parser = argparse.ArgumentParser( description = "Test nearest neighbor model partitition function for RNA sequence" )
    parser.add_argument("-v","--verbose", action='store_true', default=False, help='output dynamic programming matrices')
    parser.add_argument("--simple", action='store_true', default=False, help='Use simple recursions (fast!)')
    args     = parser.parse_args()
    test_alphafold( verbose = args.verbose, use_simple_recursions = args.simple )

