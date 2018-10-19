#!/usr/bin/python
from alphafold.output_helpers import *
from alphafold.partition import *

def test_alphafold():
    ( Kd_BP, C_init, l, l_BP, C_eff_stacked_pair, C_std, min_loop_length ) = AlphaFoldParams().get_variables()

    # test of sequences where we know the final partition function.
    sequence = 'CAAAGAA'
    (Z, bpp, dZ) = partition( sequence, circle = True )
    output_test( Z, C_init  * (l**7) * (1 + (C_init * l_BP**2) / Kd_BP ) / C_std, \
                 bpp, [0,4], (C_init * l_BP**2/ Kd_BP) / ( 1 + C_init * l_BP**2/ Kd_BP) )

    sequence = 'CAG'
    (Z, bpp, dZ) = partition( sequence )
    output_test( Z, 1 + C_init * l**2 * l_BP/ Kd_BP, \
                 bpp, [0,2], (C_init * l**2 * l_BP/Kd_BP)/( 1 + C_init * l**2 * l_BP/Kd_BP ) )

    sequences = ['C','G']
    (Z, bpp, dZ) = partition( sequences ) # note that Z sums over only base pair (not dissociated strands!)
    output_test( Z, C_std/ Kd_BP, \
                 bpp, [0,1], 1.0 )

    sequences = ['GC','GC']
    (Z, bpp, dZ) = partition( sequences )
    output_test( Z, (C_std/Kd_BP)*(2 + l**2 * l_BP**2 *C_init/Kd_BP + C_eff_stacked_pair/Kd_BP ), \
                 bpp, [0,3], (1 + l**2 * l_BP**2 * C_init/Kd_BP + C_eff_stacked_pair/Kd_BP )/(2 + l**2 * l_BP**2 *C_init/Kd_BP + C_eff_stacked_pair/Kd_BP ) )

    sequence = 'CAGGC'
    (Z, bpp, dZ) = partition( sequence )
    output_test( Z, 1 + C_init * l**2 *l_BP/Kd_BP * ( 2 + l ), \
                 bpp, [0,2], C_init*l**2*l_BP/Kd_BP /(  1+C_init*l**2*l_BP/Kd_BP * ( 2 + l )) )

    sequence = 'CGACG'
    (Z, bpp, dZ) = partition( sequence )
    output_test( Z, 1 + C_init*l**2*l_BP/Kd_BP +
                 C_init*l**4*l_BP/Kd_BP  +
                 C_init**2 * (l_BP**3) * l**4 /Kd_BP /Kd_BP +
                 C_init * l_BP * l**2 * C_eff_stacked_pair/Kd_BP /Kd_BP , \
                 bpp, [0,4], ( C_init*l**4*l_BP/Kd_BP  + C_init**2 * (l_BP**3) * l**4 /Kd_BP /Kd_BP  + C_init * l_BP * l**2 * C_eff_stacked_pair/Kd_BP /Kd_BP) / ( 1 + C_init*l**2*l_BP/Kd_BP + C_init*l**4*l_BP/Kd_BP  + C_init**2 * (l_BP**3) * l**4 /Kd_BP /Kd_BP + C_init * l_BP * l**2 * C_eff_stacked_pair/Kd_BP /Kd_BP )  )

    sequence = 'CAGCAG'
    (Z, bpp, dZ) = partition( sequence, verbose = True )
    output_test( Z, (1 + C_init * l**2 *l_BP/Kd_BP)**2  + C_init * l**5 * l_BP/Kd_BP, \
                 bpp, [0,2], C_init * l**2 *l_BP/Kd_BP*(1 + C_init * l**2 *l_BP/Kd_BP)/((1 + C_init * l**2 *l_BP/Kd_BP)**2  + C_init * l**5 * l_BP/Kd_BP) )

    #################################################
    # let's do a numerical vs. analytic deriv test
    #################################################
    params_perturb = AlphaFoldParams()
    delta = 1.0e-12
    params_perturb.Kd_BP += delta
    (Z_perturb, bpp_perturb, dZ_perturb) = partition( sequence, params_perturb ) # note that Z sums over only base pair (not dissociated strands!)
    dZ_numerical = (Z_perturb-Z)/delta
    print "dZ_dKd (numerical) =",dZ_numerical, ";  dZ_dKd (analytic) =",dZ
    assert( abs( dZ_numerical - dZ )/ abs( dZ ) < 1.0e-5 )

if __name__=='__main__':
    test_alphafold()

