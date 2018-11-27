#!/usr/bin/env python
import argparse
from alphafold.partition import *
from tests_alphafold import test_alphafold

if __name__ =='__main__':

    parser = argparse.ArgumentParser( description = "Compute nearest neighbor model partitition function for RNA sequence" )
    parser.add_argument( "-s","-seq","--sequences",help="RNA sequences (separate by space)",nargs='*')
    parser.add_argument("-struct","--structure",type=str, default=None, help='force specific structure in dot-parens notation')
    parser.add_argument("--force_base_pairs",type=str, default=None, help='force base pairs (but allow any others) in dot-parens notation')
    parser.add_argument("-c","-circ","--circle", action='store_true', default=False, help='Sequence is a circle')
    parser.add_argument("-v","--verbose", action='store_true', default=False, help='output dynamic programming matrices')
    parser.add_argument("--simple", action='store_true', default=False, help='Use simple recursions (slow!)')
    parser.add_argument("--mfe", action='store_true', default=False, help='Get minimal free energy structure')
    parser.add_argument("--bpp", action='store_true', default=False, help='Get base pairing probability')
    parser.add_argument("--stochastic", type=int, default=0, help='Number of Boltzman-weighted stochastic structures to retrieve')
    parser.add_argument("--enumerate",action='store_true', default=False, help='Backtrack to get all structures and their Boltzmann weights')
    parser.add_argument("-params","--parameters",type=str, default='', help='Parameters to use [default: '']')
    parser.add_argument("--no_coax", action='store_true', default=False, help='Turn off coaxial stacking')
    parser.add_argument("--calc_deriv", action='store_true', default=False, help='Calculate derivative with respect to Kd_BP')
    args     = parser.parse_args()

    if args.sequences != None: # run tests
        p = partition( args.sequences, circle = args.circle, params = args.parameters, verbose = args.verbose, mfe = args.mfe, calc_deriv = args.calc_deriv, calc_bpp = args.bpp, n_stochastic = int(args.stochastic), do_enumeration = args.enumerate, structure = args.structure, force_base_pairs = args.force_base_pairs, no_coax = args.no_coax, use_simple_recursions = args.simple )
    else:
        test_alphafold( verbose = args.verbose, use_simple_recursions = args.simple )
