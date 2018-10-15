#!/usr/bin/python
import argparse
from alphafold.partition import *
from tests_alphafold import test_alphafold

if __name__=='__main__':

    parser = argparse.ArgumentParser( description = "Compute nearest neighbor model partitition function for RNA sequence" )
    parser.add_argument( "-s","-seq","--sequences",help="RNA sequences (separate by space)",nargs='*')
    parser.add_argument("--circle", action='store_true', default=False, help='Sequence is a circle')
    args     = parser.parse_args()

    if args.sequences != None: # run tests
        (Z, bpp, dZ) = partition( args.sequences, circle = args.circle )
    else:
        test_alphafold()
