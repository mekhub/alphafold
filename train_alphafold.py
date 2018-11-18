#!/usr/bin/python
from scipy.optimize import minimize
from alphafold.parameters import get_params
from alphafold.partition import partition
from alphafold.score_structure import score_structure
import numpy as np

def free_energy_gap( x, sequences, structures ):
    dG_gap = 0.0
    params = get_params( suppress_all_output = True )

    (params.C_eff_stacked_pair, params.C_init ) = 10.0**x

    for sequence,structure in zip( sequences, structures ):
        p = partition( sequence, params = params, suppress_all_output = True )
        dG = p.dG
        dG_structure = score_structure( sequence, structure, params = params )
        dG_gap += dG_structure - dG # will be a positive number, best case zero.
    print x, dG_gap
    return dG_gap

sequence  = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
structure = '(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....'
sequences  = [sequence]
structures = [structure]

loss = lambda x, sequences = sequences, structures = structures : free_energy_gap( x, sequences, structures )

# simple 1-D scan
#for x0 in range(10):  loss( [x0] )

x0 = np.array( [3, 0] )
result = minimize( loss, x0, method = 'Nelder-Mead' )

final_loss = loss( result.x )
print 'Final parameters:', result.x, 'Loss:',final_loss
