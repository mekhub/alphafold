#!/usr/bin/python
from scipy.optimize import minimize
from alphafold.parameters import get_params
from alphafold.partition import partition
from alphafold.score_structure import score_structure
import numpy as np

def free_energy_gap( x, sequence_structure_pairs ):
    dG_gap = 0.0
    params = get_params( suppress_all_output = True )
    params.K_coax = 0.0

    q = 10.0**x
    params.C_eff_stacked_pair = q[0]
    params.initialize_C_eff_stack()
    params.C_init = q[1]
    params.base_pair_types[2].Kd_BP = q[2]
    params.base_pair_types[3].Kd_BP = q[2] # A-U
    params.base_pair_types[4].Kd_BP = q[3]
    params.base_pair_types[5].Kd_BP = q[3] # G-U
    #for bpt in params.base_pair_types: print bpt.nt1, bpt.nt2, bpt.Kd_BP, q
    #for bpt1 in params.base_pair_types:
    #    for bpt2 in params.base_pair_types:
    #        print bpt1, bpt2, params.C_eff_stack[ bpt1][bpt2]

    for sequence,structure in sequence_structure_pairs:
        p = partition( sequence, params = params, suppress_all_output = True, mfe = True )
        dG = p.dG
        dG_structure = score_structure( sequence, structure, params = params )
        dG_gap += dG_structure - dG # will be a positive number, best case zero.
        print p.struct_MFE, x, dG_gap
    print
    return dG_gap

tRNA_sequence  = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
tRNA_structure = '(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....'

# P4-P6
#sequence = 'GGAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCA'
#structure = '.....((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).))))((...((((...(((((((((...)))))))))..))))...))...'

# P4-P6 outer junction
#sequence  = ['GGAAUUGCGGGAAAGGGGUC','GGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCA']
#structure = '.....((((((...(((((())))))..)).))))((...((((...(((((((((...)))))))))..))))...))...'

# P4-P6 outer junction -- further minimized
sequence  = ['GGAAUUGCGGGAAAGG','CUAACCACGCAGCCAAGUCCUAAGUC','GAUAUGGAUGCAGUUCA']
#structure =  '.....((((((...(())..)).))))((...((((...((()))..))))...))...'

# P4-P6 outer junction -- easier to read input
P4P6_outerjunction_sequence  = 'GGAAUUGCGGGAAAGG CUAACCACGCAGCCAAGUCCUAAGUC GAUAUGGAUGCAGUUCA'
P4P6_outerjunction_structure = '.....((((((...(( ))..)).))))((...((((...((( )))..))))...))...'


sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure) ]

loss = lambda x, sequence_structure_pairs = sequence_structure_pairs : free_energy_gap( x, sequence_structure_pairs )

# simple 1-D scan
#for x0 in range(10):  loss( [x0] )

x0 = np.array( [5, 1, 3, 3] )


#result = minimize( loss, x0, method = 'Nelder-Mead' )
result = minimize( loss, x0, method = 'L-BFGS-B' )

final_loss = loss( result.x )
print result
print 'Final parameters:', result.x, 'Loss:',final_loss
