#!/usr/bin/python
from scipy.optimize import minimize
from alphafold.parameters import get_params
from alphafold.partition import partition
from alphafold.score_structure import score_structure
import numpy as np

def free_energy_gap( x, sequence_structure_pairs, apply_params ):
    dG_gap = 0.0
    params = get_params( suppress_all_output = True )
    apply_params( params, x )
    print x
    for sequence_structure_pair in sequence_structure_pairs:
        sequence, structure = sequence_structure_pair[0:2]
        force_base_pairs = None
        if len( sequence_structure_pair ) > 2: force_base_pairs = sequence_structure_pair[2]
        p = partition( sequence, params = params, suppress_all_output = True, mfe = True )
        dG = p.dG
        dG_structure = score_structure( sequence, structure, params = params )
        dG_gap += dG_structure - dG # will be a positive number, best case zero.
        print p.struct_MFE, dG_gap
    print
    return dG_gap

tRNA_sequence  = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
tRNA_structure = '(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....'

# P4-P6
P4P6_sequence = 'GGAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCA'
P4P6_structure = '.....((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).))))((...((((...(((((((((...)))))))))..))))...))...'

# P4-P6 outer junction
#sequence  = ['GGAAUUGCGGGAAAGGGGUC','GGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCA']
#structure = '.....((((((...(((((())))))..)).))))((...((((...(((((((((...)))))))))..))))...))...'

# P4-P6 outer junction -- further minimized
#sequence  = ['GGAAUUGCGGGAAAGG','CUAACCACGCAGCCAAGUCCUAAGUC','GAUAUGGAUGCAGUUCA']
#structure =  '.....((((((...(())..)).))))((...((((...((()))..))))...))...'

# P4-P6 outer junction -- easier to read input
P4P6_outerjunction_sequence  = 'GGAAUUGCGGGAAAGG CUAACCACGCAGCCAAGUCCUAAGUC GAUAUGGAUGCAGUUCA'
P4P6_outerjunction_structure = '.....((((((...(( ))..)).))))((...((((...((( )))..))))...))...'
P4P6_outerjunction_force_bps = '...............( )........................( )................'

def apply_params_Ceff_Cinit_KdAU_KdGU( params, x ):
    q = 10.0**x
    params.C_eff_stacked_pair = q[0]
    params.initialize_C_eff_stack()
    params.C_init = q[1]
    params.base_pair_types[2].Kd_BP = q[2]
    params.base_pair_types[3].Kd_BP = q[2] # A-U
    params.base_pair_types[4].Kd_BP = q[3]
    params.base_pair_types[5].Kd_BP = q[3] # G-U
    params.K_coax = 0.0

def apply_params_Cinit_CeffSix( params, x ):
    q = 10.0**x
    params.C_init = q[0]
    bpts_WC = params.base_pair_types[0:4]
    bpt_GU  = params.base_pair_types[4]
    bpt_UG  = params.base_pair_types[5]
    for bpt1 in bpts_WC:
        for bpt2 in bpts_WC:
            params.C_eff_stack[bpt1][bpt2] = q[1]
    for bpt in bpts_WC:  params.C_eff_stack[bpt][bpt_GU] = q[2]
    for bpt in bpts_WC:  params.C_eff_stack[bpt][bpt_UG] = q[3]
    params.C_eff_stack[bpt_GU][bpt_GU] = q[4]
    params.C_eff_stack[bpt_GU][bpt_UG] = q[5]
    params.C_eff_stack[bpt_UG][bpt_GU] = q[6]
    params.K_coax = 0.0

def apply_params_Cinit_CeffSix_Kcoax( params, x ):
    apply_params_Cinit_CeffSix( params, x[:-1] )
    params.K_coax = 10.0**x[-1]

def apply_params_Cinit_l_lBP_CeffSix( params, x ):
    apply_params_Cinit_CeffSix( params, np.array([x[i] for i in [0,3,4,5,6,7,8] ]) )
    q = 10.0**x
    params.l = q[1]
    params.l_BP = q[2]

def apply_params_Ceff_Cinit_KdAU_KdGU_Kcoax( params, x ):
    q = 10.0**x
    params.C_eff_stacked_pair = q[0]
    params.initialize_C_eff_stack()
    params.C_init = q[1]
    params.base_pair_types[2].Kd_BP = q[2]
    params.base_pair_types[3].Kd_BP = q[2] # A-U
    params.base_pair_types[4].Kd_BP = q[3]
    params.base_pair_types[5].Kd_BP = q[3] # G-U
    params.K_coax = q[4]

#x0 = np.array( [5, 1, 3, 3] )
#apply_params_func = apply_params_Ceff_Cinit_KdAU_KdGU
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure) ]

#x0 = np.array( [5, 4, 4, 3, 3, 3, 3] )
#apply_params_func = apply_params_Cinit_CeffSix
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure) ]

#x0 = np.array( [1.56, 5.4, 5, 4, 4, 4, 4] )
#apply_params_func = apply_params_Cinit_CeffSix
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_sequence, P4P6_structure) ]

#x0 = np.array( [1.56, 0, 0, 5, 5, 4, 4, 4, 4] )
#apply_params_func = apply_params_Cinit_l_lBP_CeffSix
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure) ]

#x0 = np.array( [5, 1, 3, 3, 1] )
#apply_params_func = apply_params_Ceff_Cinit_KdAU_KdGU_Kcoax
#sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure) ]
#sequence_structure_pairs  = [ (P4P6_outerjunction_sequence, P4P6_outerjunction_structure) ]


x0 = np.array( [1.56, 5.4, 5, 4, 4, 4, 4, 1] )
apply_params_func = apply_params_Cinit_CeffSix_Kcoax
sequence_structure_pairs  = [ (tRNA_sequence , tRNA_structure), (P4P6_outerjunction_sequence, P4P6_outerjunction_structure, P4P6_outerjunction_force_bps) ]

loss = lambda x : free_energy_gap( x, sequence_structure_pairs, apply_params_func )

# simple 1-D scan
#for x0 in range(10):  loss( [x0] )

#result = minimize( loss, x0, method = 'Nelder-Mead' )
result = minimize( loss, x0, method = 'L-BFGS-B' )

final_loss = loss( result.x )
print result
print 'Final parameters:', result.x, 'Loss:',final_loss
