import math
from util.constants import KT_IN_KCAL
from base_pair_types import BasePairType

# Parameters developed before extensive optimization
def get_params_v0_15( params ):
    # Starting to make use of train_alphafold.py on P4-P6.
    params.name     = 'alphafold'
    params.version  = '0.15'

    params.min_loop_length = 3

    # Seven parameter model
    dG_init = +4.09 # Turner 1999, kcal/mol
    Kd_BP_CG = 1.0 * math.exp( dG_init/ KT_IN_KCAL)
    dG_terminal_AU = 0.5 # Turner 1999, kcal/mol -- NUPACK
    #Kd_BP_AU = 10.0**5.40731537 # 255000
    Kd_BP_GU = 10.0**3.5863221 # 4000
    Kd_BP_AU = 10.0**5.40731537 # 255000
    Kd_BP_GU = 10.0**3.5863221 # 4000

    dG_CG_CG = -3.30 # Turner 1999 5'-CC-3'/5'-GG-3', kcal/mol
    params.C_eff_stacked_pair = math.exp( -dG_CG_CG / KT_IN_KCAL ) * Kd_BP_CG # about 10^5

    # From nupack rna1999.params
    #>Multiloop terms: ALPHA_1, ALPHA_2, ALPHA_3
    #>ML penalty = ALPHA_1 + s * ALPHA_2 + u *ALPHA_3
    #>s = # stems adjacent to ML, u = unpaired bases in ML
    # 340   40    0
    #>AT_PENALTY:
    #>Penalty for non GC pairs that terminate a helix
    #  50
    dG_bulge = 3.4 # bulge cost is roughly 3-4 kcal/mol
    dG_multiloop_stems = 0.40 # in kcal/mol
    dG_multiloop_unpaired = 0.0 #0.40 # in kcal/mol -- ZERO in NUPACK -- fudging here.
    #params.C_init = 1.0 * math.exp( -dG_bulge / KT_IN_KCAL )
    params.C_init = 10** 1.55034499 # radically different

    params.l = math.exp( dG_multiloop_unpaired / KT_IN_KCAL )
    params.l_BP = math.exp( dG_multiloop_stems/KT_IN_KCAL ) / params.l

    base_pair_types = []
    base_pair_types.append( BasePairType( 'C', 'G', Kd_BP_CG ) )
    base_pair_types.append( BasePairType( 'A', 'U', Kd_BP_AU ) )
    base_pair_types.append( BasePairType( 'G', 'U', Kd_BP_GU ) )

    #Kd_BP_GA = Kd_BP_AU*100000  # fudge factor to make GA weaker.
    #base_pair_types.append( BasePairType( 'G', 'A', Kd_BP_GA ) ) # totally made up
    #Kd_BP_AA = Kd_BP_AU * 40 # fudge factor to make GU weaker.
    #base_pair_types.append( BasePairType( 'A', 'A', Kd_BP_AA ) ) # totally made up

    params.base_pair_types = base_pair_types

    # turn off coax!?
    params.K_coax = 0.0
    params.l_coax = 1.0

    return params

def get_params_v0_1( params ):
    params.name     = 'alphafold'
    params.version  = '0.1'

    params.min_loop_length = 3

    # Seven parameter model
    dG_init = +4.09 # Turner 1999, kcal/mol
    Kd_BP_CG = 1.0 * math.exp( dG_init/ KT_IN_KCAL)
    dG_terminal_AU = 0.5 # Turner 1999, kcal/mol -- NUPACK
    Kd_BP_AU = Kd_BP_CG *math.exp( dG_terminal_AU/ KT_IN_KCAL )
    Kd_BP_GU = Kd_BP_AU * 10 # fudge factor to make GU weaker.

    dG_CG_CG = -3.30 # Turner 1999 5'-CC-3'/5'-GG-3', kcal/mol
    params.C_eff_stacked_pair = math.exp( -dG_CG_CG / KT_IN_KCAL ) * Kd_BP_CG

    # From nupack rna1999.params
    #>Multiloop terms: ALPHA_1, ALPHA_2, ALPHA_3
    #>ML penalty = ALPHA_1 + s * ALPHA_2 + u *ALPHA_3
    #>s = # stems adjacent to ML, u = unpaired bases in ML
    # 340   40    0
    #>AT_PENALTY:
    #>Penalty for non GC pairs that terminate a helix
    #  50
    dG_bulge = 3.4 # bulge cost is roughly 3-4 kcal/mol
    dG_multiloop_stems = 0.40 # in kcal/mol
    dG_multiloop_unpaired = 0.0 #0.40 # in kcal/mol -- ZERO in NUPACK -- fudging here.
    params.C_init = 1.0 * math.exp( -dG_bulge / KT_IN_KCAL )
    params.l = math.exp( dG_multiloop_unpaired / KT_IN_KCAL )
    params.l_BP = math.exp( dG_multiloop_stems/KT_IN_KCAL ) / params.l

    base_pair_types = []
    base_pair_types.append( BasePairType( 'C', 'G', Kd_BP_CG ) )
    base_pair_types.append( BasePairType( 'A', 'U', Kd_BP_AU ) )
    base_pair_types.append( BasePairType( 'G', 'U', Kd_BP_GU ) )

    #Kd_BP_GA = Kd_BP_AU*100000  # fudge factor to make GA weaker.
    #base_pair_types.append( BasePairType( 'G', 'A', Kd_BP_GA ) ) # totally made up
    #Kd_BP_AA = Kd_BP_AU * 40 # fudge factor to make GU weaker.
    #base_pair_types.append( BasePairType( 'A', 'A', Kd_BP_AA ) ) # totally made up

    params.base_pair_types = base_pair_types
    params.K_coax = 10
    params.l_coax = 1

    return params
