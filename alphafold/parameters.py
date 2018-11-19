import math
from base_pair_types import BasePairType
from util.constants import KT_IN_KCAL
from old_parameters import *

class AlphaFoldParams:
    '''
    Parameters that define the statistical mechanical model for RNA folding
    '''
    def __init__( self ):
        self.name = 'empty'
        self.version = 0.0

        # Seven parameter model (five without coaxial stacking)
        self.C_init = 0.0  # Effective molarity for starting each loop (units of M)
        self.l      = 0.0  # Effective molarity penalty for each linkages in loop (dimensionless)
        self.l_BP   = 0.0  # Effective molarity penalty for each base pair in loop (dimensionless)
        self.C_eff_stacked_pair = 0.0 # Effective molarity for forming stacked pair (units of M)
        self.K_coax = 0.0     # coax bonus for contiguous helices (dimensionless). Set to 0 to turn off coax (dimensionless)
        self.l_coax = 0.0     # Effective molarity bonus for each coaxial stack in loop. Initial guess: C_eff_stacked_pair / (C_init*l*K_coax)

        # default base pair types -- see below for fill in, including Kd
        self.base_pair_types = []

        self.C_std  = 1.0      # 1 M. drops out in end (up to overall scale factor).
        self.min_loop_length = 1 # Disallow apical loops smaller than this size (integer)
        self.allow_strained_3WJ = False # Prevent strained three-way-junctions with two helices coaxially stacked and no spacer nucleotides to other helix.


    def get_variables( self ):
        if self.C_init == 0.0 and self.name == 'empty': print 'WARNING! C_init not defined, and params appear empty. Look at get_minimal_params() or get_latest_params() for examples'
        return ( self.C_init, self.l, self.l_BP, self.C_eff_stacked_pair, self.K_coax, self.l_coax, self.C_std, self.min_loop_length, self.allow_strained_3WJ )

def get_params( params = None, suppress_all_output = False ):
    params_object = None

    if isinstance(params,AlphaFoldParams): return params
    elif params == None or params =='': params_object = get_latest_params()
    elif params == 'minimal':         params_object = get_minimal_params()
    elif params == 'v0.1':  params_object = get_params_v0_1( AlphaFoldParams()  )
    elif params == 'v0.15': params_object = get_params_v0_15( AlphaFoldParams() )
    elif params == 'v0.16': params_object = get_params_v0_16( AlphaFoldParams() )
    else: print 'unrecognized params requested: ', params
    if not suppress_all_output: print 'Parameters: ', params_object.name, ' version', params_object.version
    return params_object

def get_latest_params():
    return get_params_v0_16( AlphaFoldParams() )

def get_minimal_params():
    params = AlphaFoldParams()
    params.name    = 'minimal'
    params.version = '0.1'

    # Seven parameter model
    params.C_init = 1.0     # Effective molarity for starting each loop (units of M)
    params.l      = 0.5     # Effective molarity penalty for each linkages in loop (dimensionless)
    params.l_BP   = 0.2     # Effective molarity penalty for each base pair in loop (dimensionless)
    params.C_eff_stacked_pair = 1e4 # Effective molarity for forming stacked pair (units of M)
    params.K_coax = 100     # coax bonus for contiguous helices (dimensionless). Set to 0 to turn off coax (dimensionless)
    params.l_coax = 200     # Effective molarity bonus for each coaxial stack in loop. Initial guess: C_eff_stacked_pair / (C_init*l*K_coax)
    params.C_std = 1.0      # 1 M. drops out in end (up to overall scale factor).
    params.min_loop_length = 1 # Disallow apical loops smaller than this size (integer)
    params.allow_strained_3WJ = False # Prevent strained three-way-junctions with two helices coaxially stacked and no spacer nucleotides to other helix.

    params.base_pair_types = []
    Kd_BP  = 0.0002  # Kd for forming base pair (units of M )
    params.base_pair_types.append( BasePairType( '*', '*', Kd_BP, match_lowercase = True ) )
    params.base_pair_types.append( BasePairType( 'C', 'G', Kd_BP ) )

    return params

def get_params_v0_16( params ):
    # Starting to make use of train_alphafold.py on tRNA and P4-P6.
    params.name     = 'alphafold'
    params.version  = '0.16'

    params.min_loop_length = 3

    # Seven parameter model
    dG_init = +4.09 # Turner 1999, kcal/mol
    Kd_BP_CG = 1.0 * math.exp( dG_init/ KT_IN_KCAL) # 762 M
    Kd_BP_AU = 10.0**5.0 # 100000 M
    Kd_BP_GU = 10.0**5.0 # 100000 M

    dG_terminal_AU = 0.5 # Turner 1999, kcal/mol -- NUPACK

    dG_CG_CG = -3.30 # Turner 1999 5'-CC-3'/5'-GG-3', kcal/mol
    params.C_eff_stacked_pair = 10**5.425 # about 10^5

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

