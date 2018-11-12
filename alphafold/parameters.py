import math
from alphafold.base_pair_types import BasePairType
KT_IN_KCAL = 0.61633135471  # 37 Celsius

class AlphaFoldParams:
    '''
    Parameters that define the statistical mechanical model for RNA folding
    '''
    def __init__( self ):
        # Seven parameter model
        self.C_init = 1.0     # Effective molarity for starting each loop (units of M)
        self.l      = 0.5     # Effective molarity penalty for each linkages in loop (dimensionless)
        self.Kd_BP  = 0.0002  # Kd for forming base pair (units of M )
        self.l_BP   = 0.2     # Effective molarity penalty for each base pair in loop (dimensionless)
        self.C_eff_stacked_pair = 1e4 # Effective molarity for forming stacked pair (units of M)
        self.K_coax = 100     # coax bonus for contiguous helices (dimensionless). Set to 0 to turn off coax (dimensionless)
        self.l_coax = 200     # Effective molarity bonus for each coaxial stack in loop. Initial guess: C_eff_stacked_pair / (C_init*l*K_coax)
        self.C_std = 1.0      # 1 M. drops out in end (up to overall scale factor).
        self.min_loop_length = 1 # Disallow apical loops smaller than this size (integer)
        self.allow_strained_3WJ = False # Prevent strained three-way-junctions with two helices coaxially stacked and no spacer nucleotides to other helix.


        # default base pair types
        self.base_pair_types = []
        self.base_pair_types.append( BasePairType( '*', '*', self.Kd_BP, match_lowercase = True ) )
        self.base_pair_types.append( BasePairType( 'C', 'G', self.Kd_BP ) )
        self.base_pair_types.append( BasePairType( 'G', 'C', self.Kd_BP ) )

    def get_variables( self ):
        return ( self.C_init, self.l, self.Kd_BP, self.l_BP, self.C_eff_stacked_pair, self.K_coax, self.l_coax, self.C_std, self.min_loop_length, self.allow_strained_3WJ )

def get_params( params = None ):
    if params == 'devel':
        return get_devel_params()
    return AlphaFoldParams()


def get_devel_params():
    params = AlphaFoldParams()
    params.min_loop_length = 3

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
    base_pair_types.append( BasePairType( 'G', 'C', Kd_BP_CG ) )
    base_pair_types.append( BasePairType( 'A', 'U', Kd_BP_AU ) )
    base_pair_types.append( BasePairType( 'U', 'A', Kd_BP_AU ) )
    base_pair_types.append( BasePairType( 'G', 'U', Kd_BP_GU ) )
    base_pair_types.append( BasePairType( 'U', 'G', Kd_BP_GU ) )


    #Kd_BP_GA = Kd_BP_AU * 40 # fudge factor to make GU weaker.
    #base_pair_types.append( BasePairType( 'G', 'A', Kd_BP_GA ) ) # totally made up
    #base_pair_types.append( BasePairType( 'A', 'G', Kd_BP_GA ) ) # totally made up
    #Kd_BP_AA = Kd_BP_AU * 40 # fudge factor to make GU weaker.
    #base_pair_types.append( BasePairType( 'A', 'A', Kd_BP_AA ) ) # totally made up

    params.base_pair_types = base_pair_types
    params.K_coax = 10
    params.l_coax = 1

    return params
