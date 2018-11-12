import math
KT_IN_KCAL = 0.61633135471  # 37 Celsius

class AlphaFoldParams:
    '''
    Parameters that Define the statistical mechanical model for RNA folding
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
    params.Kd_BP = params.C_std * math.exp( dG_init/ KT_IN_KCAL)

    dG_CG_CG = -3.30 # Turner 1999 5'-CC-3'/5'-GG-3', kcal/mol
    params.C_eff_stacked_pair = math.exp( -dG_CG_CG / KT_IN_KCAL ) * params.Kd_BP

    params.l = exp( 0.0 )

    return params
