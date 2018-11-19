class BasePairType:
    def __init__( self, nt1, nt2, Kd_BP, match_lowercase = False ):
        '''
        Two sequence characters that get matched to base pair, e.g., 'C' and 'G';
           and the dissociation constant K_d (in M) associated with initial pair formation.
        match_lowercase means match 'x' to 'x', 'y' to 'y', etc.
        TODO: also store chemical modification info.
        '''
        self.nt1 = nt1
        self.nt2 = nt2
        self.Kd_BP = Kd_BP
        self.match_lowercase = ( nt1 == '*' and nt2 == '*' and match_lowercase )
        self.flipped = self # needs up be updated later.

    def is_match( self, s1, s2 ):
        if self.match_lowercase: return ( s1.islower() and s2.islower() and s1 == s2 )
        return ( s1 == self.nt1 and s2 == self.nt2 )

def setup_base_pair_type( params, nt1, nt2, Kd_BP, match_lowercase = False ):
    bpt1 = BasePairType( nt1, nt2, Kd_BP, match_lowercase = match_lowercase )
    params.base_pair_types.append( bpt1 )

    if not match_lowercase:
        bpt2 = BasePairType( nt2, nt1, Kd_BP, match_lowercase = match_lowercase )
        bpt1.flipped = bpt2
        bpt2.flipped = bpt1
        params.base_pair_types.append( bpt2 )

def initialize_base_pair_types( self ):
    self.base_pair_types = []

