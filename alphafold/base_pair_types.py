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

##################################################################################################
def initialize_base_pair_types( self ):
    self.base_pair_types = []
