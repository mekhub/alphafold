class WrappedArray:
    '''
    For all the various cross-checks, like equality of partition function starting at any
     i and wrapping around to N and then back 1 and then i-1, need to keep applying modulo N.
    '''
    def __init__( self, N, val = None ):
        self.data = [val] * N
        self.N = N
    def __getitem__( self, idx ):
        return self.data[idx % self.N]
    def __setitem__( self, idx, item ):
        self.data[idx % self.N] = item
    def __len__( self ):
        return self.N

##################################################################################################
def initialize_matrix( N, val = None ):
    X = WrappedArray( N )
    for i in range( N ):
        X[ i ] = WrappedArray( N )
        for j in range( N ): X[ i ][ j ] = val
    return X
