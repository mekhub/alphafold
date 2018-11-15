#
# Much simpler (less intelligent) object for dynamic programming than in dynamic_programming.py --
#  forces code to explicitly figure out updates to values, derivatives, and contributions
#
class DynamicProgrammingMatrix:
    '''
    Dynamic Programming 2-D Matrix that automatically:
      knows how to update values at i,j
    '''
    def __init__( self, N, val = 0.0, diag_val = 0.0, DPlist = None, update_func = None, options = None ):
        self.N = N

        self.Q = [None]*N
        for i in range( N ): self.Q[i] = [val]*N
        for i in range( N ): self.Q[i][i] = diag_val

        self.dQ = [None]*N
        for i in range( N ): self.dQ[i] = [0.0]*N

        self.contribs = [None]*N
        for i in range( N ):
            self.contribs[i] = []
            for j in range( N ): self.contribs[i].append( [] )

        self.contribs_updated = [None]*N
        for i in range( N ): self.contribs_updated[i] = [False]*N

        if DPlist != None: DPlist.append( self )
        self.update_func = update_func

    def val( self, i, j ): return self.Q[i][j]
    def deriv( self, i, j ): return self.dQ[i][j]

    def update( self, partition, i, j ):
        self.Q[ i ][ j ] = 0
        self.dQ[ i ][ j ] = 0
        self.contribs[ i ][ j ] = []
        self.update_func( partition, i, j )

    def get_contribs( self, partition, i, j ):
        if not self.contribs_updated[i][j]:
            partition.options.calc_contrib = True
            self.update( partition, i, j )
            partition.options.calc_contrib = False
            self.contribs_updated[i][j] = True
        return self.contribs[i][j]

    def __len__( self ):
        return len( self.Q )

class DynamicProgrammingList:
    '''
    Dynamic Programming 1-D list that automatically:
      does wrapping modulo N,
      knows how to update values at i,j
    Used for Z_final
    '''
    def __init__( self, N, val = 0.0, update_func = None, options = None ):
        self.N = N
        self.Q = [ val ]*N
        self.dQ = [ 0.0 ]*N
        self.contribs = [None] * N
        for i in range( N ): self.contribs[i] = []
        self.contribs_updated = [False]*N
        self.update_func = update_func

    def __len__( self ): return self.N

    def val( self, i ): return self.Q[i]
    def deriv( self, i ): return self.dQ[i]

    def update( self, partition, i ):
        self.Q[ i ] = 0.0
        self.dQ[ i ] = 0.0
        self.contribs[ i ] = []
        self.update_func( partition, i )

    def get_contribs( self, partition, i ):
        if not self.contribs_updated[i]:
            partition.options.calc_contrib = True
            self.update( partition, i )
            partition.options.calc_contrib = False
            self.contribs_updated[i] = True
        return self.contribs[i]
