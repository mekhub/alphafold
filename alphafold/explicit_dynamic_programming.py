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

        if DPlist != None: DPlist.append( self )
        self.update_func = update_func

    def val( self, i, j ): return self.Q[i][j]
    def deriv( self, i, j ): return self.dQ[i][j]
    def get_contribs( self, i, j ): return self.contribs[i][j]

    def update( self, partition, i, j ):
        self.Q[ i ][ j ] = 0
        self.dQ[ i ][ j ] = 0
        self.contribs[ i ][ j ] = []
        self.update_func( partition, i, j )

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
        self.update_func = update_func

    def __len__( self ): return self.N

    def val( self, i ): return self.Q[i]
    def deriv( self, i ): return self.dQ[i]
    def get_contribs( self, i ): return self.contribs[i]

    def update( self, partition, i ):
        self.Q[ i ] = 0.0
        self.dQ[ i ] = 0.0
        self.contribs[ i ] = []
        self.update_func( partition, i )

class DynamicProgrammingData:
    '''
    Dynamic programming object, with derivs and contribution accumulation.
     Q   = value
     dQ  = derivative (later will generalize to gradient w.r.t. all parameters)
     contrib = contributions
    '''
    def __init__( self, val = 0.0, options = None ):
        self.Q = val
        self.dQ = 0.0
        self.contribs = []
