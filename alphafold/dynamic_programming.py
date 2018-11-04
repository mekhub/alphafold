from copy import deepcopy

class DynamicProgrammingMatrix:
    '''
    Dynamic Programming 2-D Matrix that automatically:
      does wrapping modulo N,
      knows how to update values at i,j
    '''
    def __init__( self, N, val = 0.0, diag_val = 0.0, DPlist = None, update_func = None ):
        self.N = N
        #self.DPmatrix = WrappedArray( N ) # TODO
        self.data = [None]*N
        for i in range( N ):
            #self.data[i] = WrappedArray( N ) # TODO
            self.data[i] = [None]*N
            for j in range( N ):
                self.data[i][j] = DynamicProgrammingData( val )
                self.data[i][j].info.append( (self,i,j) )

        for i in range( N ): self.data[i][i].Q = diag_val
        if DPlist != None: DPlist.append( self )
        self.update_func = update_func

    def __getitem__( self, idx ):
        return self.data[ idx ]

    def __len__( self ): return len( self.data )

    def update( self, partition, i, j ):
        self.data[ i ][ j ].zero()
        self.update_func( partition, i, j )


class DynamicProgrammingList:
    '''
    Dynamic Programming 1-D list that automatically:
      does wrapping modulo N,
      knows how to update values at i,j
    Used for Z_final
    '''
    def __init__( self, N, val = 0.0, update_func = None ):
        self.N = N
        self.data = []
        for i in range( N ):
            self.data.append( DynamicProgrammingData( val ) )

        self.update_func = update_func

    def __getitem__( self, idx ):
        return self.data[ idx ]

    def __len__( self ): return len( self.data )

    def update( self, partition, i ):
        self.data[ i ].zero()
        self.update_func( partition, i )

class DynamicProgrammingData:
    '''
    Dynamic programming object, with derivs and contribution accumulation.
     Q   = value
     dQ  = derivative (later will generalize to gradient w.r.t. all parameters)
     contrib = contributions
    '''
    def __init__( self, val = 0.0 ):
        self.Q = val
        self.dQ = 0.0
        self.contribs = []
        self.info = []

    def zero( self ):
        self.Q = 0.0
        self.dQ = 0.0
        self.contribs = []

    def __iadd__(self, other):
        self.Q += other.Q
        if len( other.info ) > 0: self.contribs.append( [other.Q, other.info] )
        return self

    def __mul__(self, other):
        prod = DynamicProgrammingData()
        if isinstance( other, DynamicProgrammingData ):
            prod.Q       = self.Q * other.Q
            info = self.info + other.info
            if len( info ) > 0:
                prod.contribs = [ [ prod.Q, info ] ]
                prod.info = info
        else:
            prod.Q = self.Q * other
            for contrib in self.contribs:
                prod.contribs.append( [contrib[0]*other, contrib[1] ] )
            prod.info = self.info
        return prod

    def __truediv__( self, other ):
        quot = deepcopy( self )
        quot.Q /= other
        return quot

    __rmul__ = __mul__
    __floordiv__ = __truediv__
    __div__ = __truediv__


class WrappedArray:
    '''
    For all the various cross-checks, like equality of partition function starting at any
     i and wrapping around to N and then back 1 and then i-1, need to keep applying modulo N.
    '''
    def __init__( self, N, val = 0.0 ):
        self.data = [val] * N
        self.N = N
    def __getitem__( self, idx ):
        return self.data[idx % self.N]
    def __setitem__( self, idx, item ):
        self.data[idx % self.N] = item
    def __len__( self ):
        return self.N
