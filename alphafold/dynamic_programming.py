from copy import deepcopy
from wrapped_array import WrappedArray

class DynamicProgrammingMatrix:
    '''
    Dynamic Programming 2-D Matrix that automatically:
      does wrapping modulo N,
      knows how to update values at i,j
    '''
    def __init__( self, N, val = 0.0, diag_val = 0.0, DPlist = None, update_func = None, options = None ):
        self.N = N
        self.data = WrappedArray(N)
        for i in range( N ):
            self.data[i] = WrappedArray( N )
            for j in range( N ):
                self.data[i][j] = DynamicProgrammingData( val, options = options )
                self.data[i][j].info.append( (self,i,j) )

        for i in range( N ): self.data[i][i].Q = diag_val
        if DPlist != None: DPlist.append( self )
        self.options = options
        self.update_func = update_func

        self.contribs_updated = [None]*N
        for i in range( N ): self.contribs_updated[i] = [False]*N

    def __getitem__( self, idx ):
        return self.data[ idx ]

    def __len__( self ): return len( self.data )

    def val( self, i, j ): return self.data[i][j].Q
    def set_val( self, i, j, val ): self.data[i][j].Q = val
    def deriv( self, i, j ): return self.data[i][j].dQ

    def update( self, partition, i, j ):
        self.data[ i ][ j ].zero()
        self.update_func( partition, i, j )

    def get_contribs( self, partition, i, j ):
        if not self.contribs_updated[i][j]:
            partition.options.calc_contrib = True
            self.update( partition, i, j )
            partition.options.calc_contrib = False
            self.contribs_updated[i][j] = True
        return self.data[i][j].contribs

class DynamicProgrammingList:
    '''
    Dynamic Programming 1-D list that automatically:
      does wrapping modulo N,
      knows how to update values at i,j
    Used for Z_final
    '''
    def __init__( self, N, val = 0.0, update_func = None, options = None ):
        self.N = N
        self.data = WrappedArray( N )
        for i in range( N ):
            self.data[i] = DynamicProgrammingData( val, options = options )
        self.options = options
        self.update_func = update_func
        self.contribs_updated = [False]*N

    def __getitem__( self, idx ):
        return self.data[ idx ]

    def __setitem__( self, idx, val ):
        self.data[ idx ] = val

    def __len__( self ): return len( self.data )

    def val( self, i ): return self.data[i].Q
    def deriv( self, i ): return self.data[i].dQ

    def get_contribs( self, partition, i ):
        if not self.contribs_updated[i]:
            partition.options.calc_contrib = True
            self.update( partition, i )
            partition.options.calc_contrib = False
            self.contribs_updated[i] = True
        return self.data[i].contribs

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
    def __init__( self, val = 0.0, options = None ):
        self.Q = val
        self.dQ = 0.0
        self.contribs = []
        self.info = []
        self.options = options

    def zero( self ):
        self.Q = 0.0
        self.dQ = 0.0
        self.contribs = []

    def __iadd__(self, other):
        self.Q  += other.Q
        if self.options and self.options.calc_deriv:
            self.dQ += other.dQ
        if self.options and self.options.calc_contrib:
            if len( other.info ) > 0: self.contribs.append( [other.Q, other.info] )
        return self

    def __mul__(self, other):
        prod = DynamicProgrammingData()
        prod.options = self.options
        if not prod.options: prod.options = other.options
        if isinstance( other, DynamicProgrammingData ):
            prod.Q  = self.Q * other.Q
            if self.options and self.options.calc_deriv:
                prod.dQ = self.Q * other.dQ + self.dQ * other.Q
            if self.options and self.options.calc_contrib:
                info = self.info + other.info
                if len( info ) > 0:
                    prod.contribs = [ [ prod.Q, info ] ]
                    prod.info = info
        else:
            prod.Q  = self.Q * other
            if self.options and self.options.calc_deriv:
                prod.dQ = self.dQ * other
            if self.options and self.options.calc_contrib:
                for contrib in self.contribs:
                    prod.contribs.append( [contrib[0]*other, contrib[1] ] )
                prod.info = self.info
        return prod

    def __truediv__( self, other ):
        return self * ( 1.0/other )

    __rmul__ = __mul__
    __floordiv__ = __truediv__
    __div__ = __truediv__


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

