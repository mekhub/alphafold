from copy import deepcopy

class DynamicProgrammingData:
    '''
    Dynamic programming object, with derivs and contribution accumulation.
     X   = values (N x N)
     dQ  = derivatives (N X N)
     X_contrib = contributions (coming soon)
    '''
    def __init__( self, N ):
        self.Q = []
        for i in range( N ): self.Q.append( [0.0]*N )

        self.dQ = deepcopy( self.Q ) # another zero matrix.

        self.contrib = []
        for i in range( N ): self.contrib.append( [[]]*N )

    def __getitem__( self, idx ):
        # overloaded []. warning: overhead! directly access object.Q[ idx ] in inner loops.
        return self.Q[ idx ]

    def __len__( self ): return len( self.Q )

    def add( self, i, j, b ):
        #  trying out a function that might make code more readable,
        #  (could hide all contribution accumulation for backtracking -- and
        #   perhaps even derivatives -- inside class!)
        # but this kind of thing appears to take up too much overhead.
        self.Q[i][j]  += b
        self.dQ[i][j] += 0
        self.Q_contrib[i][j].append( [i,j,b] )

