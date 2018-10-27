class PeriodicArray:
    '''
    For all the various cross-checks, like equality of partition function starting at any
     i and wrapping around to N and then back 1 and then i-1, need to keep applying modulo N.
    '''
    def __init__( self, N, value = 0.0 ):
        self.array = [ value ] * N
        self.N = N
    def __getitem__( self, idx ):
        return self.array[idx % self.N]
    def __setitem__(  self, idx, item ):
        self.array[idx % self.N] = item
    def __len__(self):
         return len(self.array)


