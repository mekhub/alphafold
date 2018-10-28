class PeriodicArray:
    def __init__( self, N ):
        self.array = [0] * N
        self.N = N
    def __getitem__( self, idx ):
        return self.array[idx % self.N]
    def __setitem__(  self, idx, item ):
        self.array[idx % self.N] = item


