#!/usr/bin/python
import time
import sys
import os
from copy import deepcopy
sys.path.append(os.path.join(os.getcwd(), '..'))
#from alphafold.partition import DynamicProgrammingData as DP


###################################################################################################################33
class DP:
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
        self.contrib[i][j].append( [i,j,b] )

def getval( DP, idx ):
    return DP.Q[ idx ][ idx ]

x = [[]]*500
for i in range( 500 ): x[i] = [0.0]*500
dx = deepcopy( x )
xcontrib = [[]]*500
for i in range( 500 ): xcontrib[i] = [[]]*500

xDP = DP( 500 )  # 500x500 object with other stuff in it.

N = 500000
print 'Try for ', N, 'cycles each:'
# Time getting
print 'GETTING'
t0 = time.time()
for i in range( N ): y = x[56][56]
t1 = time.time()
print t1 - t0, 'y = x[56][56]'


t0 = time.time()
for i in range( N ): y = xDP.Q[56][56]
t1 = time.time()
print t1 - t0,'y = xDP.Q[56][56]'

t0 = time.time()
for i in range( N ): y = getval(xDP,56)
t1 = time.time()
print t1 - t0,  'y = getval(xDP,56)'

t0 = time.time()
for i in range( N ): y = xDP[56][56]
t1 = time.time()
print t1 - t0,  'y = xDP[56][56]'


# Time setting
print 'SETTING'
t0 = time.time()
for i in range( N ): x[56][56] = 20 * 1.5 * 1.2 * x[1][1] * x[5][5]
t1 = time.time()
print t1 - t0, 'x[56][56] = 20 * 1.5 * 1.2 * x[1][1] * x[5][5]'

t0 = time.time()
for i in range( N ): xDP.Q[56][56] = 20 * 1.5 * 1.2 * x[1][1] * x[5][5]
t1 = time.time()
print t1 - t0,'xDP.Q[56][56] = 20 * 1.5 * 1.2 * x[1][1] * x[5][5]'

t0 = time.time()
for i in range( N ):
    val = 20  * 1.5 * 1.2 * x[1][1] * x[5][5]
    xDP.Q[56][56] = val
t1 = time.time()
print t1 - t0,'val = 20 * 1.5 * 1.2 * x[1][1] * x[5][5]; xDP.Q[56][56] = val'

t0 = time.time()
for i in range( N ):
    xDP.Q[56][56] = val = 20  * 1.5 * 1.2 * x[1][1] * x[5][5]
t1 = time.time()
print t1 - t0,'xDP.Q[56][56] = val = 20 * 1.5 * 1.2 * x[1][1] * x[5][5]'

t0 = time.time()
for i in range( N ):
    if x[1][1] > 0.0: xDP.Q[56][56] = val = 20  * 1.5 * 1.2 * x[1][1] * x[5][5]
t1 = time.time()
print t1 - t0,'if x[1][1] > 0.0: xDP.Q[56][56] = val = 20 * 1.5 * 1.2 * x[1][1] * x[5][5]'

t0 = time.time()
for i in range( N ): xDP[56][56] = 20 * 1.5 * 1.2 * x[1][1] * x[5][5]
t1 = time.time()
print t1 - t0,'xDP[56][56] = 20 * 1.5 * 1.2 * x[1][1] * x[5][5]'


# Time setting, including derivs
print 'SETTING INCLUDE DERIVS'
t0 = time.time()
for i in range( N ):
    x[56][56] = 20 * x[1][1] * x[5][5]
    dx[56][56] = 0
t1 = time.time()
print t1 - t0, 'x[56][56] = 20, dx[56][56] = 20'

t0 = time.time()
for i in range( N ):
    x[56][56] = (20,0)
t1 = time.time()
print t1 - t0, 'x[56][56] = (20,0)'

t0 = time.time()
for i in range( N ):
    xDP.Q[56][56] = 20
    xDP.dQ[56][56] = 0
t1 = time.time()
print t1 - t0,'xDP.Q[56][56] = 20,  xDP.dQ[56][56]'

t0 = time.time()
for i in range( N ):
    xDP.add(56,56,20)
t1 = time.time()
print t1 - t0,'xDP += 20'


# Time setting, including derivs and contribs
print 'SETTING INCLUDE DERIVS AND CONTRIBS'
t0 = time.time()
for i in range( N ):
    x[56][56] = 20
    dx[56][56] = 0
    xcontrib[56][56].append( [x,56,56,20] )

t1 = time.time()
print t1 - t0, 'x[56][56] = 20'

t0 = time.time()
for i in range( N ):
    xDP.Q[56][56] = 20
    xDP.dQ[56][56] = 0
    xDP.contrib[56][56].append( [x,56,56,20] )
t1 = time.time()
print t1 - t0,'xDP.Q[56][56] = 20'

t0 = time.time()
for i in range( N ):
    xDP.add(56,56,20)
t1 = time.time()
print t1 - t0,'xDP += 20'
