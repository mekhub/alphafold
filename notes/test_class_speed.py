#!/usr/bin/python
import time
import sys
import os
from copy import deepcopy
sys.path.append(os.path.join(os.getcwd(), '..'))
from alphafold.partition import DynamicProgrammingData as DP

def getval( DP, idx ):
    return DP.X[ idx ][ idx ]

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
for i in range( N ): y = xDP.X[56][56]
t1 = time.time()
print t1 - t0,'y = xDP.X[56][56]'

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
for i in range( N ): x[56][56] = 20
t1 = time.time()
print t1 - t0, 'x[56][56] = 20'

t0 = time.time()
for i in range( N ): xDP.X[56][56] = 20
t1 = time.time()
print t1 - t0,'xDP.X[56][56] = 20'

t0 = time.time()
for i in range( N ):
    val = 20
    xDP.X[56][56] = val
t1 = time.time()
print t1 - t0,'val = 20; xDP.X[56][56] = val'

t0 = time.time()
for i in range( N ): xDP[56][56] = 20
t1 = time.time()
print t1 - t0,'xDP[56][56] = 20'


# Time setting, including derivs
print 'SETTING INCLUDE DERIVS'
t0 = time.time()
for i in range( N ):
    x[56][56] = 20
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
    xDP.X[56][56] = 20
    xDP.dX[56][56] = 0
t1 = time.time()
print t1 - t0,'xDP.X[56][56] = 20,  xDP.dX[56][56]'

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
    xDP.X[56][56] = 20
    xDP.dX[56][56] = 0
    xDP.X_contrib[56][56].append( [x,56,56,20] )
t1 = time.time()
print t1 - t0,'xDP.X[56][56] = 20'

t0 = time.time()
for i in range( N ):
    xDP.add(56,56,20)
t1 = time.time()
print t1 - t0,'xDP += 20'
