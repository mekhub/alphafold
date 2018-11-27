from __future__ import print_function
import math
from .constants import KT_IN_KCAL

def _show_results( self ):
    print('sequence =', self.sequence)
    if self.structure != None:
        print('structure=',self.structure)
    cutpoint = ''
    for i in range( self.N ):
        if not self.ligated[ i ]: cutpoint += 'X'
        else: cutpoint += '-'
    print('cutpoint =', cutpoint)
    print('circle   = ', self.circle)
    print('Z =',self.Z)
    print('dG =',-KT_IN_KCAL * math.log( self.Z ))

def _show_matrices( self ):
    output_DP( "Z_BP", self.Z_BP )
    output_DP( "C_eff_basic", self.C_eff_basic )
    output_DP( "C_eff", self.C_eff, self.Z_final )
    #output_DP( "dC_eff", self.dC_eff, self.dZ_final )
    output_DP( "Z_coax", self.Z_coax )
    output_DP( "Z_linear", self.Z_linear )
    output_square( "BPP", self.bpp );

def output_DP( tag, X, X_final = []):
    N = len( X )
    print()
    print("-----", tag, "-----")
    for i in range( N ):
        for q in range( i ): print('          ', end='')# padding to line up)
        for j in range( N ):
            print(' %9.3f' % X[i][(i+j) % N], end='')
        if len( X_final ) > 0: print('==> %9.3f' % X_final.val(i), end='')
        print()

def output_square( tag, X ):
    N = len( X )
    print()
    print("-----", tag, "-----")
    for i in range( N ):
        for j in range( N ):
            print(' %9.3f' % X[i][j],)
        print()

def output_test( Z, Z_ref = 0, bpp = [], bpp_idx= [], bpp_expected = 0):
    print('Z =',Z_ref,' [expected]')
    assert( abs( (Z - Z_ref)/Z_ref )  < 1e-5 )
    print('bpp[%d,%d] = ' % (bpp_idx[0],bpp_idx[1]),bpp[ bpp_idx[0] ][ bpp_idx[1] ])
    print('bpp[%d,%d] = ' % (bpp_idx[0],bpp_idx[1]),bpp_expected,' [expected]')
    if bpp_expected > 0.0:
        assert( abs( (bpp[ bpp_idx[0] ][ bpp_idx[1] ] - bpp_expected)/bpp_expected )  < 1e-5 )
    else:
        assert( bpp[ bpp_idx[0] ][ bpp_idx[1] ] == 0.0 )
    print()

