from __future__ import print_function
import math
from .constants import KT_IN_KCAL
from .assert_equal import assert_equal

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
            print(' %9.3f' % X.val(i,i+j), end='')
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

def output_test( p, Z_ref = 0, bpp_idx= [], bpp_expected = 0,  deriv_parameters = None, log_derivs_ref = None, ):
    print('Z =',p.Z,' [calc]' )
    print('Z =',Z_ref,' [expected]')
    assert_equal( p.Z, Z_ref )
    print('bpp[%d,%d] = ' % (bpp_idx[0],bpp_idx[1]),p.bpp[ bpp_idx[0] ][ bpp_idx[1] ], ' [calc]')
    print('bpp[%d,%d] = ' % (bpp_idx[0],bpp_idx[1]),bpp_expected,' [expected]')
    assert_equal( p.bpp[ bpp_idx[0] ][ bpp_idx[1] ], bpp_expected )

    if deriv_parameters != None:
        print()
        print( 'd(logZ)/d(log parameter)' )
        log_derivs = p.get_log_derivs( deriv_parameters )
        for i,parameter in enumerate(deriv_parameters): print( parameter,':', log_derivs[i],'[calc] ',log_derivs_ref[i],'[expected]' )
        for log_deriv,log_deriv_ref in zip(log_derivs,log_derivs_ref):  assert_equal( log_deriv, log_deriv_ref )

    print()
