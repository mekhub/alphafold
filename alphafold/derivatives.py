def _get_log_derivs( self, parameters ):
    '''
    Output

       d( log Z )/ d( log parameter )

    using simple expressions that require O( N^2 ) time or less after
    the original O( N^3 ) dynamic programming calculations
    '''

    derivs = [None]*len(parameters)

    # Derivatives with respect to each Kd
    for n,parameter in enumerate(parameters):
        if parameter == 'Kd':
            Kd_derivs = []
            for base_pair_type in self.params.base_pair_types:  Kd_derivs.append( get_Kd_deriv( self, base_pair_type ) )
            derivs[ n ] = sum( Kd_derivs ) / 2.0
        else:
            Kd_tag = parameter[3:]
            for base_pair_type in self.params.base_pair_types:
                if Kd_tag == 'match_lowercase' and base_pair_type.match_lowercase:
                    derivs[ n ] = get_Kd_deriv( self, base_pair_type )
                    break
                if Kd_tag == base_pair_type.nt1 + base_pair_type.nt2:
                    derivs[ n ] = get_Kd_deriv( self, base_pair_type )
                    break

    # Derivatives with respect to loop closure parameters

    # Derivatives with respect to motifs (stacked pairs first)

    return derivs


def get_Kd_deriv( self, base_pair_type ):
    assert( self.calc_all_elements )
    bpp = 0.0
    for i in range( self.N ):
        for j in range( self.N ):
            bpp += self.Z_BPq[base_pair_type].val(i,j) * self.Z_BPq[base_pair_type.flipped].val(j,i) * base_pair_type.Kd / self.Z_final.val(0)
    return bpp

