def _get_log_derivs( self, parameters ):
    '''
    Output

       d( log Z )/ d( log parameter )

    using simple expressions that require O( N^2 ) time or less after
    the original O( N^3 ) dynamic programming calculations
    '''

    derivs = [None]*len(parameters)

    # Derivatives with respect to each Kd
    N = self.N
    for n,parameter in enumerate(parameters):
        if parameter == 'l':
            # Derivatives with respect to loop closure parameters
            num_internal_linkages = 0.0
            for i in range( N ): num_internal_linkages += self.params.l * self.C_eff.val( i+1, i )
            num_internal_linkages /= self.Z
            derivs[ n ] = num_internal_linkages
        elif parameter == 'l_BP':
            num_base_pairs_closed_by_loops = 0.0
            for i in range( N ):
                for j in range( N ):
                    num_base_pairs_closed_by_loops += self.params.l**2 * self.params.l_BP * self.C_eff.val(i+1,j-1) * self.Z_BP.val(j,i)
            num_base_pairs_closed_by_loops /= self.Z
            derivs[ n ] = num_base_pairs_closed_by_loops
        elif parameter == 'C_init':
            num_closed_loops = get_bpp_tot( self ) - self.num_strand_connections()
            derivs[ n ] = num_closed_loops
        elif len(parameter)>=2 and  parameter[:2] == 'Kd':
            if parameter == 'Kd':
                derivs[ n ] = get_bpp_tot( self )
            else:
                Kd_tag = parameter[3:]
                for base_pair_type in self.params.base_pair_types:
                    if (Kd_tag == 'match_lowercase' and base_pair_type.match_lowercase) or \
                       (Kd_tag == base_pair_type.nt1 + base_pair_type.nt2 ):
                        derivs[ n ] = get_bpp_tot_for_base_pair_type( self, base_pair_type )
                        break
        elif len(parameter)>=4 and parameter[:4] == 'C_eff':
            # Derivatives with respect to motifs (stacked pairs first)
            pass
        else:
            pass

    return derivs


def get_bpp_tot_for_base_pair_type( self, base_pair_type ):
    assert( self.calc_all_elements )
    bpp = 0.0
    N = self.N
    for i in range( N ):
        for j in range( N ):
            bpp += self.Z_BPq[base_pair_type].val(i,j) * self.Z_BPq[base_pair_type.flipped].val(j,i) * base_pair_type.Kd / self.Z_final.val(0)
    return bpp


def get_bpp_tot( self ):
    bpp_tot = []
    for base_pair_type in self.params.base_pair_types: bpp_tot.append( get_bpp_tot_for_base_pair_type( self, base_pair_type ) )
    return sum( bpp_tot ) / 2.0

