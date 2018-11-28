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
            for i in range( N ):
                if not self.ligated[i]: continue
                num_internal_linkages += self.params.l * self.C_eff_no_coax_singlet.val( i+1, i ) / self.Z
            derivs[ n ] = num_internal_linkages
        elif parameter == 'l_BP':
            num_base_pairs_closed_by_loops = 0.0
            for i in range( N ):
                for j in range( N ):
                    if ( j - i ) % N < 2: continue
                    num_base_pairs_closed_by_loops += self.params.l**2 * self.params.l_BP * self.C_eff.val(i+1,j-1) * self.Z_BP.val(j,i) / self.Z
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
        elif len(parameter)>=5 and parameter[:5] == 'C_eff':
            # Derivatives with respect to motifs (stacked pairs first)
            if parameter == 'C_eff_stacked_pair':
                print( 'YOOOO!' )
                motif_prob = 0.0
                for i in range( N ):
                    for j in range( N ):
                        if self.Z_BP.val(i,j) > 0.0 and self.Z_BP.val(j-1,i+1)>0:
                            for base_pair_type in self.params.base_pair_types:
                                if self.Z_BPq[base_pair_type].val(i,j) == 0.0: continue
                                for base_pair_type2 in self.params.base_pair_types:
                                    if self.Z_BPq[base_pair_type2].val(j-1,i+1) == 0.0: continue
                                    Z_BPq1 = self.Z_BPq[base_pair_type]
                                    Z_BPq2 = self.Z_BPq[base_pair_type2]
                                    motif_prob += self.params.C_eff_stack[base_pair_type][base_pair_type2.flipped] * Z_BPq1[i][j] * Z_BPq2[j-1][i+1] / self.Z
                derivs[n] = motif_prob
            else:
                #TODO
                pass
        else:
            #TODO some kind of informative error message that parameters is not a 'legitimate' one
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

