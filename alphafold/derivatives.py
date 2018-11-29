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
                derivs[ n ] = get_bpp_tot_for_base_pair_type( self, get_base_pair_type_for_tag( self, Kd_tag ) )

        elif len(parameter)>=11 and parameter[:11] == 'C_eff_stack':
            # Derivatives with respect to motifs (stacked pairs first)
            if parameter == 'C_eff_stacked_pair':
                motif_prob = 0.0
                for base_pair_type in self.params.base_pair_types:
                    for base_pair_type2 in self.params.base_pair_types:
                        motif_prob += get_motif_prob( self, base_pair_type, base_pair_type2 )
                derivs[n] = motif_prob
            else:
                assert( len(parameter) > 11 )
                tags = parameter[12:].split('_')
                assert( len( tags ) == 2 )
                derivs[ n ] = get_motif_prob( self, get_base_pair_type_for_tag( self, tags[0] ), get_base_pair_type_for_tag( self, tags[1] ) )
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

def get_base_pair_type_for_tag( self, tag ):
    for base_pair_type in self.params.base_pair_types:
        if (tag == 'matchlowercase' and base_pair_type.match_lowercase) or \
           (tag == base_pair_type.nt1 + base_pair_type.nt2 ):
            return base_pair_type
    print( 'Could not figure out base_pair_type for ', tag )
    return None

def get_motif_prob( self, base_pair_type, base_pair_type2 ):
    # base pair forms a stacked pair with previous pair
    #
    #      bp2
    #  i+1 ... j-1
    #    |     |
    #    i ... j
    #      bp1
    #
    motif_prob = 0.0
    Z_BPq1 = self.Z_BPq[base_pair_type.flipped]
    Z_BPq2 = self.Z_BPq[base_pair_type2]
    N = self.N
    for i in range( N ):
        for j in range( N ):
            if ( j - i ) % N < 3: continue
            if Z_BPq1.val(j  ,  i) == 0: continue
            if Z_BPq2.val(i+1,j-1) == 0: continue
            motif_prob += self.params.C_eff_stack[base_pair_type][base_pair_type2] * Z_BPq1.val(j,i) * Z_BPq2.val(i+1,j-1) / self.Z / 2.0
    return motif_prob
