def secstruct( bps, N ):
    '''
    Convert list of base pairs to dot-paren string. N is length of RNA.
    '''
    secstruct_string = ['.']*N
    for bp in bps:
        i = min( bp )
        j = max( bp )
        secstruct_string[i] = '('
        secstruct_string[j] = ')'
    return ''.join(secstruct_string)

def bps( secstruct ):
    '''
    Convert dot-paren secstruct into sorted list of base pairs
    '''
    bps_list = []
    N = len ( secstruct )
    leftbrackets = []
    for i,char in enumerate( secstruct ):
        if char == ')':
            bps_list.append( (leftbrackets[-1],i) )
            leftbrackets = leftbrackets[:-1]
        elif char == '(':
            leftbrackets.append( i )
    assert( len( leftbrackets ) == 0 )
    bps_list.sort()
    return bps_list

def get_structure_string( structure ):
    if structure == None: return None
    if isinstance( structure, list ): structure = ''.join( structure )
    return structure.replace( ' ','' ).replace('+','').replace(',','')

def parse_motifs( secstruct, N = 0 ):
    '''
    Parse secstruct into its structural motifs:
      hairpins, interior loops, multiway junctions, exterior strands
    '''
    if isinstance( secstruct, list ):
        bps_list = secstruct
        assert( N > 0 ) # must provide N if secstruct is entered as a list.
    else:
        assert( isinstance( secstruct, str ) )
        bps_list = bps( secstruct )
        N = len( secstruct )
    pair_map = {}
    for i,j in bps_list:
        pair_map[i] = j
        pair_map[j] = i

    motifs = []
    motif = []
    strand = [N-1]
    linkage_in_motif = [False]*N
    for i in range( N ):
        if i > 0 and linkage_in_motif[ i-1 ]:
            # skip if we already put this nt in a previous motif.
            strand = [i]
            continue
        # continue down a strand until we hit a base pair
        strand.append( i )
        if i != strand[ 0 ] and pair_map.has_key( i ):
            # OK found a base pair, first strand is defined
            motif.append( strand )
            motif_start = strand[0]
            j = pair_map[i]
            strand = []
            # now follow it around until either we come back (cycle) or hit N (exterior loop)
            while j != motif_start:
                strand.append( j )
                j = (j + 1) % N
                strand.append( j )
                while not pair_map.has_key( j ) and j != motif_start:
                    j = ( j + 1 ) % N
                    strand.append( j )
                motif.append( strand )
                strand = []
                if j == motif_start: break
                j = pair_map[ j ]
            if motif[-1][-1] == motif[0][0]:
                # merge last and first strand (they form a cycle that goes across the circle from N-1 back around to 0)
                motif[0] = motif[-1][:-1] + motif[0]
                motif = motif[:-1]
            for strand in motif:
                assert( len( strand ) >= 2 )
                for m in strand[:-1]: linkage_in_motif[ m ] = True
            motifs.append( motif )
            motif = []
            strand = [i]
    motifs.sort()
    return motifs

