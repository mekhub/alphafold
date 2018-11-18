#!/usr/bin/python
import argparse
import sequence_util
import secstruct
from partition import partition

parser = argparse.ArgumentParser( description = "Compute nearest neighbor model partitition function for RNA sequence" )
parser.add_argument( "-s","-seq","--sequences",help="RNA sequences (separate by space)",nargs='*')
parser.add_argument("-struct","--structure",type=str, default=None, help='force structure in dot-parens notation')
parser.add_argument("-c","-circ","--circle", action='store_true', default=False, help='Sequence is a circle')
parser.add_argument("-params","--parameters",type=str, default='', help='Parameters to use [default: '']')

args     = parser.parse_args()

# What we get if we parse out motifs
# pre-compute cost of connecting the strands into a complex
Z = 1

# Now go through each motif parsed out of the target structure
structure = args.structure
bps_list  = secstruct.bps( structure )
motifs = secstruct.parse_motifs( structure )
sequence, ligated = sequence_util.initialize_sequence_and_ligated( args.sequences, args.circle )

for motif in motifs:
    motif_res = []
    motif_sequences = []
    for strand in motif:
        strand_sequence = ''
        for i in strand:
            motif_res.append( i )
            strand_sequence += sequence[i]
            if not ligated[i]:
                motif_sequences.append( strand_sequence )
                strand_sequence = ''
        if len( strand_sequence ) > 0: motif_sequences.append( strand_sequence )
    motif_circle = ligated[ motif_res[-1] ] and ( (motif_res[0] - motif_res[-1]) % len(sequence) == 1 )
    motif_sequence = ''.join( motif_sequences )

    # each motif res better show up only once
    assert( len( set( motif_res ) ) == len( motif_res ) )

    motif_bps_list = []
    for i,j in bps_list:
        if motif_res.count( i ) == 0: continue
        if motif_res.count( j ) == 0: continue
        motif_bps_list.append( (motif_res.index(i), motif_res.index(j)) )
    motif_structure = secstruct.secstruct( motif_bps_list, len( motif_res ) )

    p = partition( motif_sequences, circle = motif_circle, structure = motif_structure, params = args.parameters, suppress_all_output = True )

    Z_motif = p.Z

    # Need to 'correct' for half-terminal penalties (a la Turner rules) and also remove extra costs
    # for connecting strands together.

    Kd_ref = p.params.base_pair_types[0].Kd_BP # Kd[G-C], a la Turner rule convention
    C_std  = p.params.C_std

    Z_motif *= ( Kd_ref / C_std ) ** sequence_util.get_num_strand_connections( motif_sequences, motif_circle )
    for i_motif, j_motif in motif_bps_list:
        # what kind of base pair is this?
        for base_pair_type in p.params.base_pair_types:
            if base_pair_type.is_match( motif_sequence[ i_motif ], motif_sequence[ j_motif ] ):
                Z_motif *= ( base_pair_type.Kd_BP / Kd_ref )**(0.5)
                break

    print "Motif: ", Z_motif, motif

    Z *= Z_motif

Z_connect = ( C_std / Kd_ref ) ** sequence_util.get_num_strand_connections( args.sequences, args.circle )
Z *= Z_connect
print "Connect strands: ", Z_connect
print 'Product of motif Z      :', Z

# Reference value from 'hacked' dynamic programming, which takes a while.
p = partition( args.sequences, circle = args.circle, structure = args.structure, params = args.parameters, suppress_all_output = True )
print 'From dynamic programming:', p.Z


