from partition_helpers import *

class AlphaFoldParams:
    def __init__( self ):
        # default
        # Four parameter model
        self.Kd_BP  = 0.001;
        self.C_init = 1          # a bit like exp(a) in multiloop
        self.l      = 0.5        # a bit like exp(b) in multiloop
        self.l_BP   = 0.1        # a bit like exp(c) in multiloop
        self.C_std = 1; # 1 M. drops out in end (up to overall scale factor).
        self.C_init_BP = self.C_init * (self.l_BP/self.l) # 0.2
        self.min_loop_length = 1

    def get_variables( self ):
        return ( self.Kd_BP, self.C_init, self.l, self.l_BP, self.C_std, self.C_init_BP, self.min_loop_length )

def partition( sequences, params = AlphaFoldParams(), circle = False, verbose = False ):
    p = Partition( sequences, params )
    p.circle  = circle
    p.verbose = verbose
    p.run()
    p.show_results()
    return ( p.Z_final[0], p.bpp, p.dZ_final[0] )

class Partition:
    def __init__( self, sequences, params ):
        self.sequences = sequences
        self.params = params
        self.verbose = False
        self.circle = False

    def run( self ):

        initialize_sequence_information( self ) # N, sequence, is_cutpoint, any_cutpoint
        initialize_dynamic_programming_matrices( self ) # ( Z_BP, C_eff, Z_linear, dZ_BP, dC_eff, dZ_linear )

        (Kd_BP, C_init, l, l_BP, C_std, C_init_BP, min_loop_length, N, sequence, is_cutpoint, any_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear ) = unpack_variables( self )

        # do the dynamic programming
        # deriv calculations are making this long and messy; this should be simplifiable
        for offset in range( 1, N ): #length of subfragment
            for i in range( N ): #index of subfragment
                j = (i + offset) % N;  # N cyclizes
                update_Z_BP( self, i, j )
                update_C_eff( self, i, j )
                update_Z_linear( self, i, j )

        # get the answer (in N ways!)
        Z_final  = []
        dZ_final = []
        for i in range( N ):
            Z_final.append( 0 )
            dZ_final.append( 0 )

            if self.is_cutpoint[(i + N - 1) % N]:
                Z_final[i]  += Z_linear[i][(i-1) % N]
                dZ_final[i] += dZ_linear[i][(i-1) % N]
            else:
                # Scaling Z_final by Kd_lig/C_std to match previous literature conventions
                Z_final[i]  += C_eff[i][(i - 1) % N] * l / C_std
                dZ_final[i] += dC_eff[i][(i - 1) % N] * l / C_std

                for c in range( i, i + N - 1):
                    if self.is_cutpoint[c % N]:
                        #any split segments, combined independently
                        Z_final[i]  += Z_linear[i][c % N] * Z_linear[(c+1) % N][ i - 1 ]
                        dZ_final[i] += ( dZ_linear[i][c % N] * Z_linear[(c+1) % N][ i - 1 ] + Z_linear[i][c % N] * dZ_linear[(c+1) % N][ i - 1 ] )


        # base pair probability matrix
        bpp = initialize_zero_matrix( N );
        bpp_tot = 0.0
        for i in range( N ):
            for j in range( N ):
                bpp[i][j] = Z_BP[i][j] * Z_BP[j][i] * Kd_BP * (l_BP / l) / Z_final[0]
                bpp_tot += bpp[i][j]/2.0 # to avoid double counting (i,j) and (j,i)

        if self.verbose:
            output_DP( "Z_BP", Z_BP )
            output_DP( "C_eff", C_eff, Z_final )
            output_DP( "Z_linear", Z_linear )
            output_square( "BPP", bpp );

        # stringent test that partition function is correct -- all the Z(i,i) agree.
        for i in range( N ):
            assert( abs( ( Z_final[i] - Z_final[0] ) / Z_final[0] ) < 1.0e-5 )
            assert( abs( ( dZ_final[i] - dZ_final[0] ) / dZ_final[0] ) < 1.0e-5 )

        # calculate bpp_tot = -dlog Z_final /dlog Kd_BP in two ways! wow cool test
        bpp_tot_based_on_deriv = -dZ_final[0] * Kd_BP / Z_final[0]
        assert( abs( ( bpp_tot - bpp_tot_based_on_deriv )/bpp_tot ) < 1.0e-5 )

        self.Z_final = Z_final
        self.bpp = bpp
        self.dZ_final = dZ_final

    def show_results( self ):
        print 'sequence =', self.sequence
        cutpoint = ''
        for i in range( self.N ):
            if self.is_cutpoint[ i ]: cutpoint += 'X'
            else: cutpoint += '-'
        print 'cutpoint =', cutpoint
        print 'circle   = ', self.circle
        print 'Z =',self.Z_final[0]

def initialize_sequence_information( self ):
    # initialize sequence
    if isinstance( self.sequences, str ): self.sequence = self.sequences
    else:
        self.sequence = ''
        for i in range( len( self.sequences ) ): self.sequence += self.sequences[i]
    self.N = len( self.sequence )

    # initialize cutpoint information
    self.is_cutpoint = [False] * self.N
    if isinstance( self.sequences, list ):
        L = 0
        for i in range( len(self.sequences)-1 ):
            L = L + len( self.sequences[i] )
            self.is_cutpoint[ L-1 ] = True
    if not self.circle: self.is_cutpoint[ self.N-1 ] = True

    self.any_cutpoint = initialize_any_cutpoint( self.is_cutpoint )

def initialize_dynamic_programming_matrices( self ):
    N = self.N
    # initialize dynamic programming matrices
    self.C_eff    = initialize_zero_matrix( N );
    self.Z_BP     = initialize_zero_matrix( N );
    self.Z_linear = initialize_zero_matrix( N );
    for i in range( N ): #length of fragment
        self.C_eff[ i ][ i ] = self.params.C_init
        self.Z_linear[ i ][ i ] = 1

    # first calculate derivatives with respect to Kd_BP
    self.dC_eff    = initialize_zero_matrix( N );
    self.dZ_BP     = initialize_zero_matrix( N );
    self.dZ_linear = initialize_zero_matrix( N );


def update_Z_BP( self, i, j ):
    (Kd_BP, C_init, l, l_BP, C_std, C_init_BP, min_loop_length, N, sequence, is_cutpoint, any_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear ) = unpack_variables( self )
    offset = ( j - i ) % N

    if (( sequence[i] == 'C' and sequence[j] == 'G' ) or ( sequence[i] == 'G' and sequence[j] == 'C' )) and \
          ( any_cutpoint[i][j] or ( ((j-i-1) % N)) >= min_loop_length  ) and \
          ( any_cutpoint[j][i] or ( ((i-j-1) % N)) >= min_loop_length  ):
        if (not is_cutpoint[ i ]) and (not is_cutpoint[ (j-1) % N]):
            Z_BP[i][j]  += (1.0/Kd_BP ) * ( C_eff[(i+1) % N][(j-1) % N] * l * l)
            dZ_BP[i][j] += (1.0/Kd_BP ) * ( dC_eff[(i+1) % N][(j-1) % N] * l * l)
        for c in range( i, i+offset ):
            if is_cutpoint[c % N]:
                Z_comp1 = 1
                Z_comp2 = 1
                dZ_comp1 = 0
                dZ_comp2 = 0
                if c != i :
                    Z_comp1  = Z_linear[i+1][c % N]
                    dZ_comp1 = dZ_linear[i+1][c % N]
                if (c+1)%N != j:
                    Z_comp2 = Z_linear[(c+1) % N][j-1]
                    dZ_comp2 = dZ_linear[(c+1) % N][j-1]
                Z_product  = Z_comp1 * Z_comp2
                dZ_product = dZ_comp1 * Z_comp2 + Z_comp1 * dZ_comp2

                Z_BP[i][j]  += (C_std/Kd_BP) * (l/l_BP) * Z_product
                dZ_BP[i][j] += (C_std/Kd_BP) * (l/l_BP) * dZ_product
    return

def update_C_eff( self, i, j ):
    (Kd_BP, C_init, l, l_BP, C_std, C_init_BP, min_loop_length, N, sequence, is_cutpoint, any_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear ) = unpack_variables( self )
    offset = ( j - i ) % N

    # key 'special sauce' for derivative w.r.t. Kd_BP
    dZ_BP[i][j] += -(1.0/Kd_BP) * Z_BP[i][j]

    if not is_cutpoint[(j-1) % N]:
        C_eff[i][j]  += C_eff[i][(j-1) % N] * l
        dC_eff[i][j] += dC_eff[i][(j-1) % N] * l

    C_eff[i][j]  += C_init_BP * Z_BP[i][j]
    dC_eff[i][j] += C_init_BP * dZ_BP[i][j]

    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            C_eff[i][j]  += C_eff[i][(k-1) % N] * Z_BP[k % N][j] * l_BP
            dC_eff[i][j] += ( dC_eff[i][(k-1) % N] * Z_BP[k % N][j] + C_eff[i][(k-1) % N] * dZ_BP[k % N][j] ) * l_BP

def update_Z_linear( self, i, j ):
    (Kd_BP, C_init, l, l_BP, C_std, C_init_BP, min_loop_length, N, sequence, is_cutpoint, any_cutpoint, Z_BP, dZ_BP, C_eff, dC_eff, Z_linear, dZ_linear ) = unpack_variables( self )
    offset = ( j - i ) % N

    if not is_cutpoint[(j-1) % N]:
        Z_linear[i][j]  += Z_linear[i][(j - 1) % N]
        dZ_linear[i][j] += dZ_linear[i][(j - 1) % N]

    Z_linear[i][j]  += Z_BP[i][j]
    dZ_linear[i][j] += dZ_BP[i][j]

    for k in range( i+1, i+offset):
        if not is_cutpoint[ (k-1) % N]:
            Z_linear[i][j]  += Z_linear[i][(k-1) % N] * Z_BP[k % N][j]
            dZ_linear[i][j] += ( dZ_linear[i][(k-1) % N] * Z_BP[k % N][j] + Z_linear[i][(k-1) % N] * dZ_BP[k % N][j] )

def unpack_variables( self ):
    '''
    This helper function just lets me write out equations without
    using "self" which obscures connection to my handwritten equations
    '''
    return (self.params.Kd_BP, self.params.C_init, self.params.l, self.params.l_BP, self.params.C_std, self.params.C_init_BP, \
            self.params.min_loop_length, \
            self.N, self.sequence, self.is_cutpoint, self.any_cutpoint,  \
            self.Z_BP, self.dZ_BP, self.C_eff, self.dC_eff, self.Z_linear, self.dZ_linear )

