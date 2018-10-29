def secstruct( bps, N ):
    secstruct_string = ['.']*N
    for bp in bps:
        i = min( bp )
        j = max( bp )
        secstruct_string[i] = '('
        secstruct_string[j] = ')'
    return ''.join(secstruct_string)

