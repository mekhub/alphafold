def secstruct( bps, N ):
    secstruct_string = ['.']*N
    for bp in bps:
        i = min( bp )
        j = max( bp )
        secstruct_string[i] = '('
        secstruct_string[j] = ')'
    return ''.join(secstruct_string)

def bps( secstruct ):
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

