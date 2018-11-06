#!/usr/bin/python
with open('recursions.py') as f:
    lines = f.readlines()


not_data_objects = ['self.Z_BPq','sequence']
not_2D_dynamic_programming_objects = ['all_ligated']
dynamic_programming_lists = ['Z_final']
dynamic_programming_data = ['Z_seg1','Z_seg2']

def find_substring(substring, string):
    """
    From stackoverflow...
    Returns list of indices where substring begins in string

    >>> find_substring('me', "The cat says meow, meow")
    [13, 19]
    """
    indices = []
    index = -1  # Begin at -1 so index + 1 is 0
    while True:
        # Find next index of substring, by starting search from index + 1
        index = string.find(substring, index + 1)
        if index == -1:
            break  # All occurrences have been found
        indices.append(index)
    return indices

lines_new = []
for line in lines:
    line_new = ''

    if line.count( '.dQ' ) or line.count( '.Q') : # if explicitly defining Q, dQ already, don't try to override
        lines_new.append( line )
        continue

    # most important thing -- need to look for get/set of Z_BP[i][j] (DynamicProgrammingMatrix)
    in_bracket = False
    in_second_bracket = False
    just_finished_first_bracket = False
    all_args = []
    args = []
    words = []
    word = ''
    bracket_word = ''
    at_beginning = True
    num_indent = 0
    for char in line:
        if (char == ' ' or char == '\n') and not in_bracket:
            if at_beginning: num_indent += 1
            if word in dynamic_programming_data: line_new += '.Q'
            line_new += char
            if len( word ) > 0:
                words.append( word )
                word = ''
            continue
        else:
            at_beginning = False
        if in_bracket:
            bracket_word += char
            arg += char
        if char == '[':
            bracket_word += char
            arg = ''
            in_bracket = True
            if just_finished_first_bracket:
                in_second_bracket = True
            else:
                args = []
                if len( word ) > 0:
                    words.append( word )
                    word = ''
        elif char == ']':

            # OMG this is so hacky
            if not words[-1].replace('(','') in not_data_objects:
                line_new += '.data'
                if len(arg[:-1]) == 1:
                    line_new += '['+arg[:-1]+'%N]'
                else:
                    line_new += '[('+arg[:-1]+')%N]'
            else:
                line_new += bracket_word
            args.append( arg[:-1] )
            if in_second_bracket:
                if not (words[-1] in not_2D_dynamic_programming_objects ):
                    assert( len( args ) == 2 )
                    all_args.append( (len(line_new),words[-1],args[0],args[1]) )
                    args = []
                    line_new += '.Q'
            else:
                just_finished_first_bracket = True
                if words[-1] in dynamic_programming_lists:
                    line_new += '.Q'
            in_bracket = False
            in_second_bracket = False
            bracket_word = ''
        else:
            if not in_bracket: line_new += char
            word += char
            just_finished_first_bracket = False

    # temporary hack -- this is for Z_seg1, Z_seg2 assignment...
    line_new = line_new.replace( '.Q  = DynamicProgrammingData',' = DynamicProgrammingData' )

    lines_new.append( line_new )

    if line == line_new: continue
    print line,
    print line_new,

    # is this an assignment? then need to create derivative and contribution lines
    assign_pos = line_new.find('+= ')
    if ( assign_pos < 0 ): assign_pos = line_new.find(' = ' )
    if assign_pos > -1:
        Qpos = find_substring( '.Q', line_new )
        if len( Qpos ) > 0 and Qpos[0] < assign_pos:
            assert( len( Qpos ) > 1 )

            lines_new.append( ' '*num_indent + 'if self.options.calc_deriv:\n' )
            print lines_new[-1],
            line_beginning = ' '*4 + line_new[:Qpos[0]] + '.dQ' # extra indent
            # deriv line
            for i in range( 1, len( Qpos ) ):
                line_deriv = line_beginning
                assert( Qpos[i] > assign_pos )
                for j in range(1, len( Qpos) ):
                    line_deriv += line_new[Qpos[j-1]+2 : Qpos[j]]
                    if i == j: line_deriv += '.dQ'
                    else: line_deriv += '.Q'
                line_deriv += line_new[ Qpos[-1]+2:]
                print line_deriv,
                lines_new.append(line_deriv)

            # contrib line
            lines_new.append( ' '*num_indent + 'if self.options.calc_contrib:\n' )
            print lines_new[-1],
            line_contrib = ' '*4 + line_new[:Qpos[0]] +  '.contribs.append( ( ' # extra indent
            line_contrib += line_new[assign_pos+3:-1] + ', ['
            for (n,info) in enumerate(all_args):
                if info[ 0 ] <= assign_pos: continue
                line_contrib += '(%s,' % info[1]
                if len(info[2])> 1:  line_contrib += '(%s)%%N' % info[2]
                else: line_contrib += '%s%%N' % info[2]
                line_contrib+=','
                if len(info[3])>1:   line_contrib += '(%s)%%N' % info[3]
                else: line_contrib += '%s%%N' % info[3]
                line_contrib+=')'
                if n < len( all_args )-1: line_contrib += ', '
            line_contrib += '] ) )\n'
            print line_contrib,
            lines_new.append( line_contrib)
    print


with open('explicit_recursions.py','w') as f:
    f.writelines( lines_new )

