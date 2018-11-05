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

    if line.count( '.dQ' ) > 0: # if i explicitly defined dQ already
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
    for char in line:
        if (char == ' ' or char == '\n') and not in_bracket:
            if word in dynamic_programming_data: line_new += '.Q'
            line_new += char
            if len( word ) > 0:
                words.append( word )
                word = ''
            continue
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
                    all_args.append( (words[-1],args) )
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

    # is this an assignment? then need to create a derivative line
    assign_pos = line_new.find('+=')
    if ( assign_pos < 0 ): assign_pos = line_new.find(' = ' )
    if assign_pos > -1:
        Qpos = find_substring( '.Q', line_new )
        if len( Qpos ) > 0 and Qpos[0] < assign_pos:
            assert( len( Qpos ) > 1 )
            # deriv line
            for i in range( 1, len( Qpos ) ):
                line_deriv = line_new[:Qpos[0]] + '.dQ'
                assert( Qpos[i] > assign_pos )
                for j in range(1, len( Qpos) ):
                    line_deriv += line_new[Qpos[j-1]+2 : Qpos[j]]
                    if i == j: line_deriv += '.dQ'
                    else: line_deriv += '.Q'
                line_deriv += line_new[ Qpos[-1]+2:]
                print line_deriv,
                lines_new.append(line_deriv)

    print


with open('explicit_recursions.py','w') as f:
    f.writelines( lines_new )

