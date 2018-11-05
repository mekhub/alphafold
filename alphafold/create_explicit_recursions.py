#!/usr/bin/python
with open('recursions.py') as f:
    lines = f.readlines()

lines_new = []
for line in lines:
    line_new = ''
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
        if char == ' ' and not in_bracket:
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
            #bracket_word += char
            if not words[-1].replace('(','') in ['self.Z_BPq','sequence']:
                line_new += '.data'
                line_new += '[('+arg[:-1]+')%N]'
            else:
                line_new += bracket_word
            args.append( arg[:-1] )
            if in_second_bracket:
                if not (words[-1] in ['all_ligated'] ):
                    all_args.append( (words[-1],args) )
                    args = []
                    line_new += '.Q'
            else:
                just_finished_first_bracket = True
                if words[-1] == 'Z_final':
                    line_new += '.Q'
            in_bracket = False
            in_second_bracket = False
            bracket_word = ''
        else:
            if not in_bracket: line_new += char
            word += char
            just_finished_first_bracket = False

    if len(all_args)>0:
        print line,
        print line_new,
        print all_args

    line_new = line_new.replace( 'DynamicProgrammingData( 1.0, options = self.options )','1.0' )
    lines_new.append( line_new )

with open('explicit_recursions.py','w') as f:
    f.writelines( lines_new )

