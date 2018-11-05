#!/usr/bin/python
with open('recursions.py') as f:
    lines = f.readlines()


not_data_objects = ['self.Z_BPq','sequence']
not_2D_dynamic_programming_objects = ['all_ligated']
dynamic_programming_lists = ['Z_final']

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
        if (char == ' ' or char == '\n') and not in_bracket:
            if word in ['Z_seg1','Z_seg2']: line_new += '.Q'
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

    # temporary hack.
    line_new = line_new.replace( '.Q  = DynamicProgrammingData',' = DynamicProgrammingData' )
    lines_new.append( line_new )

    line_is_assignment =  line.count('+=') or line.count(' = ' )

    if len(all_args)>0:
        print line,
        print line_new,
        print all_args, ' line_is_assignment: ', line_is_assignment
        print


with open('explicit_recursions.py','w') as f:
    f.writelines( lines_new )

