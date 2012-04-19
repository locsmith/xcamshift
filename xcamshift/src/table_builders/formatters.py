'''
Created on 7 Apr 2012

@author: garyt
'''
import re

from functools import partial

#TODO: extract common code


def fixup_null_values(line):
    line = line.replace('~,', " ,%10s" % "")
    line = line.replace('~}', " ,%9s}" % "")
    return line

def fixup_key_spacing_callback(m,spacing=4):
    value  =  m.group(1)[:-1]
    length = len(value)
    return value + ((spacing - length) * " ") + ": "

KEY_AND_VALUE_REGEX = "^([ ]+)([A-Z0-9][A-Z0-9 :]+){([ \-+A-Z0-9:,\.]+)}$"

def fixup_data_key_spacing(line,spacing = 4):
    matches = re.match(KEY_AND_VALUE_REGEX, line)
    if matches:
        pre_space = matches.group(1)
        key = matches.group(2)
        values = matches.group(3)
        partial_fixup = partial(fixup_key_spacing_callback, spacing=spacing)
        values = re.sub("([A-Z0-9]+:)", partial_fixup, values)
        line  = "%s%s{%s}" % (pre_space,key,values)
    return line

def fixup_line_data_key_spacing(line,spacing = 6):
    matches = re.match(KEY_AND_VALUE_REGEX, line)
    if matches:
        pre_space = matches.group(1)
        key = matches.group(2)
        values = matches.group(3)
        
        key_length = len(key)-1
        if (key_length) + 1 < spacing:
            key = key.ljust(spacing)
        line  = "%s%s{%s}" % (pre_space,key,values)
    return line

def fixup_decimal_spacing_callback(m):
    reduce_spaces_by  = len(m.group(3)) - 1
    num_spaces =  len(m.group(1))
    return " " * (num_spaces - reduce_spaces_by) +  m.group(2) +  m.group(3) + "."

def fixup_decimal_spacing(line):
    line = re.sub("(\s+)([+-])([0-9]+)\.", fixup_decimal_spacing_callback, line)
    return line


def fixup_convert_H_to_HN(line):
    line = re.sub(" H  ", " HN ", line)
    line = re.sub(" H:  ", " HN: ", line)
    line = line.replace('[H, ','[HN,')
    return line


def fixup_convert_H_to_HN_not_in_key(line):
    line = re.sub("H  ", "HN ", line)
    line = re.sub("H:  ", "HN: ", line)
    line = line.replace('[H, ','[H,  ')
    return line

def fixup_replace_plus_with_space(line):
    line = re.sub("\+", " ", line)
    return line

def fixup_complex_key_question_mark(line):
    line = re.sub("\?([\s\-]+)\[", " \g<1>[",line)
    return line
#line.replace("? [", "[")

def fixup_tuple_key_spacing_callback(m,left_spacing=3,right_spacing=3):
    values  =  m.group(1).split(',')
    value_1 =  ("%s," % values[0]).ljust(left_spacing)
    value_2 = ("%s" % values[1]).rjust(right_spacing)

    return '[%s%s]' % (value_1,value_2)

def fixup_tuple_key_spacing(line,left_spacing=3,right_spacing=3):
    callback = partial(fixup_tuple_key_spacing_callback, left_spacing = left_spacing, right_spacing= right_spacing)
    return re.sub("\[([A-Z0-9\-\,\s]+)\]", callback, line)

def fixup_put_lonely_keys_on_new_line(line):
    return re.sub("^(\s+[A-Za-z0-9\-\[][A-Za-z-0-9_, \]]+:\s*)$", "\n\g<1>", line)

def global_fixup_colons_on_same_line_as_tuple_key(lines):
    return re.sub("\]\n(\s*):", "]:\n\g<1> ", lines)


def global_fixup_data_dicts_on_same_line_as_key(lines):
    return re.sub("\]:\n\s+{", "] : {", lines)

def fixup_spaces_after_colons(line):
        return line.replace(':',':  ')
    
#def fixup_lonely_key_spacing(line, length=2):
#    match = re.match("(\s+)([0-9\-A-Z]+):$",line)
#    
#    if match != None:
#        line_spacing =  match.group(1)
#        key = match.group(2)
#        key = key.rjust(length)
#        line  = "%s%s:"  % (line_spacing,key)
#    return line

def global_fixup_multi_line_tuple_keys(lines):
    lines = re.sub(" - (\[[A-Z\,\s\-0-9]+\])", "\g<1>", lines)
    lines = lines.replace("]\n", "], ")
    lines = lines.replace("] :", "]] :")
    lines = re.sub("(\s+)\[", "\g<1>[[", lines)
    lines = re.sub("\]\,\s+\[", "], [", lines)
    lines = re.sub("\,\s\[\[",", [",lines)
    return lines