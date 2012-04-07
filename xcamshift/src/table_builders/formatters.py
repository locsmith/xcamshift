'''
Created on 7 Apr 2012

@author: garyt
'''
import re




def fixup_null_values(line):
    line = line.replace('~,', " ,%10s" % "")
    line = line.replace('~}', " ,%9s}" % "")
    return line

def fixup_key_spacing_callback(m):
    value  =  m.group(1)[:-1]
    length = len(value)
    return value + ((4 - length) * " ") + ": "

def fixup_key_spacing(line):
    line = re.sub("([^{][A-Z]+:)", fixup_key_spacing_callback, line)
    return line

def fixup_decimal_spacing_callback(m):
    reduce_spaces_by  = len(m.group(3)) - 1
    num_spaces =  len(m.group(1))
    return " " * (num_spaces - reduce_spaces_by) +  m.group(2) +  m.group(3) + "."

def fixup_decimal_spacing(line):
    line = re.sub("(\s+)([+-])([0-9]+)\.", fixup_decimal_spacing_callback, line)
    return line


def fixup_convert_H_to_HN(line):
    line = re.sub("H ", "HN", line)
    return line


def fixup_replace_plus_with_space(line):
    line = re.sub("\+", " ", line)
    return line

def fixup_complex_key_question_mark(line):
    return line.replace("? [", "[")

def fixup_tuple_key_spacing_callback(m):
    values  =  m.group(1).split(',')
    value_1 =  ("%s," % values[0]).ljust(3)
    value_2 = ("%s" % values[1]).rjust(3)

    return '[%s%s]' % (value_1,value_2)

def fixup_tuple_key_spacing(line):
    return re.sub("\[(.*)\]", fixup_tuple_key_spacing_callback, line)


def fixup_put_lonely_keys_on_new_line(line):
    return re.sub("^(\s+\[.*\]:)$", "\n\g<1>", line)

def global_fixup_colons_on_same_line_as_tuple_key(lines):
    return re.sub("\]\n(\s*):", "]:\n\g<1> ", lines)


def global_fixup_data_dicts_on_same_line_as_key(lines):
    return re.sub("\]:\n\s+{", "] : {", lines)
