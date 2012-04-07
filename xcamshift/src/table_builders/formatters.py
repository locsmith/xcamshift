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


def convert_H_to_HN(line):
    line = re.sub("H ", "HN", line)
    return line


def replace_plus_with_space(line):
    line = re.sub("\+", " ", line)
    return line
