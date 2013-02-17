#-------------------------------------------------------------------------------
# Copyright (c) 2013 gary thompson.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 24 Apr 2012

@author: garyt












'''
import sys
import re
from collection_backport import OrderedDict
from itertools import ifilter
from functools import partial
import yaml

class Camshift_output_parser(object):
    
    def __init__(self, source):
        self._source  = source

        
    def _get_num_residues_or_raise(self,lines):
        residues = None
        for line in lines:
            line = line.strip()
            if re.match('>> Number of fragments', line):
                elems = line.split()
                residues = int(elems[4])
                
        if residues ==  None:
            raise Exception("can't find a line of the form '>> Number of fragments:   3'")

            
        return residues
    
    def _open(self):
        return open(self._source)
    
    def _close(self,file_h):
        file_h.close()
        
    def parse(self):
        
        file_h  =  self._open()
        
        result = None
        try:
            lines = file_h.readlines()
            
            number_residues = self._get_num_residues_or_raise(lines)
            
            shifts = lines[-number_residues:]
            header = lines[-number_residues - 2]
            atom_types = header.strip().split()[2:]
            
            result = OrderedDict()
            for shift in shifts:
                elems = shift.strip().split()
                segid = ''
                residue_number = elems[0]
                for atom, shift_string in zip(atom_types, elems[2:]):
                    key = segid, int(residue_number), atom
                    result[key] = float(shift_string)
        finally:
            self._close(file_h)
            
        return result

def load_constants():
    constants_h = open('data/cams_1_35_0_constants_base.yaml')
    return yaml.load(constants_h)

def get_first_and_last(data):
    residue_numbers = set()
    
    for segid,residue_number,atom in data:
        residue_numbers.add(residue_number)
    residue_numbers = list(residue_numbers)
    residue_numbers.sort()
    
    return residue_numbers[0],residue_numbers[-1]

def remove_first_and_last_by_residue(data):
    keys = data.keys()
    keys.sort()
    first_and_last = get_first_and_last(data)
    
    residue_filter =  partial(ifilter, lambda x:not x[1] in first_and_last)
    
    return OrderedDict([(key, data[key]) for key in residue_filter(keys)])

def format_as_atom_selection_float_dict(data, name = 'test', float_format = '9.5f'):
    data_keys = []
    values = []
    for data_key in data:
        data_keys.append("    ('%s',%i,'%s')" % data_key)
        values.append(data[data_key])
    
    max_width = reduce(get_longest_length, data_keys, 0)
    string_format = "%%-%is  :  %%%s," % (max_width, float_format)
    
    result = []
    
    result.append( '%s  = {' % name)
    for data_key_value in zip(data_keys, values):
        result.append(string_format % data_key_value)
    result.append('}')
    
    return '\n'.join(result)

def get_longest_length(accumulator ,value):
    result = accumulator
    
    value_length = len(value)
    
    if value_length > accumulator: 
        result = value_length
    
    return result 
    

def offset_by_well(data, constants, scale = 1.0, epsilon = 0.1):
    for key in data:
        segid, residue_num, atom_name = key
        
        offset = constants['flat_bottom_limit'][atom_name] * scale + epsilon
        data[key] += offset
        
    return data


if __name__ == '__main__':
    
    camshift_parser = Camshift_output_parser(sys.argv[1])

    data  = camshift_parser.parse()
    data  = remove_first_and_last_by_residue(data)
    
    constants = load_constants()
    
#    data  = offset_by_well(data,constants)
    
    print format_as_atom_selection_float_dict(data)
            
        
