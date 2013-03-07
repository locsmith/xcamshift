#-------------------------------------------------------------------------------
# Copyright (c) 2013 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 13 Feb 2013

@author: garyt
'''

from common_constants import DATA, h_keys#, HA,CA,H,N,C,CB,O, h_keys,  XTRA
from collection_backport import OrderedDict
from table_builders.formatters import fixup_null_values,\
    fixup_decimal_spacing, fixup_convert_H_to_HN, fixup_replace_plus_with_space,\
    fixup_tuple_key_spacing, fixup_put_lonely_keys_on_new_line,\
    fixup_complex_key_question_mark,\
    global_fixup_colons_on_same_line_as_tuple_key,\
    global_fixup_data_dicts_on_same_line_as_key, fixup_spaces_after_colons
from table_builders.table_extractor import Table_extractor



XTRA_DATA = 'XTRADISTS'
DISU = 'DISU'
class DISU_table_extractor(Table_extractor):
    
    def __init__(self,data):
        super(DISU_table_extractor, self).__init__(data)
        
    
    @classmethod
    def get_name(self):
        return DISU

    def is_table_required(self,table_residue_type):
        result = False
        
        if table_residue_type == '':
            result = True
        
        return result
    
    def serialize(self,data):
        out_data = OrderedDict()
        out_data[DATA] =  OrderedDict()
        
        out_line = out_data[DATA].setdefault(('DISU','CYS'),OrderedDict())
                
        for h_key in h_keys:
            out_line[h_key] = out_line[h_key] = data[XTRA_DATA][h_key][-1]

        return out_data
    
    def format_lines(self,lines):
        result = []
        lines = global_fixup_colons_on_same_line_as_tuple_key(lines)
        lines = global_fixup_data_dicts_on_same_line_as_key(lines)
        for line in lines.split('\n'):
            line = fixup_spaces_after_colons(line)
            line = fixup_complex_key_question_mark(line)
            line = fixup_null_values(line)
            line = fixup_tuple_key_spacing(line)
            line = fixup_decimal_spacing(line)
            line = fixup_replace_plus_with_space(line)
            line = fixup_convert_H_to_HN(line)
            result.append(line)
        return result
    
