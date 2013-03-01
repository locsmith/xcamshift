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
Created on 6 Apr 2012

@author: garyt
'''
import re
from common_constants import HA,CA,H,N,C,O, DATA, HBOND
from common_constants import h_keys
#
from collection_backport import OrderedDict

    
from table_builders.formatters import fixup_null_values,\
    fixup_decimal_spacing, fixup_convert_H_to_HN, fixup_replace_plus_with_space,\
    fixup_tuple_key_spacing, fixup_put_lonely_keys_on_new_line,\
    fixup_complex_key_question_mark,\
    global_fixup_colons_on_same_line_as_tuple_key,\
    global_fixup_data_dicts_on_same_line_as_key, fixup_spaces_after_colons
from ..yaml_patches import apply_ordered_dict_patch, apply_patch_float_format_with_nulls, \
    apply_tuple_patch

from table_builders.table_extractor import Table_extractor

HBOND_V_KEYS_1 = (-1, 'O', 'HN'), (0, 'HN', 'O'), (0, 'O', 'HN'), (1, 'HN', 'O')
HBOND_V_KEYS_2 = 'DIST', 'ANG1', 'ANG2'        


class HBOND_table_extractor(Table_extractor):
    
    
    def __init__(self,data):

        self._data = data
        
        apply_patch_float_format_with_nulls()
        apply_tuple_patch()
        apply_ordered_dict_patch()

    @classmethod
    def get_name(self):
        return HBOND

    def extract(self, file_type  = ''):
        
        data =  self._data[file_type]
        
        serialized_data = self.serialize(data)
        
        lines = self.build_output_lines(serialized_data)
        
        lines =  self.format_lines(lines)
        
        return "\n".join(lines)

    def serialize(self,data):
        out_data = OrderedDict()
        out_data[DATA] = OrderedDict()
                
        for i,v_key in enumerate(HBOND_V_KEYS_1):
            out_data[DATA][v_key] = OrderedDict()
            for j,parameter in enumerate(HBOND_V_KEYS_2):
                out_data[DATA][v_key][parameter] = OrderedDict()
                for h_key in h_keys:
                    out_data[DATA][v_key][parameter][h_key] = data['HBONDS'][h_key][i*3+j]
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