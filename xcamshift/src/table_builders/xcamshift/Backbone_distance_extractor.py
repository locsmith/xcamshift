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
Created on 6 Apr 2012

@author: garyt
'''
import re
from common_constants import HA,CA,H,N,C,O, DATA, BACK_BONE
from common_constants import h_keys

from collection_backport import OrderedDict

from table_builders.formatters import fixup_convert_H_to_HN,fixup_decimal_spacing, \
     fixup_data_key_spacing, fixup_null_values, fixup_replace_plus_with_space,\
    fixup_put_lonely_keys_on_new_line, fixup_line_data_key_spacing
from ..yaml_patches import apply_ordered_dict_patch, apply_patch_float_format_with_nulls

from table_builders.table_extractor import Table_extractor



BACK_BONE_DATA = 'COBB2'

missing_v_keys = (H, -1),(O, 1)
bb_v_keys = ((N,-1),
          (H, -1),
          (CA,-1),
          (HA,-1),
          (C,-1),
          (O,-1),
          
          (N,0),
          (H,0),
          (CA,0),
          (HA,0),
          (C,0),
          (O,0),
          
          (N,1),
          (H,1),
          (CA,1),
          (HA,1),
          (C,1),
          (O,1)
          )



class BB_table_extractor(Table_extractor):
    
    
    def __init__(self,data):

        self._data = data
        
        apply_patch_float_format_with_nulls()
        apply_ordered_dict_patch()
        
    @classmethod
    def get_name(self):
        return BACK_BONE
    
        
    def extract(self, file_type  = ''):
        
        data =  self._data[file_type]
        
        serialized_data = self.serialize(data)
        
        lines = self.build_output_lines(serialized_data)
        
        lines =  self.format_lines(lines)
        
        return "\n".join(lines)
        

    def add_missing_keys_to_copy_of_data(self, data):
        data = dict(data)
        for i, v_key in enumerate(bb_v_keys):
            if v_key in missing_v_keys:
                for h_key in h_keys:
                    data[BACK_BONE_DATA][h_key].insert(i, 0.0)
        
        return data

    def serialize(self,data):
        
        data = self.add_missing_keys_to_copy_of_data(data) 
                
        out_data = OrderedDict()
        out_data[DATA] = OrderedDict()
                    
        for i,v_key in enumerate(bb_v_keys):
    
    
            out_line = out_data[DATA].setdefault(v_key[1],OrderedDict()).setdefault(v_key[0],OrderedDict())

            for h_key in h_keys:
                out_line[h_key] = data[BACK_BONE_DATA][h_key][i]
        
        return out_data

            
    def format_lines(self,lines):
        result = []
        for line in lines.split('\n'):
            line = fixup_null_values(line)
            line = fixup_data_key_spacing(line)
            line = fixup_line_data_key_spacing(line)
            line = fixup_decimal_spacing(line)
            line = re.sub("^(\s+[0-9])","\n\g<0>",line)
            line = fixup_replace_plus_with_space(line)
            line = fixup_convert_H_to_HN(line)
#            line = fixup_lonely_key_spacing(line)
            line = fixup_put_lonely_keys_on_new_line(line)
            result.append(line)
        return result
