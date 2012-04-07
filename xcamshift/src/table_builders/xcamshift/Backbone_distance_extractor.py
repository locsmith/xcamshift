'''
Created on 6 Apr 2012

@author: garyt
'''
import re
from table_builders.common_constants import HA,CA,H,N,C,CB
import yaml
from collection_backport import OrderedDict
from abc import abstractmethod

from table_builders.formatters import convert_H_to_HN,fixup_decimal_spacing, \
     fixup_key_spacing, fixup_null_values, replace_plus_with_space
from table_builders.xcamshift.yaml_patches import apply_patch_float_format_with_nulls,\
     apply_ordered_dict_patch
from table_builders.table_extractor import Table_extractor



BB = 'COBB2'

O = 'O'

h_keys = (HA,CA,H,N,C,CB)
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

DEFAULT_INDENT = 6


class BB_table_extractor(Table_extractor):
    
    
    def __init__(self,data):

        self._data = data
        
        apply_patch_float_format_with_nulls()
        apply_ordered_dict_patch()
    
    def _get_data(self,file_type=''):
        return super(BB_table_extractor, self)._get_data(file_type)
        
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
                    data[BB][h_key].insert(i, 0.0)
        
        return data

    def serialize(self,data):
        
        data = self.add_missing_keys_to_copy_of_data(data) 
                
        out_data = OrderedDict()
        out_data['data'] = OrderedDict()
                    
        for i,v_key in enumerate(bb_v_keys):
    
    
            out_line = out_data['data'].setdefault(v_key[1],OrderedDict()).setdefault(v_key[0],OrderedDict())

            for h_key in h_keys:
                out_line[h_key] = data[BB][h_key][i]
        
        return out_data
    
    def build_output_lines(self,serialized_data):
        return  yaml.dump(serialized_data,default_flow_style=None,width=1000,indent=6)
   

    def format_lines(self,lines):
        result = []
        for line in lines.split('\n'):
            line = fixup_null_values(line)
            line = fixup_key_spacing(line)
            line = fixup_decimal_spacing(line)
            line = re.sub("^(\s+[0-9])","\n\g<0>",line)
            line = replace_plus_with_space(line)
            line = convert_H_to_HN(line)
            result.append(line)
        return result
