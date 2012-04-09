'''
Created on 7 Apr 2012

@author: garyt
'''

from table_builders.common_constants import HA,CA,H,N,C,CB,O, h_keys
from ..yaml_patches import apply_ordered_dict_patch, apply_patch_float_format_with_nulls, \
     apply_tuple_patch
from collection_backport import OrderedDict
from table_builders.formatters import fixup_null_values,\
    fixup_decimal_spacing, fixup_convert_H_to_HN, fixup_replace_plus_with_space,\
    fixup_tuple_key_spacing, fixup_put_lonely_keys_on_new_line,\
    fixup_complex_key_question_mark,\
    global_fixup_colons_on_same_line_as_tuple_key,\
    global_fixup_data_dicts_on_same_line_as_key
from table_builders.table_extractor import Table_extractor

CG = 'CG'

xtra_v_keys = (   (H,0,HA,0),
                  (H,0,C,0),
                  (H,0,CB,0),
                  
                  (C,-1,HA,0),
                  (C,-1,C,0),
                  (C,-1,CB,0),
                  
                  (O,0,HA,0),
                  (O,0,N,0),
                  (O,0,CB,0),
                  
                  (N,1,HA,0),
                  (N,1,HA,0),
                  (N,1,HA,0),
                  
                  (O,-1,HA,-1),
                  (O,-1,N,-1),
                  (O,-1,CB,-1),
                  
                  (N,0,HA,-1),
                  (N,0,N,-1),
                  (N,0,CB,-1),
                  
                  (CG,0,HA,0),
                  (CG,0,N,0),
                  (CG,0,C,0),
                  
                  (CG,0,C,-1),
                  (CG,0,N,+1),
                  
                  (CG,-1,CA,0),
                  (CG,1,CA,0),
                  
                  (CA,-1,CA,1)
)

XTRA = 'XTRADISTS'
class XTRA_table_extractor(Table_extractor):
    
    def __init__(self,data):
        self._data =  data
        
        apply_patch_float_format_with_nulls()
        apply_tuple_patch()
        apply_ordered_dict_patch()

    def serialize(self,data):
        out_data = OrderedDict()
        out_data['data'] =  OrderedDict()
                
        for i,v_key in enumerate(xtra_v_keys):
    
            out_line = out_data['data'].setdefault(v_key[0:2],OrderedDict()).setdefault(v_key[2:4],OrderedDict())
    
            for h_key in h_keys:
                out_line[h_key] = data[XTRA][h_key][i]

        return out_data
    

    def _get_data(self, file_type=''):
        return super(XTRA_table_extractor, self)._get_data(self, file_type=file_type)

    def format_lines(self,lines):
        result = []
        lines = global_fixup_colons_on_same_line_as_tuple_key(lines)
        lines = global_fixup_data_dicts_on_same_line_as_key(lines)
        for line in lines.split('\n'):
            line = fixup_complex_key_question_mark(line)
            line = fixup_null_values(line)
            line = fixup_tuple_key_spacing(line)
            line = fixup_decimal_spacing(line)
            line = fixup_replace_plus_with_space(line)
            line = fixup_convert_H_to_HN(line)
            line = fixup_put_lonely_keys_on_new_line(line)
            result.append(line)
        return result
    
