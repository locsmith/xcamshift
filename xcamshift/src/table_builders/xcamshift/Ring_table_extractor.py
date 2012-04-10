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


TYR='TYR'
PHE='PHE'
TRP='TRP'
HIS='HIS'

xtra_v_keys = (  
    (PHE,6),
    (TYR,6),
    (TRP,5),
    (TRP,6),
    (HIS,6)
)

RINGS = 'RINGS'
class RING_table_extractor(Table_extractor):
    
    def __init__(self,data):
        self._data =  data
        
        apply_patch_float_format_with_nulls()
        apply_tuple_patch()
        apply_ordered_dict_patch()

    def serialize(self,data):
        out_data = OrderedDict()
        out_data['data'] =  OrderedDict()
                
        for i,v_key in enumerate(xtra_v_keys):
    
            out_line = out_data['data'].setdefault(v_key,OrderedDict())
    
            for h_key in h_keys:
                out_line[h_key] = data[RINGS][h_key][i]

        return out_data
    

    def _get_data(self, file_type=''):
        return super(RING_table_extractor, self)._get_data(self, file_type=file_type)

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
    
