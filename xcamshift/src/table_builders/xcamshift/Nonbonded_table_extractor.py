'''
Created on 7 Apr 2012

@author: garyt
'''

from table_builders.common_constants import H,N,C,O, h_keys
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
import re
from table_builders.yaml_patches import apply_no_aliases_patch

S='S'

SP3='SP3'
SP2='SP2'
nonbonded_v_keys = (  
               
    
    (C, SP3),
    (H,'NONE'),
    (N, SP3),
    (O, SP3),
    (S, SP3),
    
    (C, SP2),
    (N, SP2),
    (O, SP2),
)

nonbonded_v_keys_required_order = (
    (H,'NONE'),
    
    (C, SP3),
    (N, SP3),
    (O, SP3),
    (S, SP3),
    
    (C, SP2),
    (N, SP2),
    (O, SP2),                                   
)

NONBONDED = 'nonbonded'
SPHERE1 = 'SPHERE1'
SPHERE2 = 'SPHERE2'

DATA = 'data'
EXPONENT_ID = 'exponent'
SPHERE_ID_1 = 'sphere_1'
SPHERE_ID_2 = 'sphere_2'

SPHERE_1_EXPONENT = -3.0
SPHERE_2_EXPONENT =  1.0
        
def ignore_aliases(self,_data):
    return True

class Nonbonded_table_extractor(Table_extractor):
    
    def __init__(self,data):
        self._data =  data
        
        apply_patch_float_format_with_nulls()
        apply_tuple_patch()
        apply_ordered_dict_patch()
        apply_no_aliases_patch()


    def serialize(self,data):
        
        raw_out_data = self._basic_serialise(data)
        
        return self._tidy_key_order(raw_out_data)
        

    def _build_outline(self, data, i,sphere):
        out_line = OrderedDict()
        for h_key in h_keys:
            out_line[h_key] = data[sphere][h_key][i]
            
        return out_line

    def _basic_serialise(self, data):
        raw_out_data = {}
        raw_out_data[DATA] = {}
        
        for i, v_key in enumerate(nonbonded_v_keys):
            
            out_line = self._build_outline(data, i, SPHERE1)
            raw_out_data[DATA].setdefault(SPHERE_ID_1, OrderedDict())[v_key] = out_line
            
            out_line = self._build_outline(data, i, SPHERE2)
            raw_out_data[DATA].setdefault(SPHERE_ID_2, OrderedDict())[v_key] = out_line

        
        return raw_out_data


    def _tidy_key_order(self, raw_out_data):
        out_data = OrderedDict()
        out_data[DATA] = OrderedDict()
        out_data[DATA][SPHERE_ID_1] = OrderedDict()
        out_data[DATA][SPHERE_ID_2] = OrderedDict()
        
        out_data[DATA][SPHERE_ID_1][EXPONENT_ID] = SPHERE_1_EXPONENT
        out_data[DATA][SPHERE_ID_2][EXPONENT_ID] = SPHERE_2_EXPONENT
        for v_key in nonbonded_v_keys_required_order:
            out_data[DATA].setdefault(SPHERE_ID_1, OrderedDict())[v_key] = raw_out_data[DATA][SPHERE_ID_1][v_key]
            out_data[DATA].setdefault(SPHERE_ID_2, OrderedDict())[v_key] = raw_out_data[DATA][SPHERE_ID_2][v_key]
        
        return out_data

    

    def _get_data(self, file_type=''):
        return super(Nonbonded_table_extractor, self)._get_data(self, file_type=file_type)

    def fixup_CA_values(self,line):
        return re.sub("CA:([ \-0-9]{3,3}?)\.","CA: \g<1>.",line)
    

    def split_catergories(self, line):
        tests = ['[H,  NONE]','[C,   SP3]','[C,   SP2]']
        for test in tests:
            if line.strip().startswith(test):
                line = "\n" + line
        return line
    

    def indent_exponent(self, line):
        if line.strip().startswith('exponent'):
            line = '  '  + line 
        return line
    
    
    def fixup_multiple_zeros(self, line):
        return re.sub('\.[0]+$','.0',line)   
    
    def format_lines(self,lines):
        result = []
        lines = global_fixup_colons_on_same_line_as_tuple_key(lines)
        lines = global_fixup_data_dicts_on_same_line_as_key(lines)
        for line in lines.split('\n'):
            line = fixup_complex_key_question_mark(line)
            line = fixup_null_values(line)
            line = fixup_tuple_key_spacing(line,3,5)
            line = fixup_decimal_spacing(line)
            line = fixup_replace_plus_with_space(line)
            line = fixup_convert_H_to_HN(line)
            line = fixup_put_lonely_keys_on_new_line(line)
            line = self.fixup_CA_values(line)
            line = self.split_catergories(line)
            line = self.indent_exponent(line)
            line = self.fixup_multiple_zeros(line)

            result.append(line)
        return result
    