#-------------------------------------------------------------------------------
# Copyright (c) 2013-2015 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 7 Apr 2012

@author: garyt
'''

from common_constants import HA,CA,H,N,C,CB,O, h_keys, DATA, XTRA
from collection_backport import OrderedDict
from table_builders.formatters import fixup_null_values,\
    fixup_decimal_spacing, fixup_convert_H_to_HN, fixup_replace_plus_with_space,\
    fixup_tuple_key_spacing, fixup_put_lonely_keys_on_new_line,\
    fixup_complex_key_question_mark,\
    global_fixup_colons_on_same_line_as_tuple_key,\
    global_fixup_data_dicts_on_same_line_as_key, fixup_spaces_after_colons
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
                  (N,1,N,0),
                  (N,1,CB,0),
                  
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

XTRA_DATA = 'XTRADISTS'
class XTRA_table_extractor(Table_extractor):
    
    def __init__(self,data):
        super(XTRA_table_extractor, self).__init__(data)
        
    
    @classmethod
    def get_name(self):
        return XTRA

    def serialize(self,data):
        out_data = OrderedDict()
        out_data[DATA] =  OrderedDict()
                
        for i,v_key in enumerate(xtra_v_keys):
    
            out_line = out_data[DATA].setdefault(v_key[0:2],OrderedDict()).setdefault(v_key[2:4],OrderedDict())
    
            for h_key in h_keys:
                out_line[h_key] = data[XTRA_DATA][h_key][i]

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
            line = fixup_put_lonely_keys_on_new_line(line)
            result.append(line)
        return result
    
