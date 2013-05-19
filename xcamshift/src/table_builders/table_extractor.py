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
from table_builders.Table_modifier import Table_modifier
'''
Created on 7 Apr 2012

@author: garyt
'''
from abc import abstractmethod, ABCMeta
import yaml
from yaml_patches import apply_ordered_dict_patch, apply_patch_float_format_with_nulls, \
                         apply_tuple_patch,apply_no_aliases_patch

DEFAULT_INDENT = 6

class Table_extractor(object):
    
    __metaclass__ = ABCMeta
    
    def __init__(self,data):
        self._data = data
        self._modifier = Table_modifier({})
        
                
        apply_patch_float_format_with_nulls()
        apply_tuple_patch()
        apply_ordered_dict_patch()
        apply_no_aliases_patch()
        
    
    def is_table_required(self,table_residue_type):
        return True
    
    def extract(self,file_type=''):
        data =  self._data[file_type]
        
        serialized_data = self.serialize(data)
        
        serialized_data = self.modify(serialized_data, self._modifier)
        
        lines = self.build_output_lines(serialized_data)
        
        lines =  self.format_lines(lines)
        
        return "\n".join(lines)
    
    def set_modifier(self,modifier):
        self._modifier = modifier
    
    def modify(self,serialized_data,modifier):
        return modifier.run(serialized_data)
    
    def build_output_lines(self,serialized_data):
        return  yaml.dump(serialized_data,default_flow_style=None,width=1000,indent=DEFAULT_INDENT)
    
    @abstractmethod
    def serialize(self,data):
        pass
    
    @abstractmethod
    def format_lines(self,lines):
        return lines
    
    @classmethod
    def get_name(self):
        raise Exception("name of extractor must be defined!")
