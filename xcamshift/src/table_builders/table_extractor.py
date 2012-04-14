'''
Created on 7 Apr 2012

@author: garyt
'''
from abc import abstractmethod, ABCMeta
import yaml

DEFAULT_INDENT = 6

class Table_extractor(object):
    
    __metaclass__ = ABCMeta
         
    def extract(self,file_type=''):
        data =  self._data[file_type]
        
        serialized_data = self.serialize(data)
        
        lines = self.build_output_lines(serialized_data)
        
        lines =  self.format_lines(lines)
        
        return "\n".join(lines)
    
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
