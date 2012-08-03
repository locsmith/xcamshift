'''
Created on 11 May 2012

@author: garyt
'''
from abc import abstractmethod, ABCMeta

TRANSLATIONS = 'translations'
TABLE_TYPE = 'table_type'
RESIDUE_TYPE = 'residue_type'
TABLE_INDEX = 'index'

class Table_base(object):
    __metaclass__ = ABCMeta
    
    def __init__(self,table):
        self._table = table
        
        self._translation_to_table = self._build_from_translation_table(table)
        self._translation_from_table = self._build_from_translation_table(table)

        
    def _build_to_translation_table(self,data_table):
        result  = {}
        if TRANSLATIONS in data_table:
            result  =  data_table[TRANSLATIONS]
        return result
    
    def _build_from_translation_table(self,data_table):
        result  = {}
        if TRANSLATIONS in data_table:
            for key,to_atom in data_table[TRANSLATIONS].items():
                residue, from_atom = key
                new_key = (residue,to_atom)
                result[new_key]= from_atom
        return result
    
    def get_translation_to_table(self,residue,atom):
        result =  atom
        
        key = residue, atom
        if key in self._table[TRANSLATIONS]:
            result =  self._table[TRANSLATIONS][key]
        return result

    def get_translation_from_table(self, residue, atom):
        result = atom
        
        key = residue,atom
        if key in self._translation_from_table:
            result = self._translation_from_table[key]
        return result
    
    @abstractmethod
    def get_target_atoms(self):
        pass
    
    def get_table_residue_type(self):
        return self._table[RESIDUE_TYPE]
    
    def get_table_type(self):
        return self._table[TABLE_TYPE]
    
    def get_table_index(self):
        return self._table[TABLE_INDEX]

