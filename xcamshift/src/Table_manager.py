'''
Created on 30 Dec 2011

@author: garyt
'''
import os
from yaml import load
from test.distance_table import Distance_table

class Table_manager(object):
    '''
    class to load an store a set of yaml tables for the camshift forcefield, 
    results are looked up in order favouring specific residue types 
    and then a table name
    '''

    BASE_TABLE =  'base'
    TEMPLATE = '%s_%s_%s_%s.yaml'
    TYPE = 'cams'
    VERSION = '1_35_0'
    DEFAULT_DIRECTORY = 'data'
    
    BACKBONE = "bb"
    
    def __init__(self):
        '''
        Constructor
        '''
        
        self.search_paths = ['.',self.DEFAULT_DIRECTORY]
        self.tables ={}
        
        
    
    def __get_table_name(self, table_type, residue_type):
        return self.TEMPLATE % (self.TYPE,self.VERSION,table_type,residue_type)
    
    def add_search_path(self,path):
        self.search_paths.insert(0, path)
        

    def __raise_table_load_error(self, table_name, detail):
        raise IOError("ERROR: couldn't load %s because %s" % (table_name, detail))

    def __load_table(self, table_type, residue_type):
        
        for residue_type in (residue_type,self.BASE_TABLE):
            table_name = self.__get_table_name(table_type,residue_type)
            
            new_table = None
            for search_path in self.search_paths:
                path = os.path.join(search_path, table_name)
                
                if os.path.exists(path):
                    key = (table_type,residue_type)
                    
                    try:
                        table_fp = open(path)
                        new_table = load(table_fp)
                        if new_table != None:
                            break
                    except Exception as detail:
                        self.__raise_table_load_error(table_name, detail)
                        
        if new_table == None:
            self.__raise_table_load_error(table_name, "the table couldn't be found in %s" % ", ".join(self.search_paths))
        self.tables[key]=new_table
    

    def __get_seach_keys(self,table_type, residue_type):
        result = ((table_type,residue_type),(table_type,self.BASE_TABLE))
        return result
    
    def __search_for_table(self, table_type, residue_type):
        keys= self.__get_seach_keys(table_type, residue_type)
        tables = self.tables
        
        result = None
        for key in keys:
            if key in  tables:
                result = tables[key]
                break
        
        return result

    def __get_table(self,table_type,residue_type):

        
        result = self.__search_for_table(table_type,residue_type)
        
        if result == None:
            self.__load_table(table_type,residue_type)
            
        result = self.__search_for_table(table_type,residue_type)
        return result
            
    def get_BB_Distance_Table(self,residue_type):
        return Distance_table(self.__get_table(self.BACKBONE, residue_type))
            