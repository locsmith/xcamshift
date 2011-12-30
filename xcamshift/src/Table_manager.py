'''
Created on 30 Dec 2011

@author: garyt
'''
import os
from yaml import load

class Table_manager(object):
    '''
    classdocs
    '''

    BASE_TABLE =  'base'
    TEMPLATE = '%s_%s_%s_%s.json'
    TYPE = 'cams'
    VERSION = '1_35_0'
    DEFAULT_DIRECTORY = 'data'
    
    def __init__(self):
        '''
        Constructor
        '''
        
        self.search_paths = ['.',self.DEFAULT_DIRECTORY]
        self.tables ={}
        
        
    
    def get_table_name(self, table_type, residue_type):
        return self.TEMPLATE % (self.TYPE,self.VERSION,table_type,residue_type)
    
    def add_search_path(self,path):
        self.search_paths.insert(0, path)
        
    def load_table(self, table_type, residue_type):
        
        for residue_type in (residue_type,self.BASE_TABLE):
            table_name = self.get_table_name(table_type,residue_type)
            
            for search_path in self.search_paths:
                print os.getcwd()
                path = os.path.join(search_path, table_name)
                print path
                if os.path.exists(path):
                    key = (table_type,residue_type)
                    
                    try:
                        table_fp = open(path)
                        self.tables[key] = load(table_fp)
                        print self.tables
                    except Exception as detail:
                        raise IOError("ERROR: couldn't load %s because %s"  % (table_name, detail))
    
    
    def get_table(self,table_type,residue_type):
        key = (table_type,residue_type)
        base_key = (table_type,self.BASE_TABLE)
        
        tables = self.tables
        result = None
        if key in tables:
            result = tables[key]
        elif base_key in tables:
            result = tables[key]
        else:
            result = self.load_table(table_type,residue_type)
        
        return result
            
            