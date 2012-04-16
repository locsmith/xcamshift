'''
Created on 30 Dec 2011

@author: garyt
'''
import os
from yaml import load
from distance_table import Distance_table
from random_coil_table import Random_coil_table
import yaml
import extra_table
from dihedral_table import Dihedral_table, Dihedral_parameter_table,\
    Composite_dihedral_table
from python_utils import tupleit
from sidechain_table import Sidechain_table
from constants_table import Constants_table
from ring_table import Ring_table
from non_bonded_table import Non_bonded_table
from table_builders.yaml_patches import add_access_to_yaml_list_based_keys


class Table_manager(object):

    NON_BONDED = "nb"
    

        
    __default = None
    
    @staticmethod
    def get_default_table_manager():
        if Table_manager.__default == None:
            Table_manager.__default = Table_manager()
        return Table_manager.__default
    
    '''
    class to load and store a set of yaml tables for the camshift forcefield, 
    results are looked up in order favouring specific residue types 
    and then a table name
    '''

    BASE_TABLE =  'base'
    TEMPLATE_3 = '%s_%s_%s_%s.yaml'
    TYPE = 'cams'
    VERSION = '1_35_0'
    DEFAULT_DIRECTORY = 'data'
    
    #TODO use common constants
    BACKBONE = "bb"
    RANDOM_COIL = "rc"
    EXTRA="extra"
    DIHEDRAL="dihedral"
    DIHEDRAL_PARS="dihedral_pars"
    SIDECHAIN="sidechain"
    CONSTANTS="constants"
    RING="ring"
    
    def __init__(self,paths=[]):
        '''
        Constructor
        '''
        
        self.search_paths = paths + ['.',self.DEFAULT_DIRECTORY]
        self.tables ={}
        
        add_access_to_yaml_list_based_keys()
                
    def __get_table_name(self, table_type, residue_type):
        return self.TEMPLATE_3 % (self.TYPE,self.VERSION,table_type,residue_type)
    
    def add_search_path(self,path):
        self.search_paths.insert(0, path)
        

    def __raise_table_load_error(self, table_name, detail):
        raise IOError("ERROR: couldn't load %s because %s" % (table_name, detail))


    def __get_residue_types(self, residue_type):
        if residue_type != None:
            residue_types = residue_type, self.BASE_TABLE
        else:
            residue_types = (self.BASE_TABLE,)
        return residue_types

    def __load_table(self, table_type, residue_type=None):

        residue_types = self.__get_residue_types(residue_type)
            
        for residue_type in residue_types:
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

    def _get_table(self,table_type,residue_type=None):

        result = self.__search_for_table(table_type,residue_type)
        
        if result == None:
            self.__load_table(table_type,residue_type)
            
        result = self.__search_for_table(table_type,residue_type)
        return result
            
    def get_BB_Distance_Table(self,residue_type):
        return Distance_table(self._get_table(self.BACKBONE, residue_type))
    
    def get_random_coil_table(self, residue_type):
        return Random_coil_table(self._get_table(self.RANDOM_COIL, residue_type))
    
    def get_extra_table(self,residue_type):
        return extra_table.Extra_table(self._get_table(self.EXTRA,residue_type))
    
    def _get_dihedral_table(self,residue_type):
        return Dihedral_table(self._get_table(self.DIHEDRAL,residue_type))
    
    def _get_dihedral_parameter_table(self,residue_type):
        return Dihedral_parameter_table(self._get_table(self.DIHEDRAL_PARS,residue_type))
    
    def get_dihedral_table(self,residue_type):
        dihedral_table = self._get_dihedral_table(residue_type)
        dihedral_parameter_table =  self._get_dihedral_parameter_table(residue_type)
        return Composite_dihedral_table(dihedral_table, dihedral_parameter_table)
    
    def get_sidechain_table(self,residue_type):
        return Sidechain_table(self._get_table(self.SIDECHAIN,residue_type))
    
    def get_constants_table(self,residue_type):
        return Constants_table(self._get_table(self.CONSTANTS,residue_type))
    
    def get_ring_table(self,residue_type):
        return Ring_table(self._get_table(self.RING,residue_type))
    
    def get_non_bonded_table(self,residue_type):
        return Non_bonded_table(self._get_table(self.NON_BONDED,residue_type))
            