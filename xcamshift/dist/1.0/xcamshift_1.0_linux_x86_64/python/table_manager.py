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
Created on 30 Dec 2011

@author: garyt
'''
import os,sys
from yaml import load
from distance_table import Distance_table
from random_coil_table import Random_coil_table
import yaml
import extra_table
from dihedral_table import Dihedral_table, Dihedral_parameter_table,\
    Composite_dihedral_table
from hydrogen_bond_table import Hydrogen_bond_table
from python_utils import tupleit, Hierarchical_dict
from sidechain_table import Sidechain_table
from constants_table import Constants_table
from ring_table import Ring_table
from non_bonded_table import Non_bonded_table
from disulphide_table import Disulphide_table
from table_builders.yaml_patches import add_access_to_yaml_list_based_keys
import utils
from table_base import RESIDUE_TYPE, TABLE_TYPE, TABLE_INDEX, TABLE_LOADED
from cache_manager import Cache_manager

BASE_TABLE = 'base'

TABLE_MANAGER = 'TABLE_MANAGER'
def _build_table_manager():
    return Table_manager()

Cache_manager.get_cache_manager().add_cache_builder(TABLE_MANAGER,_build_table_manager)

class Database_not_found_exception(Exception):
    def __init__(self,*args,**kwargs):
        super(Database_not_found_exception, self).__init__(*args,**kwargs)

def _get_residue_type_cache():
    return  Cache_manager.get_cache(RESIDUE_TYPE)
#TODO: cleanup internal structure, caching needs a better implementation
#TODO: remove specific functions to load tables?
#TODO: needs composite table support and hierachical item integration
#TODO: caching should work on the wrapper object not the internal _table object
class Table_manager(object):



    @staticmethod
    def get_default_table_manager():
        return Cache_manager.get_cache_manager().get_cache(TABLE_MANAGER)


    '''
    class to load and store a set of yaml tables for the camshift forcefield,
    results are looked up in order favouring specific residue types
    and then a table name
    '''

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
    NON_BONDED = "nb"
    DISULPHIDE = 'disu'
    HBOND = 'hbond'

    def __init__(self,paths=[]):
        '''
        Constructor
        '''


        self.search_paths = paths + ['.',
                                     self.DEFAULT_DIRECTORY,
                                     self._get_program_database__path(),
                                     self._get_distribution_database_path()]

        self.tables ={}
        self._known_table_types = set()
        self._table_index = {}

        add_access_to_yaml_list_based_keys()

    def _get_program_database__path(self):
        return os.path.join(os.path.dirname(os.path.realpath(sys.executable)),'..','databases','XCamShift')

    def _get_distribution_database_path(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data')

    def __str__(self):
        result = []
        result.append('table manager')
        result.append('')
        result.append('search paths: %s' % ', '.join(self.search_paths))
        result.append('')
        key_list = ', '.join(['%s.%s' % (table,residue) for table,residue in self.searched_for_tables])
#        print key_list
        result.append('searched for tables: %s' %  key_list)
        result.append('')
        # add table indices in output
        key_list = ', '.join(['%s.%s' % (table,residue) for table,residue in self.tables.keys()])
        result.append('tables: %s' % key_list)
        result.append('')


        return '\n'.join(result)

    def get_table_types(self):

        types = set([key[0]for key in self.tables.keys()])
        types = list(types)
        types.sort()
        return tuple(types)


    def _make_presorted_tuple(self, result):
        result = list(result)
        result.sort()
        result = tuple(result)

        return result


    def _load_tables(self, table_type):
        for residue_type in utils.iter_residue_types():
            self._get_table(table_type, residue_type)


    def get_residue_types_for_table(self,table_type):
        table_type = table_type.lower()
        self._load_tables(table_type)
        residues = set([key[1] for key in self.tables.keys() if key[0] == table_type])

        self._make_presorted_tuple(residues)
        return residues

    #TODO: note this only works for loaded potentials
    #TODO: note the above statement is wrong and maybe should be corrected?
    #TODO: in actual fact it could just be replaced by a call to the segment manager!
    def get_all_known_residue_types(self):

        result = set()
        for table_type in self.get_table_types():
            for residue_type in self.get_residue_types_for_table(table_type):
                result.add(residue_type)


        return self._make_presorted_tuple(result)

    def get_all_known_table_residue_types(self):

        result  =  set(self.tables[key][RESIDUE_TYPE] for key in self.tables)

        return self._make_presorted_tuple(result)


    def __get_table_name(self, table_type, residue_type):
        return self.TEMPLATE_3 % (self.TYPE,self.VERSION,table_type,residue_type)

    def add_search_path(self,path):
        self.search_paths.insert(0, path)


    def __raise_table_load_error(self, table_name, detail):
        raise IOError("ERROR: couldn't load %s because %s" % (table_name, detail))


    #TODO: this could come from the file itself...
    def _find_parent_table(self, table_type, residue_type):
        result = None
        if residue_type != BASE_TABLE:
            result = self._get_table(table_type,None)
        return result
#



    def _register_new_table(self, new_table, key):
        table_type =  key[0]
        if new_table != None:
            self.tables[key] = new_table
        else:
            if (table_type,BASE_TABLE) in self.tables:
                self.tables[key]= self.tables[table_type,BASE_TABLE]
            else:
                table_name  = self.__get_table_name(table_type, BASE_TABLE)
                msg =  'ERROR: table %s is missing, the path to the database is bad!\n' % table_name
                msg += '       search path (length = %i) is: %s\n' % (len(self.search_paths),'\n'.join(self.search_paths))
                raise Database_not_found_exception(msg)

        self._known_table_types.add(table_type)



    def _ornament_table(self, new_table, table_type, residue_type):
        table_index = self._table_index.setdefault(table_type,0)
        if not TABLE_LOADED in new_table:
            ornaments  = {
                            RESIDUE_TYPE : residue_type,
                            TABLE_TYPE : table_type,
                            TABLE_INDEX : table_index,
                            TABLE_LOADED : True,
                          }
            new_table.update(ornaments)
            self._table_index[table_type]+=1

    #TODO: make this take a key instead
    def __load_table(self, table_type, residue_type=BASE_TABLE):

#        residue_types = self.__get_residue_types(residue_type)

        new_table = None
#        for residue_type in residue_types:
        table_search_key = table_type,residue_type
#            print 'search for ', table_search_key
#            print residue_type
        table_name = self.__get_table_name(table_type,residue_type)

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
#
#        if new_table != None:
#            break


        if new_table != None:
            parent  = self._find_parent_table(table_type,residue_type)
            self._ornament_table(new_table,table_type,residue_type)
            new_table = Hierarchical_dict(new_table, parent=parent)


#        print self.tables.keys()
        return new_table


    def __get_seach_keys(self,table_type, residue_type):
        result = ((table_type,residue_type),(table_type,self.BASE_TABLE))
        return result

    #TODO: remove __'s
    def _seach_for_loaded_table(self,table_type, residue_type):

        result  = None
        key = table_type, residue_type
        if key in self.tables:
            result  = self.tables[key]
        return result

    #TODO:  combine with search for loaded table and and rename
    def __search_for_table(self, table_type, residue_type):
        keys= self.__get_seach_keys(table_type, residue_type)
        tables = self.tables

        result = None
        for key in keys:
            if key in  tables:
                result = tables[key]
                break

        if result == None:
            search_paths = ", ".join(self.search_paths)
            working_directory =  os.getcwd()
            args =  search_paths, working_directory
            self.__raise_table_load_error(table_type, "the table couldn't be found in %s working directory is: %s" % args)

        return result


    def _force_residue_type_lowercase(self, residue_type):
        if residue_type != None:
            residue_type = residue_type.lower()
        else:
            residue_type =  BASE_TABLE.lower()
        return residue_type

    def iter_residue_types(self):
        yield(BASE_TABLE)

        for residue_type in utils.iter_residue_types():
            yield residue_type

    def load_tables_for_known_residues(self, table_type):
        table_type =  table_type.lower()
        # add sequence lookup delegate to allow testing (currently we need a real molecule)
        for residue_type in self.iter_residue_types():

            residue_type = self._force_residue_type_lowercase(residue_type)
            key =  table_type,residue_type
            new_table  = self.__load_table(table_type, residue_type)

            self._register_new_table(new_table, key)

        self._known_table_types.add(table_type)


    #TODO: this is a really inefficient method it will do lots of disk accesses, why not load all residue type tables on first call
    def _get_table(self,table_type,residue_type=None):
        residue_type = self._force_residue_type_lowercase(residue_type)

        if table_type not in self._known_table_types:
            self.load_tables_for_known_residues(table_type)



        key= table_type,residue_type
        if key in self.tables:
            result = self.tables[key]
        else:
            result = self._seach_for_loaded_table(table_type,residue_type)
            if result == None:
                self.__raise_table_load_error(table_type, "couldn't find table for residue type %s" % residue_type)
            self.tables[key]=result

#        self.tables.keys()



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

    def get_disulphide_table(self,residue_type):
        return Disulphide_table(self._get_table(self.DISULPHIDE,residue_type))

    def get_hydrogen_bond_table(self,residue_type):
        return Hydrogen_bond_table(self._get_table(self.HBOND,residue_type))


