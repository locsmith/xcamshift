'''
Created on Apr 10, 2012

@author: garyt
'''
import unittest2
from table_manager import Table_manager
from table_builders.table_extractor import Table_extractor
from table_builders.xcamshift import Backbone_distance_extractor
from table_builders.xcamshift.camshift_table_builder import extract
from common_constants import BACK_BONE, SIDE_CHAIN
from yaml import Loader, load
import common_constants
from python_utils import IsMappingType, Dict_walker, value_from_key_path,\
    filter_dict
from numbers import Number



#make this a global default               
DEFAULT_DECIMAL_PLACES = 5

class Test_table_importers(unittest2.TestCase):


    def assertDictAlmostEqual(self, dict_1, dict_2, places = DEFAULT_DECIMAL_PLACES):
        dict_1_aggregator = Dict_walker.Aggregator()
        dict_2_aggregator = Dict_walker.Aggregator()
        
        dict_walker_1 = Dict_walker()
        dict_walker_2 = Dict_walker()
        
        dict_walker_1.walk_dict(dict_1, dict_1_aggregator)
        dict_walker_2.walk_dict(dict_2, dict_2_aggregator)
        
        dict_1_keys = dict_1_aggregator.get_result().keys()
        dict_1_keys.sort()
        
        dict_2_keys = dict_2_aggregator.get_result().keys()
        dict_2_keys.sort()
        

        self.assertSequenceEqual(dict_1_keys, dict_2_keys, "keys for the two dictionarys are not equal")
        
        for key_path_1,key_path_2 in zip(dict_1_keys,dict_2_keys):
            data_1_value = value_from_key_path(dict_1,key_path_1)
            data_2_value = value_from_key_path(dict_2,key_path_2)
            
            if isinstance(data_1_value, Number):
                self.assertAlmostEqual(data_1_value, data_2_value, places)
            else:
                self.assertEqual(data_1_value, data_2_value)
                
                
                
    
    
    def test_Tables(self):
        for sub_potential in common_constants.CAMSHIFT_SUB_POTENTIALS:

            sub_potential_file_id = sub_potential.lower().strip()
            
            table_manager  = Table_manager()
            expected_table_data = table_manager._get_table(sub_potential_file_id, 'ALA')
            
            read_table = extract("data/camshift",sub_potential,'')
            
            read_table_data = load(read_table[sub_potential, ''])
            

    
            self.maxDiff = None
            data_entry = read_table_data['data']
            
            # filter data that are not in the hand built tables
            if sub_potential == SIDE_CHAIN:
                filter_dict(data_entry, lambda k,v : k != 'ALA')
                filter_dict(data_entry['ALA'], lambda k,v : k != 'CB')

            self.assertDictAlmostEqual(data_entry, expected_table_data['data'])

        

        


        
        
#if __name__ == "__main__":
#    #import sys;sys.argv = ['', 'Test.testName']
#    unittest2.main()