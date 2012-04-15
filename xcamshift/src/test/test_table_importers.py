'''
Created on Apr 10, 2012

@author: garyt
'''
import unittest2
from table_manager import Table_manager
from table_builders.table_extractor import Table_extractor
from table_builders.xcamshift import Backbone_distance_extractor
from table_builders.xcamshift.camshift_table_builder import extract
from common_constants import BACK_BONE, SIDE_CHAIN, BASE, GLY, PRO, C, CA, DATA,\
    N, HN, XTRA, DIHEDRAL, RING, NON_BONDED, CB, SPHERE_1, SP3, HA, SP2, O,\
    SPHERE_2, CG, CAMSHIFT_RESIDUE_TYPES
from yaml import Loader, load
import common_constants
from python_utils import IsMappingType, Dict_walker, value_from_key_path,\
    filter_dict
from numbers import Number
from table_builders.yaml_patches import add_access_to_yaml_list_based_keys



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
                
                
                
    
    
    def test_tables_against_original(self):
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
            
    def test_table_entries(self):
        add_access_to_yaml_list_based_keys()
        
        # list of top left and botttom right significant entries in the paper for testing
        # ammended to update with the values in the 1.35.0 data files from almost 1.04 (which are not the same as those in camshift 1.35.0?)
        # not most changes are for CB values and it has been checked that they only differ in detail (i.e. typically the changes are after ~3 dp)
        # this test does not in clude values for 
        #
        #  1. random coil values 
        #  2. hydrogen bonds
        #  3. disulphide cystines
        
        expected = {
                    BACK_BONE : {
                          BASE : {(DATA,-1,N,HA) :  0.1049538,
                                  (DATA,1,C,CB)  :  0.11095702}, # note this is different from the paper (0.1112132)
                          GLY :  {(DATA,-1,N,HA) :  0.0663175,
                                  (DATA,1,C,C)   : -0.9356821},
                          PRO :  {(DATA,-1,N,HA) :  0.0297782,
                                  (DATA,1,C,CB)  :  0.09684666 } # note this is different from the paper (0.0970634)
                    },
                    
                    SIDE_CHAIN : {
                          BASE : {(DATA,'ALA',CB,CA)     : -3.4713212,
                                  (DATA,'VAL','HB',C)    : -0.5444149},
                          
                          GLY : {(DATA,'GLY','HA2',HA)   : 0.059929,
                                 (DATA,'GLY','HA2',C)    : 0.4653395},
                                  
                          PRO : {(DATA,'PRO',CB,HA)      : 0.7663159,
                                 (DATA,'PRO','HD2',CB)   :-2.86174391}  # note this is different from the paper (-2.8604453)
                    },
                    
                    NON_BONDED  : {
                          BASE : {(DATA,SPHERE_1,(C,SP3),HA)      : -2.9677854,
                                  (DATA,SPHERE_1,(O,SP2),CB)      : -9.21507952586192,      # note this is different from the paper (-9.2063201)
                                  (DATA,SPHERE_2,(C,SP3),HA)      :  0.0026981,
                                  (DATA,SPHERE_2,(O,SP2),CB)      :  0.03094620239855296 }, # note this is different from the paper (0.0308797)
                                 
                          GLY : {(DATA,SPHERE_1,(C,SP3),HA)       : -3.0168132,
                                 (DATA,SPHERE_1,(O,SP2),C)        : 24.239118,  
                                 (DATA,SPHERE_2,(C,SP3),HA)       :  0.0019265,
                                 (DATA,SPHERE_2,(O,SP2),C)        : -0.1020419 }, 
                                  
                          PRO : {(DATA,SPHERE_1,(C,SP3),HA)       : -2.8890318,
                                 (DATA,SPHERE_1,(O,SP2),CB)       : -8.63293246,  # note this is different from the paper (-8.62477719)
                                 (DATA,SPHERE_2,(C,SP3),HA)       :  0.0023439,
                                 (DATA,SPHERE_2,(O,SP2),CB)       :  0.03006555}  # note this is different from the paper (0.0300053)
                    },
                    
                    RING  : {
                          BASE : {(DATA,("PHE","6"),HA) : 0.0197106,
                                  (DATA,("HIS","5"),CB) : 0.013111927804067058}, # note this is different from the paper (0.0131036)
                                  
                          GLY :  {(DATA,("PHE","6"),HA) : 0.0187097,
                                  (DATA,("HIS","5"),C)  : 0.0169792},
                             
                          PRO :  {(DATA,("PHE","6"),HA) : 0.0186948,
                                  (DATA,("HIS","5"),CB) :  0.0144458}            #  note this is different from the paper (0.0144385)
                    },
                    
                    DIHEDRAL  : {
                          BASE : {(DATA,((C, -1),  (N, 0), (CA, 0),  (C, 0)),HA) : 0.3960648,
                                  (DATA,((N, 0),  (CA, 0), (CB, 0), (CG, 0)),CB) : 0.1438778 },  # note this is different from the paper (0.1438652)
                                  
                          GLY  : {(DATA,((C, -1),  (N, 0), (CA, 0),  (C, 0)),HA) : 0.4864764,
                                  (DATA,((N, 0),  (CA, 0), (CB, 0), (CG, 0)),C)  : 0.4117375},
                                 
                          PRO  : {(DATA,((C, -1),  (N, 0), (CA, 0),  (C, 0)),HA) : 0.4106709,
                                  (DATA,((N, 0),  (CA, 0), (CB, 0), (CG, 0)),CB) : 0.1205188 }   # note this is different from the paper (0.1204433)
                    },
                    
                    XTRA : {
                          BASE : {(DATA, (HN, 0),  (HA, 0), HA) : -0.2327538,
                                  (DATA, (CA, -1), (CA, 1), CB) :  1.09990498 }, # note this is different from the bpaper value (1.0992811)
                                  
                          GLY  : {(DATA, (HN, 0),  (HA, 0), HA) : -0.1629661,
                                  (DATA, (CA, -1), (CA, 1), C)  :  0.6711163},
                            
                          PRO  : {(DATA, (C , -1), (HA, 0), HA) :  0.02395468,
                                  (DATA, (CA, -1), (CA, 1), CB) :  0.3278771 }   # note this is different from the bpaper value (0.3274878)
                    },
                    
        }
        for sub_potential in common_constants.CAMSHIFT_SUB_POTENTIALS:
            for residue_type in CAMSHIFT_RESIDUE_TYPES:
                read_table = extract("data/camshift",sub_potential,residue_type)
                
                read_table_data = load(read_table[sub_potential, residue_type])
                
                for key in value_from_key_path(expected,(sub_potential,residue_type)):
                    value = value_from_key_path(read_table_data, key)
                    expected_key = sub_potential,residue_type,key
                    expected_value  = value_from_key_path(expected,expected_key)
                    self.assertAlmostEqual(value, expected_value, places=DEFAULT_DECIMAL_PLACES,msg=`expected_key`)
        

        


        
        
#if __name__ == "__main__":
#    #import sys;sys.argv = ['', 'Test.testName']
#    unittest2.main()