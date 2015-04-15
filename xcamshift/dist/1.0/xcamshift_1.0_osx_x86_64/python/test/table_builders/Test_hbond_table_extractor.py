'''
Created on 19 Feb 2013

@author: garyt
'''
import unittest
from table_builders.xcamshift.hbond_data_extractor import HBOND_table_extractor
from common_constants import DATA, h_keys
from collection_backport import OrderedDict


class Test(unittest.TestCase):

    def assertDictHasKey(self,dict,key,msg=''):
        msg = '%s not in dict' % `key`
        self.assertTrue(key in dict, msg = msg)
      
    def setUp(self):
            

        HBOND = 'HBOND' 
        self.test_data_N =  (-0.0000752745753528108, -0.018041072764128266,   0.0018754630721177568,
                             -0.013428139650030523,   0.0003982251798295468, -0.0001003438625950712,
                             -0.00008440091204673355, 0.006562726645723401,  -0.0006374109394798102,
                             -0.015572931717349981,   0.01177279379780734,   -0.0012193090771868422)
        
        self.expected_key = (((-1,'O','HN'),'DIST'), ((-1,'O','HN'),'ANG1'), ((-1,'O','HN'),'ANG2'),
                             (( 0,'HN','O'),'DIST'), (( 0,'HN','O'),'ANG1'), (( 0,'HN','O'),'ANG2'),
                             (( 0,'O','HN'),'DIST'), (( 0,'O','HN'),'ANG1'), (( 0,'O','HN'),'ANG2'),
                             (( 1,'HN','O'),'DIST'), (( 1,'HN','O'),'ANG1'), (( 1,'HN','O'),'ANG2'))
        
        self.test_data = OrderedDict()
        self.test_data['HBONDS'] = OrderedDict()
        
        for h_key in h_keys:
            self.test_data['HBONDS'][h_key] = [0.0]*12
        
        self.test_data['HBONDS']['N'] = self.test_data_N
            
        self.data_extractor = HBOND_table_extractor(self.test_data)
    
    def test_serialization_has_data_key(self):
        result  = self.data_extractor.serialize(self.test_data)
        self.assertDictHasKey(result,DATA)
        

    def test_serialization_has_prior_current_and_post_hbonds_keys(self):
        
        result  = self.data_extractor.serialize(self.test_data)
        
        for key in ((-1,'O','HN'),(0,'HN','O'),(0,'O','HN'),(1,'HN','O')):
            self.assertDictHasKey(result[DATA],key)
                    
    def test_serialization_has_distance_and_angle_parameters(self):
        result  = self.data_extractor.serialize(self.test_data)
        
        for key in ((-1,'O','HN'),(0,'HN','O'),(0,'O','HN'),(1,'HN','O')):
            for param in ['DIST','ANG1','ANG2']:
                self.assertDictHasKey(result[DATA][key],param)

    def test_serialization_has_param_keys(self):

        result  = self.data_extractor.serialize(self.test_data)
        
        for key in ((-1,'O','HN'),(0,'HN','O'),(0,'O','HN'),(1,'HN','O')):
            for param in ['DIST','ANG1','ANG2']:
                self.assertDictHasKey(result[DATA][key][param],'N')            
    
    def test_serialization_of_N(self):
    
        result  = self.data_extractor.serialize(self.test_data)
        
        for (h_bond_key,param_key), value in zip(self.expected_key,self.test_data_N):
            self.assertAlmostEqual(result[DATA][h_bond_key][param_key]['N'],value) 
        
if __name__ == "__main__":
    unittest.main()