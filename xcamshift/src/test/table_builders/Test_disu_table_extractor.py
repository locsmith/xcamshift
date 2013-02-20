'''
Created on 19 Feb 2013

@author: garyt
'''
import unittest
from table_builders.xcamshift.Disulphide_data_extractor import DISU_table_extractor
from common_constants import DATA, h_keys


class Test(unittest.TestCase):

    def assertDictHasKey(self,dict,key,msg=''):
        msg = '%s not in dict' % `key`
        self.assertTrue(key in dict, msg = msg)
        
    def setUp(self):
        XTRADISTS = 'XTRADISTS' 
        test_values=  {'HA' : 0.0115552, 'CA' : -2.851700383944911, 'H' : -0.106055, 'N': -1.41314, 'C' :  -0.370881, 'CB' : 11.3616}

        self.test_data= {XTRADISTS : {}}
        line_data = self.test_data[XTRADISTS]
        for h_key in h_keys:
            line_data[h_key]  = [0.0] * 27
            line_data[h_key][-1] = test_values[h_key]
            
        self.data_extractor = DISU_table_extractor(self.test_data)

    def tearDown(self):
        pass

    
    def test_serialization_has_data_key(self):
        result  = self.data_extractor.serialize(self.test_data)
        self.assertDictHasKey(result,DATA)
    
    def test_serialization_has_data_disu_cys_key(self):
        
        result  = self.data_extractor.serialize(self.test_data)
        self.assertDictHasKey(result[DATA],('DISU','CYS'))

    def test_serialization_has_data_disu_cys_ca_key(self):
        
        result  = self.data_extractor.serialize(self.test_data)
        self.assertDictHasKey(result[DATA][('DISU','CYS')],'CA') 
        
    def test_serialization_has_data_disu_cys_ca_value(self):
        
        result  = self.data_extractor.serialize(self.test_data) 
        self.assertAlmostEqual(result[DATA][('DISU','CYS')]['CA'], -2.851700383944911)   


if __name__ == "__main__":
    unittest.main()