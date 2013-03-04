'''
Created on 1 Mar 2013

@author: garyt
'''
import unittest
from collection_backport import OrderedDict
from copy import deepcopy
from table_builders.Table_modifier import Table_modifier


TEST_DATA=None

class Test(unittest.TestCase):

    
    def setUp(self):
        self._reset_test_data()
        self.test_target_data= deepcopy(TEST_DATA)
    
    def _reset_test_data(self):
        global TEST_DATA
        TEST_DATA = OrderedDict()
        TEST_DATA['A'] = OrderedDict()
        TEST_DATA['A']['B']='1'
        TEST_DATA['C']='2'
            
    def test_append_operation(self):
        program = (
                   ('append', ('D','3')),
                   )
        
        modifier = Table_modifier(program)
        
        result = modifier.run(self.test_target_data)
        
        EXPECTED_RESULT = deepcopy(TEST_DATA)
        EXPECTED_RESULT['D']='3'
        self.assertEqual(result,EXPECTED_RESULT)




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()