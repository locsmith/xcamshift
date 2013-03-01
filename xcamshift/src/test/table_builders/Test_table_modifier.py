'''
Created on 1 Mar 2013

@author: garyt
'''
import unittest
from collection_backport import OrderedDict
from copy import deepcopy
from table_builders.Table_modifier import Table_modifier

TEST_DATA = OrderedDict()
TEST_DATA['A'] = OrderedDict()
TEST_DATA['A']['B']='1'
TEST_DATA['C']='2'

class Test(unittest.TestCase):


    def setUp(self):
        self.test_data= deepcopy(TEST_DATA)
        
    def test_append_operation(self):
        program = (
                   ('append', ('D','3')),
                   )
        
        modifier = Table_modifier(program)
        
        result = modifier.run(self.test_data)
        
        EXPECTED_RESULT = deepcopy(self.test_data)
        EXPECTED_RESULT['D']='3'
        self.assertEqual(result,EXPECTED_RESULT)




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()