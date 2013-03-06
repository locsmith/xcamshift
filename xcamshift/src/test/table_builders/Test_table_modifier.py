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
                   ('append', ('D','3'), ('E','4')),
                   )
        
        modifier = Table_modifier(program)
        
        result = modifier.run(self.test_target_data)
        
        EXPECTED_RESULT = deepcopy(TEST_DATA)
        EXPECTED_RESULT['D']='3'
        EXPECTED_RESULT['E']='4'
        self.assertEqual(result,EXPECTED_RESULT)
        
                    
    def test_append_multiple_operation(self):
        program = (
                   ('append', ('D','3')),
                   )
        
        modifier = Table_modifier(program)
        
        result = modifier.run(self.test_target_data)
        
        EXPECTED_RESULT = deepcopy(TEST_DATA)
        EXPECTED_RESULT['D']='3'
        self.assertEqual(result,EXPECTED_RESULT)
        
    def test_append_to_operation(self):
        program = (
                   ('append_to', ('A',),('D','3')),
                   )
        
        modifier = Table_modifier(program)
        
        result = modifier.run(self.test_target_data)
        
        EXPECTED_RESULT = deepcopy(TEST_DATA)
        EXPECTED_RESULT['A']['D']='3'
        
    def test_append_to_multiple_operation(self):
        program = (
                   ('append_to', ('A',),('D','3'),('E','4')),
                   )
        
        modifier = Table_modifier(program)
        
        result = modifier.run(self.test_target_data)
        
        EXPECTED_RESULT = deepcopy(TEST_DATA)
        EXPECTED_RESULT['A']['D']='3'
        EXPECTED_RESULT['A']['E']='4'
        
    def test_prepend_to_operation(self):
        program = (
                   ('prepend_to', ('A',),('D','3')),
                   )
        
        modifier = Table_modifier(program)
        
        result = modifier.run(self.test_target_data)
        
        EXPECTED_RESULT = deepcopy(TEST_DATA)
        EXPECTED_RESULT['A']=OrderedDict()
        EXPECTED_RESULT['A']['D']='3'
        EXPECTED_RESULT['A']['B']=TEST_DATA['A']['B']
        
        self.assertEqual(result,EXPECTED_RESULT)        

    def test_multiple_prepend_to_operation(self):
        program = (
                   ('prepend_to', ('A',),('D','3'),('E','4')),
                   )
        
        modifier = Table_modifier(program)
        
        result = modifier.run(self.test_target_data)
        
        EXPECTED_RESULT = deepcopy(TEST_DATA)
        EXPECTED_RESULT['A']=OrderedDict()
        EXPECTED_RESULT['A']['D']='3'
        EXPECTED_RESULT['A']['E']='4'
        EXPECTED_RESULT['A']['B']=TEST_DATA['A']['B']
        
        self.assertEqual(result,EXPECTED_RESULT) 

    def test_prepend_operation(self):
        program = (
                   ('prepend', ('D','3')),
                   )
        
        modifier = Table_modifier(program)
        result = modifier.run(self.test_target_data)
        
        EXPECTED_RESULT = OrderedDict()
        EXPECTED_RESULT['D'] = '3'
        EXPECTED_RESULT['A'] = TEST_DATA['A']
        EXPECTED_RESULT['C'] = TEST_DATA['C']
        
        self.assertEqual(result,EXPECTED_RESULT) 

    def test_prepend_multiple_operation(self):
        program = (
                   ('prepend', ('D','3'),('E','4')),
                   )
        
        modifier = Table_modifier(program)
        result = modifier.run(self.test_target_data)
        
        EXPECTED_RESULT = OrderedDict()
        EXPECTED_RESULT['D'] = '3'
        EXPECTED_RESULT['E'] = '4'
        EXPECTED_RESULT['A'] = TEST_DATA['A']
        EXPECTED_RESULT['C'] = TEST_DATA['C']
        
        self.assertEqual(result,EXPECTED_RESULT) 
        
    def test_replace(self):
        program = (
                   ('replace', ('A','B'), '3'),
                   )

        modifier = Table_modifier(program)
        result = modifier.run(self.test_target_data)
        
        EXPECTED_RESULT = deepcopy(TEST_DATA)
        EXPECTED_RESULT['A']['B'] =  '3'
        
        self.assertEqual(result, EXPECTED_RESULT)
        
    def test_add(self):
        program = (
                   ('add', ('A','B'), '4'),
                   )

        modifier = Table_modifier(program)
        result = modifier.run(self.test_target_data)
        
        EXPECTED_RESULT = deepcopy(TEST_DATA)
        EXPECTED_RESULT['A']['B'] =  '14'
        
        self.assertEqual(result, EXPECTED_RESULT)        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()