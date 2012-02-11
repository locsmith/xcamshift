'''
Created on 11 Feb 2012

@author: garyt
'''
import unittest2
from component_list import Component_list

TEST_DATA_1 = ((1,2),
                       (2,4),
                       (1,3))
TEST_DATA_1 = ((1,2),
               (2,4),
               (1,3))

EXPECTED_1 =  (((1,2),(1,3)),
               ((2,4),))

EXPECTED_COMPONENTS_1 = (1,2)

TEST_DATA_2 = ((1,4),(0,1))

EXPECTED_ALL_1 = ((1,2),(1,3),(2,4))

EXPECTED_2 =  (((0,1),),
               ((1,2),(1,3),(1,4)),
               ((2,4),))

EXPECTED_COMPONENTS_2 = (0,1,2)

EXPECTED_ALL_2 = ((0,1),(1,2),(1,3),(1,4),(2,4))

class Test_component_list(unittest2.TestCase):


    def setUp(self):
        self._component_list =  Component_list()

    def testEmptyComponentList(self):
        with self.assertRaises(Exception):
            self._component_list.get_components_for__atom_id(0)
    
    def testExpectedRanges(self):
        
        self._component_list.add_components(TEST_DATA_1)
        for atom_id,expected in zip (EXPECTED_COMPONENTS_1,EXPECTED_1):
            result = self._component_list.get_components_for__atom_id(atom_id)
            self.assertSequenceEqual(result, expected)
    
    def testExpectedRaangesAfterAdd(self):
        self._component_list.add_components(TEST_DATA_1)
        self._component_list.add_components(TEST_DATA_2)
        for atom_id,expected in zip (EXPECTED_COMPONENTS_2,EXPECTED_2):
            result = self._component_list.get_components_for__atom_id(atom_id)
            self.assertSequenceEqual(result, expected)
            
    def testComponentIds(self):
        self._component_list.add_components(TEST_DATA_1)
        result = self._component_list.get_component_atom_ids()
        self.assertSequenceEqual(EXPECTED_COMPONENTS_1,result)
        
        self._component_list.add_components(TEST_DATA_2)
        result = self._component_list.get_component_atom_ids()
        self.assertSequenceEqual(EXPECTED_COMPONENTS_2,result)
        
    def testAddComponent(self):
        for component in TEST_DATA_1:
            self._component_list.add_component(component)
        result = self._component_list.get_all_components()
        
        self.assertSequenceEqual(result,EXPECTED_ALL_1)
        
        for component in TEST_DATA_2:
            self._component_list.add_component(component)
        result = self._component_list.get_all_components()
        
        self.assertSequenceEqual(result,EXPECTED_ALL_2)

#        TODO add missing tests
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testEmptyComponentList']

    unittest2.main()