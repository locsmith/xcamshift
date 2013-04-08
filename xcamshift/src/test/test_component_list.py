#-------------------------------------------------------------------------------
# Copyright (c) 2013 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 11 Feb 2012

@author: garyt
'''
import unittest2
from component_list import Component_list
from struct import  Struct
from array import array
from  cython.shift_calculators import test_dump_dist_comp
import ctypes
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

TEST_DATA_3 =  ((1,2,3.0,4.0),(5,6,7.0,8.0))

TEST_DATA_4 =  ((1,2,3.0,'bad'),(5,6,7.0,8.0))

TEST_DATA_5  = ((1,2,3,4.0,5.0),(6,7,8,9.0,10.0))

EXPECTED_3 = (1, 2, 3, 4.0, 5.0), (2, 4, 6, 8.0, 10.0)
        

class Test_component_list(unittest2.TestCase):


    def setUp(self):
        self._component_list =  Component_list()

    def testEmptyComponentList(self):
        with self.assertRaises(Exception):
            self._component_list.get_components_for__atom_id(0)
    
    def testExpectedRanges(self):
        
        self._component_list.add_components(TEST_DATA_1)
        for atom_id,expected in zip (EXPECTED_COMPONENTS_1,EXPECTED_1):
            result = self._component_list.get_components_for_atom_id(atom_id)
            self.assertSequenceEqual(result, expected)
    
    def testExpectedRaangesAfterAdd(self):
        self._component_list.add_components(TEST_DATA_1)
        self._component_list.add_components(TEST_DATA_2)
        for atom_id,expected in zip (EXPECTED_COMPONENTS_2,EXPECTED_2):
            result = self._component_list.get_components_for_atom_id(atom_id)
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
    
    def testIterComponents(self):
        self._component_list.add_components(TEST_DATA_1)
        
        for expected,component in zip(EXPECTED_ALL_1,self._component_list):
            self.assertSequenceEqual(expected, component)
            
    def testGetItems(self):
        self._component_list.add_components(TEST_DATA_1)
        
        for i,expected in enumerate (EXPECTED_ALL_1):
            result  = self._component_list[i]
            self.assertSequenceEqual(expected, result)

    def testGetItemsBeyondRange(self):
        self._component_list.add_components(TEST_DATA_1)
        
        with self.assertRaises(IndexError):
            self._component_list[len(TEST_DATA_1)]
            
    def testLen(self):
        
        self.assertEqual(len(self._component_list),0)
        self._component_list.add_components(TEST_DATA_1)
        self.assertEqual(len(self._component_list),len(TEST_DATA_1))
        
    def test_struct_translations(self):

        distance_component_struct = Struct('iiiff')
        struct_size = distance_component_struct.size
        bytes = ctypes.create_string_buffer(struct_size*2) 
        for j in range(2):
            i = j+1
            distance_component_struct.pack_into(bytes,struct_size*j, 1*i,2*i,3*i,4*i,5*i)
        result = test_dump_dist_comp(bytes)
        expected = EXPECTED_3
        self.assertEqual(result, expected)
    
    def test_get_data_type(self):
        self._component_list.add_components(TEST_DATA_3)
        result = self._component_list._get_struct_type()
        
        self.assertEqual(result, 'iiff')
    
    def test_get_bad_data_type(self):
        self._component_list.add_components(TEST_DATA_4)
        with self.assertRaises(Exception):
            result = self._component_list._get_struct_type()
        
    def test_get_data_type_empty_list(self):
        with self.assertRaises(Exception):
            result = self._component_list._get_struct_type()    
    
    def test_get_native_component(self): 
        self._component_list.add_components(TEST_DATA_5)
        result = self._component_list.get_native_components()
        
        result = test_dump_dist_comp(result)
        expected = TEST_DATA_5
        self.assertEqual(result, expected)
        
               
        
if __name__ == "__main__":
    unittest2.main()
