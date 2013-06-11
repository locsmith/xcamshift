'''
Created on 20 Apr 2013

@author: garyt
'''
import unittest
from cython.shift_calculators import  Non_bonded_interaction_list, test_dump_component_index_pair

class Test_cython_non_boned_list(unittest.TestCase):

    def setUp(self):
        self.test_list = Non_bonded_interaction_list()
        
    def test_create(self):
        self.assertEqual(len(self.test_list),0)
        
    def test_append(self):
        self.test_list.test_append(0,1,2,3)
        self.assertEqual(self.test_list[0], (0,1,2,3))
        self.assertEqual(len(self.test_list), 1)
        
        self.test_list.test_append(4,5,6,7)
        self.assertEqual(self.test_list[1], (4,5,6,7))
        self.assertEqual(len(self.test_list), 2)


    def _calc_expected_allocation(self, test_list):
        return 6 * test_list.get_size_increment() * 4
    

    def _get_test_data_elem(self, offset, i):
        return tuple([((i + offset) * 4) + j for j in range(4)])

    def _build_101(self,offset=0):
        for i in range(101):
            test_data = self._get_test_data_elem(offset, i)
            self.test_list.test_append(*test_data)

    
    def _test_101(self, test_list, offset=0):
        for i in range(101):
            test_data = self._get_test_data_elem(offset, i)
            self.assertEqual(test_list[i],test_data)
        
        self.assertEqual(len(test_list), 101)
        if hasattr(test_list,'get_allocation'):
            self.assertEqual(test_list.get_allocation(), self._calc_expected_allocation(self.test_list))

    def test_allocation(self):
        
        self._build_101() 
        self._test_101(self.test_list)
    
    def test_clear(self):
        self._build_101() 
        self._test_101(self.test_list)
        self.assertEqual(self.test_list.get_allocation(), self._calc_expected_allocation(self.test_list))
        
        self.test_list.clear() 
        self.assertEqual(self.test_list.get_allocation(), self._calc_expected_allocation(self.test_list))
                
        self._build_101(1000) 
        self._test_101(self.test_list,1000)
        self.assertEqual(self.test_list.get_allocation(), self._calc_expected_allocation(self.test_list))
            

    def test_get_native(self):
        self._build_101() 
        
        for i in range(100):
            result  = test_dump_component_index_pair(self.test_list,i)
            self.assertEqual(result, self._get_test_data_elem(0, i))
#         
    def test_get_all_components(self):
        self._build_101()
        self._test_101(self.test_list.get_all_components())

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()