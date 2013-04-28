'''
Created on 20 Apr 2013

@author: garyt
'''
import unittest
from cython.shift_calculators import  Non_bonded_list, test_dump_component_index_pair

class Test_cython_non_boned_list(unittest.TestCase):

    def setUp(self):
        self.test_list = Non_bonded_list()
        
    def test_create(self):
        self.assertEqual(len(self.test_list),0)
        
    def test_append(self):
        self.test_list.test_append(0,1,2)
        self.assertEqual(self.test_list[0], (0,1,2))
        self.assertEqual(len(self.test_list), 1)
        
        self.test_list.test_append(3,4,5)
        self.assertEqual(self.test_list[1], (3,4,5))
        self.assertEqual(len(self.test_list), 2)


    def _calc_expected_allocation(self, test_list):
        return 6 * test_list.get_size_increment() * 3
    
    def _build_101(self,offset=0):
        for i in range(101):
            self.test_list.test_append((i+offset) * 2, ((i+offset) * 2) + 1, ((i+offset) * 2) + 2)


    def _test_101(self,offset=0):
        for i in range(101):
            self.assertEqual(self.test_list[i], ((i+offset) * 2, ((i+offset) * 2) + 1, ((i+offset) * 2) + 2))
        
        self.assertEqual(len(self.test_list), 101)
        self.assertEqual(self.test_list.get_allocation(), self._calc_expected_allocation(self.test_list))

    def test_allocation(self):
        
        self._build_101() 
        self._test_101()
    
    def test_reset(self):
        self._build_101() 
        self._test_101()
        self.assertEqual(self.test_list.get_allocation(), self._calc_expected_allocation(self.test_list))
        
        self.test_list.reset() 
        self.assertEqual(self.test_list.get_allocation(), self._calc_expected_allocation(self.test_list))
                
        self._build_101(1000) 
        self._test_101(1000)
        self.assertEqual(self.test_list.get_allocation(), self._calc_expected_allocation(self.test_list))
            

    def test_get_native(self):
        self._build_101() 
        
        for i in range(100):
            target,remote  = test_dump_component_index_pair(self.test_list,i)
            self.assertEqual(target, i*2+1)
            self.assertEqual(remote, i*2+2)
        
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()