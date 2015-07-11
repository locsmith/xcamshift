'''
Created on 21 Feb 2013

@author: garyt
'''
import unittest2
from cache_manager import Cache_manager


class Test(unittest2.TestCase):

    def setUp(self):
         Cache_manager._cache_manager=None
       

    def assertSameInstance(self, instance_1, instance_2):
        return self.assertEqual(instance_1, instance_2)

    def test_get_cache_manager_returns_a_singleton(self):
        cache_manager_1 = Cache_manager.get_cache_manager()
        cache_manager_2 = Cache_manager.get_cache_manager()
        self.assertSameInstance(cache_manager_1, cache_manager_2)
    
    
    def test_requesting_a_non_exitent_item_raises_an_exception(self):
        with self.assertRaises(Exception):
            Cache_manager.get_cache_manager().get_cache('TEST') 
    
    def test_builder_creates_cache(self):
        EXPECTED = {'1':2} 
        TEST = 'TEST'
        def test_builder():
            return EXPECTED
        
        cache_manager = Cache_manager.get_cache_manager()
        cache_manager.add_cache_builder(TEST, test_builder)
        self.assertEqual(EXPECTED,cache_manager.get_cache(TEST))
        
    
    def test_cache_manager_clear_cache_calls_builder(self):
        TEST = 'TEST'
        EXPECTED_0 = {0:0}
        EXPECTED_1 = {1:1}
        class Test_builder:
            def __init__(self):
                self._calls = 0
                
            def __call__(self):
                result = {self._calls:self._calls}
                self._calls += 1
                return result
            
        cache_manager = Cache_manager.get_cache_manager()
        cache_manager.add_cache_builder(TEST, Test_builder())
        
        self.assertEqual(EXPECTED_0,cache_manager.get_cache(TEST))
        self.assertEqual(EXPECTED_0,cache_manager.get_cache(TEST))   
         
        cache_manager.clear_cache()
        
        self.assertEqual(EXPECTED_1,cache_manager.get_cache(TEST))
        self.assertEqual(EXPECTED_1,cache_manager.get_cache(TEST))   

    def test_adding_a_builder_twice_raises_an_exception(self):
        TEST = 'TEST'
        def test_builder():
            return EXPECTED
        
        cache_manager = Cache_manager.get_cache_manager()
        cache_manager.add_cache_builder(TEST, test_builder)
        
        with self.assertRaises(Exception):
            cache_manager.add_cache_builder(TEST, test_builder)      

if __name__ == "__main__":
    unittest2.main()