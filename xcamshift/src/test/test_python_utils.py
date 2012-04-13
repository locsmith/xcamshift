'''
Created on 17 Jan 2012

@author: garyt
'''
import unittest2
from python_utils import tupleit, IsMappingType, Dict_walker
from UserDict import UserDict
from numpy.distutils.misc_util import dict_append


class Test_python_utils(unittest2.TestCase):


    def testTupleIt(self):
        test = [[1, 2], [3, 4]]
       
        result = tupleit(test)
       
        self.assertTrue(isinstance(result, tuple))
        self.assertTrue(isinstance(result[0], tuple))
        self.assertTrue(isinstance(result[1], tuple))
        
        
    def testMappingChecker(self):
        
        self.assertTrue(IsMappingType({}))
        self.assertFalse(IsMappingType([]))
        self.assertTrue(IsMappingType(UserDict()))
        
    def testDictWalker(self):
        

            
        test = {1: {2:{3:4,5:6}},7:8}
        expected = {(1,2,3) :4, (1,2,5) : 6, (7,) : 8}
        
        walker = Dict_walker()
        
        aggregator = Dict_walker.Aggregator()
        walker.walk_dict(test, aggregator)
        
        
        self.assertDictEqual(expected, aggregator.get_result())
        
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test_python_utils.testName']
    unittest2.main()