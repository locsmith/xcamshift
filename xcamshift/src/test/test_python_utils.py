'''
Created on 17 Jan 2012

@author: garyt
'''
import unittest2
from python_utils import tupleit, IsMappingType, Dict_walker,\
    value_from_key_path, filter_dict, Hierarchical_dict
from UserDict import UserDict
from numpy.distutils.misc_util import dict_append
from copy import deepcopy
import collections
import abc
from unittest2.util import safe_repr
import difflib
import pprint
import copy

TEST_DICT_DATA_1 = {1: {2:{3:4,5:6}},7:8}
        
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
        
        expected = {(1,2,3) :4, (1,2,5) : 6, (7,) : 8}
        
        walker = Dict_walker()
        
        aggregator = Dict_walker.Aggregator()
        walker.walk_dict(TEST_DICT_DATA_1, aggregator)
        
        
        self.assertDictEqual(expected, aggregator.get_result())
        
        
    def test_dict_key_path(self):
        
        self.assertEqual(value_from_key_path(TEST_DICT_DATA_1, (1,2,3)), 4)
        
        with(self.assertRaises(KeyError)):
            self.assertEqual(value_from_key_path(TEST_DICT_DATA_1, (1,2,2)), 4)
        
    def test_filter_dict(self):
        test_data  = deepcopy(TEST_DICT_DATA_1)
        
        filter_dict(test_data, lambda k,v: k==1)
        
        self.assertDictEqual(test_data, {7:8})
        
        test_data  = deepcopy(TEST_DICT_DATA_1)
        
        filter_dict(test_data, lambda k,v: k==1, invert=True)
        
        self.assertDictEqual(test_data,{1: {2: {3: 4, 5: 6}}})



    def assertDictEqual(self, d1, d2, msg=None):
        self.assertTrue(IsMappingType(d1), 'First argument is not a dictionary')
        self.assertTrue(IsMappingType(d2), 'Second argument is not a dictionary')

        if d1 != d2:
            standardMsg = '%s != %s' % (safe_repr(d1, True), safe_repr(d2, True))
            diff = ('\n' + '\n'.join(difflib.ndiff(
                           pprint.pformat(d1).splitlines(),
                           pprint.pformat(d2).splitlines())))
            standardMsg = self._truncateMessage(standardMsg, diff)
            self.fail(self._formatMessage(msg, standardMsg))
            
                

    def build_test_hier_dict(self):
        dict_0 = {1:-1, 3:-30}
        hdict_1 = Hierarchical_dict({2:-2}, parent=dict_0)
        hdict_2 = Hierarchical_dict({3:-3}, parent=hdict_1)
        return hdict_2, hdict_1, dict_0

    def test_hier_dict_lookup(self):
        hdict_2, hdict_1, dict_0 = self.build_test_hier_dict() #@UnusedVariable
        
        
        self.assertDictEqual(hdict_2, {1:-1,2:-2,3:-3})
        self.assertDictEqual(hdict_1, {1:-1,2:-2,3:-30})
        
        hdict_2, hdict_1, dict_0 = self.build_test_hier_dict()
        hdict_2[3]=3
        
        self.assertEqual(hdict_2[3], 3)
        self.assertEqual(dict_0[3], -30)

        hdict_2, hdict_1, dict_0 = self.build_test_hier_dict()
        
        del hdict_2[3]
        self.assertEqual(hdict_2[3], -30)
        
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test_python_utils.testName']
    unittest2.main()