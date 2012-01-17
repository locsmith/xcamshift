'''
Created on 17 Jan 2012

@author: garyt
'''
import unittest2
from utils import tupleit


class Test(unittest2.TestCase):


    def testTupleIt(self):
        test = [[1, 2], [3, 4]]
       
        result = tupleit(test)
       
        self.assertTrue(isinstance(result, tuple))
        self.assertTrue(isinstance(result[0], tuple))
        self.assertTrue(isinstance(result[1], tuple))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest2.main()