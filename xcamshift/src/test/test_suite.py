'''
Created on 25 Jan 2012

@author: garyt
'''
import unittest2
from test.test_segment_manager import Test_segment_manager
from test.test_observed_chemical_shifts import TestObservedShiftTable
from test.test_xcamshift import TestXcamshift
from test.Test_table_manager import Test_table_manager


def load_tests(loader, tests, pattern):
    suite = unittest2.TestSuite()
    suite.addTests(tests = loader.loadTestsFromTestCase(TestObservedShiftTable))
    suite.addTests(tests = loader.loadTestsFromTestCase(Test_segment_manager))
    suite.addTests(tests = loader.loadTestsFromTestCase(TestXcamshift))
    suite.addTests(tests = loader.loadTestsFromTestCase(Test_table_manager))
    return suite

if __name__=="__main__":
    
    unittest2.main(verbosity=10)
    
