'''
Created on 25 Jan 2012

@author: garyt
'''
import unittest2
from test.test_segment_manager import Test_segment_manager
from test.test_observed_chemical_shifts import TestObservedShiftTable
from test.test_xcamshift import TestXcamshift
from test.Test_table_manager import Test_table_manager
from test.test_component_list import Test_component_list
from test.test_xcamshift_afa import TestXcamshiftAFA
from test.test_xcamshift_a4 import TestXcamshiftA4
from test.test_python_utils import Test_python_utils
from test.test_table_importers import Test_table_importers


def load_tests(loader, tests, pattern):
    suite = unittest2.TestSuite()
    suite.addTests(tests = loader.loadTestsFromTestCase(TestObservedShiftTable))
    suite.addTests(tests = loader.loadTestsFromTestCase(Test_segment_manager))
    suite.addTests(tests = loader.loadTestsFromTestCase(TestXcamshift))
    suite.addTests(tests = loader.loadTestsFromTestCase(Test_table_manager))
    suite.addTests(tests = loader.loadTestsFromTestCase(Test_component_list))
    suite.addTests(tests = loader.loadTestsFromTestCase(TestXcamshiftAFA))
    suite.addTests(tests = loader.loadTestsFromTestCase(TestXcamshiftA4))
    suite.addTests(tests = loader.loadTestsFromTestCase(Test_python_utils))
    suite.addTests(tests = loader.loadTestsFromTestCase(Test_table_importers))
    return suite

if __name__=="__main__":
    
    unittest2.main(verbosity=10)
    

