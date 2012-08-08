'''
Created on 1 Aug 2012

@author: garyt
'''
from test import test_suite
from test.test_suite import load_tests as target_load_tests

load_tests = target_load_tests

if __name__ == '__main__':
    test_suite.fast=True
    test_suite.run_tests()