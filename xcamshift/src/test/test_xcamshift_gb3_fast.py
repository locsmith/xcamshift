'''
Created on 1 Aug 2012

@author: garyt
'''
from test import test_xcamshift_gb3
import cProfile


if __name__ == '__main__':
    test_xcamshift_gb3.fast=True
    test_xcamshift_gb3.run_tests()
#    cProfile.run('test_xcamshift_gb3.run_tests()')
