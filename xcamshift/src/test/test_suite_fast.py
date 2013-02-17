#-------------------------------------------------------------------------------
# Copyright (c) 2013 gary thompson.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v2.1
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
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
