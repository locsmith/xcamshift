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
from test import test_xcamshift_gb3
import cProfile


if __name__ == '__main__':
    test_xcamshift_gb3.fast=True
    test_xcamshift_gb3.run_tests()
#    cProfile.run('test_xcamshift_gb3.run_tests()')
