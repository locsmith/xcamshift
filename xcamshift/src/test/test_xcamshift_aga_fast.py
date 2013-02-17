#-------------------------------------------------------------------------------
# Copyright (c) 2013 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 8 Aug 2012

@author: garyt
'''
from test import test_xcamshift_aga

if __name__ == '__main__':
    test_xcamshift_aga.fast=True
    test_xcamshift_aga.run_tests()
