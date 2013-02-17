#-------------------------------------------------------------------------------
# Copyright (c) 2013 gary thompson.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
from cython.test_distance import get_distance, test_vec3_pointer
import unittest2


class Test_cython(unittest2.TestCase):
    
    # TODO add extra places
    DEFAULT_DECIMAL_PLACES = 5
    
    def test_distance(self):
        self.assertAlmostEqual(1.03922182425, get_distance(), self.DEFAULT_DECIMAL_PLACES)
        
    def test_vec3_pointer(self):
        self.assertSequenceEqual(test_vec3_pointer(), (1.0,2.0,3.0))

if __name__ ==  '__main__':
    unittest2.main()
     





