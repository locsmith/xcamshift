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
Created on 27 Feb 2013

@author: garyt
'''
import unittest2
from protocol import initStruct
from pdbTool import PDBTool
from cython.fast_segment_manager import Segment_Manager
from utils import Atom_utils
from xcamshift import Disulphide_shift_calculator
from test.acaggaca import acaggaca_ss_shifts

class  TestXcamshiftAcaggaca(unittest2.TestCase):


    def setUp(self):
        initStruct("test_data/acaggaca/acaggaca.psf")
        PDBTool("test_data/acaggaca/acaggaca.pdb").read()
        Segment_Manager.reset_segment_manager()
        Atom_utils.clear_cache()

    def testDisulphideComponents(self):
        
        disulphide_calculator = Disulphide_shift_calculator()
        components = disulphide_calculator._get_component_list('ATOM')
        
        expected = [(Atom_utils.find_atom_ids(*elem)[0],acaggaca_ss_shifts[elem]) for elem in acaggaca_ss_shifts]
        expected.sort()  
        self.assertSequenceEqual([elem[0] for elem in components], [elem[0] for elem in expected])
       



if __name__ == "__main__":
    unittest2.main()