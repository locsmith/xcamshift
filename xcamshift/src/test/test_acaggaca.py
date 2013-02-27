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

    # TODO add extra places
    DEFAULT_DECIMAL_PLACES = 5

    def setUp(self):
        initStruct("test_data/acaggaca/acaggaca.psf")
        PDBTool("test_data/acaggaca/acaggaca.pdb").read()
        Segment_Manager.reset_segment_manager()
        Atom_utils.clear_cache()

    def check_almost_equal(self, list_1, list_2, delta = 1e-7):
        difference_offset = -1
        for i, (elem_1, elem_2) in enumerate(zip(list_1, list_2)):
            diff = abs(elem_1 - elem_2)
            if diff > delta:
                difference_offset = i
        
        return difference_offset

    def assertSequenceAlmostEqual(self,result,expected, delta = 1e-7, msg=""):
        len_result = len(result)
        len_expected = len(expected)
        if len_result != len_expected:
            raise AssertionError("the two lists are of different length %i and %i" % (len_result,len_expected))
        
        difference_offset = self.check_almost_equal(result, expected, delta)
        
            
        if difference_offset > 0:
            if msg != "":
                msg = msg + " "
                
            template = "%slists differ at item %i: %s - %s > %s"
            elem_1 = result[difference_offset]
            elem_2 = expected[difference_offset]
            message = template % (msg,difference_offset, `elem_1`,`elem_2`,delta)
            raise AssertionError(message)

    
    def make_result_array(self):
        num_atoms = len(AtomSel('(all)').indices())
        result = [0.0] * num_atoms
        return result
    
    def testDisulphideComponents(self):
        
        disulphide_calculator = Disulphide_shift_calculator()
        components = disulphide_calculator._get_component_list('ATOM')
        
        expected = [(Atom_utils.find_atom_ids(*elem)[0],acaggaca_ss_shifts[elem]) for elem in acaggaca_ss_shifts]
        expected.sort()  
        self.assertSequenceEqual([elem[0] for elem in components], [elem[0] for elem in expected])
       

    def testDisulphidePotential(self):
        initStruct("test_data/acaggaca/acaggaca.psf") 
        
        disulphide_calculator = Disulphide_shift_calculator()
        result = self.make_result_array() 
        expected  = self.make_result_array()
        for elem in acaggaca_ss_shifts:
            index = Atom_utils.find_atom_ids(*elem)[0]
            expected[index] += acaggaca_ss_shifts[elem]

        disulphide_calculator.set_shifts(result)
        self.assertSequenceAlmostEqual(expected, result, self.DEFAULT_DECIMAL_PLACES)


if __name__ == "__main__":
    unittest2.main()