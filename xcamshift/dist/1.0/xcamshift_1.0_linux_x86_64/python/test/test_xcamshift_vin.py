#-------------------------------------------------------------------------------
# Copyright (c) 2013-2015 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on Apr 23, 2012

@author: garyt
'''
from common_constants import BACK_BONE
from pdbTool import PDBTool
from protocol import initStruct
from test.vin import vin_shifts, vin_subpotential_shifts
from utils import Atom_utils
from xcamshift import Xcamshift
import sys
import unittest2
from cython.fast_segment_manager import Segment_Manager


def almostEqual(first, second, places = 7):
    result  = False
    if round(abs(second-first), places) == 0:
        result=True
    return result

class TestXcamshiftVIN(unittest2.TestCase):
    DEFAULT_DECIMAL_PLACES = 5
    DEFAULT_ERROR = 10 ** -DEFAULT_DECIMAL_PLACES

    def remove_zero_valued_keys(self, expected_force_factors):
        for key, value in expected_force_factors.items():
            if almostEqual(value, 0.0, self.DEFAULT_DECIMAL_PLACES):
                 del expected_force_factors[key]
    
    def assertEmpty(self, expected_force_factors,msg = None):
        return self.assertEqual(len(expected_force_factors), 0)

                
    def setUp(self):
        initStruct("test_data/vin/vin.psf")
        PDBTool("test_data/vin/vin.pdb").read()
        Segment_Manager.reset_segment_manager()
        Atom_utils.clear_cache()
        

    def _get_xcamshift(self):
        xcamshift = Xcamshift()
        xcamshift.remove_named_sub_potential('HBOND')
        return xcamshift

    def test_vin_shifts(self):
        xcamshift = self._get_xcamshift()
        
        vin_shifts_copy = dict(vin_shifts)
        for key in vin_shifts:
            segment, residue_number, atom = key
            
            atom_ids = Atom_utils.find_atom_ids(segment, residue_number, atom)
            
            shift = xcamshift._calc_single_atom_shift(atom_ids[0])
            expected_shift = vin_shifts_copy[key]
            
            self.assertAlmostEqual(expected_shift, shift, places=self.DEFAULT_DECIMAL_PLACES - 2, msg=`key`)
            del vin_shifts_copy[key]
        self.assertEmpty(vin_shifts_copy)
        
    def test_vin_component_shifts(self):
        xcamshift = self._get_xcamshift()
        
        sub_potential_shifts = dict(vin_subpotential_shifts)
        for key in vin_subpotential_shifts:
            segment, residue_number, atom, sub_potential_name = key

            atom_ids = Atom_utils.find_atom_ids(segment, residue_number, atom)
            
            sub_potential = xcamshift.get_named_sub_potential(sub_potential_name)
            
            shift = sub_potential._calc_single_atom_shift(atom_ids[0])
            expected_shift = vin_subpotential_shifts[key]

            self.assertAlmostEqual(expected_shift, shift, places=self.DEFAULT_DECIMAL_PLACES - 1, msg=`key`)
            
            del sub_potential_shifts[key]
            
        self.assertEmpty(sub_potential_shifts)

            
#    def test_component_shifts_bb(self):
#        
#        xcamshift = Xcamshift()
#        bb_subpotential = xcamshift.get_named_sub_potential(BACK_BONE)
#        
#        expected_bb_shifts = dict(agaga_component_shifts_bb)
#        expected_component_keys = expected_bb_shifts.keys()
#        for component_index, component in enumerate(bb_subpotential._get_component_list()):
#            from_atom_id, to_atom_id = component[0:2]
#            from_atom_key = Atom_utils._get_atom_info_from_index(from_atom_id)
#            to_atom_key = Atom_utils._get_atom_info_from_index(to_atom_id)
#            
#            expected_key = from_atom_key, to_atom_key
#
#            self.assertIn(expected_key, expected_component_keys, `expected_key` + " exists")
#            
#            shift = bb_subpotential._calc_component_shift(component_index)
#            
#            self.assertAlmostEqual(expected_bb_shifts[expected_key], shift, places=self.DEFAULT_DECIMAL_PLACES - 1)
#            
#            del expected_bb_shifts[expected_key]
#            
#        self.remove_zero_valued_keys(expected_bb_shifts)
#        self.assertEmpty(expected_bb_shifts)
#        
#        
#    def test_test_shifts_add_up(self):
#        
#        test_data_sum = {}
#        for key in agaga_subpotential_shifts:
#            sub_key = key[:3]
#
#            test_data_sum.setdefault(sub_key, 0.0)
#            test_data_sum[sub_key] += agaga_subpotential_shifts[key]
#            
#        for key in test_data_sum:
#            #TODO: correct this its due to limited output formats?
#            self.assertAlmostEqual(test_data_sum[key], agaga_shifts[key], places=self.DEFAULT_DECIMAL_PLACES - 2, msg=key)
            
def run_tests():
    unittest2.main(module='test.test_xcamshift_vin')
#    unittest2.main(module='test.test_xcamshift_afa',defaultTest='TestXcamshiftAFA.test_calc_component_shift')
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestXcamshifAGA.testName']
    run_tests()
    
