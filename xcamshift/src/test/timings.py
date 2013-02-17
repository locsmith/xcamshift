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
Created on 29 Jul 2012

@author: garyt
'''
import unittest
import timeit
from common_constants import NON_BONDED, RING
from test import gb3
from utils import Atom_utils
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from xcamshift import Xcamshift
from observed_chemical_shifts import Observed_shift_table
from atomSel import AtomSel
import common_constants


test_function = None


    
class Test(unittest2.TestCase):


    def setUp(self):
        initStruct("test_data/gb3/gb3.psf")
        PDBTool("test_data/gb3/gb3_refined_II.pdb").read()
        Segment_Manager.reset_segment_manager()
        Atom_utils.clear_cache()


    def tearDown(self):
        pass
    
    def _setup_xcamshift_with_shifts_table(self, test_shifts):
        xcamshift = Xcamshift()
        observed_shifts = Observed_shift_table(test_shifts)
        xcamshift.set_observed_shifts(observed_shifts)
        return xcamshift
    
    def make_xcamshift(self, shifts):
        xcamshift = self._setup_xcamshift_with_shifts_table(shifts)

        
        non_bonded = xcamshift.get_named_sub_potential(NON_BONDED)
        non_bonded.update_non_bonded_list()
        

        ring_potential = xcamshift.get_named_sub_potential(RING)
        ring_potential._build_ring_data_cache()
        
        return xcamshift
    
    translations = {('GLY','HA')  :'HA1',
                    ('ILE','CD')  : 'CD1',
                    ('ILE','HD1') : 'HD11',
                    ('ILE','HD2') : 'HD12',
                    ('ILE','HD3') : 'HD13'}
    
    def get_atom_index(self, key):
        residue_type = Atom_utils._get_residue_type(key[0], key[1])

        atom_name = key[2]
        
        residue_key = list(key)
        residue_key[2] =  self.translations.setdefault((residue_type,atom_name),atom_name)
        residue_key = tuple(residue_key)
        
        target_atom_index = Atom_utils.find_atom(*residue_key)[0].index()
        return target_atom_index


    def set_cache_shifts(self,xcamshift, shift_list):
        shift_cache =  xcamshift._shift_cache
        
        for key in shift_list:
            index = self.get_atom_index(key)
            shift_cache[index] = shift_list[key]
            
    def make_result_array_forces(self):
#        TODO: use segment manager
        num_atoms = len(AtomSel('(all)').indices())
        result = [None] * num_atoms
        return result
    
    def test_timings(self):

        xcamshift = self.make_xcamshift(gb3.gb3_zero_shifts)

        self.set_cache_shifts(xcamshift, gb3.gb3_shifts)
        
        non_bonded = xcamshift.get_named_sub_potential(NON_BONDED)
#        timings.test_function = self.doiit
#        test_function()
#        Test.doiit = self.done
        class Test_nonbonded:
            def __init__(self, non_bonded):
                self.non_bonded =  non_bonded
              
            def __call__(self):
                self.non_bonded.update_non_bonded_list()
                
        class Test_xcamshift_forces:
            def __init__(self,xcamshift,derivs):
                self.derivs = derivs
                self.xcamshift = xcamshift
            def __call__(self):
                self.xcamshift.calcEnergyAndDerivs(self.derivs)
                
        class Test_xcamshift_subpotential():
            def __init__(self,xcamshift,potential_list,derivs):
                self.derivs = derivs
                self.xcamshift = xcamshift
                self.potential_list = potential_list
                
            def __call__(self):
                self.xcamshift._calc_derivs(self.derivs,self.potential_list)
                

            
#        print 'non_bonded',timeit.Timer(Test_nonbonded(non_bonded)).timeit(number=10).0
#        print 'all',timeit.Timer(Test_xcamshift_forces(xcamshift, self.make_result_array_forces())).timeit(number=10)/10.0
        for sub_potential in common_constants.CAMSHIFT_SUB_POTENTIALS:
            potential_list =  [xcamshift.get_named_sub_potential(sub_potential)]
            print sub_potential,timeit.Timer(Test_xcamshift_subpotential(xcamshift, potential_list, self.make_result_array_forces())).timeit(number=10)/10.0
if __name__ == "__main__":
#    import sys;sys.argv = ['', 'Test.testName']
    unittest2.main()
