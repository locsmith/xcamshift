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
Created on 31 Dec 2011
 
@author: garyt
'''
from vec3 import Vec3
from simulation import currentSimulation,makeCurrent
from derivList import DerivList
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from xcamshift import  Xcamshift 
from shift_calculators import CDSSharedVectorFloat
from atomSel import AtomSel
from test import gb3, gb3_10_steps
from test import  util_for_testing
from observed_chemical_shifts import Observed_shift_table
from common_constants import  SIDE_CHAIN, NON_BONDED, RING,\
    CAMSHIFT_SUB_POTENTIALS
from common_constants import TARGET_ATOM_IDS_CHANGED, ROUND_CHANGED
from test.gb3 import gb3_component_shifts_sc, gb3_component_shifts_ring
from utils import Atom_utils
import common_constants
import sys
from cython.fast_segment_manager import Segment_Manager
from test.util_for_testing import  _check_shift_results, _shift_cache_as_result,\
    get_atom_index, get_key_for_atom_index, keyed_atoms_to_list
TOTAL_ENERGY = 'total'
from cython.shift_calculators import Out_array
from array import array
from ensembleSimulation import EnsembleSimulation
from cython.shift_calculators import CDSVectorFloat
 
def almostEqual(first, second, places = 7):
    result  = False
    if round(abs(second-first), places) == 0:
        result=True
    return result
 
ZEROS_3 = [0.0] * 3
class TestXcamshiftEnsembleGB3(unittest2.TestCase):

    def __init__(self,*args,**kwargs):
        super(TestXcamshiftEnsembleGB3, self).__init__(*args,**kwargs)
        self._esim = None
        
    
    def get_single_member_ensemble_simulation(self):
        if self._esim.__class__ ==  None.__class__:
            #TODO note EnsembleSimulation can't have a single member that causes a crash!
            # therefore a hack
            self._esim =  Xcamshift().ensembleSimulation()
        return self._esim


    # TODO add extra places
    DEFAULT_DECIMAL_PLACES = 5
    DEFAULT_ERROR = 10**-DEFAULT_DECIMAL_PLACES

#     #TODO: this should be triggered on structure change events
    def _clear_caches(self):
        Atom_utils.clear_cache()
        Segment_Manager.reset_segment_manager()
  
    def setUp(self):
        initStruct("test_data/gb3/gb3.psf")
        PDBTool("test_data/gb3/gb3_refined_II.pdb").read()
        self._clear_caches()
        print "In method", self._testMethodName
        self._ensembleSize=2
        self._esim = EnsembleSimulation('ensemble',self._ensembleSize)
          
    def _setup_xcamshift_with_shifts_table(self, test_shifts):
        xcamshift = Xcamshift()
        xcamshift.remove_named_sub_potential('HBOND')
        observed_shifts = Observed_shift_table(test_shifts)
        xcamshift.set_observed_shifts(observed_shifts)
        return xcamshift
   
    def test_shift_averaging_identical_structures(self):
        xcamshift = self._setup_xcamshift_with_shifts_table(gb3.gb3_zero_shifts)
        energy = xcamshift.calcEnergy()
         
        self.assertAlmostEqual(energy, gb3_10_steps.gb3_energies[0][TOTAL_ENERGY], places = 0)
        active_target_atoms_ids =  xcamshift._get_active_target_atom_ids()
        measured = xcamshift._shift_cache
        expected = gb3.gb3_shifts
        _check_shift_results(active_target_atoms_ids, measured,expected)
     
 
    def _get_sim_index(self):
        member_ids = [self._esim.members(i).id() for i in range(self._esim.size())]
        index = member_ids.index(self._esim.member().id())
        return index
 
    def set_ensemble_member_coords(self, coordinates):
         
        orig_simulation =  currentSimulation()
        index = self._get_sim_index()
 
        makeCurrent(self._esim.member())
        PDBTool("test_data/gb3_10_steps/%s" % coordinates[index]).read()
        makeCurrent(orig_simulation)
 
    def _select_shift_subset(self,esim,ensemble_shift_cache,index):
        data = [elem for elem in ensemble_shift_cache]
        result = data[index::esim.size()]
        return result
 
    def test_shift_averaging_two_structures(self):
        
        index = self._get_sim_index()
 
        xcamshift = self._setup_xcamshift_with_shifts_table(gb3.gb3_zero_shifts)
        print xcamshift.ensembleSimulation().size(),index
         
        files = gb3_10_steps.gb3_files[:self._ensembleSize]
        self.set_ensemble_member_coords(files)
         
        energy = xcamshift.calcEnergy()
         
        active_target_atoms_ids =  xcamshift._get_active_target_atom_ids()
        measured = self._select_shift_subset(self._esim, xcamshift._ensemble_shift_cache, index)
        expected_dict = gb3_10_steps.gb3_shifts[index]
 
        _check_shift_results(active_target_atoms_ids, measured,expected_dict)
         
        measured =  xcamshift._shift_cache
        expected_1 = keyed_atoms_to_list(active_target_atoms_ids, gb3_10_steps.gb3_shifts[0])
        expected_2 = keyed_atoms_to_list(active_target_atoms_ids, gb3_10_steps.gb3_shifts[1])
 
        for atom_id,value_1,value_2 in zip(active_target_atoms_ids,expected_1,expected_2):
            key =  get_key_for_atom_index(atom_id)
            expected_dict[key]=(value_1+value_2)/2.0
         
        _check_shift_results(active_target_atoms_ids, measured,expected_dict)
    
if __name__ == "__main__":
#     TODO: add a way to run the complete test suite
      unittest2.main(module='test.test_xcamshift_ensemble_gb3')