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
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from xcamshift import Ring_Potential,Xcamshift
from atomSel import AtomSel
from test.xdists import xdists_ala_3
from test.dihedrals import dihedrals_ala_3
from cython.fast_segment_manager import Segment_Manager
from test.sidechains import sidechains_ala_3
from test.ala_3 import ala_3_total_shifts
from test import sidechains, ala_3, AFA
from observed_chemical_shifts import Observed_shift_table
from utils import Atom_utils
import sys
from table_manager import Table_manager
from vec3 import Vec3
from common_constants import TARGET_ATOM_IDS_CHANGED, ROUND_CHANGED
TOTAL_ENERGY = 'total'
from cython.shift_calculators import Out_array


def text_keys_to_atom_ids(keys, segment = '*'):
    result = []
    for key in keys:
        atom_ids = single_text_key_to_atom_ids(key, segment)
        result.extend(atom_ids)
    return result

def single_text_key_to_atom_ids(key,segment = "*"):
    residue_number, atom_name  =  key
    atoms = Atom_utils.find_atom(segment, residue_number, atom_name)
    return [atom.index() for atom in atoms]

#TODO: test the tests and isolate
def almostEqual(first, second, places = 7):
    result  = False
    if round(abs(second-first), places) == 0:
        result=True
    return result

#class testSegmentManager(object):
class TestXcamshiftAFA(unittest2.TestCase):

    def __init__(self,*args,**kwargs):
        super(TestXcamshiftAFA, self).__init__(*args,**kwargs)
        self._esim = None

    # TODO add extra places
    DEFAULT_DECIMAL_PLACES = 5
    DEFAULT_ERROR = 10**-DEFAULT_DECIMAL_PLACES


    def assertEmpty(self, expected_force_factors,msg = None):
        return self.assertEqual(len(expected_force_factors), 0)

    def check_almost_equal(self, list_1, list_2, delta = 1e-7):
        difference_offset = -1
        for i, (elem_1, elem_2) in enumerate(zip(list_1, list_2)):
            diff = abs(elem_1 - elem_2)
            if diff > delta:
                difference_offset = i
                break
        return difference_offset

    def are_almost_equal_sequences(self, list_1, list_2, delta =  1e-7):
        result = True
        if self.check_almost_equal(list_1, list_2, delta) > 0:
            result = False
        return result

    def assertSequenceAlmostEqual(self,result,expected,places = 7):
        delta  = 10**-places
        len_result = len(result)
        len_expected = len(expected)
        if len_result != len_expected:
            raise AssertionError("the two lists are of different length %i and %i" % (len_result,len_expected))

        difference_offset = self.check_almost_equal(result, expected, delta)

        if difference_offset > -1:
            template = "lists differ at item %i: %s - %s > %s"
            elem_1 = result[difference_offset]
            elem_2 = expected[difference_offset]
            message = template % (difference_offset, `elem_1`,`elem_2`,delta)
            raise AssertionError(message)

    def get_single_member_ensemble_simulation(self):
        if self._esim.__class__ ==  None.__class__:
            #TODO note EnsembleSimulation can't have a single member that causes a crash!
            # therefore a hack
            self._esim =  Xcamshift().ensembleSimulation()
        return self._esim

    def setUp(self):
        initStruct("test_data/afa/AFA.psf")
        PDBTool("test_data/afa/AFA.pdb").read()
        Atom_utils.clear_cache()
        Segment_Manager.reset_segment_manager()
#         print "In method", self._testMethodName

#TODO: shoulf be private
    def make_out_array(self):
#        TODO: use segment manager
        num_atoms = len(AtomSel('(all)').indices())
        result = Out_array(num_atoms,self.get_single_member_ensemble_simulation())
        return result

    def make_result_array(self):
        num_atoms = len(AtomSel('(all)').indices())
        result = [0.0] * num_atoms
        return result

    def test_ring_table_names(self):
        ring_potential = self.make_ring_potential()

        expected_table_names = ('ATOM','COEF','RING')
        table_names = ring_potential.get_component_table_names()
        self.assertSequenceEqual(expected_table_names,table_names)

    def test_ring__bb_atoms(self):
        ring_potential = self.make_ring_potential()

        ring_table = Table_manager.get_default_table_manager().get_ring_table('PHE')
        all_components = ring_potential._get_component_list('COEF').get_all_components()
        for atom_id,ring_id,coefficient in all_components:
            self.assertEqual(ring_id,0)
            target_atoms = ring_table.get_target_atoms()
            atom_name = target_atoms[atom_id]
            expected_coefficient = ring_table.get_ring_coefficient(atom_name,'PHE','6')
            self.assertAlmostEqual(coefficient, expected_coefficient, self.DEFAULT_DECIMAL_PLACES)

    def test_ring_sidechain_atoms(self):
        ring_potential = self.make_ring_potential()

        ATOM_NAMES_PHE_RING = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]

#        ring_table = Table_manager.get_default_table_manager().get_ring_table('PHE')
        all_components = ring_potential._get_component_list('RING').get_all_components()
        for ring_id,ring_atom_ids in all_components:
            self.assertEqual(ring_id,0)
            for ring_atom_id in ring_atom_ids:
                segment, residue_number,atom_name = Atom_utils._get_atom_info_from_index(ring_atom_id)
                self.assertEqual(segment,'')
                self.assertEqual(residue_number,2)
                self.assertIn(atom_name, ATOM_NAMES_PHE_RING)
                del ATOM_NAMES_PHE_RING[ATOM_NAMES_PHE_RING.index(atom_name)]
            self.assertEmpty(ATOM_NAMES_PHE_RING)

    def test_ring_bb_atoms(self):
        ring_potential = self.make_ring_potential()



        ring_table = Table_manager.get_default_table_manager().get_ring_table('PHE')
        target_atoms = ring_table.get_target_atoms()

        all_component_names = ring_potential._get_component_list('ATOM').get_all_components()
        all_component_names_index = dict((i,atom_name) for i,atom_name in enumerate(target_atoms))
        for atom_id,atom_type_id in all_component_names:
            segment,residue_number,atom_name = Atom_utils._get_atom_info_from_index(atom_id)
            self.assertEqual(segment, '')
            self.assertEqual(residue_number, 2)
            atom_name  = Atom_utils._get_atom_name_from_index(atom_id)
            self.assertEqual(all_component_names_index[atom_type_id], atom_name)

            del all_component_names_index[atom_type_id]
        self.assertEmpty(all_component_names_index)


    def assertListVec3AlmostEqual(self, ring_centres, expected, places = DEFAULT_DECIMAL_PLACES):
        self.assertEqual(len(ring_centres), 1)
        for result, expected in zip(ring_centres, expected):
            for item1, item2 in zip(result, expected):
                self.assertAlmostEqual(item1, item2, places)

    def test_ring_calculate_centre(self):
        ring_potential = self.make_ring_potential()

        ring_potential. _build_ring_data_cache()
        ring_centres = [elem[1] for elem in ring_potential._get_cache_list('CENT')]
        expected = [Vec3(*AFA.expected_ring_centre),]

        self.assertListVec3AlmostEqual(ring_centres, expected)

    def test_ring_calculate_normals(self):
        ring_potential = self.make_ring_potential()

        ring_potential._build_ring_data_cache()
        ring_normals = [elem[1] for elem in ring_potential._get_cache_list('NORM')]
        expected = [Vec3(*AFA.expected_ring_normals),]

        #TODO: check difference why only 3 places consistent, different atom selections for normals?
        self.assertListVec3AlmostEqual(ring_normals, expected, self.DEFAULT_DECIMAL_PLACES)

    def get_target_atom_ids(self, components):
        return  []



    def _get_potential_target_atom_ids(self, ring_potential):
        target_components = ring_potential._get_component_list('ATOM')
        target_atom_ids = [component[0] for component in target_components]
        return target_atom_ids

    def _prepare_potential(self, ring_potential):
        target_atom_ids = self._get_potential_target_atom_ids(ring_potential)
        ring_potential._prepare(TARGET_ATOM_IDS_CHANGED, target_atom_ids)
        ring_potential._prepare(ROUND_CHANGED, None)


    def _get_potential_target_components(self, potential):
        return potential._get_component_list('ATOM')

    def test_calc_component_shift(self):
        ring_potential = self.make_ring_potential()
        self._prepare_potential(ring_potential)
        target_components = self._get_potential_target_components(ring_potential)


        for i, component in enumerate(target_components):
            component_shift = ring_potential._calc_component_shift(i)
            #TODO: this looks inside the object too much
            atom_id = component[0]
            atom_key = Atom_utils._get_atom_info_from_index(atom_id)
            expected_shift = AFA.expected_ring_shifts[atom_key]
            #TODO: is the difference in error down to coordinates?
            self.assertAlmostEqual(component_shift, expected_shift, places=self.DEFAULT_DECIMAL_PLACES)


    def make_ring_potential(self):
        ring_potential = Ring_Potential(self.get_single_member_ensemble_simulation())

        return ring_potential

    def test_calc_component_forces(self):
        ring_potential = self.make_ring_potential()
        ring_potential._setup_ring_calculator(ring_potential._force_calculator)
        ring_potential._force_calculator._set_components(ring_potential._get_components())
        for target_atom_id,atom_type_id in ring_potential._get_component_list('ATOM'):
            for ring_id,ring_atom_ids in ring_potential._get_component_list('RING'):

                target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
                force_factor_key = target_atom_key
                force_factor = AFA.force_factors_harmonic[force_factor_key]

                #TODO: these should not be needed
                ring_potential._get_component_list('RING')
                ring_potential._build_ring_data_cache()

                out_array = self.make_out_array()
                force_terms = ring_potential._force_calculator._build_force_terms(target_atom_id, ring_id)


                ring_potential._force_calculator._calc_target_atom_forces(target_atom_id, force_factor, force_terms, out_array)
                forces =out_array.add_forces_to_result()

                target_atom_forces = forces[target_atom_id]
                expected_forces = AFA.target_forces_harmonic[force_factor_key]
                self.assertSequenceAlmostEqual(target_atom_forces, expected_forces, self.DEFAULT_DECIMAL_PLACES)

                out_array = self.make_out_array()
                ring_potential._force_calculator._calculate_ring_forces(atom_type_id, ring_id, force_factor, force_terms, out_array)
                forces =out_array.add_forces_to_result()

                for ring_atom_id in ring_atom_ids:
                    ring_atom_key  = Atom_utils._get_atom_info_from_index(ring_atom_id)
                    ring_force_key =target_atom_key,ring_atom_key
                    expected_ring_forces = AFA.ring_forces_harmonic[ring_force_key]
                    ring_atom_forces = forces[ring_atom_id]
                    self.assertSequenceAlmostEqual(ring_atom_forces, expected_ring_forces, self.DEFAULT_DECIMAL_PLACES)


def run_tests():
    unittest2.main(module='test.test_xcamshift_afa')
#    unittest2.main(module='test.test_xcamshift_afa',defaultTest='TestXcamshiftAFA.test_calc_component_shift')

if __name__ == "__main__":
    run_tests()
#    TestXcamshift.list_test_shifts()
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testCalcForceSetTanh')
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testSingleFactorHarmonic')
