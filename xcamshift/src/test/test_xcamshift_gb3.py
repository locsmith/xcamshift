'''
Created on 31 Dec 2011

@author: garyt
'''
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from xcamshift import Ring_Potential, Xcamshift
from atomSel import AtomSel
from test.xdists import xdists_ala_3
from test.dihedrals import dihedrals_ala_3
from segment_manager import Segment_Manager
from test.sidechains import sidechains_ala_3
from test.ala_3 import ala_3_total_shifts
from test import sidechains, ala_3, AFA, gb3
from observed_chemical_shifts import Observed_shift_table
from utils import Atom_utils
import sys
from table_manager import Table_manager
from vec3 import Vec3
from common_constants import BACK_BONE
TOTAL_ENERGY = 'total'
#def text_keys_to_atom_ids(keys, segment = '*'):
#    result = []
#    for key in keys:
#        atom_ids = single_text_key_to_atom_ids(key, segment)
#        result.extend(atom_ids)
#    return result
#
#def single_text_key_to_atom_ids(key,segment = "*"):
#    residue_number, atom_name  =  key
#    atoms = Atom_utils.find_atom(segment, residue_number, atom_name)
#    return [atom.index() for atom in atoms]
#    
##TODO: test the tests and isolate
#def almostEqual(first, second, places = 7):
#    result  = False
#    if round(abs(second-first), places) == 0:
#        result=True
#    return result
#
##class testSegmentManager(object):
class TestXcamshiftAFA(unittest2.TestCase):

    # TODO add extra places
    DEFAULT_DECIMAL_PLACES = 5
    DEFAULT_ERROR = 10**-DEFAULT_DECIMAL_PLACES
#
#    
#    def assertEmpty(self, expected_force_factors,msg = None):
#        return self.assertEqual(len(expected_force_factors), 0)
#    
#    def check_almost_equal(self, list_1, list_2, delta = 1e-7):
#        difference_offset = -1
#        for i, (elem_1, elem_2) in enumerate(zip(list_1, list_2)):
#            diff = abs(elem_1 - elem_2)
#            if diff > delta:
#                difference_offset = i
#                break
#        return difference_offset
#    
#    def are_almost_equal_sequences(self, list_1, list_2, delta =  1e-7):
#        result = True
#        if self.check_almost_equal(list_1, list_2, delta) > 0:
#            result = False
#        return result
#        
#    def assertSequenceAlmostEqual(self,result,expected,places = 7):
#        delta  = 10**-places
#        len_result = len(result)
#        len_expected = len(expected)
#        if len_result != len_expected:
#            raise AssertionError("the two lists are of different length %i and %i" % (len_result,len_expected))
#        
#        difference_offset = self.check_almost_equal(result, expected, delta)
#        
#        if difference_offset > -1:
#            template = "lists differ at item %i: %s - %s > %s"
#            elem_1 = result[difference_offset]
#            elem_2 = expected[difference_offset]
#            message = template % (difference_offset, `elem_1`,`elem_2`,delta)
#            raise AssertionError(message)
#            
    def setUp(self):
        initStruct("test_data/gb3/gb3.psf")
        PDBTool("test_data/gb3/gb3_refined_II_almost.pdb").read()
#
##TODO: shoulf be private
#    def make_result_array_forces(self):
##        TODO: use segment manager
#        num_atoms = len(AtomSel('(all)').indices())
#        result = [None] * num_atoms
#        return result
#    
#    def make_result_array(self):
#        num_atoms = len(AtomSel('(all)').indices())
#        result = [0.0] * num_atoms
#        return result
#
#    def test_ring_table_names(self):
#        ring_potential = Ring_Potential()
#        
#        expected_table_names = ('ATOM','COEF','RING')
#        table_names = ring_potential.get_component_table_names()
#        self.assertSequenceEqual(expected_table_names,table_names)
#    
#    def test_ring__bb_atoms(self):
#        ring_potential = Ring_Potential()
#        
#        ring_table = Table_manager.get_default_table_manager().get_ring_table('PHE')
#        all_components = ring_potential._get_component_list('COEF').get_all_components()
#        for atom_id,ring_id,coefficient in all_components:
#            self.assertEqual(ring_id,0)
#            target_atoms = ring_table.get_target_atoms()
#            atom_name = target_atoms[atom_id]
#            expected_coefficient = ring_table.get_ring_coefficient(atom_name,'PHE','6')
#            self.assertAlmostEqual(coefficient, expected_coefficient, self.DEFAULT_DECIMAL_PLACES)
#        
#    def test_ring_sidechain_atoms(self):
#        ring_potential = Ring_Potential()
#        
#        ATOM_NAMES_PHE_RING = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
#
##        ring_table = Table_manager.get_default_table_manager().get_ring_table('PHE')
#        all_components = ring_potential._get_component_list('RING').get_all_components()
#        for ring_id,ring_atom_ids in all_components:
#            self.assertEqual(ring_id,0)
#            for ring_atom_id in ring_atom_ids:
#                segment, residue_number,atom_name = Atom_utils._get_atom_info_from_index(ring_atom_id)
#                self.assertEqual(segment,'')
#                self.assertEqual(residue_number,2)
#                self.assertIn(atom_name, ATOM_NAMES_PHE_RING)
#                del ATOM_NAMES_PHE_RING[ATOM_NAMES_PHE_RING.index(atom_name)]
#            self.assertEmpty(ATOM_NAMES_PHE_RING)
#            
#    def test_ring_bb_atoms(self):
#        ring_potential = Ring_Potential()
#        
#        
#
#        ring_table = Table_manager.get_default_table_manager().get_ring_table('PHE')
#        target_atoms = ring_table.get_target_atoms()
#        
#        all_component_names = ring_potential._get_component_list('ATOM').get_all_components()
#        all_component_names_index = dict((i,atom_name) for i,atom_name in enumerate(target_atoms))
#        for atom_id,atom_type_id in all_component_names:
#            segment,residue_number,atom_name = Atom_utils._get_atom_info_from_index(atom_id)
#            self.assertEqual(segment, '')
#            self.assertEqual(residue_number, 2)
#            atom_name  = Atom_utils._get_atom_name_from_index(atom_id)
#            self.assertEqual(all_component_names_index[atom_type_id], atom_name)
#            
#            del all_component_names_index[atom_type_id]
#        self.assertEmpty(all_component_names_index)
#
#
#    def assertListVec3AlmostEqual(self, ring_centres, expected, places = DEFAULT_DECIMAL_PLACES):
#        self.assertEqual(len(ring_centres), 1)
#        for result, expected in zip(ring_centres, expected):
#            for item1, item2 in zip(result, expected):
#                self.assertAlmostEqual(item1, item2, places)
#
#    def test_ring_calculate_centre(self):
#        ring_potential = Ring_Potential()
#        
#        ring_centres = ring_potential._calculate_ring_centres()
#        expected = [Vec3(*AFA.expected_ring_centre),]
#        
#        self.assertListVec3AlmostEqual(ring_centres, expected)
#
#    def test_ring_calculate_normals(self):
#        ring_potential = Ring_Potential()
#        
#        ring_normals = ring_potential._calculate_ring_normals()
#        expected = [Vec3(*AFA.expected_ring_normals),]
#        
#        #TODO: check difference why only 3 places consistent, different atom selections for normals?
#        self.assertListVec3AlmostEqual(ring_normals, expected, self.DEFAULT_DECIMAL_PLACES)
#    
#    def test_calc_component_shift(self):
#        ring_potential = Ring_Potential()
#        
#        for i in range(len(ring_potential._get_component_list('ATOM'))):
#            component_shift = ring_potential._calc_component_shift(i)
#            #TODO: this looks inside the object too much
#            atom_id = ring_potential._get_component_list('ATOM')[i][0]
#            atom_key = Atom_utils._get_atom_info_from_index(atom_id)
#            expected_shift = AFA.expected_ring_shifts[atom_key]
#            #TODO: is the difference in error down to coordinates?
#            self.assertAlmostEqual(component_shift, expected_shift, places=self.DEFAULT_DECIMAL_PLACES)
#            
#    def test_calc_component_forces(self):
#        ring_potential = Ring_Potential()
#        
#        for target_atom_id,atom_type_id in ring_potential._get_component_list('ATOM'):
#            for ring_id,ring_atom_ids in ring_potential._get_component_list('RING'):
#                
#                target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
#                force_factor_key = target_atom_key
#                force_factor = AFA.force_factors_harmonic[force_factor_key]
#                
#                #TODO: these should not be needed
#                ring_potential._get_component_list('RING')
#                ring_potential._build_ring_data_cache()
#                
#                forces = self.make_result_array_forces()
#                force_terms = ring_potential._build_force_terms(target_atom_id, ring_id)
#                ring_potential._calc_target_atom_forces(target_atom_id, ring_id, force_factor, force_terms, forces)
#                
#                target_atom_forces = forces[target_atom_id]
#                expected_forces = AFA.target_forces_harmonic[force_factor_key]
#                self.assertSequenceAlmostEqual(target_atom_forces, expected_forces, self.DEFAULT_DECIMAL_PLACES)
#                
#                forces = self.make_result_array_forces()
#                ring_potential._calculate_ring_forces(atom_type_id, ring_id, force_factor, force_terms, forces)
#                
#                for ring_atom_id in ring_atom_ids:
#                    ring_atom_key  = Atom_utils._get_atom_info_from_index(ring_atom_id)
#                    ring_force_key =target_atom_key,ring_atom_key
#                    expected_ring_forces = AFA.ring_forces_harmonic[ring_force_key]
#                    ring_atom_forces = forces[ring_atom_id]
#                    self.assertSequenceAlmostEqual(ring_atom_forces, expected_ring_forces, self.DEFAULT_DECIMAL_PLACES)
    def test_chemical_shifts(self):
        xcamshift  = Xcamshift()
#        print xcamshift.print_shifts()
        bad_residues =  set()
        for key in gb3.gb3_component_shifts:
            segment, residue_number,atom,sub_potential = key
            sub_potential = xcamshift.get_named_sub_potential(sub_potential)
            if sub_potential.get_abbreviated_name()  ==  BACK_BONE:
                atom_ids  =  Atom_utils.find_atom_ids(segment, residue_number, atom)
                if len(atom_ids) > 0:
    #            find_atom(segment, residue_number, atom)
                    shift  = sub_potential.calc_single_atom_shift(atom_ids[0])
                    expected_shift = gb3.gb3_component_shifts[key]
                    residue_type = Atom_utils._get_residue_type_from_atom_id(atom_ids[0])
                    if abs(shift-expected_shift) > 0.0001:
                        print ` key` +' ' + residue_type, expected_shift, shift
                        bad_residues.add(key[1])
        print bad_residues
#        print xcamshift.print_shifts()
        
if __name__ == "__main__":
    unittest2.main()
#    TestXcamshift.list_test_shifts()
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testCalcForceSetTanh')
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testSingleFactorHarmonic')
