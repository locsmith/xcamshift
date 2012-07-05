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
from common_constants import BACK_BONE, RANDOM_COIL, XTRA, DIHEDRAL, SIDE_CHAIN,\
    NON_BONDED, RING
from test.gb3 import gb3_component_shifts_sc, gb3_component_shifts_ring
import os
import time
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

def almostEqual(first, second, places = 7):
    result  = False
    if round(abs(second-first), places) == 0:
        result=True
    return result

class TestXcamshiftGB3(unittest2.TestCase):

    # TODO add extra places
    DEFAULT_DECIMAL_PLACES = 5
    DEFAULT_ERROR = 10**-DEFAULT_DECIMAL_PLACES
    
    def assertEmpty(self, expected_force_factors,msg = None):
        return self.assertEqual(len(expected_force_factors), 0)
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
        PDBTool("test_data/gb3/gb3_refined_II.pdb").read()
        Atom_utils.clear_cache()
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

        bad_residues =  set()
        component_shifts_keys = gb3.gb3_subpotential_shifts.keys()
        component_shifts_keys.sort()
        do_update = True
        for i,key in enumerate(component_shifts_keys):
            segment, residue_number,atom,sub_potential = key
            sub_potential = xcamshift.get_named_sub_potential(sub_potential)
            if  do_update and key[3] == NON_BONDED:
                sub_potential.update_non_bonded_list()
                do_update = False

            atom_ids  =  Atom_utils.find_atom_ids(segment, residue_number, atom)
            if len(atom_ids) > 0:
                shift  = sub_potential.calc_single_atom_shift(atom_ids[0])
                expected_shift = gb3.gb3_subpotential_shifts[key]
                residue_type = Atom_utils._get_residue_type_from_atom_id(atom_ids[0])
#               self.assertAlmostEqual(shift, expected_shift, self.DEFAULT_DECIMAL_PLACES-1, msg=`key` + " " + residue_type)

#        print xcamshift.print_shifts()
    
#    def test_non_bonded_components(self):
#        xcamshift =  Xcamshift()
#        sub_potential = xcamshift.get_named_sub_potential(NON_BONDED)
#        sub_potential.update_non_bonded_list()
#        
#        non_bonded_components =  dict(gb3.gb3_component_shifts_non_bonded)
#        
#        exponent_keys = {-3 : 0, 1 : 1}
#        for target_atom_id, remote_atom_id, coefficient, exponent in sub_potential._get_component_list('NBLT'):
#            target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
#            remote_atom_key =  Atom_utils._get_atom_info_from_index(remote_atom_id)
#            
#            if target_atom_key[2]== 'HA1':
#                target_atom_key =  list(target_atom_key)
#                target_atom_key[2]= 'HA'
#                target_atom_key =  tuple(target_atom_key)
#                
#            exponent_key = exponent_keys[int(exponent)]
#            non_bonded_component_key = (target_atom_key, remote_atom_key, exponent_key)
#            
#            
#            if Atom_utils._calculate_distance(target_atom_id, remote_atom_id) > 5.0:
#                continue
#            
#            self.assertIn(non_bonded_component_key, non_bonded_components, `non_bonded_component_key`)
#            del non_bonded_components[non_bonded_component_key]
#            
#        self.assertEmpty(non_bonded_components)
        
    
    def test_non_bonded_components(self):
        xcamshift =  Xcamshift()
        sub_potential = xcamshift.get_named_sub_potential(NON_BONDED)
        sub_potential.update_non_bonded_list()
        
        non_bonded_components =  dict(gb3.gb3_component_shifts_non_bonded)
        
        exponent_keys = {-3 : 0, 1 : 1}
        for target_atom_id, remote_atom_id, coefficient, exponent in sub_potential._get_component_list('NBLT'):
            target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
            remote_atom_key =  Atom_utils._get_atom_info_from_index(remote_atom_id)
            
            if target_atom_key[2]== 'HA1':
                target_atom_key =  list(target_atom_key)
                target_atom_key[2]= 'HA'
                target_atom_key =  tuple(target_atom_key)
                
            exponent_key = exponent_keys[int(exponent)]
            non_bonded_component_key = (target_atom_key, remote_atom_key, exponent_key)
            
            
            if Atom_utils._calculate_distance(target_atom_id, remote_atom_id) > 5.0:
                continue
            
            self.assertIn(non_bonded_component_key, non_bonded_components, `non_bonded_component_key`)
            del non_bonded_components[non_bonded_component_key]
            
        self.assertEmpty(non_bonded_components)
        
    def test_non_bonded_component_shifts(self):
        xcamshift =  Xcamshift()
        sub_potential = xcamshift.get_named_sub_potential(NON_BONDED)
        sub_potential.update_non_bonded_list()
        
        non_bonded_components =  dict(gb3.gb3_component_shifts_non_bonded)
        
        exponent_keys = {-3 : 0, 1 : 1}
        for i, (target_atom_id, remote_atom_id, coefficient, exponent) in enumerate(sub_potential._get_component_list('NBLT')):
            target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
            remote_atom_key =  Atom_utils._get_atom_info_from_index(remote_atom_id)
            
            if target_atom_key[2]== 'HA1':
                target_atom_key =  list(target_atom_key)
                target_atom_key[2]= 'HA'
                target_atom_key =  tuple(target_atom_key)
                
            exponent_key = exponent_keys[int(exponent)]
            non_bonded_component_key = (target_atom_key, remote_atom_key, exponent_key)
            
            
            if Atom_utils._calculate_distance(target_atom_id, remote_atom_id) > 5.0:
                continue
            shift = sub_potential._calc_component_shift(i)
            self.assertAlmostEqual(shift,non_bonded_components[non_bonded_component_key],self.DEFAULT_DECIMAL_PLACES-2,`i` + ' '  + `non_bonded_component_key`)
            del non_bonded_components[non_bonded_component_key]
            
        self.assertEmpty(non_bonded_components)
    def remove_zero_valued_keys(self, expected_force_factors):
        for key, value in expected_force_factors.items():
            if almostEqual(value, 0.0, self.DEFAULT_DECIMAL_PLACES):
                del expected_force_factors[key]

#    def test_component_shifts_sidechain(self):
#        return
#        
#        xcamshift = Xcamshift()
#        sidechain_subpotential = xcamshift.get_named_sub_potential(SIDE_CHAIN)
#        
#        expected_sidechain_shifts = dict(gb3_component_shifts_sc)
#        expected_component_keys = expected_sidechain_shifts.keys()
#        for component_index, component in enumerate(sidechain_subpotential._get_component_list()):
#            from_atom_id, to_atom_id = component[0:2]
#            from_atom_key = Atom_utils._get_atom_info_from_index(from_atom_id)
#            to_atom_key = Atom_utils._get_atom_info_from_index(to_atom_id)
#            
##            print from_atom_key, to_atom_key
#            if from_atom_key[2] == 'HA1':
#                from_atom_key = from_atom_key[0],from_atom_key[1],'HA'
#                
#            expected_key = from_atom_key, to_atom_key
#            
#            self.assertIn(expected_key, expected_component_keys, `expected_key` + " exists")
#            
#            shift = sidechain_subpotential._calc_component_shift(component_index)
#            
#            self.assertAlmostEqual(expected_sidechain_shifts[expected_key], shift, places=self.DEFAULT_DECIMAL_PLACES - 1, msg=`expected_key`)
#            
#            del expected_sidechain_shifts[expected_key]
#            
#        self.remove_zero_valued_keys(expected_sidechain_shifts)
#        self.assertEmpty(expected_sidechain_shifts)
#        

#    def test_component_shifts_ring(self):
#        
#        xcamshift = Xcamshift()
#        ring_subpotential = xcamshift.get_named_sub_potential(RING)
#        
#        expected_ring_shifts = dict(gb3_component_shifts_ring)
#        expected_component_keys = expected_ring_shifts.keys()
#        for component_index, component in enumerate(ring_subpotential._get_component_list()):
#            from_atom_id, atom_type_id = component
#            from_atom_info_list = ring_subpotential._get_component_list('COEF').get_components_for_atom_id(atom_type_id)
#            
##            print Atom_utils._get_atom_info_from_index(from_atom_id), from_atom_info_list
#            from_atom_key = list(Atom_utils._get_atom_info_from_index(from_atom_id))
#            
#            if from_atom_key[2] == 'HA1':
#                from_atom_key[2] =  'HA'
#            from_atom_key = tuple(from_atom_key)
#            
#            for sub_component_index, from_atom_info in enumerate(from_atom_info_list):
#                from_atom_id, ring_id, coefficent = from_atom_info
##                print 'debug ', (from_atom_id, atom_type_id), from_atom_info
#                ring_info = ring_subpotential._get_component_list('RING').get_components_for_atom_id(ring_id)
##                print ring_info[0][1]
#                
#                ring_atoms =  ring_info[0][1]
#                ring_residue_type = Atom_utils._get_residue_type_from_atom_id(ring_atoms[1])
#                ring_residue = Atom_utils._get_atom_info_from_index(ring_atoms[0])[1]
#                expected_key =   from_atom_key,(ring_residue,ring_residue_type,len(ring_atoms))
##                print expected_key
#                self.assertIn(expected_key, expected_component_keys, `expected_key`)
#            
#
#
#                shift = ring_subpotential._calc_sub_component_shift(component_index, sub_component_index)
#                self.assertAlmostEqual(expected_ring_shifts[expected_key], shift, places=self.DEFAULT_DECIMAL_PLACES - 2, msg=`expected_key`)
#                if abs(expected_ring_shifts[expected_key] - shift) > 0.001:
#                    print 'fail', expected_key, expected_ring_shifts[expected_key], shift, Atom_utils._get_residue_type_from_atom_id(from_atom_id)
#                    print
#                
#                del expected_ring_shifts[expected_key]
#                print
#
#            
#        self.remove_zero_valued_keys(expected_ring_shifts)
#        print expected_ring_shifts
#        self.assertEmpty(expected_ring_shifts)
import cProfile



if __name__ == "__main__":
    unittest2.main()
#    cProfile.run('unittest2.main()')
    
#    TestXcamshift.list_test_shifts()
#    unittest2.main(module='test.test_xcamshift_gb3',defaultTest='TestXcamshiftGB3.test_non_bonded_component_shifts')
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testSingleFactorHarmonic')
