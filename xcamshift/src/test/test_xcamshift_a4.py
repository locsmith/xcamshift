'''
Created on 31 Dec 2011

@author: garyt
'''
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
#from xcamshift import RandomCoilShifts, Distance_potential,  Extra_potential,\
#    Dihedral_potential, Base_potential, Sidechain_potential, Xcamshift,\
from xcamshift import    Non_bonded_potential, Non_bonded_list
from atomSel import AtomSel
#from test.xdists import xdists_ala_3
#from test.dihedrals import dihedrals_ala_3
#from segment_manager import Segment_Manager
#from test.sidechains import sidechains_ala_3
#from test.ala_3 import ala_3_total_shifts
#from test import sidechains, ala_3
from observed_chemical_shifts import Observed_shift_table
from utils import Atom_utils
import sys
from table_manager import Table_manager
from component_list import Component_list
from test import ala_4

TOTAL_ENERGY = 'total'
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
    

def almostEqual(first, second, places = 7):
    result  = False
    if round(abs(second-first), places) == 0:
        result=True
    return result

#class testSegmentManager(object):
class TestXcamshiftA4(unittest2.TestCase):

    def assertLengthIs(self, components_0, length):
        return self.assertEqual(len(components_0), length)


    # TODO add extra places
    DEFAULT_DECIMAL_PLACES = 5
    

    def check_almost_equal(self, list_1, list_2, delta = 1e-7):
        difference_offset = -1
        for i, (elem_1, elem_2) in enumerate(zip(list_1, list_2)):
            diff = abs(elem_1 - elem_2)
            if diff > delta:
                difference_offset = i
        
        return difference_offset
    
    def are_almost_equal_sequences(self, list_1, list_2, delta =  1e-7):
        result = True
        if self.check_almost_equal(list_1, list_2, delta) > 0:
            result = False
        return result
        
    def assertSequenceAlmostEqual(self,result,expected, delta = 1e-7):
        len_result = len(result)
        len_expected = len(expected)
        if len_result != len_expected:
            raise AssertionError("the two lists are of different length %i and %i" % (len_result,len_expected))
        
        difference_offset = self.check_almost_equal(result, expected, delta)
                
        if difference_offset > 0:
            template = "lists differ at item %i: %s - %s > %s"
            elem_1 = result[difference_offset]
            elem_2 = expected[difference_offset]
            message = template % (difference_offset, `elem_1`,`elem_2`,delta)
            raise AssertionError(message)
            
    def setUp(self):
        initStruct("test_data/4_ala/A4.psf")
        PDBTool("test_data/4_ala/A4.pdb").read()


    def make_result_array_forces(self):
#        TODO: use segment manager
        num_atoms = len(AtomSel('(all)').indices())
        result = [None] * num_atoms
        return result
    
    def make_result_array(self):
        num_atoms = len(AtomSel('(all)').indices())
        result = [0.0] * num_atoms
        return result

    
    def assertIsSorted(self, test_elems):
        test_elems_copy = list(test_elems)
        test_elems_copy.sort()
        self.assertEqual(test_elems,test_elems_copy)
    

    @staticmethod
    def assertElemeInSet( elem, xdist_set):
        if not elem in xdist_set:
            raise AssertionError("%s not found in set" % `elem`)

    def dihedral_key_to_atom_ids(self, dihedral_element):
        dihedral_atoms = []
        for vector in dihedral_element[0][1]:
            dihedral_atoms.extend(vector)
        
        return dihedral_atoms


    def remove_zero_valued_keys(self, expected_force_factors):
        for key, value in expected_force_factors.items():
            if almostEqual(value, 0.0, self.DEFAULT_DECIMAL_PLACES):
                del expected_force_factors[key]


    def assertEmpty(self, expected_force_factors,msg = None):
        return self.assertEqual(len(expected_force_factors), 0)

    def get_force_triplet(self, target_atom_index, forces):
        target_atom_forces = forces[target_atom_index]
        return target_atom_forces


    def remove_almost_zero_force_elems(self, expected_forces_dict, delta = 1e-7):
        ZEROS_3 = 0.0, 0.0, 0.0
        for key, value in expected_forces_dict.items():
            if self.check_almost_equal(ZEROS_3, value):
                del expected_forces_dict[key]


    # TODO this is a shim and will disapear once we unify keys formats
    def _convert_dump_dihedral_key(self, dump_selection_data):
        DIHEDRAL_ATOMS_INDEX = 1
        DIHEDRAL_PAIR_1 = 0
        DIHEDRAL_PAIR_2 = 1
        dihedral_angle_key = []
        dihedral_angle_key.extend(dump_selection_data[DIHEDRAL_ATOMS_INDEX][DIHEDRAL_PAIR_1])
        dihedral_angle_key.extend(dump_selection_data[DIHEDRAL_ATOMS_INDEX][DIHEDRAL_PAIR_2])
        dihedral_angle_key = tuple(dihedral_angle_key)
        return dihedral_angle_key
    
    def _extract_dihedral_forces(self, dihedral_atom_ids, forces):
        result = []
        
        for atom_id in dihedral_atom_ids:
            result.append(self.get_force_triplet(atom_id,forces))
            
        return result
                
    


#    def _setup_xcamshift_with_shifts_table(self, test_shifts):
#        xcamshift = Xcamshift()
#        observed_shifts = Observed_shift_table(test_shifts)
#        xcamshift.set_observed_shifts(observed_shifts)
#        return xcamshift
#
#    def _test_force_sets(self, xcamshift, expected_energy, expected_forces):
#        expected_forces = dict(expected_forces)
#        number_atoms = Segment_Manager().get_number_atoms()
#        derivs = [None] * number_atoms
#        energy = xcamshift.calcEnergyAndDerivs(derivs)
#        
#        for atom_id in range(number_atoms):
#            
#            atom_key  =  Atom_utils._get_atom_info_from_index(atom_id)
#            if atom_key in expected_forces:
#                force_triplet = derivs[atom_id]
#                expected_force_triplet = expected_forces[atom_key]
#                
#                self.assertSequenceAlmostEqual(force_triplet, expected_force_triplet, self.DEFAULT_DECIMAL_PLACES)
#                del expected_forces[atom_key]
#                
#        self.remove_almost_zero_force_elems(expected_forces, self.DEFAULT_DECIMAL_PLACES)
#        self.assertEmpty(expected_forces)
#                
#        self.assertAlmostEqual(energy, expected_energy, self.DEFAULT_DECIMAL_PLACES)
        

    def _split_tuple_pair_to_column_sets(self, expected_non_bonded):
        non_bonded_1 = set()
        non_bonded_2 = set()
        for non_bonded_elem_1, non_bonded_elem_2 in expected_non_bonded:
            non_bonded_1.add(non_bonded_elem_1)
            non_bonded_2.add(non_bonded_elem_2)
        
        return non_bonded_1, non_bonded_2

#    def testNonBondedComponents(self):
#        non_bonded_potential = Non_bonded_potential()
#        
#        non_bonded_table = Table_manager().get_non_bonded_table("ALA")
#        target_atom_types = list(non_bonded_table.get_target_atoms())
#        
##        spheres = non_bonded_table.get_spheres()
#        
#        for component in non_bonded_potential._get_all_components('ATOM'):
#            atom_info = Atom_utils._get_atom_info_from_index(component[0])
#            atom_name  =  atom_info[2]
#            
#            self.assertIn(atom_name, target_atom_types)
#            
#            del target_atom_types[target_atom_types.index(atom_name)]
#            
#        self.assertEmpty(target_atom_types)
#        
#        remote_component_spheres_and_exponents = {}
#        for component in non_bonded_potential._get_all_components('NBRM'):
#            atom_info = Atom_utils._get_atom_info_from_index(component[0])
#            
#            remote_component_spheres_and_exponents.setdefault(atom_info,[]).append(component[1:3])
#            
#        
#        for sphere_exponent in remote_component_spheres_and_exponents.values():
#            sphere_exponent.sort()
#            self.assertSequenceEqual(sphere_exponent, [(0,1.0),(1,-3.0)])
#            
#            
#        
#        non_bonded_1, non_bonded_2 = self._split_tuple_pair_to_column_sets(ala_4.ala4_putative_non_bonded_pairs)
#        self.maxDiff = None
#        test_non_bonded_set_1 = set()
#        for component in non_bonded_potential._get_all_components('ATOM'):
#            atom_info = Atom_utils._get_atom_info_from_index(component[0])
#            self.assertElemeInSet(atom_info, non_bonded_1)
#            test_non_bonded_set_1.add(atom_info)
#        self.assertSequenceEqual(test_non_bonded_set_1, non_bonded_1)
#        
#        test_non_bonded_set_2 = set()
#        for component in non_bonded_potential._get_all_components('NBRM'):
#            atom_info = Atom_utils._get_atom_info_from_index(component[0])
#            self.assertElemeInSet(atom_info, non_bonded_2)
#            test_non_bonded_set_2.add(atom_info)
#        self.assertSequenceEqual(test_non_bonded_set_2, non_bonded_2)
#        
#        components_0 = non_bonded_potential._get_components_for_id('NBRM',0)
#        self.assertLengthIs(components_0,2)
#        self.assertEqual(components_0[0][:3], (0,0,1.0))
#        self.assertEqual(components_0[1][:3], (0,1,-3.0))
#        
#        expected_coefficients = (31.32384112711744, -12.5982466223797, -14.86605302072972, -1.6214566324137984, 2.505391454472044,  -88.71419716149445)
#        self.assertSequenceAlmostEqual(components_0[0][3:], expected_coefficients, self.DEFAULT_DECIMAL_PLACES)
#    
#    def testDistances(self):
#        non_bonded_list = Non_bonded_list(min_residue_separation=1)
##        sel = AtomSel("((residue 2 and name CA) around 5.2)")
##        
##        atom_indices = [atom.index() for atom in sel]
##        atom_indices.sort()
#        
#        expected_box_atoms_1 = set(ala_4.ala4_expected_non_bonded_pairs)
#        expected_box_atoms_3 = set(ala_4.ala4_expected_non_bonded_pairs)
#        non_bonded_potential = Non_bonded_potential()
#        target_atoms = non_bonded_potential._get_all_components('ATOM')
#        remote_atoms  = non_bonded_potential._get_all_components('NBRM')
#        
#        component_list =  Component_list()
#        
#        boxes = non_bonded_list.get_boxes(target_atoms, remote_atoms,component_list)
#        for box_component in boxes:
#            target_atom_id, distant_atom_id,coefficient,exponent  = box_component
#            
#            target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
#            box_atom_key = Atom_utils._get_atom_info_from_index(distant_atom_id)
#            atom_key = target_atom_key,box_atom_key
#                
#            if exponent == 1.0:
#                self.assertElemeInSet(atom_key, expected_box_atoms_1)
#                expected_box_atoms_1.remove(atom_key)
#            elif exponent == -3.0:
#                self.assertElemeInSet(atom_key, expected_box_atoms_3)
#                expected_box_atoms_3.remove(atom_key)
#                
#                
#        self.assertEmpty(expected_box_atoms_1)
#        self.assertEmpty(expected_box_atoms_3)    
#        
##        print target_atoms



    def testNonBondedComponents(self):
            non_bonded_potential = Non_bonded_potential()
            non_bonded_potential.update_non_bonded_list()
    
            expected_components = dict(ala_4.ala_components_non_bond)
            for component in non_bonded_potential._get_all_components():
                target_atom_id,remote_atom_id,coefficient,exponent = component
                
                
                target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
                remote_atom_key = Atom_utils._get_atom_info_from_index(remote_atom_id)
            
                expected_shift_key = target_atom_key,remote_atom_key,int(exponent)
                self.assertElemeInSet(expected_shift_key, expected_components)
                self.assertAlmostEqual(expected_components[expected_shift_key], coefficient, self.DEFAULT_DECIMAL_PLACES)
                
                del expected_components[expected_shift_key]
            
            self.assertEmpty(expected_components)
                
                

    def _test_non_bonded_shifts(self, non_bonded_potential, non_bonded_shifts):
        non_bonded_potential.update_non_bonded_list()
        for i, component in enumerate(non_bonded_potential._get_all_components()):
            target_atom_id, remote_atom_id = component[:2]
            exponent = component[3]
            target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
            remote_atom_key = Atom_utils._get_atom_info_from_index(remote_atom_id)
            expected_shift_key = target_atom_key, remote_atom_key, int(exponent)
            if expected_shift_key in non_bonded_shifts:
                expected_shift = non_bonded_shifts[expected_shift_key]
                calculated_shift = non_bonded_potential._calc_component_shift(i)
    #                print expected_shift_key,calculated_shift
                self.assertAlmostEqual(expected_shift, calculated_shift, self.DEFAULT_DECIMAL_PLACES, expected_shift_key)
                del non_bonded_shifts[expected_shift_key]
            else:
                distance = Atom_utils._calculate_distance(target_atom_id, remote_atom_id)
                self.assertTrue(distance >= 5.0, expected_shift_key)
        
        self.assertEmpty(non_bonded_shifts)

    def testNonBondedShiftsNoSmoothing(self):
        non_bonded_potential = Non_bonded_potential(smoothed=False)
        
        non_bonded_shifts = dict(ala_4.ala4_predicted_shifts_non_smoothed)
        
        
        self._test_non_bonded_shifts(non_bonded_potential, non_bonded_shifts)

    def testNonBondedShiftsSmoothed(self):
        non_bonded_potential = Non_bonded_potential()
        
        non_bonded_shifts = dict(ala_4.ala4_predicted_shifts_non_bond)
        
        
        self._test_non_bonded_shifts(non_bonded_potential, non_bonded_shifts)
#        for (target_atom_key,remote_atom_key),forces in non_bonded_forces.items():
#            target_atom_id = Atom_utils.find_atom(*target_atom_key)[0].index()
#            remote_atom_id = Atom_utils.find_atom(*remote_atom_key)[0].index()
#            
#            force_factor =  ala_4.ala4_force_factors_non_bonded[target_atom_key]
#            shift = non_bonded_potential.calc_single_atom_shift(target_atom_id)
#            print 'shift',(target_atom_key,remote_atom_key), target_atom_key,shift

    def _test_non_bonded_force_factors(self, non_bonded_potential, non_bonded_force_factors, active_shifts, factors):
        
        for i, component in enumerate(non_bonded_potential._get_all_components()):
            target_atom_id, remote_atom_id, coefficient, exponent = component
            target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
            remote_atom_key = Atom_utils._get_atom_info_from_index(remote_atom_id)
            
            if active_shifts[target_atom_key] == 0:
                continue
    #
            distance = Atom_utils._calculate_distance(target_atom_id, remote_atom_id)
    #            print target_atom_id,remote_atom_id,distance,target_atom_key, remote_atom_key
            if distance  >= 5.0:
                continue
            
            expected_key = target_atom_key, remote_atom_key, int(exponent)
            factor = factors[target_atom_key]
    #            print expected_key,factor
            force_factor = non_bonded_potential._calc_single_force_factor(i, factor)
            self.assertAlmostEqual(force_factor, non_bonded_force_factors[expected_key])
            del non_bonded_force_factors[expected_key]
        self.assertEmpty(non_bonded_force_factors)

    def testNonBondedForcesNotSmoothed(self):
        non_bonded_potential = Non_bonded_potential()
        non_bonded_potential.set_observed_shifts(ala_4.ala_4_expected_shifts)
        non_bonded_potential.update_non_bonded_list()
        
        non_bonded_force_factors = dict(ala_4.ala_4_force_factors_not_smoothed)
        
        self._test_non_bonded_force_factors(non_bonded_potential, non_bonded_force_factors, 
                                            ala_4.active_shifts, ala_4.ala4_factors_non_bonded)
        
#            forces = self.make_result_array_forces()
#            non_bonded_potential._calc_single_force_set(i, factor, forces)
#                
#            
#            expected_forces = non_bonded_forces[expected_key]
#                print forces,expected_forces
#                    ring_atom_key  = Atom_utils._get_atom_info_from_index(ring_atom_id)
#                    ring_force_key =target_atom_key,ring_atom_key
#                    expected_ring_forces = AFA.ring_forces_harmonic[ring_force_key]
#                    ring_atom_forces = forces[ring_atom_id]
#                    self.assertSequenceAlmostEqual(ring_atom_forces, expected_ring_forces, self.DEFAULT_DECIMAL_PLACES)
if __name__ == "__main__":
    unittest2.main()
#    TestXcamshift.list_test_shifts()
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testNonBondedComponents')
#    unittest2.main(module='test.test_xcamshift_a4',defaultTest='TestXcamshiftA4.testNonBondedForces')
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testSingleFactorHarmonic')
