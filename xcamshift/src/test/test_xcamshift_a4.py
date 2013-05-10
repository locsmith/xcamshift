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
#from xcamshift import RandomCoilShifts, Distance_potential,  Extra_potential,\
#    Dihedral_potential, Base_potential, Sidechain_potential, Xcamshift,\
from xcamshift import    Non_bonded_potential, Non_bonded_list, Xcamshift
from atomSel import AtomSel
from cython.shift_calculators import allocate_array, zero_array
#from test.xdists import xdists_ala_3
#from test.dihedrals import dihedrals_ala_3
#from cython.fast_segment_manager import Segment_Manager
#from test.sidechains import sidechains_ala_3
#from test.ala_3 import ala_3_total_shifts
#from test import sidechains, ala_3
from observed_chemical_shifts import Observed_shift_table
from utils import Atom_utils
import sys
from table_manager import Table_manager
from component_list import Component_list
from test import ala_4
from cython.fast_segment_manager import Segment_Manager
from cython.shift_calculators import Out_array



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

    def _setup_xcamshift_with_shifts_table(self, test_shifts):
        xcamshift = self._get_xcamshift()
        observed_shifts = Observed_shift_table(test_shifts)
        xcamshift.set_observed_shifts(observed_shifts)
        return xcamshift
    
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
        Segment_Manager.reset_segment_manager()
        Atom_utils.clear_cache()

    def make_out_array(self):
#        TODO: use segment manager
        num_atoms = len(AtomSel('(all)').indices())
        result = Out_array(num_atoms)
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

#TODO: replace use of delta with decimals or sgignificant figures (better)
    def remove_almost_zero_force_elems(self, expected_forces_dict, delta = 1e-7):
        ZEROS_3 = 0.0, 0.0, 0.0
        for key, value in expected_forces_dict.items():
            if self.check_almost_equal(ZEROS_3, value):
                del expected_forces_dict[key]

    def remove_almost_zero_force_elems_from_list(self, forces, delta = 1e-7):
        ZEROS_3 = 0.0, 0.0, 0.0
        
        to_delete = []
        for i,value in enumerate(forces):
            if value == None or not self.check_almost_equal(ZEROS_3, value):
                to_delete.append(i)
        to_delete.reverse()
        
        for i in to_delete:
            del forces[i]
        return forces

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




    def _get_non_bonded_potential(self, smoothed = True):
        non_bonded_potential = Non_bonded_potential(smoothed)
        return non_bonded_potential

    def testNonBondedComponents(self):
            non_bonded_potential = self._get_non_bonded_potential()
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
        data =  non_bonded_potential._get_all_components()
        for i, components in enumerate(zip(data[0::2], data[1::2])):
             
            
            expected_shift = 0.0
            key_good = False
            for j, component in enumerate(components):
                target_atom_id, remote_atom_id = components[0][:2]
                exponent = component[3]
                target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
                remote_atom_key = Atom_utils._get_atom_info_from_index(remote_atom_id)
                expected_shift_key = target_atom_key, remote_atom_key, int(exponent)
                if expected_shift_key in non_bonded_shifts:
                    expected_shift += non_bonded_shifts[expected_shift_key]
                    key_good=True
                elif j == 0:
                    distance = Atom_utils._calculate_distance(target_atom_id, remote_atom_id)
                    self.assertTrue(distance >= 5.0, expected_shift_key)
                
            if key_good:
                calculated_shift = non_bonded_potential._calc_component_shift(i)
                self.assertAlmostEqual(expected_shift, calculated_shift, self.DEFAULT_DECIMAL_PLACES, expected_shift_key)
                


    def testNonBondedShiftsNoSmoothing(self):
        non_bonded_potential = self._get_non_bonded_potential(smoothed=False)
        
        non_bonded_shifts = dict(ala_4.ala4_predicted_shifts_non_smoothed)
        
        
        self._test_non_bonded_shifts(non_bonded_potential, non_bonded_shifts)

    def testNonBondedShiftsSmoothed(self):
        non_bonded_potential = self._get_non_bonded_potential()
        
        non_bonded_shifts = dict(ala_4.ala4_predicted_shifts_non_bond)
        
        
        self._test_non_bonded_shifts(non_bonded_potential, non_bonded_shifts)
#        for (target_atom_key,remote_atom_key),forces in non_bonded_forces.items():
#            target_atom_id = Atom_utils.find_atom(*target_atom_key)[0].index()
#            remote_atom_id = Atom_utils.find_atom(*remote_atom_key)[0].index()
#            
#            force_factor =  ala_4.ala4_force_factors_non_bonded[target_atom_key]
#            shift = non_bonded_potential._calc_single_atom_shift(target_atom_id)
#            print 'shift',(target_atom_key,remote_atom_key), target_atom_key,shift

    def _test_non_bonded_force_factors(self, non_bonded_potential, non_bonded_force_factors, active_shifts, factors):
        
        components = dict(non_bonded_potential._get_components())
        nb_interaction_list = components['NBLT']
        active_array = allocate_array(1,'i')
        components['ACTI'] =  active_array
        non_bonded_potential._force_calculator._set_components(components)
        
        for i, component in enumerate(nb_interaction_list):
            zero_array(active_array)
            active_array[0] = i
            
            target_atom_id =  non_bonded_potential._get_component_list('ATOM')[component[1]][0]
            remote_atom_id =  non_bonded_potential._get_component_list('NBRM')[component[2]][0]

            target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
            remote_atom_key = Atom_utils._get_atom_info_from_index(remote_atom_id)
            
            if active_shifts[target_atom_key] == 0:
                continue
    #
            distance = Atom_utils._calculate_distance(target_atom_id, remote_atom_id)

            if distance  < 5.0:
               
                for exponent_index, exponent in enumerate((1,-3)):
                    expected_key = target_atom_key, remote_atom_key, int(exponent)
                    factor = factors[target_atom_key]
    
                    non_bonded_potential._force_calculator._build_component(0,exponent_index)
    
                    force_factor = non_bonded_potential._force_calculator._calc_single_force_factor(0, factor)
                    #TODO: check change from 7 to 5 dp is ok
                    self.assertAlmostEqual(force_factor, non_bonded_force_factors[expected_key],places=self.DEFAULT_DECIMAL_PLACES)
                    del non_bonded_force_factors[expected_key]
        self.assertEmpty(non_bonded_force_factors)

    def testNonBondedForceFactorsNotSmoothed(self):
        non_bonded_potential = self._get_non_bonded_potential(smoothed=False)
        non_bonded_potential.set_observed_shifts(ala_4.ala_4_expected_shifts)
        non_bonded_potential.update_non_bonded_list()
        
        non_bonded_force_factors = dict(ala_4.ala_4_force_factors_not_smoothed)
        
        self._test_non_bonded_force_factors(non_bonded_potential, non_bonded_force_factors, 
                                            ala_4.active_shifts, ala_4.ala4_factors_non_bonded)
        
    def testNonBondedForceFactorsSmoothed(self):
        non_bonded_potential = self._get_non_bonded_potential()
        non_bonded_potential.set_observed_shifts(ala_4.ala_4_expected_shifts)
        non_bonded_potential.update_non_bonded_list()
        
        non_bonded_force_factors = dict(ala_4.ala_4_force_factors_smoothed)
        
        self._test_non_bonded_force_factors(non_bonded_potential, non_bonded_force_factors, 
                                            ala_4.active_shifts, ala_4.ala4_factors_non_bonded)

    def _test_non_bonded_forces(self, non_bonded_potential, non_bonded_forces, active_shifts, factors):

        components = dict(non_bonded_potential._get_components())
        nb_interaction_list = components['NBLT']
        active_array = allocate_array(1,'i')
        components['ACTI'] =  active_array
        non_bonded_potential._force_calculator._set_components(components)
        
        for i, component in enumerate(nb_interaction_list):
            zero_array(active_array)
            active_array[0] = i
            
            target_atom_id =  non_bonded_potential._get_component_list('ATOM')[component[1]][0]
            remote_atom_id =  non_bonded_potential._get_component_list('NBRM')[component[2]][0]

            target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
            remote_atom_key = Atom_utils._get_atom_info_from_index(remote_atom_id)
             
            if active_shifts[target_atom_key] == 0:
                continue
    #
            distance = Atom_utils._calculate_distance(target_atom_id, remote_atom_id)

            if distance  < 5.0:
               
                for exponent_index, exponent in enumerate((1,-3)):
                    expected_key = target_atom_key, remote_atom_key, int(exponent)
                    factor = factors[target_atom_key]
    
                    non_bonded_potential._force_calculator._build_component(0,exponent_index)
                    
                    factor = factors[target_atom_key]

                    out_array = self.make_out_array()
                    non_bonded_potential._force_calculator._calc_single_force_set(0, factor,out_array)
                     
                    forces = out_array.add_forces_to_result()
                    target_force_triplet = forces[target_atom_id]
                    remote_force_triplet  = forces[remote_atom_id]
                     
                    expected_target_forces = non_bonded_forces[expected_key]
                    expected_remote_forces  = [-elem for elem in expected_target_forces]
                     
                    #TODO: check change from 7 to 5 dp is ok also improve assertSequenceAlmostEqual ro take a places argument
                    self.assertSequenceAlmostEqual(target_force_triplet, expected_target_forces,delta= 1e-1**self.DEFAULT_DECIMAL_PLACES)
                    self.assertSequenceAlmostEqual(remote_force_triplet, expected_remote_forces,delta= 1e-1**self.DEFAULT_DECIMAL_PLACES)

                    atom_ids =  [target_atom_id,remote_atom_id]
                    atom_ids.sort(reverse=True)
                    for atom_id in atom_ids:
                        del forces[atom_id]
                    forces = self.remove_almost_zero_force_elems_from_list(forces)
                     
                    self.assertEmpty(forces)
                    
                    del non_bonded_forces[expected_key]
        self.assertEmpty(non_bonded_forces)

    def testNonBondedForces(self):
        non_bonded_potential = self._get_non_bonded_potential(smoothed=True)
        non_bonded_potential.set_observed_shifts(ala_4.ala_4_expected_shifts)
        non_bonded_potential.update_non_bonded_list()
        
        non_bonded_forces = dict(ala_4.ala4_non_bonded_forces)
        
        self._test_non_bonded_forces(non_bonded_potential, non_bonded_forces, 
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


    def _get_xcamshift(self):
        xcamshift_potential = Xcamshift()
        return xcamshift_potential

    def  test_overall_shifts_a4(self):
        xcamshift_potential = self._get_xcamshift()
        
        shifts = self.make_result_array()
        shifts = xcamshift_potential.set_shifts(shifts)
        
        expected  = [0.0] * len(shifts)
        
        for target_atom_key in ala_4.ala_4_expected_shifts:
            target_atom_id  = Atom_utils.find_atom_ids(*target_atom_key)[0]
            expected[target_atom_id] = ala_4.ala_4_expected_shifts[target_atom_key]
            
        self.assertSequenceAlmostEqual(expected, shifts,self.DEFAULT_DECIMAL_PLACES)
        

    def testComponentShiftsA4(self):
        xcamshift_potential =  self._get_xcamshift()
        
        expected_shift_components = dict(ala_4.ala_4_expected_shift_components)
        for sub_potential_name in xcamshift_potential.get_sub_potential_names():
            sub_potential = xcamshift_potential.get_named_sub_potential(sub_potential_name)
            
            for target_atom_id in sub_potential.get_target_atom_ids():
                expected_key = Atom_utils._get_atom_info_from_index(target_atom_id),sub_potential_name
                single_atom_shift = sub_potential._calc_single_atom_shift(target_atom_id)

                #TODO: remove limitation to 4decimal points due to poor float reading in almost
                self.assertAlmostEqual(expected_shift_components[expected_key],single_atom_shift,self.DEFAULT_DECIMAL_PLACES-1, `expected_key`)
                
                del expected_shift_components[expected_key]
        self.remove_zero_valued_keys(expected_shift_components)
        self.assertEmpty(expected_shift_components)
        
    def _test_force_sets(self, xcamshift, expected_energy, expected_forces):
        expected_forces = dict(expected_forces)
        number_atoms = Segment_Manager().get_number_atoms()
        derivs = [None] * number_atoms
        energy = xcamshift.calcEnergyAndDerivs(derivs)
        
        for atom_id in range(number_atoms):
            
            atom_key  =  Atom_utils._get_atom_info_from_index(atom_id)
            if atom_key in expected_forces:
                force_triplet = derivs[atom_id]
                if force_triplet == None:
                    force_triplet = (0.0,0.0,0.0)
                expected_force_triplet = expected_forces[atom_key]
                
                self.assertSequenceAlmostEqual(force_triplet, expected_force_triplet, self.DEFAULT_DECIMAL_PLACES)
                del expected_forces[atom_key]
                
        self.remove_almost_zero_force_elems(expected_forces, self.DEFAULT_DECIMAL_PLACES)
        self.assertEmpty(expected_forces)
        
        #TODO: improve accuracy of test energy
        self.assertAlmostEqual(energy, expected_energy, self.DEFAULT_DECIMAL_PLACES-1)
        
    
    def testCalcForceSetHarmonic(self):
        xcamshift = self._setup_xcamshift_with_shifts_table(ala_4.ala_4_test_shifts_harmonic)
        expected_energy = ala_4.ala_4_energies_harmonic[TOTAL_ENERGY]
        expected_forces =ala_4.ala_4_total_forces_harmonic
        
        self._test_force_sets(xcamshift, expected_energy, expected_forces)
                
def run_tests():
    unittest2.main(module='test.test_xcamshift_a4')
#     unittest2.main(module='test.test_xcamshift_a4',defaultTest='TestXcamshiftA4.testNonBondedForceFactorsSmoothed')
    
if __name__ == "__main__":
    run_tests()
#    TestXcamshift.list_test_shifts()
#    unittest2.main(module='test.test_xcamshift_a4',defaultTest='TestXcamshiftA4.testNonBondedForces')
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testSingleFactorHarmonic')
