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
from simulation import currentSimulation
from derivList import DerivList
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from xcamshift import RandomCoilShifts, Distance_potential,  Extra_potential,\
    Dihedral_potential, Base_potential, Sidechain_potential, Xcamshift,\
    Non_bonded_potential, Non_bonded_list
from atomSel import AtomSel
from test.xdists import xdists_ala_3
from test.dihedrals import dihedrals_ala_3
from cython.fast_segment_manager import Segment_Manager
from test.sidechains import sidechains_ala_3
from test.ala_3 import ala_3_total_shifts
from test import sidechains, ala_3
from observed_chemical_shifts import Observed_shift_table
from utils import Atom_utils
import sys
from table_manager import Table_manager
from component_list import Component_list
from common_constants import TARGET_ATOM_IDS_CHANGED
from cython.shift_calculators import Out_array
from ensembleSimulation import   EnsembleSimulation

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
class TestXcamshift(unittest2.TestCase):
    
    def __init__(self,*args,**kwargs):
        super(TestXcamshift, self).__init__(*args,**kwargs)
        self._esim = None
        
    def _make_xcamshift(self, shifts={}):
        xcamshift = self._setup_xcamshift_with_shifts_table(shifts)
        return xcamshift
    
    def get_single_member_ensemble_simulation(self):
        if self._esim.__class__ ==  None.__class__:
            #TODO note EnsembleSimulation can't have a single member that causes a crash!
            # therefore a hack
            self._esim =  Xcamshift().ensembleSimulation()
        return self._esim
     
    def _prepare_xcamshift(self, xcamshift):
#            
        xcamshift._prepare(TARGET_ATOM_IDS_CHANGED, xcamshift._get_active_target_atom_ids())
        xcamshift._calc_shift_cache(xcamshift._get_active_target_atom_ids())
        xcamshift.update_force_factor_calculator()
        
       


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
            
    def setUp(self):
        initStruct("test_data/3_ala/3ala.psf")
        PDBTool("test_data/3_ala/3ala.pdb").read()
        Segment_Manager.reset_segment_manager()
        Atom_utils.clear_cache()

    def make_out_array(self):
#        TODO: use segment manager
        num_atoms = len(AtomSel('(all)').indices())
        result = Out_array( num_atoms)
        return result
    
    def make_result_array(self):
        num_atoms = len(AtomSel('(all)').indices())
        result = [0.0] * num_atoms
        return result

    def testRandomCoilPotential(self):
        random_coil_potential = RandomCoilShifts(self.get_single_member_ensemble_simulation())
        result = self.make_result_array() 
        expected  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                     123.6, 8.2400000000000002, 52.266440000000003, 4.4328469999999998, 
                     19.0, 0.0, 0.0, 0.0, 177.09999999999999, 0.0, 0.0, 0.0, 0.0, 0.0, 
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        random_coil_potential.set_shifts(result)
        self.assertSequenceAlmostEqual(expected, result)
        
    def testDistancePotential(self):
        distance_potential =  Distance_potential(self.get_single_member_ensemble_simulation())
        result=self.make_result_array()
        distance_potential.set_shifts(result)
        
        expected = self.make_result_array()
        expected[15] =   0.82775
        expected[14] =  41.5303
        expected[13] =  -5.56269
        expected[12] = -12.185
        expected[20] = -12.584
        expected[16] =  -3.6896
        

        self.assertSequenceAlmostEqual(result, expected,delta=0.001)
        
    
    def assertIsSorted(self, test_elems):
        test_elems_copy = list(test_elems)
        test_elems_copy.sort()
        self.assertEqual(test_elems,test_elems_copy)
    
    
    def testIndexComponents(self):
        distance_potential =  Distance_potential(self.get_single_member_ensemble_simulation())
        test_elems = []
        for elem in distance_potential._get_all_components():
            test_elems.append(elem[0])
        self.assertIsSorted(test_elems)
        
    @staticmethod
    def assertElemeInSet( elem, xdist_set):
        if not elem in xdist_set:
            raise AssertionError("%s not found in set" % `elem`)


    def testExtraPotentialComponentsCorrect(self):
        extra_potential = Extra_potential(self.get_single_member_ensemble_simulation())
        
        xdists_ala_3_copy = dict(xdists_ala_3)
        for extra_elem in extra_potential.dump():
            elem_key = extra_elem[:-1]
            self.assertElemeInSet(elem_key, xdists_ala_3)
            del xdists_ala_3_copy[elem_key]
        
        self.assertEqual(0, len(xdists_ala_3_copy))


    def testExtraPotentialCoefficientsCorrect(self):
        extra_potential = Extra_potential(self.get_single_member_ensemble_simulation())
        for extra_elem in extra_potential.dump():
            elem_key = extra_elem[:-1]
            coefficient = xdists_ala_3[elem_key][0]

            self.assertAlmostEqual(coefficient, extra_elem[-1], places=self.DEFAULT_DECIMAL_PLACES)
        

    def testExtraPotentialComponentShiftsCorrect(self):
        extra_potential = Extra_potential(self.get_single_member_ensemble_simulation())
        
        result=self.make_result_array()
        extra_potential.set_shifts(result)
        
        for i,extra_elem in enumerate(extra_potential.dump()):
            
            elem_key = extra_elem[:-1]
            shift = extra_potential._calc_component_shift(i)
            
            expected_shift = xdists_ala_3[elem_key][2]
            
            self.assertAlmostEqual(shift, expected_shift,places=4)

    def testDihdedralPotentialComponentsCorrect(self):
        dihedral_portential = Dihedral_potential(self.get_single_member_ensemble_simulation())
        
        dihedrals_ala_3_copy = dict(dihedrals_ala_3)
        for dihedral_elem in dihedral_portential.dump():
            self.assertElemeInSet(dihedral_elem[0], dihedrals_ala_3_copy)
            del dihedrals_ala_3_copy[dihedral_elem[0]]

        self.assertEqual(0, len(dihedrals_ala_3_copy))
        
    def testDihedralPotentialCoefficientsCorrect(self):
        dihedral_potential = Dihedral_potential(self.get_single_member_ensemble_simulation())
        for dihedral_elem_key, coeff, param_0,param_1,param_2,param_3,param_4, exponent in  dihedral_potential.dump():
            expected_values = [dihedrals_ala_3[dihedral_elem_key][0]]
            expected_values.extend(dihedrals_ala_3[dihedral_elem_key][1])

            actual_values =(coeff, param_0,param_1,param_2,param_3,param_4)
            
            self.assertEqual(len(expected_values), len(actual_values))
            for expected,actual in zip(expected_values,actual_values):
                self.assertAlmostEqual(expected, actual, self.DEFAULT_DECIMAL_PLACES)

            
            self.assertAlmostEqual(exponent, 1.0)


    def dihedral_key_to_atom_ids(self, dihedral_element):
        dihedral_atoms = []
        for vector in dihedral_element[0][1]:
            dihedral_atoms.extend(vector)
        
        return dihedral_atoms

    def testDihedralPotentialAngleCorrect(self):
        dihedral_potential = Dihedral_potential(self.get_single_member_ensemble_simulation())
        
        for dihedral_element in dihedral_potential.dump():
            
            dihedral_atoms = self.dihedral_key_to_atom_ids(dihedral_element)
            atom_ids  = text_keys_to_atom_ids(dihedral_atoms)
            
            expected  =  dihedrals_ala_3[dihedral_element[0]][2]
            angle =  dihedral_potential._get_dihedral_angle(*atom_ids)
            
            self.assertAlmostEqual(expected, angle,self.DEFAULT_DECIMAL_PLACES)

            
    def testDihedralPotentialComponentShiftsCorrect(self):
        dihedral_potential = Dihedral_potential(self.get_single_member_ensemble_simulation())
        
        for i,dihedral_element in enumerate(dihedral_potential.dump()):
            
            expected  =  dihedrals_ala_3[dihedral_element[0]][3]
            angle =  dihedral_potential._calc_component_shift(i)
            
            self.assertAlmostEqual(expected, angle,self.DEFAULT_DECIMAL_PLACES)
#
####
    def testSidechainPotentialComponentsCorrect(self):
        sidechain_potential = Sidechain_potential(self.get_single_member_ensemble_simulation())
        
        sidechains_ala_3_copy = dict(sidechains_ala_3)
        for sidechain_elem  in sidechain_potential.dump():
            expected_key = tuple(sidechain_elem[:2])
            self.assertElemeInSet(expected_key, sidechains_ala_3_copy)
            del sidechains_ala_3_copy[expected_key]

        self.assertEqual(0, len(sidechains_ala_3_copy))
        
    def testSidechainPotentialCoefficientsCorrect(self):
        sidechain_potential = Sidechain_potential(self.get_single_member_ensemble_simulation())
        for target_atom,distant_atom, coeff, exponent in  sidechain_potential.dump():
            sidechain_potential_key=target_atom,distant_atom
            expected_values = sidechains_ala_3[sidechain_potential_key][:2]
            
            actual_values  = (coeff,exponent)
            
            self.assertEqual(len(expected_values), len(actual_values))
            for expected,actual in zip(expected_values,actual_values):
                self.assertAlmostEqual(expected, actual, self.DEFAULT_DECIMAL_PLACES)

            
            self.assertAlmostEqual(exponent, 1.0)
         
    def testSidechainlPotentialComponentShiftsCorrect(self):
        sidechain_potential = Sidechain_potential(self.get_single_member_ensemble_simulation())
        
        for i,elem in enumerate(sidechain_potential.dump()):

            key = tuple(elem[:2])
            
            shift =  sidechain_potential._calc_component_shift(i)
            expected_shift  = sidechains.sidechains_ala_3[key][2]
            self.assertAlmostEqual(expected_shift, shift, self.DEFAULT_DECIMAL_PLACES)

    def testXcamshift(self):
        xcamshift_potential =  self._make_xcamshift()
        
        shifts = self.make_result_array()
        shifts = xcamshift_potential.set_shifts(shifts)
        
        expected  = [0.0] * len(shifts)
        expected[12] = ala_3_total_shifts[2]['N']
        expected[13] = ala_3_total_shifts[2]['HN']
        expected[14] = ala_3_total_shifts[2]['CA']
        expected[15] = ala_3_total_shifts[2]['HA']
        expected[16] = ala_3_total_shifts[2]['CB']
        expected[20] = ala_3_total_shifts[2]['C']
        
        self.assertSequenceAlmostEqual(expected, shifts, delta=0.0001)
        
#        xcamshift_potential.print_shifts()
        

    def _test_single_energies_ala_3(self, test_shifts, expected_energys):
        xcamshift = self._make_xcamshift()
        shift_table = Observed_shift_table(test_shifts)
        xcamshift.set_observed_shifts(shift_table)

        for atom_index in shift_table.get_atom_indices():
            key = Atom_utils._get_atom_info_from_index(atom_index)[1:]
            expected_energy = expected_energys[key]
            energy = xcamshift._calc_single_atom_energy(atom_index)
            self.assertAlmostEqual(energy, expected_energy, self.DEFAULT_DECIMAL_PLACES-1)
        

    def testSingleEnergiesHarmonic(self):
        test_shifts = ala_3.ala_3_test_shifts_harmonic
        expected_energys = ala_3.ala_3_energies_harmonic
        
        self._test_single_energies_ala_3(test_shifts, expected_energys)
        
    def testSingleEnergiesWell(self):
        test_shifts = ala_3.ala_3_test_shifts_well
        expected_energys = ala_3.ala_3_energies_well
        
        self._test_single_energies_ala_3(test_shifts, expected_energys)
        
    def testSingleEnergiesTanh(self):
        test_shifts = ala_3.ala_3_test_shifts_tanh
        expected_energys = ala_3.ala_3_energies_tanh
        
        self._test_single_energies_ala_3(test_shifts, expected_energys)

    def _test_single_factor_set(self, test_shifts, expected_factors):
        xcamshift = self._make_xcamshift(test_shifts)
        self._prepare_xcamshift(xcamshift)
        shift_table = Observed_shift_table(test_shifts)
        for atom_index in shift_table.get_atom_indices():
            key = Atom_utils._get_atom_info_from_index(atom_index)[1:]
            factor = xcamshift._calc_single_factor(atom_index)
            expected_factor = expected_factors[key]
            self.assertAlmostEqual(factor, expected_factor, self.DEFAULT_DECIMAL_PLACES-1)

    def testSingleFactorHarmonic(self):
        test_shifts = ala_3.ala_3_test_shifts_harmonic
        expected_factors = ala_3.ala_3_factors_harmonic
        
        self._test_single_factor_set(test_shifts, expected_factors)

    def testSingleFactorWell(self):
        test_shifts = ala_3.ala_3_test_shifts_well
        expected_factors = ala_3.ala_3_factors_well
        
        self._test_single_factor_set(test_shifts, expected_factors)

    def testSingleFactorTanh(self):
        test_shifts = ala_3.ala_3_test_shifts_tanh
        expected_factors = ala_3.ala_3_factors_tanh
        
        self._test_single_factor_set(test_shifts, expected_factors)


    def remove_zero_valued_keys(self, expected_force_factors):
        for key, value in expected_force_factors.items():
            if almostEqual(value, 0.0, self.DEFAULT_DECIMAL_PLACES):
                del expected_force_factors[key]


    def assertEmpty(self, expected_force_factors,msg = None):
        return self.assertEqual(len(expected_force_factors), 0)

    def testDistanceBasedPotentialSingleForceFactor(self):
        # TODO make this a dummy distance based potential, 
        # _calc_single_force_factor is in distance based potential
        #TODO: test the force calculator not the potential
        distance_potential = Distance_potential(self.get_single_member_ensemble_simulation())
        distance_potential._force_calculator._set_components(distance_potential._get_component_list().get_native_components())
            
        expected_force_factors = dict(ala_3.ala_3_distance_forces_well)
        test_force_factors = ala_3.ala_3_factors_harmonic
        
        for i,data in enumerate(distance_potential.dump()):
            target_atom_key = data[0][1:]
            distant_atom_key = data[1][1:]
            expected_key = (target_atom_key,distant_atom_key)
            test_factor  =  test_force_factors[target_atom_key]
            
            force = distance_potential._force_calculator._calc_single_force_factor(i, test_factor)
            expected_force_factor = expected_force_factors[expected_key]
            del expected_force_factors[expected_key]

            self.assertAlmostEqual(force, expected_force_factor, self.DEFAULT_DECIMAL_PLACES-2)
        
        self.remove_zero_valued_keys(expected_force_factors)
        self.assertEmpty(expected_force_factors)


    def get_force_triplet(self, target_atom_index, forces):
        target_atom_forces = forces[target_atom_index]
        return target_atom_forces


    def remove_almost_zero_force_elems(self, expected_forces_dict, delta = 1e-7):
        ZEROS_3 = 0.0, 0.0, 0.0
        for key, value in expected_forces_dict.items():
            if self.check_almost_equal(ZEROS_3, value):
                del expected_forces_dict[key]


    def _test_distance_forces(self, test_factors, distance_potential,expected_forces):
        expected_forces_dict = dict(expected_forces)
        test_factors = ala_3.ala_3_factors_harmonic
        
        distance_potential._force_calculator._set_components(distance_potential._get_component_list().get_native_components())
        indices = distance_potential._get_indices()
        
        #TODO don't like using dump here
        for i, data in enumerate(distance_potential.dump()):
            target_atom_key = data[indices.target_atom_index][1:]
            distant_atom_key_1 = data[indices.distance_atom_index_1][1:]
           
            distant_atom_key_2 = data[indices.distance_atom_index_2][1:]
            distant_atom_index_1 = single_text_key_to_atom_ids(distant_atom_key_1)[0]
            distant_atom_index_2 = single_text_key_to_atom_ids(distant_atom_key_2)[0]
            
            expected_key = target_atom_key,distant_atom_key_1,distant_atom_key_2
            expected_forces = expected_forces_dict[expected_key]
            negative_expected_forces = [elem * -1 for elem in expected_forces]
            
            test_factor = test_factors[target_atom_key]
            
            out_array = self.make_out_array()
            
            #TODO: test the force calculator not the potential
            distance_potential._force_calculator._calc_single_force_set(i, test_factor, out_array)
            forces = out_array.add_forces_to_result()
            
            distant_atom_forces_1 = self.get_force_triplet(distant_atom_index_1, forces)
            distant_atom_forces_2 = self.get_force_triplet(distant_atom_index_2, forces)
            
            self.assertSequenceAlmostEqual(distant_atom_forces_1, expected_forces, self.DEFAULT_DECIMAL_PLACES)
            self.assertSequenceAlmostEqual(distant_atom_forces_2, negative_expected_forces, self.DEFAULT_DECIMAL_PLACES)

            del expected_forces_dict[expected_key]
        del expected_forces_dict['name']
        self.remove_almost_zero_force_elems(expected_forces_dict)
        self.assertEmpty(expected_forces_dict)

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

    def testDihedralPotentialSingleForceFactor(self):
        expected_force_factors = dict(ala_3.ala_3_dihedral_force_factors_tanh)
        potential = Dihedral_potential(self.get_single_member_ensemble_simulation())
        
        potential._force_calculator._set_components(potential._get_component_list().get_native_components())
        for i, data in enumerate(potential.dump()):
            SELECTION_INDEX = 0
            TARGET_ATOM_INDEX = 0
            dump_selection_data = data[SELECTION_INDEX]
            target_atom_key = dump_selection_data[TARGET_ATOM_INDEX]
            dihedral_angle_key = self._convert_dump_dihedral_key(dump_selection_data)
            
            force_factor_key = target_atom_key,dihedral_angle_key
            expected_force_factor =  expected_force_factors[force_factor_key]
            force_factor = potential._force_calculator._calc_single_force_factor(i)
            self.assertAlmostEqual(expected_force_factor, force_factor,self.DEFAULT_DECIMAL_PLACES)

            del expected_force_factors[force_factor_key]
        self.assertEmpty(expected_force_factors)

        
#        potential._calc_single_force_factor(index)
        

    def _extract_dihedral_forces(self, dihedral_atom_ids, forces):
        result = []
        
        for atom_id in dihedral_atom_ids:
            result.append(self.get_force_triplet(atom_id,forces))
            
        return result
                
    
    
    def _test_dihedral_forces(self, test_factors, potential,expected_forces):
        expected_forces_dict = dict(expected_forces)
        
        potential._force_calculator._set_components(potential._get_component_list().get_native_components())
        for i, data in enumerate(potential.dump()):
            
            
            out_array = self.make_out_array()
            #TODO imporve names and dihdedral atom lookup
            SELECTION_INDEX = 0
            TARGET_ATOM_INDEX = 0
            dump_selection_data = data[SELECTION_INDEX]
            target_atom_key = dump_selection_data[TARGET_ATOM_INDEX]
            dihedral_atoms_key = self._convert_dump_dihedral_key(dump_selection_data)
            expected_key = target_atom_key,dihedral_atoms_key
            dihedral_atom_ids =  []
            for elem in dihedral_atoms_key:
                dihedral_atom_ids.extend(single_text_key_to_atom_ids(elem))
            test_factor = test_factors[target_atom_key]
            potential._force_calculator._calc_single_force_set(i,test_factor,out_array)
            forces  =  out_array.add_forces_to_result()
            dihedral_forces = self._extract_dihedral_forces(dihedral_atom_ids,forces)

            
            for forces,expected in zip(dihedral_forces,expected_forces_dict[expected_key]):
                self.assertSequenceAlmostEqual(forces, expected, self.DEFAULT_DECIMAL_PLACES,msg=`expected_key`)

            del expected_forces_dict[expected_key]
        del expected_forces_dict['name']
        self.remove_almost_zero_force_elems(expected_forces_dict)
        self.assertEmpty(expected_forces_dict)
        
    #TODO for completeness there ought to be tests that the well forces are zero here
    def testDistancePotentialSingleForceHarmonic(self):
        distance_potential = Distance_potential(self.get_single_member_ensemble_simulation())
        expected_forces = ala_3.ala_3_distance_real_forces_harmonic
        factors_harmonic = ala_3.ala_3_factors_harmonic
        
        self._test_distance_forces(factors_harmonic, distance_potential,expected_forces)

    def testExtraPotentialSingleForceHarmonic(self):
        extra_potential = Extra_potential(self.get_single_member_ensemble_simulation())
        expected_forces = ala_3.ala_3_extra_real_forces_harmonic
        factors_harmonic = ala_3.ala_3_factors_harmonic
        
        self._test_distance_forces(factors_harmonic, extra_potential ,expected_forces)
        
    def testSidechainPotentialSingleForceHarmonic(self):
        sidechain_potential = Sidechain_potential(self.get_single_member_ensemble_simulation())
        expected_forces = ala_3.ala_3_sidechain_real_forces_harmonic
        factors_harmonic = ala_3.ala_3_factors_harmonic

        self._test_distance_forces(factors_harmonic, sidechain_potential,expected_forces)
        
    def testDihedralPotentialSingleForceTanh(self):
        dihedral_potential = Dihedral_potential(self.get_single_member_ensemble_simulation())
        expected_forces = ala_3.ala_3_dihedral_forces_tanh
        factors_tanh = ala_3.ala_3_factors_tanh

        self._test_dihedral_forces(factors_tanh, dihedral_potential,expected_forces)
    

    def _test_total_energy(self, xcamshift, expected_energy):
        calculated_energy = xcamshift.calcEnergy()
        self.assertAlmostEqual(calculated_energy, expected_energy, self.DEFAULT_DECIMAL_PLACES-1)


    def _setup_xcamshift_with_shifts_table(self, test_shifts):
        xcamshift = Xcamshift()
        observed_shifts = Observed_shift_table(test_shifts)
        xcamshift.set_observed_shifts(observed_shifts)
        return xcamshift

    def testCalcTotalEnergyWell(self):
        xcamshift = self._make_xcamshift(ala_3.ala_3_test_shifts_well)
        expected_energy = ala_3.ala_3_energies_well[TOTAL_ENERGY]
        
        self._test_total_energy(xcamshift, expected_energy)
        
    def testCalcTotalEnergyHarmonic(self):
        xcamshift = self._make_xcamshift(ala_3.ala_3_test_shifts_harmonic)
        expected_energy = ala_3.ala_3_energies_harmonic[TOTAL_ENERGY]
        
        self._test_total_energy(xcamshift, expected_energy)

    def testCalcTotalEnergyTanh(self):
        xcamshift = self._make_xcamshift(ala_3.ala_3_test_shifts_tanh)
        expected_energy = ala_3.ala_3_energies_tanh[TOTAL_ENERGY]
        
        self._test_total_energy(xcamshift, expected_energy)
    

    def _test_force_sets(self, xcamshift, expected_energy, expected_forces):
        expected_forces = dict(expected_forces)
        number_atoms = Segment_Manager().get_number_atoms()
        derivs = DerivList()
        derivs.init(currentSimulation())
        energy = xcamshift.calcEnergyAndDerivs(derivs)
        derivs_array =  derivs.get(currentSimulation())
        
        self.assertAlmostEqual(energy, expected_energy, self.DEFAULT_DECIMAL_PLACES-1)
        for atom_id in range(number_atoms):
            
            atom_key  =  Atom_utils._get_atom_info_from_index(atom_id)
            if atom_key in expected_forces:
                force_triplet = derivs_array[atom_id]
                expected_force_triplet = expected_forces[atom_key]
                
                self.assertSequenceAlmostEqual(force_triplet, expected_force_triplet, self.DEFAULT_DECIMAL_PLACES, `atom_key`)
                del expected_forces[atom_key]
                
        self.remove_almost_zero_force_elems(expected_forces, self.DEFAULT_DECIMAL_PLACES)
        self.assertEmpty(expected_forces)
                
        self.assertAlmostEqual(energy, expected_energy, self.DEFAULT_DECIMAL_PLACES-1)
        
    
    def testCalcForceSetWell(self):
        xcamshift = self._make_xcamshift(ala_3.ala_3_test_shifts_well)
        expected_energy = 0.0
        expected_forces = {}
        
        self._test_force_sets(xcamshift, expected_energy, expected_forces)

    
    def testCalcForceSetTanh(self):
        xcamshift = self._make_xcamshift(ala_3.ala_3_test_shifts_tanh)

#        xcamshift._prepare(self.)
        xcamshift.update_force_factor_calculator()
        expected_energy = ala_3.ala_3_energies_tanh[TOTAL_ENERGY]
        expected_forces =ala_3.ala_3_total_forces_tanh
        
        self._test_force_sets(xcamshift, expected_energy, expected_forces)
    

    def _split_tuple_pair_to_column_sets(self, expected_non_bonded):
        non_bonded_1 = set()
        non_bonded_2 = set()
        for non_bonded_elem_1, non_bonded_elem_2 in expected_non_bonded:
            non_bonded_1.add(non_bonded_elem_1)
            non_bonded_2.add(non_bonded_elem_2)
        
        return non_bonded_1, non_bonded_2

    def testNonBondedComponents(self):
        non_bonded_potential = Non_bonded_potential(self.get_single_member_ensemble_simulation())
        
        non_bonded_table = Table_manager().get_non_bonded_table("ALA")
        target_atom_types = list(non_bonded_table.get_target_atoms())
        
#        spheres = non_bonded_table.get_spheres()
        
        for component in non_bonded_potential._get_all_components('ATOM'):
            atom_info = Atom_utils._get_atom_info_from_index(component[0])
            atom_name  =  atom_info[2]
            
            self.assertIn(atom_name, target_atom_types)
            
            del target_atom_types[target_atom_types.index(atom_name)]
            
        self.assertEmpty(target_atom_types)
        
        remote_component_chem_type_ids = set()
        for component in non_bonded_potential._get_all_components('NBRM'):
            atom_info = Atom_utils._get_atom_info_from_index(component[0])
            remote_component_chem_type_ids.update(component[1:])
            
        remote_component_chem_type_ids =  list(remote_component_chem_type_ids)
        remote_component_chem_type_ids.sort()
        
        num_chem_types = len(remote_component_chem_type_ids) / 2
        remote_chem_type_ids_first_half = remote_component_chem_type_ids[:num_chem_types]
        for i, elem in enumerate(remote_chem_type_ids_first_half):
            self.assertEqual(elem, remote_component_chem_type_ids[i+num_chem_types]-8)
            
        self.assertEqual(remote_chem_type_ids_first_half, [0,1,2,4,5,6])
            
        
        non_bonded_1, non_bonded_2 = self._split_tuple_pair_to_column_sets(ala_3.ala3_putative_non_bonded_pairs)
        self.maxDiff = None
        test_non_bonded_set_1 = set()
        for component in non_bonded_potential._get_all_components('ATOM'):
            atom_info = Atom_utils._get_atom_info_from_index(component[0])
            self.assertElemeInSet(atom_info, non_bonded_1)
            test_non_bonded_set_1.add(atom_info)
        self.assertSequenceEqual(test_non_bonded_set_1, non_bonded_1)
        
        test_non_bonded_set_2 = set()
        for component in non_bonded_potential._get_all_components('NBRM'):
            atom_info = Atom_utils._get_atom_info_from_index(component[0])
            self.assertElemeInSet(atom_info, non_bonded_2)
            test_non_bonded_set_2.add(atom_info)
        self.assertSequenceEqual(test_non_bonded_set_2, non_bonded_2)
        
        components_0 = non_bonded_potential._get_components_for_id('NBRM',0)
        chem_types_0 = components_0[0][1:]
        self.assertLengthIs(components_0,1)
        self.assertEqual(components_0[0], (0,12,4))
        
        
        coefficents_0 = non_bonded_potential._get_components_for_id('COEF',chem_types_0[0])[0][3:]
        coefficents_1 = non_bonded_potential._get_components_for_id('COEF',chem_types_0[1])[0][3:]
        
        expected_coefficients_0 = (-0.0069306876125918232, 0.047420617840734335, 0.043722646930490425, 0.0093801308883550549, -0.063176276536097933, 0.34500777174208075)
        expected_coefficients_1 = (31.32384112711744, -12.5982466223797, -14.86605302072972, -1.6214566324137984, 2.505391454472044,  -88.71419716149445)
        self.assertSequenceAlmostEqual(coefficents_0, expected_coefficients_0, self.DEFAULT_DECIMAL_PLACES)
        self.assertSequenceAlmostEqual(coefficents_1, expected_coefficients_1, self.DEFAULT_DECIMAL_PLACES)

    
    def test_non_bonded_distances_found(self):
        non_bonded_list = Non_bonded_list(self.get_single_member_ensemble_simulation(),min_residue_separation=1)
        
        expected_non_bonded_pairs = set(ala_3.ala3_expected_non_bonded_pairs)
        non_bonded_potential = Non_bonded_potential(self.get_single_member_ensemble_simulation())
        
        target_atoms = non_bonded_potential._create_component_list('ATOM')
        target_atoms.add_components(non_bonded_potential._get_all_components('ATOM'))
        native_target_atoms = target_atoms.get_native_components()
        
        remote_atoms  =  non_bonded_potential._create_component_list('NBRM')
        remote_atoms.add_components(non_bonded_potential._get_all_components('NBRM'))
        native_remote_atoms = remote_atoms.get_native_components()
        
        coefficient_list = non_bonded_potential._get_component_list('COEF')
        
        component_list = non_bonded_potential._get_component_list('NBLT')
        
        non_bonded_list.get_boxes(target_atoms, remote_atoms, component_list, coefficient_list)
        
        for component in component_list:
            target_component_index,remote_component_index =  component[1:3]

            target_atom_id = target_atoms[target_component_index][0]
            distant_atom_id = remote_atoms[remote_component_index][0]

            target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)
            box_atom_key = Atom_utils._get_atom_info_from_index(distant_atom_id)
            atom_key = target_atom_key,box_atom_key

            self.assertElemeInSet(atom_key, expected_non_bonded_pairs)
            expected_non_bonded_pairs.remove(atom_key)

        self.assertEmpty(expected_non_bonded_pairs)
            
        
    @staticmethod
    def list_test_shifts():
        for item in ala_3.ala_3_test_shifts_tanh.items():
            print item
    
#    def test_fast_flag(self):
#        xcamshift = self._make_xcamshift(set_fast=False)
#        
#
#        xcamshift.calcEnergy()
#        for potential_name in xcamshift.get_sub_potential_names():
#            potential  = xcamshift.get_named_sub_potential(potential_name)
#            self.assertFalse(potential._fast, potential_name)
#        xcamshift.set_fast(True)
#        
#        xcamshift.calcEnergy()
#        for potential_name in xcamshift.get_sub_potential_names():
#            potential  = xcamshift.get_named_sub_potential(potential_name)
#            self.assertTrue(potential._fast, potential_name)        

def run_tests():
    unittest2.main(module='test.test_xcamshift')
#     unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.test_non_bonded_distances_found')
    
if __name__ == "__main__":
    run_tests()
#    TestXcamshift.list_test_shifts()
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testCalcForceSetTanh')
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testDistances')
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testSingleFactorHarmonic')
