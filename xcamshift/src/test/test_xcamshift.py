'''
Created on 31 Dec 2011

@author: garyt
'''
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from xcamshift import RandomCoilShifts, Distance_potential,  Extra_potential,\
    Dihedral_potential, Base_potential, Sidechain_potential, Xcamshift
from atomSel import AtomSel
from test.xdists import xdists_ala_3
from test.dihedrals import dihedrals_ala_3
from segment_manager import Segment_Manager
from test.sidechains import sidechains_ala_3
from test.ala_3 import ala_3_total_shifts
from test import sidechains, ala_3
from observed_chemical_shifts import Observed_shift_table
from utils import Atom_utils
import sys


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
        initStruct("test_data/3_ala/3ala.psf")
        PDBTool("test_data/3_ala/3ala.pdb").read()


    def make_result_array_forces(self):
        num_atoms = len(AtomSel('(all)').indices())
        result = [0.0] * num_atoms * 3
        return result
    
    def make_result_array(self):
        num_atoms = len(AtomSel('(all)').indices())
        result = [0.0] * num_atoms
        return result

    def testRandomCoilPotential(self):
        random_coil_potential = RandomCoilShifts()
        result = self.make_result_array() 
        expected  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                     123.6, 8.2400000000000002, 52.266440000000003, 4.4328469999999998, 
                     19.0, 0.0, 0.0, 0.0, 177.09999999999999, 0.0, 0.0, 0.0, 0.0, 0.0, 
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        random_coil_potential.set_shifts(result)
        self.assertSequenceAlmostEqual(expected, result)
        
    def testDistancePotential(self):
        distance_potential =  Distance_potential()
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
        
        
    @staticmethod
    def assertElemeInSet( elem, xdist_set):
        if not elem in xdist_set:
            raise AssertionError("%s not found in set" % `elem`)


    def testExtraPotentialComponentsCorrect(self):
        extra_potential = Extra_potential()
        
        xdists_ala_3_copy = dict(xdists_ala_3)
        for extra_elem in extra_potential.dump():
            elem_key = extra_elem[:-1]
            self.assertElemeInSet(elem_key, xdists_ala_3)
            del xdists_ala_3_copy[elem_key]
        
        self.assertEqual(0, len(xdists_ala_3_copy))


    def testExtraPotentialCoefficientsCorrect(self):
        extra_potential = Extra_potential()
        for extra_elem in extra_potential.dump():
            elem_key = extra_elem[:-1]
            coefficient = xdists_ala_3[elem_key][0]

            self.assertAlmostEqual(coefficient, extra_elem[-1], places=self.DEFAULT_DECIMAL_PLACES)
        

    def testExtraPotentialComponentShiftsCorrect(self):
        extra_potential = Extra_potential()
        
        result=self.make_result_array()
        extra_potential.set_shifts(result)
        
        for i,extra_elem in enumerate(extra_potential.dump()):
            
            elem_key = extra_elem[:-1]
            shift = extra_potential._calc_single_shift(i)
            
            expected_shift = xdists_ala_3[elem_key][2]
            
            self.assertAlmostEqual(shift, expected_shift,places=4)

    def testDihdedralPotentialComponentsCorrect(self):
        dihedral_portential = Dihedral_potential()
        
        
        dihedrals_ala_3_copy = dict(dihedrals_ala_3)
        for dihedral_elem in dihedral_portential.dump():
            self.assertElemeInSet(dihedral_elem[0], dihedrals_ala_3_copy)
            del dihedrals_ala_3_copy[dihedral_elem[0]]

        self.assertEqual(0, len(dihedrals_ala_3_copy))
        
    def testDihedralPotentialCoefficientsCorrect(self):
        dihedral_potential = Dihedral_potential()
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
        dihedral_potential = Dihedral_potential()
        
        for dihedral_element in dihedral_potential.dump():
            
            dihedral_atoms = self.dihedral_key_to_atom_ids(dihedral_element)
            atom_ids  = text_keys_to_atom_ids(dihedral_atoms)
            
            expected  =  dihedrals_ala_3[dihedral_element[0]][2]
            angle =  dihedral_potential._get_dihedral_angle(*atom_ids)
            
            self.assertAlmostEqual(expected, angle,self.DEFAULT_DECIMAL_PLACES)

            
    def testDihedralPotentialComponentShiftsCorrect(self):
        dihedral_potential = Dihedral_potential()
        
        for i,dihedral_element in enumerate(dihedral_potential.dump()):
            
            expected  =  dihedrals_ala_3[dihedral_element[0]][3]
            angle =  dihedral_potential._calc_single_shift(i)
            
            self.assertAlmostEqual(expected, angle,self.DEFAULT_DECIMAL_PLACES)
#
####
    def testSidechainPotentialComponentsCorrect(self):
        sidechain_potential = Sidechain_potential()
        
        sidechains_ala_3_copy = dict(sidechains_ala_3)
        for sidechain_elem  in sidechain_potential.dump():
            expected_key = tuple(sidechain_elem[:2])
            self.assertElemeInSet(expected_key, sidechains_ala_3_copy)
            del sidechains_ala_3_copy[expected_key]

        self.assertEqual(0, len(sidechains_ala_3_copy))
        
    def testSidechainPotentialCoefficientsCorrect(self):
        sidechain_potential = Sidechain_potential()
        for target_atom,distant_atom, coeff, exponent in  sidechain_potential.dump():
            sidechain_potential_key=target_atom,distant_atom
            expected_values = sidechains_ala_3[sidechain_potential_key][:2]
            
            actual_values  = (coeff,exponent)
            
            self.assertEqual(len(expected_values), len(actual_values))
            for expected,actual in zip(expected_values,actual_values):
                self.assertAlmostEqual(expected, actual, self.DEFAULT_DECIMAL_PLACES)

            
            self.assertAlmostEqual(exponent, 1.0)
         
    def testSidechainlPotentialComponentShiftsCorrect(self):
        sidechain_potential = Sidechain_potential()
        
        for i,elem in enumerate(sidechain_potential.dump()):

            key = tuple(elem[:2])
            
            shift =  sidechain_potential._calc_single_shift(i)
            expected_shift  = sidechains.sidechains_ala_3[key][2]
            self.assertAlmostEqual(expected_shift, shift, self.DEFAULT_DECIMAL_PLACES)

    def testXcamshift(self):
        xcamshift_potential =  Xcamshift()
        
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
        
    def testDistancePotentialSingleEnergiesHarmonic(self):
        test_shifts = ala_3.ala_3_test_shifts_harmonic
        
        xcamshift = Xcamshift()
        
        shift_table = Observed_shift_table(test_shifts)
        xcamshift.set_observed_shifts(shift_table)
        
        expected_energys = ala_3.ala_3_energies
        expected_total_energy = ala_3.ala_3_energies['total']
        
        total_energy = 0.0
        for atom_index in shift_table.get_atom_indices():
            key = Atom_utils._get_atom_info_from_index(atom_index)[1:]
            expected_energy = expected_energys[key]

            energy = xcamshift._calc_single_energy(atom_index)
            self.assertAlmostEqual(energy, expected_energy,self.DEFAULT_DECIMAL_PLACES)

            total_energy += energy
        self.assertAlmostEqual(total_energy, expected_total_energy,self.DEFAULT_DECIMAL_PLACES)
        
    def testDistancePotentialSingleEnergiesWell(self):
        test_shifts = ala_3.ala_3_test_shifts_well
        
        xcamshift = Xcamshift()
        
        shift_table = Observed_shift_table(test_shifts)
        xcamshift.set_observed_shifts(shift_table)
        
        
        total_energy = 0.0
        for atom_index in shift_table.get_atom_indices():
            key = Atom_utils._get_atom_info_from_index(atom_index)[1:]
            expected_energy = 0.0

            energy = xcamshift._calc_single_energy(atom_index)
            self.assertAlmostEqual(energy, expected_energy,self.DEFAULT_DECIMAL_PLACES,msg=`key`)
        total_energy += energy

    def testSingleFactorHarmonic(self):
        test_shifts = ala_3.ala_3_test_shifts_harmonic
        
        xcamshift = Xcamshift()
        
        shift_table = Observed_shift_table(test_shifts)
        xcamshift.set_observed_shifts(shift_table)
        
        
        for atom_index in shift_table.get_atom_indices():
            key = Atom_utils._get_atom_info_from_index(atom_index)[1:]
            factor  = xcamshift._calc_single_factor(atom_index)
            expected_factor = ala_3.ala_3_factors_harmonic[key]
            
            self.assertAlmostEqual(factor, expected_factor, self.DEFAULT_DECIMAL_PLACES)

    

    def remove_zero_valued_keys(self, expected_force_factors):
        for key, value in expected_force_factors.items():
            if almostEqual(value, 0.0, self.DEFAULT_DECIMAL_PLACES):
                del expected_force_factors[key]


    def assertEmpty(self, expected_force_factors,msg = None):
        return self.assertEqual(len(expected_force_factors), 0)

    def testDistanceBasedPotentialSingleForceFactor(self):
        
        # TODO make this a dummy distance based potential, 
        # _calc_single_force_factor is in distance based potential
        distance_potential = Distance_potential()
        expected_force_factors = dict(ala_3.ala_3_distance_forces_well)
        test_force_factors = ala_3.ala_3_factors_harmonic
        
        for i,data in enumerate(distance_potential.dump()):
            target_atom_key = data[0][1:]
            distant_atom_key = data[1][1:]
            expected_key = (target_atom_key,distant_atom_key)
            test_factor  =  test_force_factors[target_atom_key]
            
            force = distance_potential._calc_single_force_factor(i, test_factor)
            expected_force_factor = expected_force_factors[expected_key]
            del expected_force_factors[expected_key]

            self.assertAlmostEqual(force, expected_force_factor, self.DEFAULT_DECIMAL_PLACES-2)
        
        self.remove_zero_valued_keys(expected_force_factors)
        self.assertEmpty(expected_force_factors)


    def get_force_triplet(self, target_atom_index, forces):
        DIM_3 = 3
        target_atom_start_index = target_atom_index * DIM_3
        target_atom_end_index = target_atom_index * DIM_3 + DIM_3
        target_atom_forces = forces[target_atom_start_index:target_atom_end_index]
        return target_atom_forces


    def remove_almost_zero_force_elems(self, expected_forces_dict, delta = 1e-7):
        ZEROS_3 = 0.0, 0.0, 0.0
        for key, value in expected_forces_dict.items():
            if self.check_almost_equal(ZEROS_3, value):
                del expected_forces_dict[key]


    def _test_forces(self, test_factors, distance_potential,expected_forces):
        expected_forces_dict = dict(expected_forces)
        test_factors = ala_3.ala_3_factors_harmonic
        
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
            
            result_array = self.make_result_array_forces()
            
            forces = distance_potential._calc_single_force_set(i, test_factor, result_array)
            
            distant_atom_forces_1 = self.get_force_triplet(distant_atom_index_1, forces)
            distant_atom_forces_2 = self.get_force_triplet(distant_atom_index_2, forces)
            
            self.assertSequenceAlmostEqual(distant_atom_forces_1, expected_forces, self.DEFAULT_DECIMAL_PLACES)
            self.assertSequenceAlmostEqual(distant_atom_forces_2, negative_expected_forces, self.DEFAULT_DECIMAL_PLACES)

            del expected_forces_dict[expected_key]
        del expected_forces_dict['name']
        self.remove_almost_zero_force_elems(expected_forces_dict)
        self.assertEmpty(expected_forces_dict)

    def testDistancePotentialSingleForceHarmonic(self):
        distance_potential = Distance_potential()
        expected_forces = ala_3.ala_3_distance_real_forces_harmonic
        factors_harmonic = ala_3.ala_3_factors_harmonic
        
        self._test_forces(factors_harmonic, distance_potential,expected_forces)

    def testExtraPotentialSingleForceHarmonic(self):
        extra_potential = Extra_potential()
        expected_forces = ala_3.ala_3_extra_real_forces_harmonic
        factors_harmonic = ala_3.ala_3_factors_harmonic
        
        self._test_forces(factors_harmonic, extra_potential ,expected_forces)
        
    def testSidechainPotentialSingleForceHarmonic(self):
        sidechain_potential = Sidechain_potential()
        expected_forces = ala_3.ala_3_sidechain_real_forces_harmonic
        factors_harmonic = ala_3.ala_3_factors_harmonic
        self._test_forces(factors_harmonic, sidechain_potential,expected_forces)
        
if __name__ == "__main__":
    unittest2.main()
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testSidechainPotentialSingleForceHarmonic')
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testSingleFactorHarmonic')
