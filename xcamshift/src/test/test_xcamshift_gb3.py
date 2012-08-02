'''
Created on 31 Dec 2011

@author: garyt
'''
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from xcamshift import  Xcamshift
from atomSel import AtomSel
from test import gb3
from observed_chemical_shifts import Observed_shift_table
from common_constants import  SIDE_CHAIN, NON_BONDED, RING
from test.gb3 import gb3_component_shifts_sc, gb3_component_shifts_ring
from utils import Atom_utils
import common_constants
TOTAL_ENERGY = 'total'


def almostEqual(first, second, places = 7):
    result  = False
    if round(abs(second-first), places) == 0:
        result=True
    return result

ZEROS_3 = [0.0] * 3
class TestXcamshiftGB3(unittest2.TestCase):

    # TODO add extra places
    DEFAULT_DECIMAL_PLACES = 5
    DEFAULT_ERROR = 10**-DEFAULT_DECIMAL_PLACES
    
    def assertEmpty(self, expected_force_factors,msg = None):
        return self.assertEqual(len(expected_force_factors), 0)
#    
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
        if self.check_almost_equal(list_1, list_2, delta) > -1:
            result = False
        return result
        
    def assertSequenceAlmostEqual(self,result,expected,places = 7,msg=''):
        delta  = 10**-places
        if len(msg) > 0:
            msg = msg + " "
        len_result = len(result)
        len_expected = len(expected)
        if len_result != len_expected:
            raise AssertionError((msg + "the two lists are of different length %i and %i") % (len_result,len_expected))
        
        difference_offset = self.check_almost_equal(result, expected, delta)

        if difference_offset > -1:
            template = "lists differ at item %i: %s - %s > %s"
            elem_1 = result[difference_offset]
            elem_2 = expected[difference_offset]
            message = template % (difference_offset, `elem_1`,`elem_2`,delta)
            raise AssertionError(msg + message)
        
    def remove_almost_zero_force_elems(self, expected_forces_dict, delta = 1e-7):
        ZEROS_3 = 0.0, 0.0, 0.0
        for key, value in expected_forces_dict.items():
            if self.check_almost_equal(ZEROS_3, value):
                del expected_forces_dict[key]         

    def setUp(self):
        initStruct("test_data/gb3/gb3.psf")
        PDBTool("test_data/gb3/gb3_refined_II.pdb").read()
        Atom_utils.clear_cache()
#
##TODO: shoulf be private
    def make_result_array_forces(self):
#        TODO: use segment manager
        num_atoms = len(AtomSel('(all)').indices())
        result = [None] * num_atoms
        return result
    
    def make_result_array(self):
        num_atoms = len(AtomSel('(all)').indices())
        result = [0.0] * num_atoms
        return result

    def assertListVec3AlmostEqual(self, ring_centres, expecteds, places = DEFAULT_DECIMAL_PLACES):
        self.assertEqual(len(ring_centres), len(expecteds))
        for i,(result, expected) in enumerate(zip(ring_centres, expecteds)):
            if (((result == None) and (expected == None)) or 
               ((result == None) and self.is_almost_zero_sequence(expected)) or 
               ((expected == None) and self.is_almost_zero_sequence(result))):
                continue
            elif result == None or expected == None:
                message = "only one item is None at %i [%s,%s]"
                raise AssertionError(message % (i,`result`,`expected`))
            else:
                self.assertSequenceAlmostEqual(result, expected, places)

    def is_almost_zero_sequence(self,sequence,places = DEFAULT_DECIMAL_PLACES):
        len_sequence = len(sequence)
        if len_sequence==3:
            zeros = ZEROS_3
        else:
            zeros = [0.0] *len_sequence
        delta  = 10**-places
        return self.are_almost_equal_sequences(sequence, zeros, delta)

    def test_component_chemical_shifts(self):
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
                self.assertAlmostEqual(shift, expected_shift, self.DEFAULT_DECIMAL_PLACES-2, msg=`key` + " " + residue_type)


        
    def _setup_xcamshift_with_shifts_table(self, test_shifts):
        xcamshift = Xcamshift()
        observed_shifts = Observed_shift_table(test_shifts)
        xcamshift.set_observed_shifts(observed_shifts)
        return xcamshift
    
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

    def test_component_shifts_sidechain(self):

        
        xcamshift = Xcamshift()
        sidechain_subpotential = xcamshift.get_named_sub_potential(SIDE_CHAIN)
        
        expected_sidechain_shifts = dict(gb3_component_shifts_sc)
        expected_component_keys = expected_sidechain_shifts.keys()
        for component_index, component in enumerate(sidechain_subpotential._get_component_list()):
            from_atom_id, to_atom_id = component[0:2]
            from_atom_key = Atom_utils._get_atom_info_from_index(from_atom_id)
            to_atom_key = Atom_utils._get_atom_info_from_index(to_atom_id)
            
#            print from_atom_key, to_atom_key
            if from_atom_key[2] == 'HA1':
                from_atom_key = from_atom_key[0],from_atom_key[1],'HA'
                
            expected_key = from_atom_key, to_atom_key
            
            self.assertIn(expected_key, expected_component_keys, `expected_key` + " exists")
            
            shift = sidechain_subpotential._calc_component_shift(component_index)
            
            self.assertAlmostEqual(expected_sidechain_shifts[expected_key], shift, places=self.DEFAULT_DECIMAL_PLACES - 1, msg=`expected_key`)
            
            del expected_sidechain_shifts[expected_key]
            
        self.remove_zero_valued_keys(expected_sidechain_shifts)
        self.assertEmpty(expected_sidechain_shifts)
        

    def test_component_shifts_ring(self):
        
        xcamshift = Xcamshift()
        ring_subpotential = xcamshift.get_named_sub_potential(RING)
        
        expected_ring_shifts = dict(gb3_component_shifts_ring)
        expected_component_keys = expected_ring_shifts.keys()
        for component_index, component in enumerate(ring_subpotential._get_component_list()):
            from_atom_id, atom_type_id = component
            from_atom_info_list = ring_subpotential._get_component_list('COEF').get_components_for_atom_id(atom_type_id)
            
#            print Atom_utils._get_atom_info_from_index(from_atom_id), from_atom_info_list
            from_atom_key = list(Atom_utils._get_atom_info_from_index(from_atom_id))
            
            if from_atom_key[2] == 'HA1':
                from_atom_key[2] =  'HA'
            from_atom_key = tuple(from_atom_key)
            
            for sub_component_index, from_atom_info in enumerate(from_atom_info_list):
                from_atom_id, ring_id, coefficent = from_atom_info
                ring_info = ring_subpotential._get_component_list('RING').get_components_for_atom_id(ring_id)
                
                ring_atoms =  ring_info[0][1]
                ring_residue_type = Atom_utils._get_residue_type_from_atom_id(ring_atoms[1])
                ring_residue = Atom_utils._get_atom_info_from_index(ring_atoms[0])[1]
                expected_key =   from_atom_key,(ring_residue,ring_residue_type,len(ring_atoms))
                self.assertIn(expected_key, expected_component_keys, `expected_key`)
            
                
                shift = ring_subpotential._shift_calculator._calc_sub_component_shift(component, from_atom_info)
                self.assertAlmostEqual(expected_ring_shifts[expected_key], shift, places=self.DEFAULT_DECIMAL_PLACES - 2, msg=`expected_key`)
                if abs(expected_ring_shifts[expected_key] - shift) > 0.001:
                    print 'fail', expected_key, expected_ring_shifts[expected_key], shift, Atom_utils._get_residue_type_from_atom_id(from_atom_id)
                    print
                
                del expected_ring_shifts[expected_key]

            
        self.remove_zero_valued_keys(expected_ring_shifts)
        self.assertEmpty(expected_ring_shifts)

    def test_force_components_add_up(self):
        summary = {}
        gb3_forces_copy  = dict(gb3.gb3_forces)
        
        
        for elem in gb3.gb3_component_forces:
            key =  self.get_force_atom_key(elem)
            summary_forces = summary.setdefault(key,[0.0,]*3)
            for i, value in enumerate(gb3.gb3_component_forces[elem]):
                summary_forces[i] += value
        
        for i,elem in enumerate(sorted(summary)):

            self.assertSequenceAlmostEqual(summary[elem], gb3_forces_copy[elem], 2, `elem`)

            del gb3_forces_copy[elem]
        
        self.remove_almost_zero_force_elems(gb3_forces_copy)
        
        self.assertEmpty(gb3_forces_copy)
        

    @staticmethod
    def filter_dict(the_dict, predicate=lambda k, v: True):
        for k, v in the_dict.iteritems():
            if predicate(k,v):
                yield k, v    
    
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

    def test_energies(self):
        xcamshift  = self._setup_xcamshift_with_shifts_table(gb3.gb3_zero_shifts)

        
        non_bonded_subpotential = xcamshift.get_named_sub_potential(NON_BONDED)
        non_bonded_subpotential.update_non_bonded_list()
        
        total_component_energies = 0.0
        for key in sorted(gb3.gb3_energies):
            if key ==  'total':
                continue

            target_atom_index = self.get_atom_index(key)
            
            expected_energy =  gb3.gb3_energies[key]
            energy = xcamshift._calc_single_atom_energy(target_atom_index)
            
            total_component_energies += energy
            
            self.assertAlmostEqual(energy,expected_energy, self.DEFAULT_DECIMAL_PLACES-3,key)
        
        total_energy = xcamshift.calcEnergy()
        expected_total_energy = gb3.gb3_energies['total']
        
        self.assertAlmostEqual(total_component_energies, total_energy, self.DEFAULT_DECIMAL_PLACES)
        #TODO: add the ability to test almost equals with significant figures instead
        self.assertAlmostEqual(total_energy,expected_total_energy, self.DEFAULT_DECIMAL_PLACES-4,)
        
    def test_shift_differences(self):
        xcamshift  = self._setup_xcamshift_with_shifts_table(gb3.gb3_zero_shifts)
        
        non_bonded  = xcamshift.get_named_sub_potential(NON_BONDED)
        non_bonded.update_non_bonded_list()        
        
        
        for key in sorted(gb3.gb3_shift_diffs):
            expected_shift_diff  = gb3.gb3_shift_diffs[key]
            
            target_atom_index = self.get_atom_index(key)
            
            shift_diff = xcamshift.calc_single_atom_shift(target_atom_index)

            self.assertAlmostEqual(shift_diff, - expected_shift_diff, self.DEFAULT_DECIMAL_PLACES-2,key)
            
    def set_cache_shifts(self,xcamshift, shift_list):
        shift_cache =  xcamshift._shift_cache
        
        for key in shift_list:
            index = self.get_atom_index(key)
            shift_cache[index] = shift_list[key]
        

    def add_forces_to_array(self, target_array, index, force_component):
        expected_forces_set = target_array[index]
        if expected_forces_set == None:
            target_array[index] = [0.0] * 3
            expected_forces_set = target_array[index]
        for i, new_force in enumerate(force_component):
            expected_forces_set[i] += new_force
        
        return target_array


    def filter_force_components_by_potential(self, potential_name):
        return dict(self.filter_dict(gb3.gb3_component_forces, lambda k, v:k[-1] == potential_name))

    def make_total_forces_from_components(self):
        expected_total_result_forces = self.make_result_array_forces()
        for key in gb3.gb3_forces:
            index = self.get_atom_index(key)
            expected_total_result_forces[index] = gb3.gb3_forces[key]
        
        return expected_total_result_forces

    def test_total_forces_and_energy(self):
        xcamshift =  self.make_xcamshift(gb3.gb3_zero_shifts)
        
        total_result_forces = self.make_result_array_forces()
        energy = xcamshift.calcEnergyAndDerivs(total_result_forces)

        expected_total_result_forces = self.make_total_forces_from_components()
        
        self.assertListVec3AlmostEqual(total_result_forces, expected_total_result_forces, self.DEFAULT_DECIMAL_PLACES - 3)
        
        self.assertAlmostEqual(energy, gb3.gb3_energies['total'],self.DEFAULT_DECIMAL_PLACES-4)


        
    def make_xcamshift(self, shifts):
        xcamshift = self._setup_xcamshift_with_shifts_table(shifts)

        
        non_bonded = xcamshift.get_named_sub_potential(NON_BONDED)
        non_bonded.update_non_bonded_list()
        

        ring_potential = xcamshift.get_named_sub_potential(RING)
        ring_potential._build_ring_data_cache()
        
        return xcamshift


    def print_force_component_diffs(self, potential_name, target_atom_key, expected_forces, result_forces):
        diffs = []
        for i, (result_force, expected_force) in enumerate(zip(result_forces, expected_forces)):
            if result_force == None and expected_force == None:
                continue
            if (result_force == None) and (max(expected_force) < 0.0000001):
                continue
            if (expected_force == None) and (max(result_force) < 0.0000001):
                continue                      
            if result_force == None or expected_force == None:
                print 'BAD!', potential_name, i, target_atom_key, Atom_utils._get_atom_info_from_index(i),result_force,expected_force
                continue
#                    print 'here', result_force, expected_force
            current_diffs = []
            for num_1,num_2 in zip(result_force,expected_force):
#                        print num_1,num_2,abs(num_1-num_2)
                diffs.append(abs(num_1-num_2))
                current_diffs.append(abs(num_1-num_2))
            if max(current_diffs) > 0.01:
                print 'BAD!!',  potential_name, i, target_atom_key, Atom_utils._get_atom_info_from_index(i),max(current_diffs),result_force,expected_force
                continue
            else:
                print 'OK!!'
                
    @staticmethod
    def get_target_atom_key(raw_backbone_force_component_key):
        return raw_backbone_force_component_key[0]
    
    @staticmethod
    def get_force_atom_key(raw_backbone_force_component_key):
        return raw_backbone_force_component_key[2]
    

        
    def test_force_components(self):

        
        def get_target_atom_keys_for_potential(raw_backbone_force_components):
            target_atom_keys = set()
            for raw_backbone_force_component_key in raw_backbone_force_components:
                target_atom_key = self.get_target_atom_key(raw_backbone_force_component_key)
                target_atom_keys.add(target_atom_key)
        
                return target_atom_keys

        def build_expected_forces_for_potential(target_atom_key, raw_backbone_force_components):
            
            expected_forces = self.make_result_array_forces()
            for raw_backbone_force_component_key in raw_backbone_force_components:
                if target_atom_key == self.get_target_atom_key(raw_backbone_force_component_key):
                    force_atom_key = self.get_force_atom_key(raw_backbone_force_component_key)
                    force_atom_index = self.get_atom_index(force_atom_key)
                    expected_force_component = raw_backbone_force_components[raw_backbone_force_component_key]
                    self.add_forces_to_array(expected_forces, force_atom_index, expected_force_component)
        
            return expected_forces
        
        
#       TODO: xcamshift should be created anew  for each passage through the loop but the non bonded list is too slow
        xcamshift = self.make_xcamshift(gb3.gb3_zero_shifts)
        self.set_cache_shifts(xcamshift, gb3.gb3_shifts)
        
        for potential_name in common_constants.CAMSHIFT_SUB_POTENTIALS:
            
            raw_backbone_force_components = self.filter_force_components_by_potential(potential_name)
            
            target_atom_keys =  get_target_atom_keys_for_potential(raw_backbone_force_components)
            
            potentials_list = [xcamshift.get_named_sub_potential(potential_name)]
            for target_atom_key in reversed(sorted(target_atom_keys)):
                
        
                expected_forces = build_expected_forces_for_potential(target_atom_key, raw_backbone_force_components)

                target_atom_id = self.get_atom_index(target_atom_key)
                result_forces = self.make_result_array_forces()

                xcamshift._calc_single_atom_force_set_with_potentials(target_atom_id, result_forces, potentials_list)
                self.assertListVec3AlmostEqual(result_forces, expected_forces, self.DEFAULT_DECIMAL_PLACES-3)

                # to list differences as well as check them uncomment this! 
                # self.print_force_component_diffs(potential_name,target_atom_key, expected_forces, result_forces)
    
#    def compare_shift_diffs_and_expected(self):
#        for elem in sorted(gb3.gb3_shifts):
#            self.assertAlmostEqual(gb3.gb3_shifts[elem], -gb3.gb3_shift_diffs[elem], self.DEFAULT_DECIMAL_PLACES-3,  elem)



if __name__ == "__main__":
    unittest2.main()
#    cProfile.run('unittest2.main()')
    
#    TestXcamshift.list_test_shifts()
#    unittest2.main(module='test.test_xcamshift_gb3',defaultTest='TestXcamshiftGB3.test_forces')
#    unittest2.main(module='test.test_xcamshift_gb3',defaultTest='TestXcamshiftGB3.test_energies')
#    unittest2.main(module='test.test_xcamshift_gb3',defaultTest='TestXcamshiftGB3.test_component_chemical_shifts')
#    unittest2.main(module='test.test_xcamshift_gb3',defaultTest='TestXcamshiftGB3.test_force_components_add_up')
#    unittest2.main(module='test.test_xcamshift_gb3',defaultTest='TestXcamshiftGB3.test_force_components2')
#    unittest2.main(module='test.test_xcamshift_gb3',defaultTest='TestXcamshiftGB3.test_total_forces_and_energy')
    
#    unittest2.main(module='test.test_xcamshift',defaultTest='TestXcamshift.testSingleFactorHarmonic')
