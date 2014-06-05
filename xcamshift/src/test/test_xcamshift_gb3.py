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
from cython.shift_calculators import CDSVectorFloat, New_fast_non_bonded_calculator, Fast_non_bonded_calculator
from cython.shift_calculators import  Non_bonded_interaction_list
from nanotime import now
from functools import partial

def almostEqual(first, second, places = 7):
    result  = False
    if round(abs(second-first), places) == 0:
        result=True
    return result

ZEROS_3 = [0.0] * 3
class TestXcamshiftGB3(unittest2.TestCase):

    def __init__(self,*args,**kwargs):
        super(TestXcamshiftGB3, self).__init__(*args,**kwargs)
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

    #TODO: this should be triggered on structure change events
    def _clear_caches(self):
        Atom_utils.clear_cache()
        Segment_Manager.reset_segment_manager()

    def setUp(self):
        initStruct("test_data/gb3/gb3.psf")
        PDBTool("test_data/gb3/gb3_refined_II.pdb").read()
        self._clear_caches()
#         print "In method", self._testMethodName

#
##TODO: shoulf be private
    def make_out_array(self):
#        TODO: use segment manager
        num_atoms = len(AtomSel('(all)').indices())
        result = Out_array( num_atoms, self.get_single_member_ensemble_simulation())
        return result

    def make_results_array(self):
        result = DerivList()
        result.init(currentSimulation())
        return result

    def assertListVec3AlmostEqual(self, ring_centres, expecteds, places = DEFAULT_DECIMAL_PLACES, msg=''):
        self.assertEqual(len(ring_centres), len(expecteds))
        for i,(result, expected) in enumerate(zip(ring_centres, expecteds)):
            if (((result == None) and (expected == None)) or
               ((result == None) and self.is_almost_zero_sequence(expected)) or
               ((expected == None) and self.is_almost_zero_sequence(result))):
                continue
            elif result == None or expected == None:
                message = "only one item is None at %i [%s,%s] %s"
                raise AssertionError(message % (i,`result`,`expected`,msg))
            else:
                self.assertSequenceAlmostEqual(result, expected, places,msg=("pos %i " % i) +msg)

    def is_almost_zero_sequence(self,sequence,places = DEFAULT_DECIMAL_PLACES):
        len_sequence = len(sequence)
        if len_sequence==3:
            zeros = ZEROS_3
        else:
            zeros = [0.0] *len_sequence
        delta  = 10**-places
        return self.are_almost_equal_sequences(sequence, zeros, delta)



    def _do_test_component_shifts(self, xcamshift, component_shifts):
        xcamshift._prepare(ROUND_CHANGED, None)
        component_shifts_keys = component_shifts.keys()
        component_shifts_keys.sort()
        for i, key in enumerate(component_shifts_keys):
            segment, residue_number, atom, sub_potential = key
            sub_potential = xcamshift.get_named_sub_potential(sub_potential)
            atom_ids = Atom_utils.find_atom_ids(segment, residue_number, atom)
            if len(atom_ids) > 0:

                xcamshift._prepare(TARGET_ATOM_IDS_CHANGED, atom_ids)
                shift = sub_potential._calc_single_atom_shift(atom_ids[0])
                expected_shift = component_shifts[key]
                residue_type = Atom_utils._get_residue_type_from_atom_id(atom_ids[0])
                self.assertAlmostEqual(shift, expected_shift, self.DEFAULT_DECIMAL_PLACES - 2, msg=`key` + " " + residue_type)

    #TODO check all tests correct

    def _get_xcamshift(self):
        xcamshift = Xcamshift()
        return xcamshift

    def _get_xcamshift_no_hbond(self):
        xcamshift =self._get_xcamshift()
        xcamshift.remove_named_sub_potential('HBOND')
        return xcamshift

    def test_component_chemical_shifts(self):
        xcamshift = self._get_xcamshift()


        component_shifts = gb3.gb3_subpotential_shifts
        self._do_test_component_shifts(xcamshift, component_shifts)



    def _setup_xcamshift_with_shifts_table(self, test_shifts):
        xcamshift = self._get_xcamshift_no_hbond()
        observed_shifts = Observed_shift_table(test_shifts)
        xcamshift.set_observed_shifts(observed_shifts)
        return xcamshift

    def test_non_bonded_components(self):
        #TODO: add common loading method for xcamshift
        xcamshift =  self._get_xcamshift()
        sub_potential = xcamshift.get_named_sub_potential(NON_BONDED)
        components = sub_potential._get_component_list('NBLT')
        target_atom_ids = [component[0] for component in components]
        xcamshift._prepare(TARGET_ATOM_IDS_CHANGED,target_atom_ids)
        xcamshift._prepare(ROUND_CHANGED,None)


        non_bonded_components =  dict(gb3.gb3_component_shifts_non_bonded)

        target_components = sub_potential._get_component_list('ATOM')
        remote_components = sub_potential._get_component_list('NBRM')

        for  target_atom_id, target_index, remote_index, component_index in components:

            target_atom_id = target_components[target_index][0]
            target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)

            remote_atom_id = remote_components[remote_index][0]
            remote_atom_key =  Atom_utils._get_atom_info_from_index(remote_atom_id)


            if target_atom_key[2]== 'HA1':
                target_atom_key =  list(target_atom_key)
                target_atom_key[2]= 'HA'
                target_atom_key =  tuple(target_atom_key)



            if Atom_utils._calculate_distance(target_atom_id, remote_atom_id) > 5.0:
                continue

            for exponent_key in 0,1:
                non_bonded_component_key = (target_atom_key, remote_atom_key, exponent_key)

                self.assertIn(non_bonded_component_key, non_bonded_components, `non_bonded_component_key`)
                del non_bonded_components[non_bonded_component_key]

        self.assertEmpty(non_bonded_components)

    def test_non_bonded_component_shifts(self):
        xcamshift =  self._get_xcamshift()
        sub_potential = xcamshift.get_named_sub_potential(NON_BONDED)
        sub_potential.update_non_bonded_list()
        components = sub_potential._get_component_list('NBLT')
        target_atom_ids = [component[0] for component in components]
        xcamshift._prepare(TARGET_ATOM_IDS_CHANGED,target_atom_ids)
        xcamshift._prepare(ROUND_CHANGED,None)


        non_bonded_components =  dict(gb3.gb3_component_shifts_non_bonded)

        target_components = sub_potential._get_component_list('ATOM')
        remote_components = sub_potential._get_component_list('NBRM')
        for  i,(target_atom_id, target_index, remote_index, component_index) in enumerate(components):

            target_atom_id = target_components[target_index][0]
            target_atom_key = Atom_utils._get_atom_info_from_index(target_atom_id)

            remote_atom_id = remote_components[remote_index][0]
            remote_atom_key =  Atom_utils._get_atom_info_from_index(remote_atom_id)


            if target_atom_key[2]== 'HA1':
                target_atom_key =  list(target_atom_key)
                target_atom_key[2]= 'HA'
                target_atom_key =  tuple(target_atom_key)



            if Atom_utils._calculate_distance(target_atom_id, remote_atom_id) > 5.0:
                continue

            shift = sub_potential._calc_component_shift(i)

            expected_shift = 0.0
            seen_keys = []
            for exponent_key in 0,1:
                non_bonded_component_key = (target_atom_key, remote_atom_key, exponent_key)
                seen_keys.append(non_bonded_component_key)
                expected_shift += non_bonded_components[non_bonded_component_key]

            self.assertAlmostEqual(shift,expected_shift,self.DEFAULT_DECIMAL_PLACES-2,`i` + ' '  + `non_bonded_component_key`)

            for seen_key in seen_keys:
                del non_bonded_components[seen_key]

        self.assertEmpty(non_bonded_components)

    def remove_zero_valued_keys(self, expected_force_factors):
        for key, value in expected_force_factors.items():
            if almostEqual(value, 0.0, self.DEFAULT_DECIMAL_PLACES):
                del expected_force_factors[key]

    def test_component_shifts_sidechain(self):


        xcamshift = self._get_xcamshift()
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

            self.assertAlmostEqual(expected_sidechain_shifts[expected_key], shift, places=self.DEFAULT_DECIMAL_PLACES - 2, msg=`expected_key`)

            del expected_sidechain_shifts[expected_key]

        self.remove_zero_valued_keys(expected_sidechain_shifts)
        self.assertEmpty(expected_sidechain_shifts)


    def test_component_shifts_ring(self):

        xcamshift = self._get_xcamshift()

        xcamshift._prepare(TARGET_ATOM_IDS_CHANGED, xcamshift._get_all_component_target_atom_ids())
        xcamshift._prepare(ROUND_CHANGED,None)

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
                from_atom_type, ring_id, coefficent = from_atom_info
                ring_info = ring_subpotential._get_component_list('RING').get_components_for_atom_id(ring_id)

                ring_atoms =  ring_info[0][1]
                ring_residue_type = Atom_utils._get_residue_type_from_atom_id(ring_atoms[1])
                ring_residue = Atom_utils._get_atom_info_from_index(ring_atoms[0])[1]
                expected_key =   from_atom_key,(ring_residue,ring_residue_type,len(ring_atoms))
                self.assertIn(expected_key, expected_component_keys, `expected_key`)

                shift = ring_subpotential._shift_calculator._calc_sub_component_shift(from_atom_id, ring_id,coefficent)
                self.assertAlmostEqual(expected_ring_shifts[expected_key], shift, places=self.DEFAULT_DECIMAL_PLACES - 2, msg=`expected_key`)
                if abs(expected_ring_shifts[expected_key] - shift) > 0.001:
                    print 'fail', expected_key, expected_ring_shifts[expected_key], shift, Atom_utils._get_residue_type_from_atom_id(from_atom_id)
                    print

                del expected_ring_shifts[expected_key]


        self.remove_zero_valued_keys(expected_ring_shifts)
        self.assertEmpty(expected_ring_shifts)



    def test_force_components_add_up(self):
        forces  =  gb3.gb3_forces
        component_forces =  gb3.gb3_component_forces
        self._do_test_force_componenets_add_up(forces, component_forces)

    def test_force_components_add_up_10_step(self):
        for i,file_name in enumerate(gb3_10_steps.gb3_files):
            print file_name
            forces  =  gb3_10_steps.gb3_forces[i]
            component_forces =  gb3_10_steps.gb3_component_forces[i]['forces']
            self._do_test_force_componenets_add_up(forces, component_forces)


    def _do_test_force_componenets_add_up(self, forces, component_forces):
        gb3_forces_copy = dict(forces)
        summary = self._build_component_force_summary(component_forces)


        for i, elem in enumerate(sorted(summary)):

            self.assertSequenceAlmostEqual(summary[elem], gb3_forces_copy[elem], 2, `elem`)

            del gb3_forces_copy[elem]


        self.remove_almost_zero_force_elems(gb3_forces_copy)

        self.assertEmpty(gb3_forces_copy)

    def _build_component_force_summary(self, component_forces):
        summary = {}
        for elem in component_forces:
            key = self.get_force_atom_key(elem)
            summary_forces = summary.setdefault(key, [0.0] * 3)
            for i, value in enumerate(component_forces[elem]):
                summary_forces[i] += value

        return summary

    @staticmethod
    def filter_dict(the_dict, predicate=lambda k, v: True):
        for k, v in the_dict.iteritems():
            if predicate(k,v):
                yield k, v

    def test_atom_energies(self):
        xcamshift  = self._setup_xcamshift_with_shifts_table(gb3.gb3_zero_shifts)


        expected_energies = gb3.gb3_energies

        self._do_test_atom_energies(xcamshift, expected_energies)

    def test_atom_energies_10_step(self):
        xcamshift  = self._setup_xcamshift_with_shifts_table(gb3.gb3_zero_shifts)
        for i,file_name in enumerate(gb3_10_steps.gb3_files):
            PDBTool("test_data/gb3_10_steps/%s" % file_name).read()
            print file_name
            if i == 0:
                xcamshift.reset()
                self._clear_caches()


            expected_energies = gb3_10_steps.gb3_energies[i]
            shifts =  gb3_10_steps.gb3_shifts[i]

            self._do_test_atom_energies(xcamshift, expected_energies)

    def _do_test_atom_energies(self, xcamshift, expected_energies):
        expected_energies = dict(expected_energies)
        xcamshift._prepare(ROUND_CHANGED, None)
        xcamshift._prepare(TARGET_ATOM_IDS_CHANGED, xcamshift._get_active_target_atom_ids())
        total_component_energies = 0.0

        for key in sorted(expected_energies):
            if key ==  'total':
                continue
            target_atom_index = get_atom_index(key)
            expected_energy = expected_energies[key]

            energy = xcamshift._calc_single_atom_energy(target_atom_index)

            self.assertAlmostEqual(energy, expected_energy, self.DEFAULT_DECIMAL_PLACES - 5, key)





    def test_total_energy(self):
        xcamshift  = self._setup_xcamshift_with_shifts_table(gb3.gb3_zero_shifts)

        expected_energies = gb3.gb3_energies

        self._do_test_total_energy(xcamshift, expected_energies)

    def test_total_energy_10_step(self):
        xcamshift  = self._setup_xcamshift_with_shifts_table(gb3.gb3_zero_shifts)
        for i,file_name in enumerate(gb3_10_steps.gb3_files):
            PDBTool("test_data/gb3_10_steps/%s" % file_name).read()
            print file_name
            if i == 0:
                xcamshift.reset()
                self._clear_caches()


            expected_energies = gb3_10_steps.gb3_energies[i]
            shifts = gb3_10_steps.gb3_shifts[i]

            self._do_test_total_energy(xcamshift, expected_energies)


    def _do_test_total_energy(self, xcamshift, expected_energies):
        expected_total_energy = expected_energies['total']
        total_energy = xcamshift.calcEnergy()
        self.assertAlmostEqual(total_energy, expected_total_energy, self.DEFAULT_DECIMAL_PLACES - 5)


    def test_batch_shift_calculator(self):
        xcamshift  = self._setup_xcamshift_with_shifts_table(gb3.gb3_zero_shifts)
        xcamshift._prepare(ROUND_CHANGED, None)
        active_target_atoms_ids = xcamshift._get_active_target_atom_ids()

        result  =  CDSVectorFloat(len(active_target_atoms_ids))

        xcamshift.calc_shifts(active_target_atoms_ids, result)

        _check_shift_results(active_target_atoms_ids, result,gb3.gb3_shifts)


    def test_shift_differences(self):
        xcamshift  = self._setup_xcamshift_with_shifts_table(gb3.gb3_zero_shifts)
        xcamshift._prepare(ROUND_CHANGED, None)


        for key in sorted(gb3.gb3_shift_diffs):
            expected_shift_diff  = gb3.gb3_shift_diffs[key]

            target_atom_index = get_atom_index(key)

            shift_diff = xcamshift._calc_single_atom_shift(target_atom_index)

            self.assertAlmostEqual(shift_diff, - expected_shift_diff, self.DEFAULT_DECIMAL_PLACES-3,key)

    def set_cache_shifts(self,xcamshift, shift_list):
        shift_cache =  xcamshift._shift_cache

        for key in shift_list:
            index = get_atom_index(key)
            shift_cache[index] = shift_list[key]


    def add_forces_to_array(self, target_array, index, force_component):
        expected_forces_set = target_array[index]
        if expected_forces_set == None:
            target_array[index] = [0.0] * 3
            expected_forces_set = target_array[index]
        for i, new_force in enumerate(force_component):
            expected_forces_set[i] += new_force

        return target_array


    def filter_force_components_by_potential(self, force_components, potential_name):
        return dict(self.filter_dict(force_components, lambda k, v:k[-1] == potential_name))

    def _make_expected_array(self):
        num_atoms = len(AtomSel('(all)').indices())
        result = [None] * num_atoms
        return result

    def make_forces_from_map(self,forces_map):
        expected_total_result_forces = self._make_expected_array()
        for key in forces_map:
            index = get_atom_index(key)
            expected_total_result_forces[index] = forces_map[key]

        return expected_total_result_forces

    def test_total_forces_and_energy_frozen(self):
        xcamshift =  self.make_xcamshift(gb3.gb3_zero_shifts)
        xcamshift.setup()

        total_result_forces = self.make_results_array()
        energy = xcamshift.calcEnergyAndDerivs(total_result_forces)
        result_forces =  self.deriv_list_as_array(total_result_forces, currentSimulation())
        expected_total_result_forces = self.make_forces_from_map(gb3.gb3_forces)

        self.assertListVec3AlmostEqual(result_forces, expected_total_result_forces, self.DEFAULT_DECIMAL_PLACES - 3)
#        for elem in zip(total_result_forces,expected_total_result_forces):
#            print elem

        self.assertAlmostEqual(energy, gb3.gb3_energies['total'],self.DEFAULT_DECIMAL_PLACES-5)


    def test_total_forces_and_energy(self):
        xcamshift =  self.make_xcamshift(gb3.gb3_zero_shifts)

        forces = gb3.gb3_forces
        energy = gb3.gb3_energies['total']

        xcamshift._prepare(TARGET_ATOM_IDS_CHANGED, None)

        self._do_test_total_forces_and_energy(xcamshift, energy, forces)

    def test_total_forces_and_energy_10_step(self):
        xcamshift =  self.make_xcamshift(gb3.gb3_zero_shifts)

        for i,file_name in enumerate(gb3_10_steps.gb3_files):
                PDBTool("test_data/gb3_10_steps/%s" % file_name).read()
                print 'coord file',file_name
                if i == 0:
                    xcamshift.reset()
                    self._clear_caches()

                energy = gb3_10_steps.gb3_energies[i]['total']
                force_components =  gb3_10_steps.gb3_forces[i]

                self._do_test_total_forces_and_energy(xcamshift, energy, force_components)


    def deriv_list_as_array(self, total_result_forces_deriv_list, simulation):
        derivs_vec3 = total_result_forces_deriv_list.get(simulation)
        #todo: report to charles  print derivs_vec3 segfaults
        result = [None]*len(derivs_vec3)
        for i,vec in enumerate(derivs_vec3):
            if vec != Vec3(0.0,0.0,0.0):
                result[i]=Vec3(vec)
        return result


    def _do_test_total_forces_and_energy(self, xcamshift, expected_energy,  force_components):

        total_result_forces_deriv_list = DerivList()
        total_result_forces_deriv_list.init(currentSimulation())

        energy = xcamshift.calcEnergyAndDerivList(total_result_forces_deriv_list)

        expected_total_result_forces = self.make_forces_from_map(force_components)
        total_result_forces = self.deriv_list_as_array(total_result_forces_deriv_list, currentSimulation())

        self.assertListVec3AlmostEqual(total_result_forces, expected_total_result_forces, self.DEFAULT_DECIMAL_PLACES - 3)

        self.assertAlmostEqual(energy, expected_energy, self.DEFAULT_DECIMAL_PLACES - 5)




    def make_xcamshift(self, shifts):
        xcamshift = self._setup_xcamshift_with_shifts_table(shifts)
        xcamshift.setup()

        #non_bonded = xcamshift.get_named_sub_potential(NON_BONDED)
        #non_bonded.update_non_bonded_list()


        #ring_potential = xcamshift.get_named_sub_potential(RING)
        #ring_potential._build_ring_data_cache()

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
        xcamshift =  self._setup_xcamshift_with_shifts_table(gb3.gb3_zero_shifts)


        expected_force_components = gb3.gb3_component_forces
        shifts =  gb3.gb3_shifts

        self._do_test_force_components(xcamshift, shifts, expected_force_components)

    def _do_test_force_components(self, xcamshift, shifts, expected_force_components):


        def get_target_atom_keys_for_potential(raw_backbone_force_components):
            target_atom_keys = set()
            for raw_backbone_force_component_key in raw_backbone_force_components:
                target_atom_key = self.get_target_atom_key(raw_backbone_force_component_key)
                target_atom_keys.add(target_atom_key)

            return target_atom_keys

        def make_result_array_forces(self):
#        TODO: use segment manager
            num_atoms = len(AtomSel('(all)').indices())
            result = [None] * num_atoms
            return result

        def build_expected_forces_for_potential(target_atom_key, raw_backbone_force_components):
            expected_forces = [None] * len(AtomSel('(all)').indices())
            for raw_backbone_force_component_key in raw_backbone_force_components:
                if target_atom_key == self.get_target_atom_key(raw_backbone_force_component_key):
                    force_atom_key = self.get_force_atom_key(raw_backbone_force_component_key)
                    force_atom_index = get_atom_index(force_atom_key)
                    expected_force_component = raw_backbone_force_components[raw_backbone_force_component_key]
                    self.add_forces_to_array(expected_forces, force_atom_index, expected_force_component)

            return expected_forces


#       TODO: xcamshift should be created anew  for each passage through the loop but the non bonded list is too slow
        xcamshift._prepare(ROUND_CHANGED, None)
        xcamshift.calc_shifts(xcamshift._get_active_target_atom_ids())
        xcamshift.update_force_factor_calculator()

        for i, potential_name in enumerate(CAMSHIFT_SUB_POTENTIALS):
            out_string = '%s %i/%i:  ' % (potential_name,i+1, len(common_constants.CAMSHIFT_SUB_POTENTIALS))
            out_string = '%-20s' % out_string
            print out_string,
            sys.stdout.flush()

            raw_backbone_force_components = self.filter_force_components_by_potential(expected_force_components, potential_name)

            target_atom_keys =  get_target_atom_keys_for_potential(raw_backbone_force_components)

            potentials_list = [xcamshift.get_named_sub_potential(potential_name)]

            threshold=10.0
            for j,target_atom_key in enumerate(reversed(sorted(target_atom_keys))):

                percent_done = float(j) / float(len(target_atom_keys)) *100.0

                if percent_done > threshold:
                    print '%i%% ' % percent_done,
                    sys.stdout.flush()
                    threshold += 10.0

                expected_forces = build_expected_forces_for_potential(target_atom_key, raw_backbone_force_components)

                target_atom_id = get_atom_index(target_atom_key)
                out_array = self.make_out_array()

                atom_ids = [target_atom_id]
                xcamshift._prepare(TARGET_ATOM_IDS_CHANGED, atom_ids)
                xcamshift._calc_single_atom_force_set_with_potentials(target_atom_id, out_array, potentials_list)
                result_forces = out_array.add_forces_to_result()

                self.assertListVec3AlmostEqual(result_forces, expected_forces, self.DEFAULT_DECIMAL_PLACES-3, msg='%s -  %s' % (potential_name,target_atom_key))

            print "100%"
        print

    def test_force_components_10_step(self):
        xcamshift =  self._setup_xcamshift_with_shifts_table(gb3.gb3_zero_shifts)
        for i,file in enumerate(gb3_10_steps.gb3_files):
            PDBTool("test_data/gb3_10_steps/%s" % file).read()
            print '%i/%i' % (i+1,len(gb3_10_steps.gb3_files)), file
            if i == 0:
                xcamshift.reset()
                self._clear_caches()
#            print gb3_10_steps.gb3_component_forces[i]
            expected_force_components = gb3_10_steps.gb3_component_forces[i]['forces']
            shifts =  gb3_10_steps.gb3_shifts[i]

            self._do_test_force_components(xcamshift, shifts, expected_force_components)

                # to list differences as well as check them uncomment this!
                # self.print_force_component_diffs(potential_name,target_atom_key, expected_forces, result_forces)

#    def compare_shift_diffs_and_expected(self):
#        for elem in sorted(gb3.gb3_shifts):
#            self.assertAlmostEqual(gb3.gb3_shifts[elem], -gb3.gb3_shift_diffs[elem], self.DEFAULT_DECIMAL_PLACES-3,  elem)

    def test_component_chemical_shifts_10_step(self):
        xcamshift  = self._get_xcamshift()
        for i,file_name in enumerate(gb3_10_steps.gb3_files):
            PDBTool("test_data/gb3_10_steps/%s" % file_name).read()
            print file_name
            if i == 0:
                xcamshift.reset()
                self._clear_caches()

            component_shifts = gb3_10_steps.gb3_subpotential_shifts[i]
            self._do_test_component_shifts(xcamshift, component_shifts)

    def test_total_chemical_shifts(self):
            xcamshift  = self._get_xcamshift_no_hbond()

            component_shifts = gb3.gb3_subpotential_shifts
            self._do_test_shifts(xcamshift, component_shifts)

    def test_total_chemical_shifts_10_step(self):
            xcamshift  = self._get_xcamshift_no_hbond()
            for i,file_name in enumerate(gb3_10_steps.gb3_files):
                PDBTool("test_data/gb3_10_steps/%s" % file_name).read()
                print 'coord file',file_name
                if i == 0:
                    xcamshift.reset()
                    self._clear_caches()
                component_shifts = gb3_10_steps.gb3_subpotential_shifts[i]
                self._do_test_shifts(xcamshift, component_shifts)


    def _do_test_shifts(self, xcamshift, expected_component_shifts):
        xcamshift._prepare(ROUND_CHANGED, None)
        component_shifts_keys = expected_component_shifts.keys()
        component_shifts_keys.sort()

        target_atom_selections, calculated_shifts = xcamshift.calc_shifts()

        self._check_shifts(target_atom_selections, expected_component_shifts, calculated_shifts)

    def _sum_sub_potential_shifts(self,component_shifts):
        result = {}
        for i, key in enumerate(component_shifts.keys()):
            segment, residue_number, atom, sub_potential = key
            new_key = (segment, residue_number, atom)
            value = result.setdefault(new_key,0.0)
            value+=component_shifts[key]
            result[new_key] = value
        return result

    def _check_shifts(self, target_atom_selections, expected_component_shifts, calculated_shifts):

        expected_shifts = self._sum_sub_potential_shifts(expected_component_shifts)
        shift_diffs = []
        for target_atom_selection, shift in zip(target_atom_selections, calculated_shifts):
            target_atom_selection = util_for_testing.translate_atom_key(target_atom_selection)
            shift_diffs.append(expected_shifts[target_atom_selection]-shift)
            self.assertAlmostEqual(expected_shifts[target_atom_selection], shift, self.DEFAULT_DECIMAL_PLACES - 2, target_atom_selection)

    def test_atom_id_stability(self):
        test = None
        for i,file in enumerate(gb3_10_steps.gb3_files):
            PDBTool("test_data/gb3_10_steps/%s" % file).read()
            print file
            latest = {}
            for atom_sel in AtomSel('all'):
                index  = atom_sel.index()
                info  = atom_sel.segmentName(), atom_sel.residueNum(), atom_sel.atomName()
                latest[index] = info

            if test != None:
                self.assertDictEqual(latest, test)
            test=latest


    def _create_naive_and_fast_non_bonded_lists(self, nb_potential):
        simulation = self.get_single_member_ensemble_simulation()

        new_nb_calculator = New_fast_non_bonded_calculator(simulation, 1)
        old_nb_calculator = Fast_non_bonded_calculator(simulation, 1)


        components = nb_potential._get_components()

        atom_list_1 = components['ATOM']
        atom_list_2 = components['NBRM']

        non_bonded_lists = Non_bonded_interaction_list()
        old_non_bonded_lists = Non_bonded_interaction_list()

        active_components = None

        new_calc = partial(new_nb_calculator,atom_list_1, atom_list_2, non_bonded_lists, active_components)
        old_calc = partial(old_nb_calculator,atom_list_1, atom_list_2, old_non_bonded_lists, active_components)

        return new_calc, old_calc, non_bonded_lists, old_non_bonded_lists

    def test_new_fast_non_bonded_list(self):

        nb_potential = self._get_xcamshift().get_named_sub_potential(NON_BONDED)

        new_calc, old_calc, non_bonded_lists, old_non_bonded_lists = self._create_naive_and_fast_non_bonded_lists(nb_potential)

        new_calc()
        old_calc()

        remote_components = nb_potential._get_component_list('NBRM')

        seen_dists = set()
        expected_dists = set()
        for elem in non_bonded_lists.dump():
            atom_id_1 = elem[0]
            atom_id_2 = remote_components[elem[2]][0]

            key = atom_id_1,atom_id_2

            seen_dists.add(key)

        for elem in old_non_bonded_lists.dump():
            atom_id_1 = elem[0]
            atom_id_2 = remote_components[elem[2]][0]

            key = atom_id_1,atom_id_2

            expected_dists.add(key)

        self.assertEqual(seen_dists,expected_dists)

    def test_new_fast_non_bonded_list_timing(self):

        nb_potential = self._get_xcamshift().get_named_sub_potential(NON_BONDED)

        new_calc, old_calc, non_bonded_lists, old_non_bonded_lists = self._create_naive_and_fast_non_bonded_lists(nb_potential)

        new_calc()
        non_bonded_lists.clear()
        start =  now()
        new_calc()
        end = now()

        new_time = (end - start).seconds()

        old_calc()
        old_non_bonded_lists.clear()
        start =  now()
        old_calc()
        end = now()

        old_time = (end - start).seconds()

        self.assertTrue(new_time < old_time)

        print 'new  %4.3f ms / cycle' % (new_time*1000.0)
        print 'new  %4.3f ms / cycle' % (old_time*1000.0)


if __name__ == "__main__":
#     TODO: add a way to run the complete test suite
      unittest2.main(module='test.test_xcamshift_gb3',failfast=True, defaultTest='TestXcamshiftGB3.test_new_fast_non_bonded_list_timing')
#     unittest2.main(module='test.test_xcamshift_gb3',defaultTest='TestXcamshiftGB3.test_total_forces_and_energy_10_step', exit=False)
#     unittest2.main(module='test.test_xcamshift_gb3',defaultTest='TestXcamshiftGB3.test_force_components')
#     unittest2.main(module='test.test_xcamshift_gb3',defaultTest='TestXcamshiftGB3.test_shift_averaging_two_structures')
#      unittest2.main(module='test.test_xcamshift_gb3',defaultTest='TestXcamshiftGB3.test_shift_averaging_identical_structures')