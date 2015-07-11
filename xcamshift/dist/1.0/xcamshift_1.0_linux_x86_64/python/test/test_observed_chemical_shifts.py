#-------------------------------------------------------------------------------
# Copyright (c) 2013-2015 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
#
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 24 Jan 2012

@author: garyt
'''
import unittest2
from test import ala_3
from observed_chemical_shifts import Observed_shift_table
from protocol import initStruct
from pdbTool import PDBTool
from test.ala_3 import ala_3_test_shifts_harmonic
from utils import Atom_utils
from cython.fast_segment_manager import Segment_Manager


class TestObservedShiftTable(unittest2.TestCase):

    EXPECTED_ALA_3_KEYS_SHORT = ala_3.ala_3_test_shifts_harmonic.keys()
    EXPECTED_ALA_3_KEYS_LONG =  [('',key[0],key[1]) for key in EXPECTED_ALA_3_KEYS_SHORT]

    def setUp(self):
        initStruct("test_data/3_ala/3ala.psf")
        PDBTool("test_data/3_ala/3ala.pdb").read()
        Atom_utils.clear_cache()
        Segment_Manager.reset_segment_manager()

    def assertEmpty(self, expected_keys, msg=""):
        return self.assertEqual(len(expected_keys), 0, msg)

    def create_ala_3_shift_table(self):
        test_shifts = ala_3.ala_3_test_shifts_harmonic
        shift_table = Observed_shift_table(test_shifts)
        return shift_table


    def _test_shift_keys(self, test_keys, expected_keys):

        for test_key in test_keys:
            self.assertIn(test_key, expected_keys)
            expected_keys.remove(test_key)

        self.assertEmpty(expected_keys)

    def testObservedShiftTableKeysShort(self):
        test_shifts = ala_3.ala_3_test_shifts_harmonic
        shift_table = Observed_shift_table(test_shifts)

        expected_keys = set(self.EXPECTED_ALA_3_KEYS_SHORT)

        test_keys = []
        for key in shift_table.dump_observed_shifts():
            short_atom_key = key[0][1:]
            test_keys.append(short_atom_key)


        self._test_shift_keys(test_keys,expected_keys)

#TODO add test for keys which are the wrong length


    def _get_keys(self, shift_table):
        return [item[0] for item in shift_table.dump_observed_shifts()]

    def testObservedShiftTableKeysLong(self):
        test_shifts = ala_3.ala_3_test_shifts_harmonic
        shift_table = Observed_shift_table(test_shifts)

        test_keys = self._get_keys(shift_table)

        expected_keys = set(self.EXPECTED_ALA_3_KEYS_LONG)
        self._test_shift_keys(test_keys,expected_keys)


    def _check_shifts(self, shift_table, test_shifts):
        for elem in shift_table.dump_observed_shifts():
            short_atom_key = elem[0][1:]
            value = elem[1]
            expected = test_shifts[short_atom_key]
            self.assertAlmostEqual(expected, value)

    def testObservedShiftTableValues(self):
        shift_table = self.create_ala_3_shift_table()

        self._check_shifts(shift_table,ala_3_test_shifts_harmonic)

    def testIndexToAtomId(self):
        shift_table = self.create_ala_3_shift_table()
        atom_id_indices = dict(shift_table.get_indices_for_atom_id())

        for i,elem in enumerate(shift_table.dump_observed_shifts()):
            atom_ids = Atom_utils.find_atom_ids(*elem[0])
            self.assertEqual(len(atom_ids), 1)
            atom_id = atom_ids[0]
            self.assertEqual(atom_id_indices[atom_id],i)
            del atom_id_indices[atom_id]
        self.assertEmpty(atom_id_indices)
#            short_atom_key = elem[0][1:]


    def test_delete_shift(self):
        shift_table = self.create_ala_3_shift_table()

        del shift_table[12]
        self.assertEqual(len(shift_table),5)
        self.assertTrue(('', 2, 'N') not in self._get_keys(shift_table))

    def _check_native_shifts(self,atom_ids,shift_array,shift_map):
        for atom_id,shift in zip(atom_ids,shift_array):
            atom_key = tuple(Atom_utils._get_atom_info_from_index(atom_id)[1:])
            self.assertTrue(atom_key in shift_map)
            self.assertAlmostEqual(shift_map[atom_key],shift,places=3)

    def test_add_shifts(self):
        shift_table = self.create_ala_3_shift_table()
        shift_table_2 = Observed_shift_table({(2, "C") : 179.7749})
        target_atom_ids = 12, 13, 14, 15, 16, 20

        self._check_shifts(shift_table,ala_3_test_shifts_harmonic)
        self._check_native_shifts(target_atom_ids, shift_table.get_native_shifts(target_atom_ids),ala_3.ala_3_test_shifts_harmonic)

        shift_table.add_shifts(shift_table_2)

        altered_shifts = dict(ala_3_test_shifts_harmonic)
        altered_shifts[(2, "C")] =  179.7749
        self._check_shifts(shift_table,altered_shifts)
        self._check_native_shifts(target_atom_ids, shift_table.get_native_shifts(target_atom_ids),altered_shifts)

        del shift_table[12]
        del altered_shifts[(2, "N")]
        altered_target_atom_ids =  target_atom_ids[1:]
        self._check_shifts(shift_table,altered_shifts)
        self._check_native_shifts(altered_target_atom_ids, shift_table.get_native_shifts(altered_target_atom_ids),altered_shifts)

        altered_shifts[(2, "N")] =  120
        shift_table_3 =  Observed_shift_table({(2, "N") : 120.00})
        shift_table.add_shifts(shift_table_3)
        self._check_shifts(shift_table,altered_shifts)
        self._check_native_shifts(target_atom_ids, shift_table.get_native_shifts(target_atom_ids),altered_shifts)



#TODO this doen't have a test of native access...

if __name__ == "__main__":
    unittest2.main()
