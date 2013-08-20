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
Created on 24 Jan 2012

@author: garyt
'''
import unittest2
from test import ala_3
from cython.observed_chemical_shifts import Observed_shift_table
from protocol import initStruct
from pdbTool import PDBTool
from test.ala_3 import ala_3_test_shifts_harmonic
from utils import Atom_utils
from cython.fast_segment_manager import Segment_Manager


class TestObservedShiftTable(unittest2.TestCase):
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

    def testObservedShiftTableKeysLong(self):
        test_shifts = ala_3.ala_3_test_shifts_harmonic
        shift_table = Observed_shift_table(test_shifts)
        
        test_keys = [item[0] for item in shift_table.dump_observed_shifts()]
        
        expected_keys = set(self.EXPECTED_ALA_3_KEYS_LONG)
        self._test_shift_keys(test_keys,expected_keys)
        
    def testObservedShiftTableValues(self):
        shift_table = self.create_ala_3_shift_table()
        
        for elem in shift_table.dump_observed_shifts():
            short_atom_key = elem[0][1:]
            value = elem[1]
            
            expected = ala_3_test_shifts_harmonic[short_atom_key]
            self.assertAlmostEqual(expected, value,places=4)
            
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

    def test_native_shifts(self):
        shift_table = self.create_ala_3_shift_table()
        native_shifts  = shift_table.py_get_native_shifts([12,13,14,15,16,20])
        keys =   (2, "N"), (2, "HN"), (2, "CA"), (2, "HA"), (2, "CB"), (2, "C")
        expected_shifts = [ala_3.ala_3_test_shifts_harmonic[key] for key in keys]
        self.assertSequenceAlmostEqual(native_shifts, expected_shifts, delta = 0.1**4)
        
    def test_native_shifts_not_cached(self):
        shift_table = self.create_ala_3_shift_table()
        
        atom_ids_1 = [12,13,14,15,16,20]
        atom_ids_2 = [12,13,14]
        
        native_shifts_1  = shift_table.py_get_native_shifts(atom_ids_1)
        native_shifts_2  = shift_table.py_get_native_shifts(atom_ids_2)
        
        keys_1 =   (2, "N"), (2, "HN"), (2, "CA"), (2, "HA"), (2, "CB"), (2, "C")
        keys_2 =  (2, "N"), (2, "HN"), (2, "CA")
        
        expected_shifts_1 = [ala_3.ala_3_test_shifts_harmonic[key] for key in keys_1]
        expected_shifts_2 = [ala_3.ala_3_test_shifts_harmonic[key] for key in keys_2]
        
        self.assertSequenceAlmostEqual(native_shifts_1, expected_shifts_1, delta = 0.1**4)
        self.assertSequenceAlmostEqual(native_shifts_2, expected_shifts_2, delta = 0.1**4)

if __name__ == "__main__":
    unittest2.main()
