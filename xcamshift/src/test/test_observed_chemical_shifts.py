'''
Created on 24 Jan 2012

@author: garyt
'''
import unittest2
from test import ala_3
from observed_chemical_shifts import Observed_shift_table
import sys
from protocol import initStruct
from pdbTool import PDBTool
from test.ala_3 import ala_3_test_shifts


class TestObservedShiftTable(unittest2.TestCase):
    
    EXPECTED_ALA_3_KEYS = ala_3.ala_3_test_shifts.keys()
    
    def setUp(self):
        initStruct("test_data/3_ala/3ala.psf")
        PDBTool("test_data/3_ala/3ala.pdb").read()
        
    def assertEmpty(self, expected_keys, msg=""):
        return self.assertEqual(len(expected_keys), 0, msg)
    
    def create_ala_3_shift_table(self):
        test_shifts = ala_3.ala_3_test_shifts
        shift_table = Observed_shift_table(test_shifts)
        return shift_table

    def testObservedShiftTableKeys(self):
        test_shifts = ala_3.ala_3_test_shifts
        
        shift_table = Observed_shift_table(test_shifts)
        
        
        expected_keys = set(self.EXPECTED_ALA_3_KEYS)
        for key in shift_table.dump_observed_shifts():
            short_atom_key = key[0][1:]
            self.assertIn(short_atom_key, expected_keys)
            expected_keys.remove(short_atom_key)
        self.assertEmpty(expected_keys)
        
    def testObservedShiftTableValues(self):
        shift_table = self.create_ala_3_shift_table()
        
        for elem in shift_table.dump_observed_shifts():
            short_atom_key = elem[0][1:]
            value = elem[1]
            
            expected = ala_3_test_shifts[short_atom_key]
            self.assertAlmostEqual(expected, value)

if __name__ == "__main__":
    unittest2.main()