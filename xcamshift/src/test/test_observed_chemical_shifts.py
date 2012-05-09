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


class TestObservedShiftTable(unittest2.TestCase):
    
    EXPECTED_ALA_3_KEYS_SHORT = ala_3.ala_3_test_shifts_harmonic.keys()
    EXPECTED_ALA_3_KEYS_LONG =  [('',key[0],key[1]) for key in EXPECTED_ALA_3_KEYS_SHORT]
    
    def setUp(self):
        initStruct("test_data/3_ala/3ala.psf")
        PDBTool("test_data/3_ala/3ala.pdb").read()
        
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
            self.assertAlmostEqual(expected, value)
            
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
            

if __name__ == "__main__":
    unittest2.main()