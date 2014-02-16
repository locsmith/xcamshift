'''
Created on 11 Feb 2014

@author: garyt
'''
import unittest
from io.xplor_reader import Xplor_reader
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from test.alvin_2L1A_test_shifts import alvin_shift_format
from utils import Atom_utils

class Test(unittest2.TestCase):


    def setUp(self):
        self.file_name = 'test_data/io/cs.dat'
        PDBTool("test_data/io/2L1A_single.pdb").read()
        initStruct("test_data/io/2L1A.psf")
         
       


    def testName(self):
        data = Xplor_reader(self.file_name).read()
        self.assertEqual(len(data),601)
        for elem in data:
            TEST_DATA_SHIFT = 0
            TEST_DATA_ERROR = 1

            key = Atom_utils._get_atom_info_from_index(elem[0])
            self.assertTrue(key in alvin_shift_format, key)
            self.assertAlmostEqual(alvin_shift_format[key][TEST_DATA_SHIFT],elem[1])
            self.assertAlmostEqual(alvin_shift_format[key][TEST_DATA_ERROR],elem[2])
            self.assertAlmostEqual(1.0,elem[3])
            self.assertEqual('DEFAULT', elem[4])
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()