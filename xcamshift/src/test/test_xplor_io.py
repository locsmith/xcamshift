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

#TODO this acts as a stream and a context manager and an Xplor_reader
class Test_xplor_reader(Xplor_reader):
    def __init__(self,strings):
        Xplor_reader.__init__(self,'no file!')
        self._strings = strings
        
    def get_file(self):
        return self
    
    def __exit__(self,type, value, traceback):
        pass
    
    def __enter__(self):
        return self
    
    def __len__(self):
        return len(self._strings)
    
    def __getitem__(self, key):
        return self._strings[key]
            
class Test(unittest2.TestCase):


    def setUp(self):
        self.file_name = 'test_data/io/cs.dat'
        PDBTool("test_data/io/2L1A_single.pdb").read()
        initStruct("test_data/io/2L1A.psf")
         
       


    def test_alvins_data(self):
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
             
    def test_bad_atom_selection(self):

            
        with self.assertRaises(Exception) as exception:       
            Test_xplor_reader(("assign (resid 20 and (name HA or name C)) 127.0 0.1",)).read()
            
        with self.assertRaises(Exception) as exception:       
            Test_xplor_reader(("assign (resid 20 and name HK) 127.0 0.1",)).read()
            
    def test_unweighted_shift(self):
        results = Test_xplor_reader(("assign (resid 20 and name HA) 127.0",)).read()
        self.assertEqual(len(results),1)
        
        self.assertAlmostEqual(results[0].error,0.1)
        self.assertAlmostEqual(results[0].shift,127.0)
        self.assertAlmostEqual(results[0].weight,1.0)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_alvins_data']
    unittest.main()