'''
Created on 11 Feb 2014

@author: garyt
'''
import unittest
from io.xplor_reader import Xplor_reader
from protocol import initStruct
from pdbTool import PDBTool
import unittest2

class Test(unittest2.TestCase):


    def setUp(self):
        self.file_name = 'test_data/io/cs.dat'
        PDBTool("test_data/io/2L1A_single.pdb").read()
        initStruct("test_data/io/2L1A.psf")
         
       


    def testName(self):
        data = Xplor_reader(self.file_name).read()
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()