'''
Created on 31 Dec 2011

@author: garyt
'''
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from xcamshift import RandomCoilShifts, Distance_potential
from atomSel import AtomSel


#class testSegmentManager(object):
class Test(unittest2.TestCase):


    def setUp(self):
        initStruct("test_data/3_ala/3ala.psf")
        PDBTool("test_data/3_ala/3ala.pdb").read()



    def make_result_array(self):
        num_atoms = len(AtomSel('(all)').indices())
        result = [0.0] * num_atoms
        return result

    def testRandomCoilPotential(self):
        random_coil_potential = RandomCoilShifts()
        result = self.make_result_array() 
        expected  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                     123.8, 8.2400000000000002, 52.266440000000003, 4.4328469999999998, 
                     19.0, 0.0, 0.0, 0.0, 177.09999999999999, 0.0, 0.0, 0.0, 0.0, 0.0, 
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        random_coil_potential.set_shifts(result)
        self.assertSequenceEqual(expected, result)
        
    def testDistancePotential(self):
        distance_potential =  Distance_potential()
        result=self.make_result_array()
        distance_potential.set_shifts(result)
        
        print result
        


if __name__ == "__main__":
    unittest2.main()
