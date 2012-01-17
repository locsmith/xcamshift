'''
Created on 31 Dec 2011

@author: garyt
'''
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from xcamshift import RandomCoilShifts, Distance_potential, Extra_potential,\
    Dihedral_potential
from atomSel import AtomSel
from test.xdists import xdists_ala_3
from test.dihedrals import dihedrals_ala_3


#class testSegmentManager(object):
class Test(unittest2.TestCase):

    def assertSequenceAlmostEqual(self,result,expected, delta = 1e-7):
        len_result = len(result)
        len_expected = len(expected)
        if len_result != len_expected:
            raise AssertionError("the two lists are of different length %i and %i" % (len_result,len_expected))
        
        for i,(elem_1, elem_2) in enumerate(zip(result,expected)):
            
            diff = abs(elem_1 - elem_2)
            if diff > delta:
                template = "lists differ at item %i: %s - %s > %s"
                message = template % (i, `elem_1`,`elem_2`,delta)
                raise AssertionError(message)
            
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
        self.assertSequenceAlmostEqual(expected, result)
        
    def testDistancePotential(self):
        distance_potential =  Distance_potential()
        result=self.make_result_array()
        distance_potential.set_shifts(result)
        
        expected = self.make_result_array()
        expected[15] =   0.82775
        expected[14] =  41.5303
        expected[13] =  -5.56269
        expected[12] = -12.185
        expected[20] = -12.584
        expected[16] =  -3.6896
        

        self.assertSequenceAlmostEqual(result, expected,delta=0.001)
        
        
    @staticmethod
    def assertElemeInSet( elem, xdist_set):
        if not elem in xdist_set:
            raise AssertionError("%s not found in set" % `elem`)


    def testExtraPotentialComponentsCorrect(self):
        extra_potential = Extra_potential()
        
        xdists_ala_3_copy = dict(xdists_ala_3)
        for extra_elem in extra_potential.dump():
            elem_key = extra_elem[:-1]
            self.assertElemeInSet(elem_key, xdists_ala_3)
            del xdists_ala_3_copy[elem_key]
        
        self.assertEqual(0, len(xdists_ala_3_copy))


    def testExtraPotentialCoefficientsCorrect(self):
        extra_potential = Extra_potential()
        for extra_elem in extra_potential.dump():
            elem_key = extra_elem[:-1]
            coefficient = xdists_ala_3[elem_key][0]

            # TODO add extra places
            self.assertAlmostEqual(coefficient, extra_elem[-1], places=5)
        

    def testExtraPotentialComponentShiftsCorrect(self):
        extra_potential = Extra_potential()
        
        result=self.make_result_array()
        extra_potential.set_shifts(result)
        
        for i,extra_elem in enumerate(extra_potential.dump()):
            
            elem_key = extra_elem[:-1]
            shift = extra_potential._calc_single_shift(i)
            
            expected_shift = xdists_ala_3[elem_key][2]
            
            self.assertAlmostEqual(shift, expected_shift,places=4)

    def testDihdedralPotentialComponentsCorrect(self):
        dihedral_portential = Dihedral_potential()
        
        
        dihedrals_ala_3_copy = dict(dihedrals_ala_3)
        for dihedral_elem_key in dihedral_portential.dump():
            self.assertElemeInSet(dihedral_elem_key[0], dihedrals_ala_3_copy)
            del dihedrals_ala_3_copy[dihedral_elem_key[0]]

        self.assertEqual(0, len(dihedrals_ala_3_copy))
        

if __name__ == "__main__":
    unittest2.main()
