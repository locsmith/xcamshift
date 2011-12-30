'''
Created on 30 Dec 2011

@author: garyt
'''
import unittest2
from Table_manager import Table_manager


class Test(unittest2.TestCase):

   

    def setUp(self):
        self.table_manager = Table_manager()
        self.table_manager.add_search_path('../../data')

    def testLoadTable(self):
        table = self.table_manager.get_BB_Distance_Table('ala')
        self.assertTrue(table != None)

    def testLoadTableValues(self):
        table = self.table_manager.get_BB_Distance_Table('glu') 
        
        exponent = table.get_exponent()
        self.assertAlmostEqual(exponent,1.0)
        
        coefficient = table.get_distance_coeeficent('N', -1, 'HA')
        self.assertAlmostEqual(coefficient,0.10495380090732243)
        
    
    def testOffsets(self):
        table = self.table_manager.get_BB_Distance_Table('glu') 
        offsets = table.get_offsets()
        
        expected = set((-1,0,1))
        self.assertItemsEqual(expected,offsets)
        
    def testFromAtomList(self):
        table = self.table_manager.get_BB_Distance_Table('glu') 
        from_atom_list = table.get_from_atoms()
        expected = set(('N','HN','CA','HA','C','O'))
        self.assertItemsEqual(expected, from_atom_list)

    def testToAtomList(self):
        table = self.table_manager.get_BB_Distance_Table('glu') 
        to_atom_list = table.get_to_atoms()
        expected = set(('N','HN','CA','HA','C','CB'))
        self.assertItemsEqual(expected, to_atom_list)
        
if __name__ == "__main__":
    unittest2.main()