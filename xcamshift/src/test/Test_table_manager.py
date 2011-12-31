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
        
        coefficient = table.get_distance_coeeficent('HA', -1, 'N')
        self.assertAlmostEqual(coefficient,0.10495380090732243)
        
        coefficient = table.get_distance_coeeficent('HA', -1, 'CA')
        self.assertAlmostEqual(coefficient,None)
        
    
    def testOffsets(self):
        table = self.table_manager.get_BB_Distance_Table('glu') 
        offsets = table.get_offsets()
        
        expected = set((-1,0,1))
        self.assertItemsEqual(expected,offsets)
        
    def testFromAtomList(self):
        table = self.table_manager.get_BB_Distance_Table('glu') 
        from_atom_list = table.get_from_atoms()
        expected = set(('N','HN','CA','HA','C','CB'))
        self.assertItemsEqual(expected, from_atom_list)

    def testToAtomList(self):
        table = self.table_manager.get_BB_Distance_Table('glu') 
        to_atom_list = table.get_to_atoms()
        expected = set(('N','HN','CA','HA','C','O'))
        self.assertItemsEqual(expected, to_atom_list)
    
    def testLoadRandomCoil(self):
        table = self.table_manager.get_random_coil_table()
        
        expected = 8.24
        shift = table.get_random_coil_shift('ALA', 'HN')
        
        self.assertAlmostEqual(expected, shift)
        
    def testGetDefaultTableManager(self):
        table_manager = Table_manager.get_default_table_manager()
        
        self.assertIsInstance(table_manager, Table_manager)
    
if __name__ == "__main__":
    unittest2.main()