'''
Created on 30 Dec 2011

@author: garyt
'''
import unittest2
from table_manager import Table_manager
from keys import Atom_key, Dihedral_key



class Test(unittest2.TestCase):
    raw_dihderal_keys = ((("C", -1), ("N", 0), ("CA", 0), ("C",  0)), 
                         (("C",  0), ("N", 0), ("CA", 0), ("C",  1)), 
                         (("N", 0), ("CA", 0), ("CB", 0), ("CG", 0)))

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
        table = self.table_manager.get_random_coil_table('ALA')
        
        expected = 8.24
        shift = table.get_random_coil_shift( 0, 'ALA','HN')
        self.assertAlmostEqual(expected, shift)
        
    def testGetDefaultTableManager(self):
        table_manager = Table_manager.get_default_table_manager()
        
        self.assertIsInstance(table_manager, Table_manager)
        
    def testLoadExtra(self):

        table = self.table_manager.get_extra_table('ALA')
        
        
#        key_1 = xcamshift.Extra_potential.Atom_key(0,"H")
#        table.get_extra_shift(0,"HN",0,"HA","HN")
        
        target_atoms = "HA","CA", "HN", "N", "C", "CB"
        
        atoms_1 =    "HN", "HN", "HN", "C", "C", "C", "O", "O", "O", "N", "N", "N", "O", "O", "O", "N", "N", "N", "CG", "CG", "CG", "CG", "CG", "CG", "CG", "CA"
        offsets_1 =  0, 0, 0, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, -1
        atoms_2 =    "HA", "C", "CB", "HA", "C", "CB", "HA", "N", "CB", "HA", "N", "CB", "HA", "N", "CB", "HA", "N", "CB", "HA", "N", "C", "C", "N", "CA", "CA", "CA"
        offets_2 =   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, 0, 0, 0, -1, 1, 0, 0, 1
            
            
            
        for atom_1,offset_1,atom_2,offset_2 in zip(atoms_1,offsets_1,atoms_2,offets_2):
            
            key_1 =  Atom_key(offset_1,atom_1)
            key_2 = Atom_key(offset_2,atom_2)
            for target_atom in target_atoms:
                extra  = table.get_extra_shift(target_atom,key_1,key_2)
                
                self.assertIsNotNone(extra)
                
    def testTupleIt(self):
        test = [[1, 2], [3, 4]]
       
        result = Table_manager.tupleit(test)
       
        self.assertTrue(isinstance(result, tuple))
        self.assertTrue(isinstance(result[0], tuple))
        self.assertTrue(isinstance(result[1], tuple))
        
    def testLoadDihedral(self):

        table = self.table_manager.get_dihedral_table('ALA')
        
        for raw_key in self.raw_dihderal_keys:
            for target_atom in ["HA",  "CA", "HN", "N", "C", "CB"]:
                keys = [Atom_key(offset,atom) for atom,offset in raw_key]
                 
                dihedral_key = Dihedral_key(*keys)
                table.get_dihedral_shift(target_atom, dihedral_key) 
                
    def testGetDihdedralKeys(self):
        table = self.table_manager.get_dihedral_table('ALA')
        
        expected = set (self.raw_dihderal_keys)
        
        for atom_key in table.get_dihedral_keys():
            raw_key = tuple([(key.atom,key.offset) for key in atom_key])
            self.assertIn(raw_key, expected)
        self.assertEqual(len(expected),len(table.get_dihedral_keys()) )
    
    def testGetDihedral(self):        
        table = self.table_manager.get_dihedral_table('ALA')
         
        dihedral_key = self.raw_key_to_dihedral_key(self.raw_dihderal_keys[0])
            
        result = table.get_dihedral_shift("HA",dihedral_key)
        self.assertAlmostEqual(result, 0.39606482662562)
    def raw_key_to_dihedral_key(self, raw_key):
        keys = [Atom_key(offset, atom) for atom, offset in raw_key]
        dihedral_key = Dihedral_key(*keys)
        return dihedral_key

    def testLoadDihedralParameters(self):

        table = self.table_manager.get_dihedral_parameter_table('ALA')
        
        
        
        for raw_key in self.raw_dihderal_keys:
            for target_atom in ["HA",  "CA", "HN", "N", "C", "CB"]:
                dihedral_key = self.raw_key_to_dihedral_key(raw_key)
                for parameter in range(5):
                 
                    table.get_dihedral_parameter(target_atom, dihedral_key, parameter)
            
    def testGetDihedralParameter(self):        
        table = self.table_manager.get_dihedral_parameter_table('ALA')
         
        dihedral_key = self.raw_key_to_dihedral_key(self.raw_dihderal_keys[0])
            
        result = table.get_dihedral_parameter("HA",dihedral_key,0)
        self.assertAlmostEqual(result, 0.3)
                
    def testGetDihdedralParameterKeys(self):
        table = self.table_manager.get_dihedral_parameter_table('ALA')
        
        expected = set (self.raw_dihderal_keys)
        
        
        for atom_key in table.get_dihedral_keys():
            raw_key = tuple([(key.atom,key.offset) for key in atom_key])
            self.assertIn(raw_key, expected)
        self.assertEqual(len(expected),len(table.get_dihedral_keys()) )
            

if __name__ == "__main__":
    unittest2.main()
