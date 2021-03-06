#-------------------------------------------------------------------------------
# Copyright (c) 2013-2015 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
#
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
from cache_manager import Cache_manager
'''
Created on 30 Dec 2011

@author: garyt
'''
import unittest2
from table_manager import Table_manager, TABLE_MANAGER, Database_not_found_exception
from keys import Atom_key, Dihedral_key
from protocol import initStruct
from pdbTool import PDBTool
from utils import Atom_utils
from cython.fast_segment_manager import Segment_Manager



class Test_table_manager(unittest2.TestCase):
    raw_dihderal_keys = ((("C", -1), ("N", 0), ("CA", 0), ("C",  0)),
                         (("N",  0), ("CA", 0), ("C", 0), ("N",  1)),
                         (("N", 0), ("CA", 0), ("CB", 0), ("CG", 0)))
    TARGET_ATOMS = ["HA", "CA", "HN", "N", "C", "CB"]

    def assertSequenceContains(self,expected,sequence):
        if sequence.index(expected) == -1:
            sequence_strings = [`elem` for elem in sequence]
            sequence_string = '\n'.join(sequence_strings)
            raise AssertionError("element %s not found in the sequence:\n%s" % sequence_string,`expected`)

    def setUp(self):
        # TODO add a delegate to get a sequence to the table manager and remove dependance on a psf file
        initStruct("test_data/gb3/gb3.psf")
        Segment_Manager.reset_segment_manager()
        cache_manager = Cache_manager.get_cache_manager()
        cache_manager.clear_cache()
        self.table_manager = cache_manager.get_cache(TABLE_MANAGER)
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
        from_atom_list = table.get_target_atoms()
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



    def testLoadDihedral(self):

        table = self.table_manager.get_dihedral_table('ALA')

        for raw_key in self.raw_dihderal_keys:
            for target_atom in ["HA",  "CA", "HN", "N", "C", "CB"]:
                keys = [Atom_key(offset,atom) for atom,offset in raw_key]

                dihedral_key = Dihedral_key(*keys)
                table.get_dihedral_shift(target_atom, dihedral_key)

    def raw_key_to_dihedral_key(self, raw_key):
        keys = [Atom_key(offset, atom) for atom, offset in raw_key]
        dihedral_key = Dihedral_key(*keys)
        return dihedral_key

    def testGetDihdedralKeys(self):
        table = self.table_manager.get_dihedral_table('ALA')

        expected = set (self.raw_dihderal_keys)

        for atom_key in table.get_dihedral_keys():
            raw_key = tuple([(key.atom,key.offset) for key in atom_key])
            self.assertIn(raw_key, expected)
        self.assertEqual(len(expected),len(table.get_dihedral_keys()) )

    def testGetDihedral(self):
        table = self.table_manager._get_dihedral_table('ALA')

        dihedral_key = self.raw_key_to_dihedral_key(self.raw_dihderal_keys[0])

        result = table.get_dihedral_shift("HA",dihedral_key)
        self.assertAlmostEqual(result, 0.39606482662562)


    def testLoadDihedralParameters(self):

        table = self.table_manager.get_dihedral_table('ALA')



        for raw_key in self.raw_dihderal_keys:
            for target_atom in self.TARGET_ATOMS:
                dihedral_key = self.raw_key_to_dihedral_key(raw_key)
                for parameter in range(5):

                    table.get_parameter(target_atom, dihedral_key, parameter)

    def testGetDihedralParameter(self):
        table = self.table_manager.get_dihedral_table('ALA')

        dihedral_key = self.raw_key_to_dihedral_key(self.raw_dihderal_keys[0])

        result = table.get_parameter("HA",dihedral_key,0)
        self.assertAlmostEqual(result, 0.3)

    def testGetSidechainTable(self):
        table = self.table_manager.get_sidechain_table('ALA')

        self.assertAlmostEqual(table.get_exponent(), 1.0)
        self.assertIn("ALA", table.get_residue_types())
        self.assertIn("CB", table.get_sidechain_atoms("ALA"))
        self.assertAlmostEqual(table.get_sidechain_coefficient("ALA", "CA", "CB"),-3.4713212)
        self.assertSequenceEqual(self.TARGET_ATOMS, table.get_target_atoms())

    def testGetConstantTable(self):
        table  =  self.table_manager.get_constants_table('ALA')
        self.assertAlmostEqual(0.4, table.get_flat_bottom_constant())
        self.assertAlmostEqual(1.151132, table.get_flat_bottom_limit("CA"))
        self.assertAlmostEqual(20.000000,table.get_end_harmonic("CA"))
        self.assertAlmostEqual(2.360000,table.get_scale_harmonic("CA"))
        self.assertAlmostEqual(1.0, table.get_weight("CA"))

    def testGetRingTable(self):
        table = self.table_manager.get_ring_table('ALA')

        result = table.get_ring_coefficient("CA", "PHE", "6")
        self.assertAlmostEqual(0.01072141*1000,result)

        EXPECTED_RING_RESIDUES = ("HIS","PHE","TRP","TYR")
        result = table.get_residue_types()
        result.sort()
        self.assertSequenceEqual(result,EXPECTED_RING_RESIDUES)

        EXPECTED_RING_TYPES = (("HIS",("5",)),("PHE",("6",)),("TRP",("5","6")),("TYR",("6",)))
        for residue_type,expected_ring_types in EXPECTED_RING_TYPES:
            ring_types = table.get_ring_types(residue_type)
            ring_types.sort()

            self.assertSequenceEqual(ring_types, expected_ring_types)

    def testGetNonBondedTable(self):
        table = self.table_manager.get_non_bonded_table('ALA')

        result = table.get_non_bonded_coefficient("CA", "sphere_1", "S", "SP3")

        self.assertAlmostEqual(result, -101.81223165473585)

        EXPECTED_SPHERE_TYPES =  ('sphere_1','sphere_2')
        spheres =   table.get_spheres()
        self.assertSequenceEqual(EXPECTED_SPHERE_TYPES, spheres)

        EXPECTED_COEFFICIENT_KEYS = [('C', 'SP3'),('H', 'None'), ('N', 'SP3'),
                                     ('O', 'SP3'), ('S', 'SP3'),('C', 'SP2'),
                                     ('N', 'SP2'),('O', 'SP2')]
        EXPECTED_COEFFICIENT_KEYS.sort()
        EXPECTED_COEFFICIENT_KEYS =  tuple(EXPECTED_COEFFICIENT_KEYS)

        coefficient_keys = table.get_remote_atom_types('sphere_1')
        self.assertSequenceEqual(EXPECTED_COEFFICIENT_KEYS,coefficient_keys)

        coefficient_keys = table.get_remote_atom_types('sphere_2')
        self.assertSequenceEqual(EXPECTED_COEFFICIENT_KEYS,coefficient_keys)

        exponent = table.get_exponent('sphere_1')
        self.assertAlmostEqual(exponent, -3.0)

        exponent = table.get_exponent('sphere_2')
        self.assertAlmostEqual(exponent, 1.0)

        chem_types = table.get_chem_types()
        chem_types.sort()
        expected_chem_types =  ['C',  'CA', 'CB', 'CC', 'CN', 'CP', 'CR', 'CT',
                                'CV', 'CW', 'CX', 'N',  'NA', 'NB', 'NC2','NH1',
                                'NH2','NH3','O',  'OC', 'OH', 'S','HA','HC','H']
        expected_chem_types.sort()
        self.assertSequenceEqual(chem_types, expected_chem_types)

        translation = table.get_chem_type_translation("CA")
        self.assertEqual(translation, ('C','SP2'))

        translation = table.get_chem_type_translation("NH1")
        self.assertEqual(translation, ('N','SP3'))

    def testGetDisulphideTable(self):
        initStruct("test_data/acaggaca/acaggaca.psf")
        table = self.table_manager.get_disulphide_table('CYS')

        self.assertAlmostEqual(table.get_disulphide_shift('HA'),0.01155515)

        self.assertSequenceEqual(table.get_residue_types(), ['CYS'])


    def test_table_hierachy(self):
        table = self.table_manager.get_random_coil_table('base')
        self.assertEqual(table._table['residue_type'], 'base',)

        table = self.table_manager.get_random_coil_table('GLY')
        self.assertEqual(table._table['residue_type'], 'base',)

        table = self.table_manager.get_non_bonded_table('base')
        self.assertEqual(table._table['residue_type'], 'base',)

        table = self.table_manager.get_non_bonded_table('GLY')
        self.assertEqual(table._table['residue_type'], 'gly',)

        expected_exponent =  -3.0
        exponent = table.get_exponent('sphere_1')
        self.assertAlmostEqual(expected_exponent, exponent)

        expected_coefficent = -90.04553541
        coefficient = table.get_non_bonded_coefficient("CA", "sphere_1", "S", "SP3")
        self.assertAlmostEqual(expected_coefficent, coefficient)

    #TODO caching shouldn't be on the _table value but on the wrapper class
    def test_caching(self):
        table_1 = self.table_manager.get_non_bonded_table('GLY')
        table_2 = self.table_manager.get_non_bonded_table('GLY')
        self.assertEqual(id(table_1._table), id(table_2._table))


    def test_hydrogen_bond_table(self):

        table = self.table_manager.get_hydrogen_bond_table('ALA')
#
        self.assertSequenceEqual(('HA', 'CA', 'HN', 'N', 'C', 'CB'), table.get_target_atoms())

        self.assertSequenceEqual(table.get_donor_types(),['HON'])
        self.assertSequenceEqual(table.get_acceptor_types(),['ON'])

        self.assertSequenceContains(['.','C', 'O'], table.get_atom_selector('ON'))
        self.assertSequenceContains(['.','N', 'HN'], table.get_atom_selector('HON'))

        with self.assertRaises(KeyError):
            table.get_atom_selector('J')

        expected_energy_terms = {
                    'p1': -0.40353099999999997,
                    'p2': 13.9148,
                    'p3': -0.42172500000000002,
                    'p4': -3.0801599999999998,
                    'p5':  0.98725099999999999,
                    'p6':  0.0,
                    's': -16.5,
                    'r': 3.0}

        self.assertDictEqual(table.get_energy_terms('HON','ON','DIST'),expected_energy_terms)

        expected_energy_term_offsets = [(-1, 'O'), (0, 'HN'), (0, 'O'), (1, 'HN')]
        self.assertSequenceEqual(table.get_energy_term_offsets(),expected_energy_term_offsets)

        self.assertAlmostEqual(-0.00000002820/table.get_energy_offset_correction((-1, 'O'),'DIST','HA'),1.0)

        table = self.table_manager.get_hydrogen_bond_table('GLY')

        self.assertDictEqual(table.get_energy_terms('HON','ON','DIST'),expected_energy_terms)

        self.assertAlmostEqual(-0.00000029080/table.get_energy_offset_correction((-1, 'O'),'DIST','HA'), 1.0)

    def test_empty_search_path(self):
        self.table_manager.search_paths = []
        with self.assertRaises(Database_not_found_exception):
            print self.table_manager.get_hydrogen_bond_table('ALA')


if __name__ == "__main__":
#     unittest2.main()
    unittest2.main(module='test.Test_table_manager')

