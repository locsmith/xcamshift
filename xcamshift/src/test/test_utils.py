#-------------------------------------------------------------------------------
# Copyright (c) 2013 gary thompson.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on Apr 23, 2012

@author: garyt
'''
import unittest2
from protocol import initStruct
from pdbTool import PDBTool
from utils import Atom_utils, iter_residue_atoms, iter_residue_atom_ids
import utils
from cython.fast_segment_manager import Segment_Manager
from test.util_for_testing import Virtual_list, Empty_loader
from test.util_for_testing import Yaml_loader
from StringIO import  StringIO

YAML_FILE_1 = '''a:
               l: {m : 'n'}
         '''
YAML_DICT_1 = {'a':{'l' :{'m':'n'}}}
YAML_FILE_2 = '''b:
               x: {y : 1}
         '''
YAML_DICT_2 = {'b':{'x' :{'y':1}}}
                             
expected_residue_atom_ids  =  (
     (1,2,3,4,5,6,7,8,9,10,11,12),
     (13,14,15,16,17,18,19),
     (20,21,22,23,24,25,26,27,28,29),
     (30,31,32,33,34,35,36),
     (37,38,39,40,41,42,43,44,45,46,47)
)

def select_gly(atom):
    result =False
    if atom.residueName() == "GLY":
        result = True
    return result

def select_name_gly(name):
    result =False
    if name == "GLY":
        result = True
    return result

def select_residue_2(residue_atoms):
    result = False
    if residue_atoms[0].residueNum() == 2:
        result = True 
    return result

def almostEqual(first, second, places = 7):
    result  = False
    if round(abs(second-first), places) == 0:
        result=True
    return result

class TestXcamshiftUtils(unittest2.TestCase):
    DEFAULT_DECIMAL_PLACES = 5
    DEFAULT_ERROR = 10 ** -DEFAULT_DECIMAL_PLACES

    def remove_zero_valued_keys(self, expected_force_factors):
        for key, value in expected_force_factors.items():
            if almostEqual(value, 0.0, self.DEFAULT_DECIMAL_PLACES):
                del expected_force_factors[key]
    
    def assertEmpty(self, expected_force_factors,msg = None):
        return self.assertEqual(len(expected_force_factors), 0)

                
    def setUp(self):
        initStruct("test_data/agaga/agaga.psf")
        PDBTool("test_data/agaga/agaga.pdb").read()
        Atom_utils.clear_cache()
        Segment_Manager.reset_segment_manager()
        
    def test_iterate_atoms(self):
        for i,atom in  enumerate(utils.iter_atoms()):
            self.assertEqual(atom.index(), i)

    def test_iterate_atom_ids(self):
        for i,atom_id in  enumerate(utils.iter_atom_ids()):
            self.assertEqual(atom_id, i)
            
    def test_iterate_atoms_type_gly_predicate(self):

        
        expected = set([13,14,15,16,17,18,19,30,31,32,33,34,35,36])
        for atom in utils.iter_atoms(select_gly):
            self.assertIn(atom.index()+1, expected)
            
    def test_iterate_residue_types(self):
        residue_types = [residue_type for residue_type in utils.iter_residue_types()]
        expected_types = ['ALA','GLY']
        
        self.assertSequenceEqual(residue_types, expected_types)
        
    def test_iterate_residue_types_gly_predicate(self):
        residue_types = [residue_type for residue_type in utils.iter_residue_types(select_name_gly)]
        expected_types = ['GLY']
        
        self.assertSequenceEqual(residue_types, expected_types)
        
    def test_iterate_residue_atoms(self):
        
        for i,residue_atoms in enumerate(iter_residue_atoms()):
            residue_atom_ids =  [atom.index()+1 for atom in residue_atoms]
            self.assertSequenceEqual(residue_atom_ids, expected_residue_atom_ids[i])
            
    def test_iterate_residue_atom_ids(self):
        
        for i,residue_atom_ids in enumerate(iter_residue_atom_ids()):
            residue_atom_ids =  [atom_id+1 for atom_id in residue_atom_ids]
            self.assertSequenceEqual(residue_atom_ids, expected_residue_atom_ids[i])

    def test_iterate_residue_atoms_predicate_residue_2(self):
        
        for residue_atoms in iter_residue_atoms(select_residue_2):
            residue_atom_ids =  [atom.index()+1 for atom in residue_atoms]
            self.assertSequenceEqual(residue_atom_ids, expected_residue_atom_ids[1]) 
            
    def test_find_all_atoms(self):
        all_atoms = Atom_utils.find_all_atoms()
        self.assertEqual(len(all_atoms),47)
        for i,atom in enumerate(all_atoms):
            self.assertEqual(atom.index(), i)
            
    def test_virtual_list_empty_loader(self):
        test_list = Virtual_list(Empty_loader())

        self.assertEmpty(test_list)
        with self.assertRaises(IndexError):
            test_list[0]
            
    class Test_yaml_loader(Yaml_loader):
        
        def __init__(self):
            super(TestXcamshiftUtils.Test_yaml_loader, self).__init__(['a','b'])
            
        def _open_stream(self,file):
            if file == "a":
                result = YAML_FILE_1
            elif file == "b":
                result = YAML_FILE_2
            return  StringIO(result)
        
        def length(self):
            return 2
                        
    def test_virtual_list_yaml_loader(self):
       
        test_list = Virtual_list(TestXcamshiftUtils.Test_yaml_loader())

        self.assertEqual(len(test_list),2)
        self.assertSequenceEqual([YAML_DICT_1, YAML_DICT_2], test_list)
        self.assertIsNot(test_list[0], test_list[0])
        
    def test_virtual_list_yaml_loader_makes_new_objects(self):
        
        test_list = Virtual_list(TestXcamshiftUtils.Test_yaml_loader())
        self.assertIsNot(test_list[0], test_list[0])
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestXcamshifAGA.testName']
    unittest2.main()
