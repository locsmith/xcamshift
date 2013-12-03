#-------------------------------------------------------------------------------
# Copyright (c) 2013 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 31 Dec 2011

@author: garyt
'''
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from xcamshift import Hbond_donor_indexer, Hbond_acceptor_indexer, Hydrogen_bond_context, Hbond_atom_type_indexer,\
     Hydrogen_bond_donor_context, Hydrogen_bond_acceptor_context, Hydrogen_bond_donor_component_factory, DONOR,ACCEPTOR,\
     Hbond_atom_type_indexer, Hydrogen_bond_acceptor_component_factory, Xcamshift, Hydrogen_bond_parameter_factory
from cython.shift_calculators import Fast_hydrogen_bond_calculator
from cython.fast_segment_manager import Segment_Manager
from utils import Atom_utils
from table_manager import Table_manager
from atomSel import AtomSel
from component_list import Native_component_list

EXPECTED_ACCEPTORS =   ((7, 'O', 'C'),
                        (8, 'O', 'C'),
                        (9, 'O', 'C'),
                        (10, 'O', 'C'),
                        (11, 'O', 'C'),
                        (12, 'O', 'C'),
                        (13, 'O', 'C'))

EXPECTED_DONORS =  ((8, 'HN', 'N'),
                    (9, 'HN', 'N'),
                    (10, 'HN', 'N'),
                    (11, 'HN', 'N'),
                    (12, 'HN', 'N'),
                    (13, 'HN', 'N'),
                    (14, 'HN', 'N'))

EXPECTED_INDIRECT_DONORS = {}
for elem in EXPECTED_DONORS:
    EXPECTED_INDIRECT_DONORS['',elem[0],elem[2]] = '',elem[0],elem[1]
    
EXPECTED_DIRECT_DONORS = [('',elem[0],elem[1]) for elem in EXPECTED_DONORS]

EXPECTED_ATOMS =  set(['O','HN'])

EXPECTED_INDIRECT_ACCEPTORS = {}
for elem in EXPECTED_ACCEPTORS:
    EXPECTED_INDIRECT_ACCEPTORS['',elem[0],elem[2]] = '',elem[0],elem[1]
    
EXPECTED_DIRECT_ACCEPTORS = [('',elem[0],elem[1]) for elem in EXPECTED_ACCEPTORS]

class TestXcamshiftHBondINGKTLKG(unittest2.TestCase):

    def __init__(self,*args,**kwargs):
        super(TestXcamshiftHBondINGKTLKG, self).__init__(*args,**kwargs)
        self._esim = None
        
    def get_single_member_ensemble_simulation(self):
        if self._esim.__class__ ==  None.__class__:
            #TODO note EnsembleSimulation can't have a single member that causes a crash!
            # therefore a hack
            self._esim =  Xcamshift().ensembleSimulation()
        return self._esim
    
    def assertEmpty(self, expected_keys, msg=""):
        return self.assertEqual(len(expected_keys), 0, msg)
    
    def assertLength(self,sequence,length,msg='bad length expected length of %i'):
        msg = msg % length
        self.assertEqual(len(sequence), length, msg)
        
    
    def check_almost_equal(self, list_1, list_2, delta = 1e-7):
        difference_offset = -1
        for i, (elem_1, elem_2) in enumerate(zip(list_1, list_2)):
            diff = abs(elem_1 - elem_2)
            if diff > delta:
                difference_offset = i
        
        return difference_offset
    
    def are_almost_equal_sequences(self, list_1, list_2, delta =  1e-7):
        result = True
        if self.check_almost_equal(list_1, list_2, delta) > 0:
            result = False
        return result
        
    def assertSequenceAlmostEqual(self,result,expected, delta = 1e-7, msg=""):
        len_result = len(result)
        len_expected = len(expected)
        if len_result != len_expected:
            raise AssertionError("the two lists are of different length %i and %i" % (len_result,len_expected))
        
        difference_offset = self.check_almost_equal(result, expected, delta)
        
            
        if difference_offset > 0:
            if msg != "":
                msg = msg + " "
                
            template = "%slists differ at item %i: %s - %s > %s"
            elem_1 = result[difference_offset]
            elem_2 = expected[difference_offset]
            message = template % (msg,difference_offset, `elem_1`,`elem_2`,delta)
            raise AssertionError(message)                 
    def setUp(self):
        initStruct("test_data/ingktlkg_hbond/INGKTLKG.psf")
        PDBTool("test_data/ingktlkg_hbond/INGKTLKG.pdb").read()
        Atom_utils.clear_cache()

        table_manager =  Table_manager.get_default_table_manager()
        self.donor_indexer  = Hbond_donor_indexer(table_manager)
        self.acceptor_indexer  = Hbond_acceptor_indexer(table_manager)
        self.atom_type_indexer =  Hbond_atom_type_indexer(table_manager)
        
        Segment_Manager.reset_segment_manager()
#         print "In method", self._testMethodName

    def test_donor_and_acceptor_indexers(self):
        

        donors = [donor for donor in self.donor_indexer.iter_keys()]
        self.assertSequenceEqual(donors, EXPECTED_DONORS)

        
        acceptors = [acceptor for acceptor in self.acceptor_indexer.iter_keys()]
        self.assertSequenceEqual(acceptors, EXPECTED_ACCEPTORS)
    
    def test_get_max_index(self):
        self.assertEqual(self.donor_indexer.get_max_index(), len(EXPECTED_DONORS))
        self.assertEqual(self.acceptor_indexer.get_max_index(), len(EXPECTED_ACCEPTORS))

    def test_get_name(self):    
        self.assertTrue('donor' in self.donor_indexer.get_name().lower())
        self.assertTrue('acceptor' in self.acceptor_indexer.get_name().lower())


    def test_get_index_for_key(self,):
        for i,acceptor in enumerate(EXPECTED_ACCEPTORS):
            self.assertEqual(i,self.acceptor_indexer.get_index_for_key(acceptor))
        self.assertEqual(i+1, self.acceptor_indexer.get_max_index())

        for j,donor in enumerate(EXPECTED_DONORS):
            self.assertEqual(j,self.donor_indexer.get_index_for_key(donor))
        self.assertEqual(j+1, self.donor_indexer.get_max_index())
                       
     
    def test_get_key_for_index(self):
        for i,acceptor in enumerate(EXPECTED_ACCEPTORS):
            self.assertEqual(acceptor,self.acceptor_indexer.get_key_for_index(i))
        self.assertEqual(i+1, self.acceptor_indexer.get_max_index())

        for j,donor in enumerate(EXPECTED_DONORS):
            self.assertEqual(donor,self.donor_indexer.get_key_for_index(j))
        self.assertEqual(j+1, self.donor_indexer.get_max_index())
    
         
    def test_hbond_context(self):
        atom = Atom_utils.find_atom('', 10, 'HN')[0]
        offset_data = (0, 'HN')
        table  =  Table_manager.get_default_table_manager().get_hydrogen_bond_table('LYS')
        hbond_context = Hydrogen_bond_context(atom,offset_data,table)
        
    def test_atom_indexer_indexers(self):
         
 
        atom_indices = [index for index in self.atom_type_indexer.iter_keys()]
        self.assertEqual(len(EXPECTED_ATOMS), len(atom_indices))
        
        self.assertEqual(sorted(EXPECTED_ATOMS), atom_indices)
        
      
    def test_get_max_index(self):
        self.assertEqual(self.atom_type_indexer.get_max_index(), len(EXPECTED_ATOMS))
 
    def test_get_name(self):    
        for elem in 'hydrogen', 'bond', 'atom','type':
            self.assertTrue(elem in self.atom_type_indexer.get_name().lower(), elem)
        
 
 
    def test_get_index_for_key(self,):
        for i,atom_name in enumerate(sorted(EXPECTED_ATOMS)):
            self.assertEqual(i,self.atom_type_indexer.get_index_for_key(atom_name))
        self.assertEqual(i+1, self.atom_type_indexer.get_max_index())
 
        
                        
      
    def test_get_key_for_index(self):
        for i,atom_name in enumerate(EXPECTED_ATOMS):
            self.assertEqual(atom_name,self.atom_type_indexer.get_key_for_index(i))
        self.assertEqual(i+1, self.atom_type_indexer.get_max_index())

     
          
    def test_hbond_context(self):
        atom = Atom_utils.find_atom('', 10, 'HN')[0]
        offset_data_0 = (0, 'O')
        table  =  Table_manager.get_default_table_manager().get_hydrogen_bond_table('LYS')
        
        hbond_context_0 = Hydrogen_bond_context(atom,offset_data_0,table)
        
        self.assertTrue(hbond_context_0.complete)        
        self.assertEqual(Atom_utils.find_atom_ids('   ',10,'HN')[0], hbond_context_0.target_atom_index)
        self.assertEqual(Atom_utils.find_atom_ids('   ',10,'O')[0], hbond_context_0.hbond_atom_index)
        
        EXPECTED_COEFFS_0 = [-0.00000547, 0.00023294, -0.00001567]
        self.assertAlmostEqual(hbond_context_0.coeffs, EXPECTED_COEFFS_0)
        self.assertSequenceAlmostEqual(EXPECTED_COEFFS_0, hbond_context_0.coeffs)
        
        offset_data_1 = (-1, 'O')
        hbond_context_1 = Hydrogen_bond_context(atom,offset_data_1,table)

        EXPECTED_COEFFS_1 = [-0.00000010, -0.00123116, 0.00012507]
        self.assertTrue(hbond_context_1.complete)
        self.assertEqual(Atom_utils.find_atom_ids('   ',10,'HN')[0], hbond_context_1.target_atom_index)
        self.assertEqual(Atom_utils.find_atom_ids('   ',9,'O')[0], hbond_context_1.hbond_atom_index)
        self.assertSequenceAlmostEqual(EXPECTED_COEFFS_1, hbond_context_1.coeffs)

        offset_data_2 = (-3, 'HN')
        hbond_context_2 = Hydrogen_bond_context(atom,offset_data_2,table)

        self.assertFalse(hbond_context_2.complete)
        
    
    def test_hydrogen_bond_donor_context(self):

        atom_0 = Atom_utils.find_atom('', 10, 'HN')[0]
        table  =  Table_manager.get_default_table_manager().get_hydrogen_bond_table('LYS')
        
        hbond_donor_context_0 = Hydrogen_bond_donor_context(atom_0,table)
        self.assertTrue(hbond_donor_context_0.complete)
        self.assertEqual(hbond_donor_context_0.atom_type_id,0)
        self.assertEqual(hbond_donor_context_0.direct_atom_id, atom_0.index())
        self.assertEqual(hbond_donor_context_0.indirect_atom_id, Atom_utils.find_atom('', 10, 'N')[0].index())

    def test_hydrogen_bond_acceptor_context(self):
        
        atom_0 = Atom_utils.find_atom('', 10, 'O')[0]
        table  =  Table_manager.get_default_table_manager().get_hydrogen_bond_table('LYS')
        
        hbond_acceptor_context_1 = Hydrogen_bond_acceptor_context(atom_0,table)
        self.assertTrue(hbond_acceptor_context_1.complete)
        self.assertEqual(hbond_acceptor_context_1.atom_type_id,1)
        self.assertEqual(hbond_acceptor_context_1.direct_atom_id, atom_0.index())
        self.assertEqual(hbond_acceptor_context_1.indirect_atom_id, Atom_utils.find_atom('', 10, 'C')[0].index())
        
    


    def _build_component_list(self, factory, format):
        component_list = Native_component_list(format)
        table_provider = Table_manager.get_default_table_manager().get_hydrogen_bond_table
        segment = '    '
        #TODO note an oddity here as this takes a residu butr builds a list for all residues...
        target_residue_number = 10
        selected_atoms = AtomSel('(all)')
        factory.create_components(component_list, table_provider, segment, target_residue_number, selected_atoms)
        return component_list

    def _do_test_donor_acceptor_components(self, factory,expected_direct_donors_or_acceptors, expected_indirect_donors_or_acceptors, donor_or_acceptor):

        component_list = self._build_component_list(factory, 'i' * 5)
        for i,component in enumerate(component_list):
            INDEX = 0
            DIRECT_ATOM_ID = 1
            INDIRECT_ATOM_ID = 2
            DONOR_OR_ACCEPTOR = 3
            ATOM_TYPE = 4
            self.assertEqual(i,component[INDEX])
            atom_key = Atom_utils._get_atom_info_from_index(component[DIRECT_ATOM_ID])
            indirect_atom_key = Atom_utils._get_atom_info_from_index(component[INDIRECT_ATOM_ID])
            self.assertIn(atom_key, expected_direct_donors_or_acceptors)
            expected_direct_donors_or_acceptors.remove(atom_key)
            self.assertEqual(atom_key,expected_indirect_donors_or_acceptors[indirect_atom_key])
            del expected_indirect_donors_or_acceptors[indirect_atom_key]
            self.assertEqual(component[DONOR_OR_ACCEPTOR], donor_or_acceptor)
            self.assertEqual(component[ATOM_TYPE], Hbond_atom_type_indexer(Table_manager.get_default_table_manager()).get_index_for_key(atom_key[INDIRECT_ATOM_ID]))
        
        self.assertEmpty(expected_direct_donors_or_acceptors)
        self.assertEmpty(expected_indirect_donors_or_acceptors)

    def test_hydrogen_bond_donor_components(self):
        self._do_test_donor_acceptor_components(Hydrogen_bond_donor_component_factory(), set(EXPECTED_DIRECT_DONORS), dict(EXPECTED_INDIRECT_DONORS), DONOR)
        
    def test_hydrogen_bond_acceptor_components(self):
        self._do_test_donor_acceptor_components(Hydrogen_bond_acceptor_component_factory(), set(EXPECTED_DIRECT_ACCEPTORS), dict(EXPECTED_INDIRECT_ACCEPTORS), ACCEPTOR)
    
    def test_fast_hydrogen_bond_calculator(self):
        test = Fast_hydrogen_bond_calculator(self.get_single_member_ensemble_simulation())
        format = 'i' * 5
        donor_components = self._build_component_list(Hydrogen_bond_donor_component_factory(),format)
        acceptor_components = self._build_component_list(Hydrogen_bond_acceptor_component_factory(),format)
        test(donor_components.get_native_components(), acceptor_components.get_native_components(), None)
    
    def test_parameter_components(self):
        factory = Hydrogen_bond_parameter_factory()
        format = ('i'*4) +('f'*8)
        component_list = self._build_component_list(factory, format)
        self.assertLength(component_list, 3)
        EXPECTED = [(0, 0, 0, 0, 12.151400000000001, -1.1225799999999999, -1.1225799999999999, 1.3095300000000001, 1.3095300000000001, 0.0, 3.0, 15.5), 
                    (0, 1, 0, 0, 15.0486, -1.7946899999999999, -1.7946899999999999, 1.4939800000000001, 1.4939800000000001, 0.0, 2.5, 15.5), 
                    (0, 2, 0, 0, -0.40353099999999997, 13.9148, -0.42172500000000002, -3.0801599999999998, 0.98725099999999999, 0.0, -16.5, 3.0)]
        
        for i in range(3):
            self.assertSequenceAlmostEqual(component_list[i], EXPECTED[i])
        component_list.get_native_components()
def run_tests():
    unittest2.main(module='test.test_xcamshift_hbond_ingktlkg')
    
if __name__ == "__main__":
    run_tests()
