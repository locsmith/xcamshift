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
from xcamshift import Hbond_backbone_donor_indexer, Hbond_backbone_acceptor_indexer, Hydrogen_bond_context, Hbond_atom_type_indexer,\
     Hydrogen_bond_donor_context, Hydrogen_bond_acceptor_context, Hydrogen_bond_donor_component_factory, DONOR,ACCEPTOR,\
     Hbond_donor_atom_type_indexer, Hbond_acceptor_atom_type_indexer, Hydrogen_bond_acceptor_component_factory, Xcamshift, Hydrogen_bond_parameter_factory, \
     Hydrogen_bond_donor_lookup_factory, Hydrogen_bond_acceptor_lookup_factory,\
    Hydrogen_bond_potential
from cython.shift_calculators import Fast_hydrogen_bond_calculator, allocate_array
from cython.fast_segment_manager import Segment_Manager
from utils import Atom_utils
from table_manager import Table_manager
from atomSel import AtomSel
from component_list import Native_component_list

BACKBONE=1
SIDE_CHAIN=0
EXPECTED_ACCEPTORS_BASE =   (((7,  'O'),   (7,  'C'),  BACKBONE),
                             ((8,  'N'),   (7,  'C'),  BACKBONE),
                             ((8,  'O'),   (8,  'C'),  BACKBONE),
                             ((9,  'N'),   (8,  'C'),  BACKBONE),
                             ((9,  'O'),   (9,  'C'),  BACKBONE),
                             ((10, 'N'),   (9,  'C'),  BACKBONE),
                             ((10, 'O'),   (10, 'C'),  BACKBONE),
                             ((10, 'NZ'),  (10, 'CE'), SIDE_CHAIN),
                             ((11, 'N'),   (10, 'C'),  BACKBONE),
                             ((11, 'O'),   (11, 'C'),  BACKBONE),
                             ((11, 'OG1'), (11, 'CB'), SIDE_CHAIN),
                             ((12, 'N'),   (11, 'C'),  BACKBONE),
                             ((12, 'O'),   (12, 'C'),  BACKBONE),
                             ((13, 'N'),   (12, 'C'),  BACKBONE),
                             ((13, 'O'),   (13, 'C'),  BACKBONE),
                             ((13, 'NZ'),  (13, 'CE'), SIDE_CHAIN),
                             ((14, 'N'),   (13, 'C'),  BACKBONE))

EXPECTED_DONORS_BASE =  (((8,  'HN'),  (8,  'N'),   BACKBONE),
                        ((9,  'HN'),  (9,  'N'),   BACKBONE),
                        ((10, 'HN'),  (10, 'N'),   BACKBONE),
                        ((10, 'HZ1'), (10, 'NZ'),  SIDE_CHAIN),
                        ((10, 'HZ2'), (10, 'NZ'),  SIDE_CHAIN),
                        ((10, 'HZ3'), (10, 'NZ'),  SIDE_CHAIN),
                        ((11, 'HN'),  (11, 'N'),   BACKBONE),
                        ((11, 'HG1'), (11, 'OG1'), SIDE_CHAIN),
                        ((12, 'HN'),  (12, 'N'),   BACKBONE),
                        ((13, 'HN'),  (13, 'N'),   BACKBONE),
                        ((13, 'HZ1'), (13, 'NZ'),  SIDE_CHAIN),
                        ((13, 'HZ2'), (13, 'NZ'),  SIDE_CHAIN),
                        ((13, 'HZ3'), (13, 'NZ'),  SIDE_CHAIN),
                        ((14, 'HN'),  (14, 'N'),   BACKBONE))

EXPECTED_DONORS = [(elem[0][0],elem[0][1],elem[1][1],elem[2]) for elem in  EXPECTED_DONORS_BASE]

DIST=0
ANG_1=1
ANG_2=2

EXPECTED_DONOR_ENERGIES = {
        ('',  9, 'HN', DIST):     3.64372,
        ('',  9, 'HN', ANG_1):  701.205,
        ('',  9, 'HN', ANG_2): 6548.52,
        ('', 11, 'HN', DIST):     4.45286,
        ('', 11, 'HN', ANG_1):  672.789,
        ('', 11, 'HN', ANG_2): 6252.28,
        ('', 12, 'HN', DIST):     4.66541,
        ('', 12, 'HN', ANG_1):  695.025,
        ('', 12, 'HN', ANG_2): 6238.2,
        ('', 14, 'HN', DIST):     3.27874313127,
        ('', 14, 'HN', ANG_1):  696.529729434,
        ('', 14, 'HN', ANG_2):  6683.4703906
}

EXPECTED_ACCEPTOR_ENERGIES = {
        ('',  7,  'O', DIST):     3.27874313127,
        ('',  7,  'O', ANG_1):  696.529729434,
        ('',  7,  'O', ANG_2):  6683.4703906,
        ('', 12,  'O', DIST):     3.64372,
        ('', 12,  'O', ANG_1):  701.205,
        ('', 12,  'O', ANG_2): 6548.52,
        ('', 10,  'N', DIST):     4.45286,
        ('', 10,  'N', ANG_1):  672.789,
        ('', 10,  'N', ANG_2): 6252.28,
        ('', 11,  'N', DIST):     4.66541,
        ('', 11,  'N', ANG_1):  695.025,
        ('', 11,  'N', ANG_2): 6238.2
}   



EXPECTED_INDIRECT_DONORS = {}
for elem in EXPECTED_DONORS:
    EXPECTED_INDIRECT_DONORS['',elem[0],elem[1]] = '',elem[0],elem[2]
    
EXPECTED_DIRECT_DONORS = [('',elem[0],elem[1]) for elem in EXPECTED_DONORS]

EXPECTED_DONOR_TYPES = sorted(['HON',])
EXPECTED_ACCEPTOR_TYPES = sorted(['ON',])

EXPECTED_INDIRECT_ACCEPTORS = {}
for elem in EXPECTED_ACCEPTORS_BASE:
    EXPECTED_INDIRECT_ACCEPTORS['',elem[0][0],elem[0][1]] = '',elem[1][0],elem[1][1]

    
EXPECTED_DIRECT_ACCEPTORS = [('',elem[0][0],elem[0][1]) for elem in EXPECTED_ACCEPTORS_BASE]


EXPECTED_BACK_BONE_DONORS = [('',elem[0],elem[1]) for elem in EXPECTED_DONORS if elem[-1] == 1]
EXPECTED_BACK_BONE_ACCEPTORS = [('',elem[0][0],elem[0][1]) for elem in EXPECTED_ACCEPTORS_BASE if elem[-1] == 1]



class TestXcamshiftHBondINGKTLKG(unittest2.TestCase):

    def __init__(self,*args,**kwargs):
        super(TestXcamshiftHBondINGKTLKG, self).__init__(*args,**kwargs)
        self._esim = None
        
        
    def assertSequenceContains(self,expected,sequence):
        if expected not in sequence:
            sequence_strings = [`elem` for elem in sequence]
            sequence_string = '\n'.join(sequence_strings)
            raise AssertionError("element %s not found in the sequence:\n%s" % sequence_string,`expected`)

        
    def assertSequenceDoesntContain(self,expected,sequence):
        if expected in sequence:
            sequence_strings = [`elem` for elem in sequence]
            sequence_string = '\n'.join(sequence_strings)
            raise AssertionError("element %s found in the sequence:\n%s" % (`expected`,sequence_string))        
        
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
        self.donor_indexer  = Hbond_backbone_donor_indexer(table_manager)
        self.acceptor_indexer  = Hbond_backbone_acceptor_indexer(table_manager)
        self.donor_atom_type_indexer =  Hbond_donor_atom_type_indexer(table_manager)
        self.acceptor_atom_type_indexer =  Hbond_acceptor_atom_type_indexer(table_manager)
        
        Segment_Manager.reset_segment_manager()
#         print "In method", self._testMethodName

    def test_backbone_donor_and_acceptor_indexers(self):
         
 
        donors = [donor for donor in self.donor_indexer.iter_keys()]
        self.assertSequenceEqual(sorted(donors), sorted(EXPECTED_BACK_BONE_DONORS)) 
         
        acceptors = [acceptor for acceptor in self.acceptor_indexer.iter_keys()]
        self.assertSequenceEqual(sorted(acceptors), sorted(EXPECTED_BACK_BONE_ACCEPTORS))
     
    def test_backbone_donor_and_acceptor_indexers_get_max_index(self):
        self.assertEqual(self.donor_indexer.get_max_index(), len(EXPECTED_BACK_BONE_DONORS))
        self.assertEqual(self.acceptor_indexer.get_max_index(), len(EXPECTED_BACK_BONE_ACCEPTORS))
 
    def test_backbone_donor_and_acceptor_indexers_get_name(self):    
        self.assertTrue('donor' in self.donor_indexer.get_name().lower())
        self.assertTrue('acceptor' in self.acceptor_indexer.get_name().lower())
 
 
    def test_backbone_donor_and_acceptor_indexers_get_index_for_key(self,):

        for i,acceptor in enumerate(EXPECTED_BACK_BONE_ACCEPTORS):
            self.assertEqual(i,self.acceptor_indexer.get_index_for_key(acceptor))
        self.assertEqual(i+1, self.acceptor_indexer.get_max_index())
 
        for j,donor in enumerate(EXPECTED_BACK_BONE_DONORS):
            self.assertEqual(j,self.donor_indexer.get_index_for_key(donor))
        self.assertEqual(j+1, self.donor_indexer.get_max_index())
                        
      
    def test_backbone_donor_and_acceptor_indexers_get_index_get_key_for_index(self):
        for i,acceptor in enumerate(EXPECTED_BACK_BONE_ACCEPTORS):
            self.assertEqual(acceptor,self.acceptor_indexer.get_key_for_index(i))
        self.assertEqual(i+1, self.acceptor_indexer.get_max_index())
 
        for j,donor in enumerate(EXPECTED_BACK_BONE_DONORS):
            self.assertEqual(donor,self.donor_indexer.get_key_for_index(j))
        self.assertEqual(j+1, self.donor_indexer.get_max_index())
    
         
    def test_hbond_context(self):
        atom = Atom_utils.find_atom('', 10, 'HN')[0]
        offset_data = (0, 'HN')
        table  =  Table_manager.get_default_table_manager().get_hydrogen_bond_table('LYS')
        hbond_context = Hydrogen_bond_context(atom,offset_data,table)
        
    def test_atom_indexer_indexers(self):
         
 
        donor_atom_indices = [index for index in self.donor_atom_type_indexer.iter_keys()]
        self.assertEqual(len(EXPECTED_DONOR_TYPES), len(donor_atom_indices))
        
        self.assertEqual(sorted(EXPECTED_DONOR_TYPES), donor_atom_indices)

        acceptor_atom_indices = [index for index in self.acceptor_atom_type_indexer.iter_keys()]
        self.assertEqual(len(EXPECTED_ACCEPTOR_TYPES), len(acceptor_atom_indices))
        
        self.assertEqual(sorted(EXPECTED_ACCEPTOR_TYPES), acceptor_atom_indices)        
      
    def test_get_max_index(self):
        self.assertEqual(self.donor_atom_type_indexer.get_max_index(), len(EXPECTED_DONOR_TYPES))
        self.assertEqual(self.acceptor_atom_type_indexer.get_max_index(), len(EXPECTED_ACCEPTOR_TYPES))
        
    def test_get_name(self):    
        for elem in 'hydrogen', 'bond', 'atom','type':
            self.assertTrue(elem in self.donor_atom_type_indexer.get_name().lower(), elem)
            self.assertTrue(elem in self.acceptor_atom_type_indexer.get_name().lower(), elem)
        
 
 
    def test_get_index_for_key(self,):
        for i,atom_name in enumerate(sorted(EXPECTED_DONOR_TYPES)):
            self.assertEqual(i,self.donor_atom_type_indexer.get_index_for_key(atom_name))
        self.assertEqual(i+1, self.donor_atom_type_indexer.get_max_index())
        
        for i,atom_name in enumerate(sorted(EXPECTED_ACCEPTOR_TYPES)):
            self.assertEqual(i,self.acceptor_atom_type_indexer.get_index_for_key(atom_name))
        self.assertEqual(i+1, self.acceptor_atom_type_indexer.get_max_index()) 
        
                        
      
    def test_get_key_for_index(self):
        for i,atom_name in enumerate(EXPECTED_DONOR_TYPES):
            self.assertEqual(atom_name,self.donor_atom_type_indexer.get_key_for_index(i))
        self.assertEqual(i+1, self.donor_atom_type_indexer.get_max_index())

        for i,atom_name in enumerate(EXPECTED_ACCEPTOR_TYPES):
            self.assertEqual(atom_name,self.acceptor_atom_type_indexer.get_key_for_index(i))
        self.assertEqual(i+1, self.acceptor_atom_type_indexer.get_max_index())
     
          
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
        self.assertEqual(hbond_acceptor_context_1.atom_type_id,0)
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

    def _do_test_donor_acceptor_components(self, factory,expected_direct_donors_or_acceptors, expected_indirect_donors_or_acceptors, donor_or_acceptor, expected_backbone):

        component_list = self._build_component_list(factory, 'i' * 5)
        for i,component in enumerate(component_list):
            INDEX = 0
            DIRECT_ATOM_ID = 1
            INDIRECT_ATOM_ID = 2
            DONOR_OR_ACCEPTOR = 3
            ATOM_TYPE = 4
            BACKBONE = 5
            
            self.assertEqual(i,component[INDEX])
            
            atom_key = Atom_utils._get_atom_info_from_index(component[DIRECT_ATOM_ID])
            indirect_atom_key = Atom_utils._get_atom_info_from_index(component[INDIRECT_ATOM_ID])
            
            if component[BACKBONE] > -1:
                self.assertSequenceContains(atom_key, expected_backbone)
                self.assertEqual(component[BACKBONE], expected_backbone.index(atom_key))
            else:
                self.assertSequenceDoesntContain(atom_key, expected_backbone)
            
            self.assertIn(atom_key, expected_direct_donors_or_acceptors)
            if atom_key in expected_direct_donors_or_acceptors:
                expected_direct_donors_or_acceptors.remove(atom_key)
            
            self.assertEqual(indirect_atom_key,expected_indirect_donors_or_acceptors[atom_key])
            del expected_indirect_donors_or_acceptors[atom_key]
            
            self.assertEqual(component[DONOR_OR_ACCEPTOR], donor_or_acceptor)
            
            if donor_or_acceptor ==  DONOR:
                self.assertEqual(component[ATOM_TYPE],0)
            elif donor_or_acceptor ==  ACCEPTOR:
                self.assertEqual(component[ATOM_TYPE],0)
        
        self.assertEmpty(expected_direct_donors_or_acceptors)
        self.assertEmpty(expected_indirect_donors_or_acceptors)

    def test_hydrogen_bond_donor_components(self):
        self._do_test_donor_acceptor_components(Hydrogen_bond_donor_component_factory(), set(EXPECTED_DIRECT_DONORS), dict(EXPECTED_INDIRECT_DONORS), DONOR, EXPECTED_BACK_BONE_DONORS)
        
    def test_hydrogen_bond_acceptor_components(self):
        self._do_test_donor_acceptor_components(Hydrogen_bond_acceptor_component_factory(), set(EXPECTED_DIRECT_ACCEPTORS), dict(EXPECTED_INDIRECT_ACCEPTORS), ACCEPTOR, EXPECTED_BACK_BONE_ACCEPTORS)
    
    def test_fast_hydrogen_bond_calculator(self):
        test = Fast_hydrogen_bond_calculator(self.get_single_member_ensemble_simulation())
        format = 'i' * 6
        donor_components = self._build_component_list(Hydrogen_bond_donor_component_factory(),format)
        acceptor_components = self._build_component_list(Hydrogen_bond_acceptor_component_factory(),format)
        parameter_format =  ('i'*4) +('f'*8)
        parameter_components = self._build_component_list(Hydrogen_bond_parameter_factory(), parameter_format)
        donor_index_format = 'i'*3
        donor_index_components = self._build_component_list(Hydrogen_bond_donor_lookup_factory(), donor_index_format)
        acceptor_index_format = 'i'*6
        acceptor_index_components = self._build_component_list(Hydrogen_bond_acceptor_lookup_factory(), acceptor_index_format)
        components = {'DONR' : donor_components.get_native_components(), 
                      'ACCP' : acceptor_components.get_native_components(), 
                      'PARA' : parameter_components.get_native_components(),
                      'DIDX' : donor_index_components.get_native_components(),
                      'AIDX' : acceptor_index_components.get_native_components()
                      }
        
        num_donors = self.donor_indexer.get_max_index()
        num_acceptors = self.acceptor_indexer.get_max_index()
        
        donor_energies = allocate_array(num_donors*3, type='f')
        acceptor_energies = allocate_array(num_acceptors*3,type='f')
        test(components, donor_energies, acceptor_energies)
        
        for donor_selector in EXPECTED_DONOR_ENERGIES.keys():
            offset =  self.donor_indexer.get_index_for_key(donor_selector[:3])*3+donor_selector[3]
            self.assertAlmostEqual(donor_energies[offset]/EXPECTED_DONOR_ENERGIES[donor_selector],1.0,places=4)
            del EXPECTED_DONOR_ENERGIES[donor_selector]
            donor_energies[offset] = 0.0

        for acceptor_selector in EXPECTED_ACCEPTOR_ENERGIES.keys():
            offset =  self.acceptor_indexer.get_index_for_key(acceptor_selector[:3])*3+acceptor_selector[3]
            self.assertAlmostEqual(acceptor_energies[offset]/EXPECTED_ACCEPTOR_ENERGIES[acceptor_selector],1.0,places=4,msg=acceptor_selector)
            del EXPECTED_ACCEPTOR_ENERGIES[acceptor_selector]
            acceptor_energies[offset] = 0.0
    
        self.assertEmpty(EXPECTED_ACCEPTOR_ENERGIES)
        self.assertEmpty(EXPECTED_DONOR_ENERGIES)
        
        self.assertSequenceAlmostEqual(donor_energies, [0.0]* len(donor_energies))
        self.assertSequenceAlmostEqual(acceptor_energies, [0.0]* len(acceptor_energies))
        
    def test_parameter_components(self):
        factory = Hydrogen_bond_parameter_factory()
        format = ('i'*4) +('f'*8)
        component_list = self._build_component_list(factory, format)
        self.assertLength(component_list, 3)
        EXPECTED = [(0, 0, 0, 0, -0.40353099999999997, 13.9148, -0.42172500000000002, -3.0801599999999998, 0.98725099999999999, 0.0, 3.0, -16.5), 
                    (1, 1, 0, 0, 12.151400000000001, -1.1225799999999999, -1.1225799999999999, 1.3095300000000001, 1.3095300000000001, 0.0, 15.5, 3.0), 
                    (2, 2, 0, 0, 15.0486, -1.7946899999999999, -1.7946899999999999, 1.4939800000000001, 1.4939800000000001, 0.0, 15.5, 2.5)]
        
        for i in range(3):
            self.assertSequenceEqual(component_list[i][:4], EXPECTED[i][:4])
            self.assertSequenceAlmostEqual(component_list[i][4:], EXPECTED[i][4:])
        component_list.get_native_components()

    def test_donor_index_components(self):
        factory = Hydrogen_bond_donor_lookup_factory()
        format = ('ii')
        component_list = self._build_component_list(factory, format)
        EXPECTED = ((0,1,0),)
        self.assertLength(component_list, 1)
        self.assertSequenceEqual(component_list[0], EXPECTED[0])

    def test_acceptor_index_components(self):
        factory = Hydrogen_bond_acceptor_lookup_factory()
        format = ('ii')
        component_list = self._build_component_list(factory, format)
        EXPECTED = ((0,0,0,0,1,2),)
        self.assertLength(component_list, 1)
        self.assertSequenceEqual(component_list[0], EXPECTED[0])
        
def run_tests():
    unittest2.main(module='test.test_xcamshift_hbond_ingktlkg')
    
if __name__ == "__main__":
    run_tests()
