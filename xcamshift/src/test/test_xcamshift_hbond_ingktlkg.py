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
from table_manager import Table_manager
'''
Created on 31 Dec 2011

@author: garyt
'''
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from xcamshift import Hbond_donor_indexer, Hbond_acceptor_indexer
from cython.fast_segment_manager import Segment_Manager
from utils import Atom_utils
from xcamshift import Hydrogen_bond_context

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

class TestXcamshiftHBondINGKTLKG(unittest2.TestCase):

                 
    def setUp(self):
        initStruct("test_data/ingktlkg_hbond/INGKTLKG.psf")
        PDBTool("test_data/ingktlkg_hbond/INGKTLKG.pdb").read()
        Atom_utils.clear_cache()

        table_manager =  Table_manager.get_default_table_manager()
        self.donor_indexer  = Hbond_donor_indexer(table_manager)
        self.acceptor_indexer  = Hbond_acceptor_indexer(table_manager)
        
        Segment_Manager.reset_segment_manager()
#         print "In method", self._testMethodName

    def test_donor_indexer(self):
        

        donors = [donor for donor in self.donor_indexer.iter()]
        self.assertSequenceEqual(donors, EXPECTED_DONORS)

        
        acceptors = [acceptor for acceptor in self.acceptor_indexer.iter()]
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
        
def run_tests():
    unittest2.main(module='test.test_xcamshift_hbond_ingktlkg')
    
if __name__ == "__main__":
    run_tests()
