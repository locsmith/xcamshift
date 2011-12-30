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
        table = self.table_manager.get_table('bb','ala')
        self.assertTrue(table != None)

    def testLoadTableValues(self):
        table = self.table_manager.get_table('bb','glu') 
        self.assertAlmostEqual(table['exponent'],1.0)
        self.assertAlmostEqual(table['data'][-1]['N']['HA'],0.10495380090732243)


if __name__ == "__main__":
    unittest2.main()