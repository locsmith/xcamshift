'''
Created on 30 Dec 2011

@author: garyt
'''
import unittest
from Table_manager import Table_manager


class Test(unittest.TestCase):


    def testLoadTable(self):
        table_manager = Table_manager()
        table_manager.add_search_path('../../data')
        table_manager.get_table('bb','ala')


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testLoadTable']
    unittest.main()