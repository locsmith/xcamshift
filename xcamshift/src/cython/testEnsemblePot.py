'''
Created on 8 Jun 2013

@author: garyt
'''
import unittest
from cython.PyEnsemblePotProxy import  PyEnsemblePotProxy
class Test(unittest.TestCase):


    def setUp(self):
        class Dummy():
            pass
        PyEnsemblePotProxy('test','test',None,Dummy())
 
    def tearDown(self):
        pass


    def testEnsemblePot(self):
        pass


if __name__ == "__main__":
    unittest.main()