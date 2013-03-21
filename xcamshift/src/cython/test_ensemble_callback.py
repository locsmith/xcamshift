import sys
import cython.cyfwk as cyfwk
import unittest

#cyfwk.AlgBase
class PyEnsemble():
    def __init__(self): 
        self.test_value =  False
        
    def run(self):
        print '[py]  .run()'


class Test_ensemble_callback(unittest.TestCase):
    
    def setUp(self):
        self._ensemble =  PyEnsemble()

    def test_virtual_call(self):
        self._ensemble.run()
        self.assertTrue(self._ensemble.test_value, msg="virtual method call did not complete")
if __name__ == '__main__':
    unittest.main()
    
 

