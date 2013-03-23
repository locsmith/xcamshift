import sys
import cython.cyfwk as cyfwk
import unittest

#
class PyEnsemble(cyfwk.AlgBase):
    def __init__(self, potential_name,instance_name,simulation): 
        super(PyEnsemble, self).__init__(potential_name,instance_name,simulation)
        self.test_value =  False
    
    def virtual_call(self):
        self.test_value = True


class Test_ensemble_callback(unittest.TestCase):
    
    def setUp(self):
        self._ensemble =  PyEnsemble('test','test2',None)

    def test_virtual_call(self):
        self._ensemble.doRun()
        self.assertTrue(self._ensemble.test_value, msg="virtual method call did not complete")
        
if __name__ == '__main__':
    unittest.main()
    
 

