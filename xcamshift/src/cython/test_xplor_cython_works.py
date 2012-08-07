

from cython.test_distance import get_distance, test_vec3_pointer
import unittest2


class Test_cython(unittest2.TestCase):
    
    # TODO add extra places
    DEFAULT_DECIMAL_PLACES = 5
    
    def test_distance(self):
        self.assertAlmostEqual(1.03922182425, get_distance(), self.DEFAULT_DECIMAL_PLACES)
        
    def test_vec3_pointer(self):
        self.assertSequenceEqual(test_vec3_pointer(), (1.0,2.0,3.0))

if __name__ ==  '__main__':
    unittest2.main()
     





