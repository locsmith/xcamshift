

from cython.test_distance import get_distance
import unittest2


class Test_cython(unittest2.TestCase):
    
    # TODO add extra places
    DEFAULT_DECIMAL_PLACES = 5
    
    def test_distance(self):
        self.assertAlmostEqual(1.03922182425, get_distance(), self.DEFAULT_DECIMAL_PLACES)

if __name__ ==  '__main__':
    unittest2.main()
     





