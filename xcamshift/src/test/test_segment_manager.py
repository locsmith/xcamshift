'''
Created on 31 Dec 2011

@author: garyt
'''
fast = True
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from utils import Atom_utils

if fast:
    from cython.fast_segment_manager import Segment_Manager
else:
    from segment_manager import Segment_Manager
    
class Test_segment_manager(unittest2.TestCase):


    def setUp(self):
        initStruct("test_data/3_ala/3ala.psf")
        PDBTool("test_data/3_ala/3ala.pdb").read()
        Segment_Manager.reset_segment_manager()
        Atom_utils.clear_cache()

    def testSegments(self):
        segments = Segment_Manager().get_segments()
        self.assertSequenceEqual(segments, ("",))
        
    def testSegmentInfo(self):
        segment_manager = Segment_Manager()
        segments = segment_manager.get_segments()
        
        segment_info = segment_manager.get_segment_info(segments[0])

        self.assertEqual(segment_info.first_residue, 1)
        self.assertEqual(segment_info.last_residue, 3)
        
        self.assertEqual(segment_info.first_atom_index, 0)
        self.assertEqual(segment_info.last_atom_index, 32)
        
        self.assertEqual(segment_info.segment_length, 3)


if __name__ == "__main__":
    unittest2.main()
#    test = testSegmentManager()
#    test.setUp()
#    test.testSegments()