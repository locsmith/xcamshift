'''
Created on 28 Jan 2014

@author: garyt
'''

from xcamshift import Xcamshift
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
import io
from test.gb3 import gb3_shifts
from utils import Atom_utils
from cython.fast_segment_manager import Segment_Manager



class TestXcamshiftPrintShifts(unittest2.TestCase):

    def _clear_caches(self):
        Atom_utils.clear_cache()
        Segment_Manager.reset_segment_manager()

    def setUp(self):
        initStruct("test_data/gb3/gb3.psf")
        PDBTool("test_data/gb3/gb3_refined_II.pdb").read()
        self._clear_caches()



    def _get_xcamshift(self):
        xcamshift = Xcamshift()
        return xcamshift

    def test_calc_shifts(self):
        class WritableObject:
            def __init__(self):
                self.content = []

            def write(self, string):
                self.content.append(string)


        output = WritableObject()

        xcamshift =  self._get_xcamshift()
        xcamshift.remove_named_sub_potential('HBOND', quiet=True)
        xcamshift.print_shifts(output)

        #print ''.join(output.content)
        contents = output.content
        contents  =  ''.join(contents)
        for line in contents.split('\n'):
            line = line.strip()
            if len(line) == 0:
                pass
            elif line.startswith('SEGID'):
                fields =  line.split()
                names =  fields [3:]
            else:

                if line[0] == '|':
                    line=line[1:]
                field_1,rest = line.split('|')
                fields = [field_1]
                fields.extend(rest.split())
                segid,resid = fields[:2]
                resid=int(resid)

                for name,shift in zip(names,fields[3:]):
                    self.assertEqual(segid, '')
                    key = ('',resid,name)
                    if shift == '.':
                        self.assertFalse(key in gb3_shifts)
                    else:
                        shift=float(shift)
                        #TODO test data is not accurate enough to test within 3dp.due to rounding errors..
                        self.assertAlmostEqual(shift,gb3_shifts[key],places=2,msg=`key`)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest2.main()