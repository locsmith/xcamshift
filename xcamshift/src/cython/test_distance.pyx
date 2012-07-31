
import xcamshift 
from protocol import initStruct
from pdbTool import PDBTool
cimport xplor_access
from xplor_access cimport norm,Vec3,currentSimulation

def get_distance():
    test_1_2 = Test_distance_1_2()
    return test_1_2.get_distance_1_2()

cdef class Test_distance_1_2:

    def setup(self):
        initStruct("test_data/3_ala/3ala.psf")
        PDBTool("test_data/3_ala/3ala.pdb").read()
         
    cpdef get_distance_1_2(self):
        self.setup()
        vec1 = currentSimulation().atomPosArr().data(0)
        vec2 = currentSimulation().atomPosArr().data(1)
        cdef Vec3 result =  vec2 - vec1
        return norm(result)
 


if __name__ == '__main__':
    print 'here'
