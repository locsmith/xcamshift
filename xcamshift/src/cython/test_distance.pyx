
import xcamshift 
from protocol import initStruct
from pdbTool import PDBTool
cimport xplor_access
from xplor_access cimport norm,Vec3,currentSimulation

def get_distance():
    test_1_2 = Test_distance_1_2()
    return test_1_2.get_distance_1_2()

cdef struct ring_atom_positions:
    int num_atoms
    Vec3 atom_pos_0
    Vec3 atom_pos_1
    Vec3 atom_pos_2
    Vec3 atom_pos_3
    Vec3 atom_pos_4
    Vec3 atom_pos_5
    
cdef Vec3* get_pos_pointer(ring_atom_positions&  positions):
    return &positions.atom_pos_0

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

cpdef test_vec3_pointer():
    cdef ring_atom_positions test
    test.atom_pos_0 =  Vec3(1,2,3)
    cdef Vec3* pos
    pos  = get_pos_pointer(test)
    return pos.x(),pos.y(),pos.z()


if __name__ == '__main__':
    print 'here'
