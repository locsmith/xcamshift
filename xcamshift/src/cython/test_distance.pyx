#-------------------------------------------------------------------------------
# Copyright (c) 2013 gary thompson.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
#import xcamshift 
from protocol import initStruct
from pdbTool import PDBTool
cimport xplor_access
from xplor_access cimport norm,Vec3,currentSimulation,Coord_holder

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
        cdef Coord_holder holder = Coord_holder(currentSimulation()[0])
        vec1 = holder.getPos(0)
        vec2 = holder.getPos(1)
#        vec1 = currentSimulation().atomPosArr().data(0)
#        vec2 = currentSimulation().atomPosArr().data(1)
        cdef Vec3 result =  vec2 - vec1
        print vec1.x()
        return norm(result)

cpdef test_vec3_pointer():
    cdef ring_atom_positions test
    test.atom_pos_0 =  Vec3(1,2,3)
    cdef Vec3* pos
    pos  = get_pos_pointer(test)
    return pos.x(),pos.y(),pos.z()


if __name__ == '__main__':
    print 'here'
