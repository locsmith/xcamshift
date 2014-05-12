#-------------------------------------------------------------------------------
# Copyright (c) 2013 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
#
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------

# cython: profile=False
# cython: boundscheck=False
# cython: wraparound=True
# cython: cdivision=True

'''
Created on 31 Jul 2012

@author: garyt

'''
cdef extern from "instantiate.hh":
    pass

from derivList import DerivList as PyDerivList
from math import ceil
cimport cython
from vec3 import Vec3 as python_vec3
from common_constants import TARGET_ATOM_IDS_CHANGED, STRUCTURE_CHANGED
from xplor_access cimport norm,Vec3, Dihedral, Atom,  dot,  cross,  Simulation, CDSVector, DerivList,\
     clearSharedVector, EnsembleSimulation, createSharedVector, getSharedVectorValue, resizeSharedVector,\
     addToSharedVectorValue, getSharedVectorValue, deleteSharedVector, setSharedVectorValue, BondAngle, \
     CDSList
from libc.math cimport cos,sin,  fabs, tanh, pow, cosh, ceil
from libc.stdlib cimport malloc, free
from libc.string cimport strcmp
from time import time
from utils import Atom_utils
from cpython cimport array
import ctypes
from component_list import  Component_list, Native_component_list
from ensembleSimulation import EnsembleSimulation as pyEnsembleSimulation
import math
cdef double PI = math.pi

#TODO: should be ensembleSimulationAsNative
cdef inline EnsembleSimulation* simulationAsNative(object simulation) except NULL:
    #TODO tell charles his __eq__ method has problems
    if id(None) == id(simulation):
        raise Exception("_simulation cannot be None!")
    if not isinstance(simulation, pyEnsembleSimulation):
        raise Exception("a %s can't be converted to a EnsembleSimulation*" % simulation.__class__.__name__)
    return <EnsembleSimulation*><size_t>int(simulation.this)

cpdef array.array allocate_array(int len, type='d'):
    result = array.array(type,[0])
    array.resize(result, len)
    array.zero(result)
    return result

cpdef zero_array(array.array in_array):
     array.zero(in_array)

cpdef resize_array(array.array in_array, int len):
    array.resize(in_array, len)

cdef struct Random_coil_component:
    int target_atom
    float shift

cdef struct Nonbonded_coefficient_component:
    int chem_type_id
    int sphere_id
    float exponent
    float[3*7] coefficients

cdef struct Distance_component:
      int target_atom
      int remote_atom_1
      int remote_atom_2
      float coefficient
      float exponent

cdef struct Dihedral_component:
      int target_atom
      int[4] dihedral_atoms
      float coefficient
      float[6] parameters
      #TODO: there seems to be an extra float here...

#TODO: same as sing target component unite and rename
cdef struct Non_bonded_target_component:
      int target_atom_id
      int atom_type_id

cdef struct Non_bonded_remote_atom_component:
      int remote_atom_id
      int chem_type[2] # one for each sphere

cdef struct Component_index_pair:
    int target_atom_id
    int target_index
    int remote_index
    int component_index

cdef struct Constant_cache:
    int     target_atom_id
    float   flat_bottom_shift_limit
    float   end_harmonic
    float   scale_harmonic
    float   weight
    float   tanh_amplitude
    float   tanh_elongation
    float   tanh_y_offset

cdef struct Hydrogen_bond_component:
    int index
    int direct_atom_id
    int indirect_atom_id
    int type
    int atom_type_id
    int backbone

cdef struct Hydrogen_bond_parameter:
    int component_index
    int param_index
    int donor_index
    int acceptor_index
    float[6] p
    float r
    float s

cdef struct Hydrogen_bond_donor_lookup:
    int donor_id
    int number_acceptors
    int[1] acceptor_offset

cdef struct Hydrogen_bond_acceptor_lookup:
    int index
    int donor_id
    int acceptor_id
    int[3] parameter_indices

cdef struct Hydrogen_bond_shift_component:
    int target_atom_index
    int hbond_index
    double dist_coef
    double ang1_coef
    double ang2_coef

def test_dump_component_index_pair(Non_bonded_interaction_list data, int index):
    cdef Component_index_pair* result =  data.get(index)

    return result[0].target_atom_id, result[0].target_index, result[0].remote_index, result[0].component_index

def test_dump_component_offsets(data):
    cdef size_t  test  = ctypes.addressof(data)
    cdef Component_Offsets* test2 = <Component_Offsets*> test
    cdef Component_Offsets[:] dummy_view
    if len(data) == 0:
        result = ()
    else:
        dummy_view = <Component_Offsets[:len(data)/ sizeof(Component_Offsets)]> &test2[0]

        result_data = []
        for i in range(len(data)/ sizeof(Component_Offsets)):
            result_data.append((dummy_view[i].id,  dummy_view[i].offset,  dummy_view[i].length))
        result = tuple(result_data)
    return result


def test_dump_dist_comp(data):
    cdef size_t  test  = ctypes.addressof(data)
    cdef Distance_component* test2 = <Distance_component*> test
    cdef Distance_component[:] dummy_view
    if len(data) == 0:
        result = ()
    else:
        dummy_view = <Distance_component[:len(data)/ sizeof(Distance_component)]> &test2[0]

        result_data = []
        for i in range(len(data)/ sizeof(Distance_component)):
            result_data.append((dummy_view[i].target_atom,  dummy_view[i].remote_atom_1,  dummy_view[i].remote_atom_2,   dummy_view[i].coefficient, dummy_view[i].exponent))
        result = tuple(result_data)
    return result

def test_dump_dihedral_comp(data):
    cdef size_t  test  = ctypes.addressof(data)

    compiled_components =  <Dihedral_component*> <size_t> ctypes.addressof(data)
    num_components =  len(data)/ sizeof(Dihedral_component)
    print 'dump num comp',num_components, test

    for i in range(num_components):

        print 'target_atom',i, compiled_components[0].target_atom

        print 'dihedral_atoms',i, compiled_components[i].dihedral_atoms[0],compiled_components[i].dihedral_atoms[1],\
                               compiled_components[i].dihedral_atoms[2], compiled_components[i].dihedral_atoms[4]
        print 'coeff',i, compiled_components[i].coefficient
        print 'param',i, compiled_components[i].parameters[0], compiled_components[i].parameters[1],\
                      compiled_components[i].parameters[2], compiled_components[i].parameters[3],\
                      compiled_components[i].parameters[4], compiled_components[i].parameters[5],\
                      compiled_components[i].parameters[6]

cdef class CDSSharedVectorFloat:
    cdef void *  data
    cdef int size
    cdef object _simulation
    cdef int _ensemble_size

    def __init__(self, int size, object ensembleSimulation):
        self._simulation = ensembleSimulation
        self._ensemble_size = ensembleSimulation.size()

        cdef EnsembleSimulation* cEnsembleSimulation = simulationAsNative(ensembleSimulation)
        self.data = <CDSVector[double]*>createSharedVector(size*self._ensemble_size, 0.0,  cEnsembleSimulation)
        self.size = size

    cpdef resize(self,int size):
        self.size = size * self._ensemble_size
        resizeSharedVector(self.data,self.size)
        return self

    cpdef clear(self):
        clearSharedVector(self.data)
        return self

    def __len__(self):
        return self.size* self._ensemble_size


    cdef inline int array_offset(self,int index, int ensemble_index) nogil:
        return (self._ensemble_size * index) + ensemble_index

    cdef inline int ensemble_member_index(self) nogil:
        cdef int member_index
        with gil:
            member_index  = self._simulation[0].member()[0].memberIndex()
        return member_index

    cdef inline int ensemble_array_offset(self,int index) nogil:
        return self.array_offset(index, self.ensemble_member_index())

    def __str__(self):
        data = []
        for i in range(len(self)):
            value  = getSharedVectorValue(self.data,i)
            data.append('%7.3f' % value)
        return 'esim size: ' +`self._ensemble_size` + ', [' + ', '.join(data) + ']'


    def average_into(self, CDSVectorFloat in_data):
        in_data.clear()
        cdef CDSVector[double] *target =  in_data.data
        in_data.resize(self.size)
        cdef int i,j
        cdef double divisor = 1.0/self._ensemble_size
        for i in range(self.size):
            target[0][i]=0.0
            for j in range(self._ensemble_size):
                target[0][i] += getSharedVectorValue(self.data,self.array_offset(i,j)) * divisor

    cdef inline void addToResult(self, int offset, double value) nogil:
        addToSharedVectorValue(self.data,offset,value)

    def __getitem__(self, int key):
        if key >= len(self):
            raise IndexError("index (%i) out of range (%i)" % (key, len(self)))
        return getSharedVectorValue(self.data,key)

    def __setitem__(self, int key,value):
        if key >= self.size:
            raise IndexError("index (%i) out of range (%i)" % (key, len(self)))
        setSharedVectorValue(self.data,key,value)

    def __del__(self):
        deleteSharedVector(self.data)
        self.data = NULL

cdef class CDSVectorFloat:
    cdef CDSVector[double] *data


    def __init__(self, int size):
        self.data = new CDSVector[double]()
        self.data[0].resize(size)

    cpdef resize(self,int size):
        self.data[0].resize(size)

    cpdef clear(self):
        self.data[0].set(0.0)

    def __len__(self):
        return self.data[0].size()

    cdef CDSVector[double]* get_data(self):
        return self.data

    def __str__(self):
        result = []
        for i in range(self.data[0].size()):
            result.append('%7.3f' % self.data[0][i])
        return '[' + ', '.join(result) + ']'

    def __getitem__(self, int key):
        if key >= self.data[0].size():
            raise IndexError("index (%i) out of range (%i)" % (key, self.data.size()))
        result = self.data[0][key]
        return result

    def __del__(self):
        del self.data





cdef  class Non_bonded_interaction_list:
    cdef CDSVector[int]  *data
    cdef int length
    cdef int size_increment
    cdef int RECORD_LENGTH
    cdef int resize_count

    def __cinit__(self, int length=0, double fill_factor=1.0):
        self.data = new CDSVector[int]()
        self.data[0].resize(<int>ceil(length*fill_factor*2))
        self.length =  0
        self.size_increment =  20
        self.RECORD_LENGTH = 4
        self.resize_count = 0

    def test_append(self,int target_atom_id, int target_id, int remote_id, component_index):
        self.append(target_atom_id,  target_id,  remote_id, component_index)

    def get_resize_count(self):
        return self.resize_count

    def clear(self):
        self.length = 0
        self.resize_count = 0

    cdef inline void append(self, int target_atom_id, int target_id, int remote_id, int component_index):

        if  self.data[0].size()*self.RECORD_LENGTH <= self.length*self.RECORD_LENGTH:
            self.resize_count+=1
            self.resize()
        self.data[0][self.length] = target_atom_id
        self.data[0][self.length+1] = target_id
        self.data[0][self.length+2] = remote_id
        self.data[0][self.length+3] = component_index

        self.length += self.RECORD_LENGTH


    cdef inline void resize(self):
        self.data[0].resize(self.data[0].size()+self.size_increment*self.RECORD_LENGTH)

    cdef inline Component_index_pair* get(self,int offset):
        return <Component_index_pair *> &self.data[0][offset*self.RECORD_LENGTH]

    def get_allocation(self):
        return self.data[0].size()

    def get_size_increment(self):
        return self.size_increment

    def __len__(self):
        return self.length/self.RECORD_LENGTH

    def __getitem__(self, int key):
        if key >= self.length/self.RECORD_LENGTH:
            raise IndexError("index (%i) out of range (%i)" % (key, self.length/self.RECORD_LENGTH))
        result  = []
        for i in range(self.RECORD_LENGTH):
            result.append(self.data[0][key*self.RECORD_LENGTH+i])
        return tuple(result)

    def build_selection_list(self, accept = lambda x: True):
        result = []
        for i,elem in enumerate(self):
            if accept(elem):
                result.append(i)
        return array.array('i',result)


    def get_all_components(self):
        result = []
        for i in range(self.length/self.RECORD_LENGTH):
            result.append(self[i])
        return result

cdef class Vec3_list:
    cdef CDSVector[Vec3] *data

    def __cinit__(self):
         self.data = new CDSVector[Vec3]()


    cdef set_length(self,int length):
        self.data.resize(length)

    @cython.profile(False)
    cdef inline Vec3*  get(self, int offset):
        return &self.data[0][offset]

    @cython.profile(False)
    cdef inline set(self, int x, Vec3& y):
        self.data[0][x] =  y

    def __len__(self):
        return self.data[0].size()

    def __dealloc__(self):

        if self.data != NULL:
            del self.data


    def __iter__(self):
        cdef Vec3 vec3



        for i in range(self.data.size()):
            vec3 = self.get(i)[0]
            yield (i,python_vec3(vec3[0],vec3[1],vec3[2]))


cdef class Out_array:
    cdef long _length
    cdef double[60000] _data
    cdef int[60000] _mask
    cdef EnsembleSimulation* _simulation

    def __init__ (self, length, simulation):
        self._length = length
        self._simulation = simulationAsNative(simulation)


    def __cinit__(self,  length, simulation):
        self._length =  60000

        if length > 60000:
            raise MemoryError("couldn't allocate array beyond length %i" % 60000)

    cpdef get_length(self):
        return self._length

    cpdef add_all(self, Out_array target_array):
        cdef Vec3 value
        for i in range(target_array._length):
            if i < self._length and target_array._mask[i]:
                value = Vec3(target_array._data[i*3+0],target_array._data[i*3+1],target_array._data[i*3+2])
                self.add(i, value)

    def __dealloc__(self):
        pass
#        free(self._data)
#        self._data = NULL
#        free(self._mask)
#        self._mask = NULL

    cdef inline bint _check_clear_length(self, long length) except -1:
        if not length <= self._length:
            raise IndexError("tried to clear to %i length is %i" % (length,self._length))
        if self._data is NULL:
            raise AttributeError("trying to access deallocated memory!")

    cdef inline bint _check_offset(self, long offset) except -1:
        if not offset < self._length:
            raise IndexError("tried to access object at %i length is %i" % (offset,self._length))
        if self._data is NULL:
            raise AttributeError("trying to access deallocated memory!")

    cpdef bint realloc(self,long length) except -1:
        if length <60000:
            self._clear(length)
        else:
            raise MemoryError("couldn't beyond of length %i" % 60000)

    cdef inline bint _check_state(self) except -1:
        if self._data is NULL:
            raise MemoryError("couldn't reallocate array of length %i" % self._length)

    cdef _clear(self, length):
        self._check_clear_length(length)
        for i in range(length*3):
            self._data[i] = 0.0
        for i in range(length):
            self._mask[i]=0

    def padd(self,offset,value):
        cdef Vec3 cvalue  = Vec3(value[0],value[1],value[2])
        self.add(offset, cvalue)

    @cython.profile(False)
    cdef inline void  add(self,long offset, Vec3& value):
        cdef long id3 = offset *3

        self._data[id3]   += value[0]
        self._data[id3+1] += value[1]
        self._data[id3+2] += value[2]

        self._mask[offset] = 1

    def clear_by_id(self,target_ids):
        cdef long id3 = 0
        for target_id in target_ids:
            self._check_offset(target_id)
            id3 = target_id *3

            self._data[id3]   = 0.0
            self._data[id3+1] = 0.0
            self._data[id3+2] = 0.0

            self._mask[target_id] = 0

    def add_forces_to_result(self, result=None, weight=1.0):
        cdef long id3
        cdef DerivList *derivList
        cdef CDSVector[Vec3]* derivs
        cdef Vec3* test

        if isinstance(result,PyDerivList):
            pointer = int(result.this)
            derivList = (<DerivList*><size_t>pointer)
            #TODO: references don't work here lots of problems -> cython mailing list
            derivs  = &derivList[0][self._simulation]
            for i in range(self._length):
                if self._mask[i] != 0:
                    id3 = i * 3

                    derivs[0][i][0]  +=  self._data[id3] * weight
                    derivs[0][i][1]  +=  self._data[id3+1] *weight
                    derivs[0][i][2]  +=  self._data[id3+2] * weight


        elif result == None:
            if result ==  None:
                result = [None] * self._length
            for i in range(self._length):
                if self._mask[i] != 0:
                    if result[i] ==  None:
                        result[i] = [0.0]*3
                    id3 = i * 3
                    result[i][0] += self._data[id3] * weight
                    result[i][1] += self._data[id3+1] *weight
                    result[i][2] += self._data[id3+2] * weight
        else:
            # raw pointer
            pointer = int(result)
            derivList = (<DerivList*><size_t>pointer)
            #TODO: references don't work here lots of problems -> cython mailing list
            derivs  = &derivList[0][self._simulation]
            for i in range(self._length):
                if self._mask[i] != 0:
                    id3 = i * 3

                    derivs[0][i][0]  +=  self._data[id3] * weight
                    derivs[0][i][1]  +=  self._data[id3+1] *weight
                    derivs[0][i][2]  +=  self._data[id3+2] * weight

        return result

    def add_forces_to_result_by_id(self, list target_atom_ids, list result=None):
        if result ==  None:
            result = [None] * self._length

        cdef double x = 0.0
        cdef double y = 0.0
        cdef double z = 0.0
        cdef long id3
        for i in target_atom_ids:
            if self._mask[i] != 0:
                if result[i] ==  None:
                    result[i] = [0.0]*3
                id3 = i * 3
                x  = self._data[id3]
                y  = self._data[id3+1]
                z  = self._data[id3+2]
                result[i][0] += x
                result[i][1] += y
                result[i][2] += z


        return result

    def __str__(self):

        cdef long id3 =0
        result =  []
        result.append('length = %i' % self._length)
        for i in range(self._length):
            if self._mask[i] ==0:
                result.append("%5i ." % i)
            else:
                id3 = i * 3
                result.append('%5i %4.7f,%4.7f,%4.7f' % (i,self._data[id3],self._data[id3+1],self._data[id3+2]))
        return '\n'.join(result)

cdef object vec3_as_tuple(Vec3& vec_3):
    return vec_3.x(), vec_3.y(), vec_3.z()




cdef struct target_distant_atom:
    int target_atom_id
    int distant_atom_id

cdef struct target_type:
    int  target_atom_id
    int  atom_type_id

cdef struct coefficient_exponent:
    float coefficient
    float exponent

cdef struct dihedral_ids:
    int atom_id_1
    int atom_id_2
    int atom_id_3
    int atom_id_4

cdef struct bond_angle_ids:
    int atom_id_1
    int atom_id_2
    int atom_id_3


cdef struct Coef_component:
    int atom_type_id
    int ring_id
    float coefficient

cdef struct Ring_force_sub_terms:
    Vec3 d
    float dLnL
    float dL3nL3
    float dn
    float dL
    float nL
    Vec3 gradUQ
    float dL3
    float u
    Vec3 gradVQ
    float dL6
    Vec3 ring_normal

cdef class Python_ring_force_sub_terms:
    cdef Ring_force_sub_terms _terms


    def __init__(self):
        pass

    cdef setup(self, Ring_force_sub_terms& terms):
        self._terms=terms

    cdef Ring_force_sub_terms get_terms(self):
        return self._terms


cdef struct Component_Offsets:
    int id
    int offset
    int length


cdef class Coef_components:
    cdef object _raw_data
    cdef Coef_component* _components
    cdef int _num_components

    cdef object _raw_component_offsets_data
    cdef Component_Offsets* _component_offsets
    cdef int num_ids

    def _bytes_to_components(self, data):

        self._raw_data =  data
        self._components =  <Coef_component*> <size_t> ctypes.addressof(data)
        self._num_components =  len(data)/ sizeof(Coef_component)

    cdef void _bytes_to_component_offsets(self, data):

        self._raw_component_offsets_data =  data
        self._component_offsets =  <Component_Offsets*> <size_t> ctypes.addressof(data)
        self.num_ids =  len(data)/ sizeof(Component_Offsets)

    def __str__(self):
        result = []

        result.append('components')
        result.append('----------')
        result.append('')
        result.append('number = %i' % self._num_components)
        for i in range(self._num_components):
            result.append('%i : %i %i %f' % (i,self._components[i].atom_type_id, self._components[i].ring_id,self._components[i].coefficient))
        result.append('')
        result.append('offsets')
        result.append('-------')
        result.append('')
        result.append('number = %i' % self.num_ids)
        for i in range(self.num_ids):
            result.append('%i : %i %i %i' % (i,self._component_offsets[i].id, self._component_offsets[i].offset,self._component_offsets[i].length))

        return '\n'.join(result)

    def __cinit__(self, object coef_components, component_offsets):


        self._bytes_to_components(coef_components)
        self._bytes_to_component_offsets(component_offsets)



    def __init__ (self, object coef_components, object components):
        pass


    cdef inline Component_Offsets* get_id_offsets(self,int id):
        return &self._component_offsets[id]

    def __dealloc__(self):
        pass


    @cython.profile(False)
    cdef inline Coef_component* get_component(self,int offset):
#        self._check_offset(offset)
        return &self._components[offset]

cdef struct ring_atom_ids:
    int num_atoms
    int[6] atom_ids






cdef struct dihedral_parameters:
    float param_0
    float param_1
    float param_2
    float param_3
    float param_4



@cython.profile(False)
cdef inline void operator_times (Vec3& vec3, float scale):
     vec3[0] =  vec3[0] * scale
     vec3[1] =  vec3[1] * scale
     vec3[2] =  vec3[2] * scale

@cython.profile(False)
cdef inline float calc_distance_simulation(Simulation* sim, int atom_index_1, int atom_index_2):

    cdef Vec3 vec1 = sim[0].atomPos(atom_index_1)
    cdef Vec3 vec2 = sim[0].atomPos(atom_index_2)
    cdef Vec3 result =  vec2 - vec1
    return norm(result)


#TODO how does currentSimulation().atomByID work inside
@cython.profile(False)
cdef inline float calc_dihedral_angle_simulation(Simulation* simulation, dihedral_ids dihedral_atom_ids):



    cdef Atom atom_1  = simulation[0].atomByID(dihedral_atom_ids.atom_id_1)
    cdef Atom atom_2  = simulation[0].atomByID(dihedral_atom_ids.atom_id_2)
    cdef Atom atom_3  = simulation[0].atomByID(dihedral_atom_ids.atom_id_3)
    cdef Atom atom_4  = simulation[0].atomByID(dihedral_atom_ids.atom_id_4)

    cdef double result = Dihedral(atom_1,atom_2,atom_3,atom_4).value()
    return result

@cython.profile(False)
cdef inline float calc_bond_angle_simulation(Simulation* simulation, bond_angle_ids bond_angle_atom_ids) nogil:


    cdef Atom atom_1  = simulation[0].atomByID(bond_angle_atom_ids.atom_id_1)
    cdef Atom atom_2  = simulation[0].atomByID(bond_angle_atom_ids.atom_id_2)
    cdef Atom atom_3  = simulation[0].atomByID(bond_angle_atom_ids.atom_id_3)

    cdef double result = BondAngle(atom_1,atom_2,atom_3).value()
    return result


cdef float DEFAULT_CUTOFF = 5.0
cdef float DEFAULT_SMOOTHING_FACTOR = 1.0
cdef float DEFAULT_NB_CUTOFF = 5.0

cdef class Base_shift_calculator:
    cdef bint _verbose
    cdef str _name
    cdef EnsembleSimulation* _simulation

    def __init__(self, simulation, str name):
            self._verbose = False
            self._name = name
            self._simulation =  simulationAsNative(simulation)



    def set_verbose(self,bint state):
        self._verbose =  state



    def _set_components(self,components):
        self._bytes_to_components(components['ATOM'])
        self._simulation = <EnsembleSimulation*><size_t>int(components['SIMU'].this)


    cdef inline int ensemble_size(self) nogil:
        return self._simulation[0].size()

    cdef inline int ensemble_member_index(self) nogil:
        return self._simulation[0].member()[0].memberIndex()

    cdef inline int ensemble_array_offset(self,int index) nogil:
        return (self.ensemble_size() * index) + self.ensemble_member_index()


cdef class Fast_random_coil_shift_calculator(Base_shift_calculator):

    cdef Random_coil_component* _compiled_components
    cdef int _num_components
    cdef object _raw_data

    def __cinit__(self, simulation,  str name = "not set"):
        self._raw_data = None
        self._compiled_components  = NULL
        self._num_components =  0

    def __init__(self, simulation, str name = "not set"):
        Base_shift_calculator.__init__(self, simulation, name)

#
#    TODO: this needs to be removed
    def _bytes_to_components(self, data):

        self._raw_data =  data
        self._compiled_components =  <Random_coil_component*> <size_t> ctypes.addressof(data)
        self._num_components =  len(data)/ sizeof(Random_coil_component)


    @cython.profile(False)
    def __call__(self, object components, CDSSharedVectorFloat shift_cache, int[:] component_to_target,  int[:] active_components):
        self._set_components(components)
        cdef double start_time = 0.0
        cdef double end_time = 0.0

        if self._verbose:
            start_time=time()

        cdef int factor_index = 0
        cdef int component_index
        cdef int offset

        #TODO: note to cython list for componnt_index in active_components produces awful code!
        for factor_index in range(len(active_components)):
            component_index = active_components[factor_index]

            offset = self.ensemble_array_offset(component_to_target[factor_index])
            shift_cache.addToResult(offset,self._compiled_components[component_index].shift)

        if self._verbose:
            end_time = time()
            print '   distance shift components ' ,self._name,len(components), 'in', "%.17g" % (end_time-start_time), "seconds"

cdef class Fast_distance_shift_calculator(Base_shift_calculator):


    cdef int _target_atom_index
    cdef int _distance_atom_index_1
    cdef int _distance_atom_index_2
    cdef int _exponent_index
    cdef int _coefficient_index
    cdef bint _smoothed
    cdef float _smoothing_factor
    cdef float _cutoff

    cdef Distance_component* _compiled_components
    cdef int _num_components
    cdef object _raw_data

    def __cinit__(self, simulation, smoothed, str name = "not set"):
        self._raw_data = None
        self._compiled_components  = NULL
        self._num_components =  0

    def __init__(self, simulation, smoothed, str name = "not set"):
        Base_shift_calculator.__init__(self, simulation, name)

        self._smoothed =  smoothed
        self._smoothing_factor =  DEFAULT_SMOOTHING_FACTOR
#
        self._cutoff =  DEFAULT_CUTOFF

#    TODO: this needs to be removed
    def _bytes_to_components(self, data):

        self._raw_data =  data
        self._compiled_components =  <Distance_component*> <size_t> ctypes.addressof(data)
        self._num_components =  len(data)/ sizeof(Distance_component)





    @cython.profile(False)
    cdef inline target_distant_atom _get_target_and_distant_atom_ids(self, int index):
        cdef target_distant_atom result

        result.target_atom_id = self._compiled_components[index].remote_atom_1
        result.distant_atom_id  = self._compiled_components[index].remote_atom_2
        return result

    @cython.profile(False)
    cdef inline coefficient_exponent _get_coefficient_and_exponent(self, int index):
        cdef coefficient_exponent result


        result.coefficient = self._compiled_components[index].coefficient
        result.exponent = self._compiled_components[index].exponent

        return result

    @cython.profile(False)
    def __call__(self, object components, CDSSharedVectorFloat shift_cache, int[:] component_to_target,  int[:] active_components):

        cdef double start_time = 0.0
        cdef double end_time = 0.0

        if self._verbose:
            start_time=time()

        self._set_components(components)
        cdef float smoothing_factor = self._smoothing_factor
        cdef float ratio
        cdef float result
        cdef coefficient_exponent coef_exp
        cdef float coefficent
        cdef float exponent
        cdef object component
        cdef int target_atom_id
        cdef int distant_atom_id
        cdef float distance



        if self._verbose:
            start_time = time()

        cdef int factor_index = 0
        cdef int component_index
        cdef int offset

        #TODO: note to cython list for componnt_index in active_components produces awful code!
        for factor_index in range(len(active_components)):
            component_index = active_components[factor_index]

            target_atom_id = self._compiled_components[component_index].remote_atom_1
            distant_atom_id  = self._compiled_components[component_index].remote_atom_2


            coef_exp = self._get_coefficient_and_exponent(component_index)
            distance =calc_distance_simulation(self._simulation, target_atom_id, distant_atom_id)

            if self._smoothed:
                ratio = distance / self._cutoff
                smoothing_factor = 1.0 - ratio ** 8

            offset = self.ensemble_array_offset(component_to_target[factor_index])
            shift_cache.addToResult(offset, smoothing_factor * pow(distance,  coef_exp.exponent) * coef_exp.coefficient)

        if self._verbose:
            end_time = time()
            print '   distance shift components ' ,self._name,len(components), 'in', "%.17g" % (end_time-start_time), "seconds"




cdef class Fast_dihedral_shift_calculator(Base_shift_calculator):

    cdef Dihedral_component *_compiled_components
    cdef int _num_components
    cdef object _raw_data

    def __cinit__(self, simulation, str name = "not set"):
        self._raw_data =  None
        self._compiled_components = NULL
        self._num_components = 0


    def __init__(self, simulation, str name = "not set"):
        Base_shift_calculator.__init__(self, simulation, name)


    def  _bytes_to_components(self, data):
        self._raw_data =  data

        self._compiled_components =  <Dihedral_component*> <size_t> ctypes.addressof(data)
        self._num_components =  len(data)/ sizeof(Dihedral_component)


    cdef inline _get_component(self,int index):
        return self._components[index]

    @cython.profile(False)
    cdef inline dihedral_ids _get_dihedral_atom_ids(self, int index):
        cdef dihedral_ids result

        result.atom_id_1 = self._compiled_components[index].dihedral_atoms[0]
        result.atom_id_2 = self._compiled_components[index].dihedral_atoms[1]
        result.atom_id_3 = self._compiled_components[index].dihedral_atoms[2]
        result.atom_id_4 = self._compiled_components[index].dihedral_atoms[3]

        return result

    cdef inline dihedral_parameters _get_parameters(self, int index):

        cdef dihedral_parameters result

        result.param_0 = self._compiled_components[index].parameters[0]
        result.param_1 = self._compiled_components[index].parameters[3]
        result.param_2 = self._compiled_components[index].parameters[1]
        result.param_3 = self._compiled_components[index].parameters[4]
        result.param_4 = self._compiled_components[index].parameters[2]




        return result

    @cython.profile(False)
    cdef inline float _get_coefficient(self, int index):
        return self._compiled_components[index].coefficient

    @cython.profile(True)
    #TODO: architecture different from force calculator no base function/class
    #TODO: still uses component to result...
    #TODO: add common force and shift base class
    def __call__(self, object components, CDSSharedVectorFloat shift_cache, int[:] component_to_target, int[:] active_components):

        cdef float angle
        cdef float angle_term
        cdef float shift
        cdef float coefficient
        cdef dihedral_parameters parameters
        cdef dihedral_ids dihedral_atom_ids
        cdef double start_time = 0.0
        cdef double end_time = 0.0



        if self._verbose:
            start_time = time()
        cdef int factor_index
        cdef int component_index
        cdef int offset

        self._set_components(components)
        for factor_index in range(len(active_components)):
            component_index = active_components[factor_index]

            dihedral_atom_ids = self._get_dihedral_atom_ids(component_index)

            coefficient = self._get_coefficient(component_index)

            parameters = self._get_parameters(component_index)

            angle = calc_dihedral_angle_simulation(self._simulation, dihedral_atom_ids)

            angle_term = parameters.param_0 * cos(3.0 * angle + parameters.param_1) + \
                         parameters.param_2 * cos(angle +  parameters.param_3) +      \
                         parameters.param_4
            shift = coefficient * angle_term

            offset = self.ensemble_array_offset(component_to_target[factor_index])
            shift_cache.addToResult(offset,shift)

        if self._verbose:
            end_time = time()
            print '   dihedral shift components ',len(components), 'in', "%.17g" %  (end_time-start_time), "seconds"

cdef class Fast_ring_shift_calculator(Base_shift_calculator):
    cdef Vec3_list _centre_cache
    cdef Vec3_list _normal_cache

    cdef Ring_target_component* _compiled_components
    cdef int _num_components
    cdef object raw_data

    cdef Ring_component* _compiled_ring_components
    cdef int _num_ring_components
    cdef object _raw_ring_component_data


    cdef Coef_components _compiled_coef_components

    def __cinit__(self, simulation, str name = "not set"):
        self.raw_data = None
        self._compiled_components = NULL
        self._num_components = 0

        self._compiled_ring_components = NULL
        self._num_ring_components = 0
        self._raw_ring_component_data = None

        #note this is not a raw array of structs it's a compiled python class
        self._compiled_coef_components = None

        self._centre_cache = None
        self._normal_cache = None

    def __init__(self, simulation, str name = "not set"):
        Base_shift_calculator.__init__(self,simulation, name)

    def _bytes_to_components(self, data):
        self.raw_data =  data
        self._compiled_components =  <Ring_target_component*> <size_t> ctypes.addressof(data)
        self._num_components =  len(data)/ sizeof(Ring_target_component)




    def _set_normal_cache(self,normals):
        self._normal_cache = <Vec3_list> normals

    def _set_centre_cache(self,centres):
        self._centre_cache = <Vec3_list> centres

    def _set_coef_components(self,coef_components, component_offsets):
        self._compiled_coef_components = Coef_components(coef_components, component_offsets)

    cdef void _bytes_to_ring_components(self, data):
        self._raw_ring_component_data =  data
        self._compiled_ring_components =  <Ring_component*> <size_t> ctypes.addressof(data)
        self._num_ring_components =  len(data)/ sizeof(Ring_component)


    #TODO: not needed ?
    def _set_ring_components(self,ring_components):
        self. _bytes_to_ring_components(ring_components)



    @cython.profile(False)
    cdef inline Vec3* _get_ring_normal(self, int ring_id):
        return  self._normal_cache.get(ring_id)

    @cython.profile(False)
    cdef inline Vec3* _get_ring_centre(self, int ring_id):
        return  self._centre_cache.get(ring_id)

    cpdef float  _calc_sub_component_shift(self, int target_atom_id, int ring_id, float coefficient):

        cdef Vec3 target_atom_pos
        cdef Vec3 ring_centre
        cdef Vec3 ring_normal
        cdef float lenght_normal
        cdef Vec3 direction_vector
        cdef float distance, distance3, angle, contrib

        target_atom_pos =  self._simulation[0].atomPos(target_atom_id)
        ring_centre = self._get_ring_centre(ring_id)[0]

        #TODO add this to a cache the same way that camshift does
        ring_normal = self._get_ring_normal(ring_id)[0]
        length_normal = norm(ring_normal)

        #correct name?
        direction_vector = target_atom_pos - ring_centre

        distance = norm(direction_vector)
        distance3 = distance ** 3

        angle = dot(direction_vector, ring_normal) / (distance * length_normal)
        contrib = (1.0 - 3.0 * angle ** 2) / distance3

        return contrib * coefficient

    @cython.profile(True)
    def __call__(self, object components, CDSSharedVectorFloat shift_cache, int[:] component_to_target, int[:] active_components):

        cdef int target_atom_id
        cdef int atom_type_id
        cdef int ring_id
        cdef float coefficient
        cdef double start_time = 0.0
        cdef double end_time =0.0
        cdef Component_Offsets* coeff_offset
        cdef Coef_component* coef_component

        if self._verbose:
            start_time = time()

        self._set_components(components)

        cdef int factor_index
        cdef int component_index
        cdef int offset

        if self._num_ring_components > 0:
            for factor_index in range(len(active_components)):
                component_index = active_components[factor_index]

                target_atom_id = self._compiled_components[component_index].target_atom_id
                atom_type_id = self._compiled_components[component_index].atom_type_id

                shift = 0.0

                coeff_offset =  self._compiled_coef_components.get_id_offsets(atom_type_id)

                for coef_offset in range(coeff_offset[0].offset,coeff_offset[0].offset+coeff_offset[0].length):
    #                TODO: remove magic numbers or add structs
                    coef_component = self._compiled_coef_components.get_component(coef_offset)

                    ring_id = coef_component[0].ring_id
                    coefficient = coef_component[0].coefficient
                    shift += self._calc_sub_component_shift(target_atom_id,  ring_id, coefficient)

                offset = self.ensemble_array_offset(component_to_target[factor_index])
                shift_cache.addToResult(offset,shift)

        if self._verbose:
            end_time = time()
            print '   ring shift components ' ,self._name,len(components), 'in', "%.17g" %  (end_time-start_time), "seconds"

#
cdef int RING_ATOM_IDS = 1

cdef class Fast_ring_data_calculator:
    cdef bint _verbose
    cdef EnsembleSimulation* _simulation

    def __init__(self,simulation):
        self._verbose = False
        self._simulation = simulationAsNative(simulation)


    def set_verbose(self,on):
        self._verbose =  on

    cdef inline void _calculate_one_ring_centre(self, ring_component, Vec3* result):

        atom_ids = ring_component[RING_ATOM_IDS]
        cdef float num_atom_ids = len(atom_ids)
        cdef Vec3 total
        cdef int atom_id

        result[0]=Vec3(0.0,0.0,0.0)
        for atom_id in atom_ids:
            result[0] += self._simulation.atomPos(atom_id)

        operator_times(result[0],1.0/float(num_atom_ids))

    cdef inline void  _calculate_normal(self, int[3] atom_id_triplet, Vec3* result):

        cdef Vec3 normal
        cdef Vec3 vec_1
        cdef Vec3 vec_2

        cdef Vec3 atom_vector_1, atom_vector_2, atom_vector_3

        atom_vector_1 = self._simulation[0].atomPos(atom_id_triplet[0])
        atom_vector_2 = self._simulation[0].atomPos(atom_id_triplet[1])
        atom_vector_3 = self._simulation[0].atomPos(atom_id_triplet[2])

        vec_1 =  atom_vector_1 -atom_vector_2
        vec_2 =  atom_vector_3- atom_vector_2

        result[0] =  cross(vec_1,vec_2)

    #TODO could try newells method http://www.opengl.org/wiki/Calculating_a_Surface_Normal
    cdef inline  void  _calculate_one_ring_normal(self, ring_component, Vec3* result):

        atom_ids = ring_component[RING_ATOM_IDS]
#        self._check_ring_size_ok(atom_ids)

        cdef int* atom_triplet_1 = [0,0,0]
        cdef int* atom_triplet_2 = [0,0,0]
        self._build_atom_triplet(atom_ids[:3], atom_triplet_1)
        self._build_atom_triplet(atom_ids[-3:], atom_triplet_2)


        cdef Vec3 normal_1
        self._calculate_normal(atom_triplet_1, &normal_1)
        cdef Vec3 normal_2
        self._calculate_normal(atom_triplet_2,  &normal_2)

        self._average_2_vec_3(normal_1, normal_2, result)

    cdef inline void  _build_atom_triplet(self,atom_ids, int[3]& result):
        cdef int i
        for i in range(3):
            result[i] = atom_ids[i]


    cdef inline void _average_2_vec_3(self, Vec3& vector_1, Vec3& vector_2, Vec3* result):
        result[0] =  vector_1 + vector_2

        operator_times(result[0], 1.0/2.0)


    @cython.profile(True)
    def __call__(self, rings, Vec3_list normals, Vec3_list centres):
        cdef double start_time = 0.0
        cdef double end_time = 0.0

        if self._verbose:
            start_time = time()

        normals.set_length(len(rings))
        centres.set_length(len(rings))
        cdef int ring_id
        for ring_id,ring_component in enumerate(rings):
            ring_id = ring_component[0]

            self._calculate_one_ring_centre(ring_component, centres.get(ring_id))


            self._calculate_one_ring_normal(ring_component, normals.get(ring_id))
#             normal_component = ring_id, normal
#             normals.add_component(normal_component)

        if self._verbose:
            end_time = time()
            print '   ring data centres: ',len(rings), ' normals ', len(normals), 'in', "%.17g" %  (end_time-start_time), "seconds"

#TODO: add common base class with Non bonded calculator
cdef class Fast_hydrogen_bond_calculator:
    cdef float _cutoff_distance
    cdef bint _verbose
    cdef EnsembleSimulation* _simulation
    cdef int _min_residue_seperation

    def __init__(self,simulation,cutoff_distance=5.0):
        self._simulation = simulationAsNative(simulation)
        self._min_residue_seperation =  1

        self._cutoff_distance =  cutoff_distance
        self._verbose =  False


    def set_verbose(self,on):
        self._verbose=on

    cdef   int _bytes_to_components(self, data, Hydrogen_bond_component** target_pointer):
        target_pointer[0] =  <Hydrogen_bond_component*> <size_t> ctypes.addressof(data)
        return len(data)/ sizeof(Hydrogen_bond_component)

    cdef  int _bytes_to_parameters(self,data, Hydrogen_bond_parameter** target_pointer):
        target_pointer[0] =  <Hydrogen_bond_parameter*> <size_t> ctypes.addressof(data)
        return len(data)/ sizeof(Hydrogen_bond_parameter)

    cdef  int _bytes_to_donor_lookup(self,data, Hydrogen_bond_donor_lookup** target_pointer):
        target_pointer[0] =  <Hydrogen_bond_donor_lookup*> <size_t> ctypes.addressof(data)
        return len(data)/ sizeof(Hydrogen_bond_donor_lookup)

    cdef  int _bytes_to_acceptor_lookup(self,data, Hydrogen_bond_acceptor_lookup** target_pointer):
        target_pointer[0] =  <Hydrogen_bond_acceptor_lookup*> <size_t> ctypes.addressof(data)
        return len(data)/ sizeof(Hydrogen_bond_acceptor_lookup)

    cdef inline bint _filter_by_residue(self, char* seg_1, int residue_1, char* seg_2, int residue_2) nogil:
        cdef int sequence_distance
        cdef bint result

        result = False

        if strcmp(seg_1, seg_2) == 0:
            sequence_distance =<int>fabs(residue_1-residue_2)
            if sequence_distance < self._min_residue_seperation:
                result =True
        return result

    cdef  inline float _calc_energy(self, float distance_or_angle, Hydrogen_bond_parameter params) nogil:


        cdef float term_1 = (params.p[1]/distance_or_angle) + params.p[3]
        cdef float term_2 = (params.p[2]/distance_or_angle) + params.p[4]

        return params.p[0] * (term_1**params.r - term_2**params.s) + params.p[5]

    @cython.profile(True)
    def __call__(self, components, float[:] energies):#,  Hydrogen_bond_energy_list energy_list):

        donor_list = components['DONR']
        acceptor_list =  components['ACCP']
        parameters =  components['PARA']
        donor_indices = components['DIDX']
        acceptor_indices = components['AIDX']

        if self._verbose:
            print '***** BUILD HYDROGEN BOND LIST ******'

        cdef  Hydrogen_bond_component* donor_components
        cdef int num_donor_components = self._bytes_to_components(donor_list,&donor_components)

        cdef Hydrogen_bond_component* acceptor_components
        cdef int num_acceptor_components = self._bytes_to_components(acceptor_list,&acceptor_components)

        cdef Hydrogen_bond_donor_lookup* donor_lookup
        cdef int num_donor_lookup_components = self._bytes_to_donor_lookup(donor_indices,&donor_lookup)

        cdef Hydrogen_bond_parameter* hydrogen_bond_parameters
        cdef int num_parameters = self._bytes_to_parameters(parameters, &hydrogen_bond_parameters)

        cdef Hydrogen_bond_acceptor_lookup* acceptor_lookup
        cdef int num_acceptor_lookup_components = self._bytes_to_acceptor_lookup(acceptor_indices,&acceptor_lookup)
        cdef double start_time = 0.0
        cdef double end_time = 0.0

        cdef Atom atom_1, atom_2
        cdef int current_acceptor
        cdef int atom_id_1,atom_id_2
        cdef float distance
        cdef float min_distance
        cdef float FLT_MAX=1000000.0 #TODO: replacce with header

        cdef bond_angle_ids bond_angle_atom_ids

        cdef char *seg_1
        cdef char *seg_2
        cdef int residue_1, residue_2

        cdef Hydrogen_bond_parameter params_dist
        cdef Hydrogen_bond_parameter params_angle_1
        cdef Hydrogen_bond_parameter params_angle_2

        cdef float angle_1,angle_2
        cdef float energy_dist
        cdef float energy_ang_1, energy_ang_2

        if self._verbose:
            start_time = time()


        cdef int i

        with nogil:
    #         hydrogen_bond_list.clear()
            for i in range(num_donor_components):
                atom_id_1 = donor_components[i].direct_atom_id

                atom_1 = self._simulation[0].atomByID(atom_id_1)
                seg_1 =  <char *>atom_1.segmentName()
                residue_1 = atom_1.residueNum()

                current_acceptor = -1
                min_distance = FLT_MAX

                for j in range(num_acceptor_components):
                    atom_id_2  = acceptor_components[j].direct_atom_id

                    atom_2 = self._simulation[0].atomByID(atom_id_2)


                    seg_2 = <char*>atom_2.segmentName()
                    residue_2 =  atom_2.residueNum()

                    distance = norm(atom_1.pos() - atom_2.pos())
                    residue_distance_ok = not self._filter_by_residue(seg_1, residue_1, seg_2, residue_2)

                    if distance < min_distance and residue_distance_ok:
                        donor_index = i
                        acceptor_index  = j
                        min_distance = distance

                if  min_distance < 3.0:
                    donor_direct_atom_id = donor_components[donor_index].direct_atom_id
                    donor_indirect_atom_id = donor_components[donor_index].indirect_atom_id
                    acceptor_direct_atom_id = acceptor_components[acceptor_index].direct_atom_id
                    acceptor_indirect_atom_id = acceptor_components[acceptor_index].indirect_atom_id

                    bond_angle_atom_ids.atom_id_1 = donor_direct_atom_id
                    bond_angle_atom_ids.atom_id_2 = acceptor_direct_atom_id
                    bond_angle_atom_ids.atom_id_3 = acceptor_indirect_atom_id

                    angle_1 = calc_bond_angle_simulation(self._simulation, bond_angle_atom_ids)

    #
                    bond_angle_atom_ids.atom_id_1 = donor_indirect_atom_id
                    bond_angle_atom_ids.atom_id_2 = donor_direct_atom_id
                    bond_angle_atom_ids.atom_id_3 = acceptor_direct_atom_id

                    angle_2 = calc_bond_angle_simulation(self._simulation, bond_angle_atom_ids)

                    if (angle_1 > PI/2.0 and angle_1 < PI *3.0/2.0) and (angle_2 > PI/2.0 and angle_2 < PI *3.0/2.0):

                        donor_atom_type_id  = donor_components[donor_index].atom_type_id
                        acceptor_atom_type_id = acceptor_components[acceptor_index].atom_type_id
                        acceptor_id = donor_lookup[donor_atom_type_id].acceptor_offset[acceptor_atom_type_id]

                        param_id_dist = acceptor_lookup[acceptor_id].parameter_indices[0]
                        param_id_angle_1 = acceptor_lookup[acceptor_id].parameter_indices[1]
                        param_id_angle_2 = acceptor_lookup[acceptor_id].parameter_indices[2]
        #
                        params_dist  = hydrogen_bond_parameters[param_id_dist]
                        params_angle_1  = hydrogen_bond_parameters[param_id_angle_1]
                        params_angle_2  = hydrogen_bond_parameters[param_id_angle_2]

                        energy_dist = self._calc_energy(min_distance,params_dist)
                        energy_ang_1 = self._calc_energy(angle_1 / PI * 180.0, params_angle_1)
                        energy_ang_2 = self._calc_energy(angle_2 / PI * 180.0, params_angle_2)

                        if donor_components[donor_index].backbone > -1:
                            energies[donor_components[donor_index].backbone*3] = energy_dist
                            energies[donor_components[donor_index].backbone*3+1] = energy_ang_1
                            energies[donor_components[donor_index].backbone*3+2] = energy_ang_2
    #
                        if acceptor_components[acceptor_index].backbone > -1:
                            energies[acceptor_components[acceptor_index].backbone*3] = energy_dist
                            energies[acceptor_components[acceptor_index].backbone*3+1] = energy_ang_1
                            energies[acceptor_components[acceptor_index].backbone*3+2] = energy_ang_2


        if self._verbose:
            end_time = time()
            print '   hydrogen donors: ',len(donor_list),' acceptors: ', len(acceptor_list),' in', "%.17g" %  (end_time-start_time), "seconds"

cdef class Fast_hydrogen_bond_shift_calculator(Base_shift_calculator):

    cdef float[:]  _hydrogen_bond_energies
    cdef int _num_hydrogen_bond_energies
    cdef Hydrogen_bond_shift_component* _compiled_components
    cdef int _num_components
    cdef object raw_data

    def __cinit__(self, simulation,  str name = "not set"):
        self._hydrogen_bond_energies = None
        self.raw_data = None
        self._num_hydrogen_bond_energies = -1

    def __init__(self, simulation, str name = "not set"):
        Base_shift_calculator.__init__(self, simulation, name)

    def _bytes_to_components(self,data):

        self.raw_data =  data
        self._compiled_components =  <Hydrogen_bond_shift_component*> <size_t> ctypes.addressof(data)
        self._num_components =  len(data)/ sizeof(Hydrogen_bond_shift_component)

    def _energies_to_components(self,float [:] energies):
        self. _hydrogen_bond_energies = energies
        self._num_hydrogen_bond_energies =energies.shape[0]/3

    def _set_components(self,components):
        Base_shift_calculator._set_components(self,components)
        self._energies_to_components(components['HBLT'])

    #TODO: add common force and shift base class
    def __call__(self, object components, CDSSharedVectorFloat shift_cache, int[:] component_to_target, int[:] active_components):
        cdef double start_time = 0.0
        cdef double end_time = 0.0


        if self._verbose:
            start_time=time()

        self._set_components(components)

        if self._verbose:
            start_time = time()

        cdef int factor_index = 0
        cdef int component_index
        cdef int offset, hbond_index
        cdef float shift

        #TODO: note to cython list for componnt_index in active_components produces awful code!
        if active_components ==  None:
            with nogil:
                for factor_index in range(self._num_components):


                    hbond_index = self._compiled_components[factor_index].hbond_index
                    target_atom_index = self._compiled_components[factor_index].target_atom_index
                    dist_coef = self._compiled_components[factor_index].dist_coef
                    ang1_coef = self._compiled_components[factor_index].ang1_coef
                    ang2_coef = self._compiled_components[factor_index].ang2_coef

                    shift = 0.0
                    if self._hydrogen_bond_energies[hbond_index*3] > 0.0:

                        shift+=self._hydrogen_bond_energies[hbond_index*3]* dist_coef
                        shift+=self._hydrogen_bond_energies[(hbond_index*3)+1]* ang1_coef
                        shift+=self._hydrogen_bond_energies[(hbond_index*3)+2]* ang2_coef

                        offset = self.ensemble_array_offset(component_to_target[factor_index])
                        shift_cache.addToResult(offset, shift)

        else:
            raise Exception("not implemented!")
        if self._verbose:
            end_time = time()
            print '   hydrogen bond shift components ' ,self._name,len(self._num_components), 'in', "%.17g" % (end_time-start_time), "seconds"

cdef class Non_bonded_bins:

    cdef int _min_residue_seperation
    cdef float _cutoff_distance
    cdef float _jitter
    cdef float _full_cutoff_distance
    cdef bint _verbose
    cdef EnsembleSimulation* _simulation

    cdef float _spacing
    cdef int x_steps, y_steps, z_steps

    cdef float _x_min,_y_min,_z_min

    cdef CDSVector[CDSVector[CDSVector[CDSList[int]]]] bins

    cdef inline int max(self,int val_1,int val_2):
        if val_1 > val_2:
            return val_1
        else:
            return val_2

    cdef inline int min(self,int val_1,int val_2):
        if val_1 < val_2:
            return val_1
        else:
            return val_2

    def __init__(self,simulation, int min_residue_seperation=2, float cutoff_distance=5.0, float jitter=0.2):
        self._simulation = simulationAsNative(simulation)

        self._min_residue_seperation =  min_residue_seperation
        self._cutoff_distance =  cutoff_distance
        self._jitter = jitter
        self._spacing =  self._cutoff_distance + self._jitter
        self._verbose =  False

    cdef inline void _clear_bins(self,):
        self.bins.resize(self.x_steps);
        for i in range(self.x_steps):
            self.bins[i].resize(self.y_steps)
            for j in range(self.y_steps):
                self.bins[i][j].resize(self.z_steps)
                for k in range(self.z_steps):
                    self.bins[i][j][k].resize(0)

    def add_to_bins(self,int[:] atom_ids):
        cdef float huge = 1e30

        cdef float x_min = huge
        cdef float y_min = huge
        cdef float z_min = huge

        cdef float x_max = -huge
        cdef float y_max = -huge
        cdef float z_max = -huge

        cdef Vec3 pos
        cdef int atom_id

        for atom_id in atom_ids:

            pos = self._simulation[0].atomByID(atom_id).pos()

            if pos.x()<x_min:
                x_min = pos.x()
            if pos.y()<y_min:
                y_min = pos.y()
            if pos.z()<z_min:
                z_min = pos.z()

            if pos.x()>x_max:
                x_max = pos.x()
            if pos.y()>y_max:
                y_max = pos.y()
            if pos.z()>z_max:
                z_max = pos.z()

        self.x_steps = self.max(1, <int>ceil( (x_max - x_min) / self._spacing))
        self.y_steps = self.max(1, <int>ceil( (y_max - y_min) / self._spacing))
        self.z_steps = self.max(1, <int>ceil( (z_max - z_min) / self._spacing))

        self._clear_bins()


        self._x_min = x_min
        self._y_min = y_min
        self._z_min = z_min

        for atom_id in atom_ids:
            self._add_to_bin(atom_id)

    cdef inline void _add_to_bin(self, int atom_id):
        cdef Vec3 pos = self._simulation[0].atomByID(atom_id).pos()
        cdef int x_bin = self._find_x_bin(pos.x())
        cdef int y_bin = self._find_y_bin(pos.y())
        cdef int z_bin = self._find_z_bin(pos.z())

        self.bins[x_bin][y_bin][z_bin].append(atom_id);

    cdef inline int _find_x_bin(self, float x):
        return self.max(0,self.min(self.x_steps-1,(int)((x-self._x_min)/self._spacing)))

    cdef inline int _find_y_bin(self, float y):
        return self.max(0,self.min(self.y_steps-1,(int)((y-self._y_min)/self._spacing)))

    cdef inline int _find_z_bin(self, float z):
        return self.max(0,self.min(self.z_steps-1,(int)((z-self._z_min)/self._spacing)))

    def dump(self):
        for i in range(self.x_steps):
            for j in range(self.y_steps):
                for k in range(self.z_steps):
                    print i,j,k,self.bins[i][j][k].size()

    def get_bin(self,x,y,z):
        result  = []
        for i in range(self.bins[x][y][z].size()):
            result.append(self.bins[x][y][z][i])
        return result
# cdef class New_fast_non_bonded_calculator:
#     cdef int _min_residue_seperation
#     cdef float _cutoff_distance
#     cdef float _jitter
#     cdef float _full_cutoff_distance
#     cdef bint _verbose
#     cdef EnsembleSimulation* _simulation
#
#     def __init__(self,simulation, min_residue_seperation,cutoff_distance=5.0,jitter=0.2):
#         self._simulation = simulationAsNative(simulation)
#
#         self._min_residue_seperation =  min_residue_seperation
#         self._cutoff_distance =  cutoff_distance
#         self._jitter = jitter
#         self._full_cutoff_distance =  self._cutoff_distance + self._jitter
#         self._verbose =  False

cdef class Fast_non_bonded_calculator:
    cdef int _min_residue_seperation
    cdef float _cutoff_distance
    cdef float _jitter
    cdef float _full_cutoff_distance
    cdef bint _verbose
    cdef EnsembleSimulation* _simulation

    def __init__(self,simulation, min_residue_seperation,cutoff_distance=5.0,jitter=0.2):
        self._simulation = simulationAsNative(simulation)

        self._min_residue_seperation =  min_residue_seperation
        self._cutoff_distance =  cutoff_distance
        self._jitter = jitter
        self._full_cutoff_distance =  self._cutoff_distance + self._jitter
        self._verbose =  False


    def set_verbose(self,on):
        self._verbose=on

    cdef inline bint _filter_by_residue(self, char* seg_1, int residue_1, char* seg_2, int residue_2):
        cdef int sequence_distance
        cdef bint result

        result = False

        if strcmp(seg_1, seg_2) == 0:
            sequence_distance =abs(residue_1-residue_2)
            if sequence_distance < self._min_residue_seperation:
                result =True
        return result

    cdef inline bint  _is_non_bonded(self, int atom_id_1, int atom_id_2):

        cdef Atom atom_1, atom_2
        cdef char *seg_1
        cdef char *seg_2
        cdef int residue_1, residue_2
        cdef bint is_non_bonded
        cdef float distance

        atom_1 = self._simulation[0].atomByID(atom_id_1)
        atom_2 = self._simulation[0].atomByID(atom_id_2)

        seg_1 =  <char *>atom_1.segmentName()
        residue_1 = atom_1.residueNum()

        seg_2 = <char*>atom_2.segmentName()
        residue_2 =  atom_2.residueNum()

        is_non_bonded = True
        if self._filter_by_residue(seg_1, residue_1, seg_2, residue_2):
            is_non_bonded = False
        else:
            distance = norm(atom_1.pos() - atom_2.pos())
            is_non_bonded =  distance < self._full_cutoff_distance
        return is_non_bonded


    cdef   int _bytes_to_target_components(self, data, Non_bonded_target_component** target_pointer):

        target_pointer[0] =  <Non_bonded_target_component*> <size_t> ctypes.addressof(data)
        return len(data)/ sizeof(Non_bonded_target_component)

    cdef   int _bytes_to_remote_components(self, data, Non_bonded_remote_atom_component** remote_pointer):

        remote_pointer[0] =  <Non_bonded_remote_atom_component*> <size_t> ctypes.addressof(data)
        return len(data)/ sizeof(Non_bonded_remote_atom_component)

    @cython.profile(True)
    def __call__(self, atom_list_1, atom_list_2,  Non_bonded_interaction_list non_bonded_lists, int[:] active_components):

        if self._verbose:
            print '***** BUILD NON BONDED ******'

        cdef  Non_bonded_target_component* target_components
        cdef int num_target_components = self._bytes_to_target_components(atom_list_1,&target_components)

        cdef Non_bonded_remote_atom_component* remote_components
        cdef int num_remote_components = self._bytes_to_remote_components(atom_list_2,&remote_components)

        cdef double start_time = 0.0
        cdef double end_time = 0.0
        if self._verbose:
            start_time = time()


        cdef int i, atom_id_1, atom_id_2

        non_bonded_lists.clear()

        if active_components == None:

            for i in range(num_target_components):
                atom_id_1 = target_components[i].target_atom_id

                for j in range(num_remote_components):
                    atom_id_2  = remote_components[j].remote_atom_id
                    if self._is_non_bonded(atom_id_1, atom_id_2):
                        non_bonded_lists.append(atom_id_1, i,j,i)
        else:
            for component_index in range(len(active_components)):
                i = active_components[component_index]
                atom_id_1 = target_components[i].target_atom_id

                for j in range(num_remote_components):
                    atom_id_2  = remote_components[j].remote_atom_id
                    if self._is_non_bonded(atom_id_1, atom_id_2):
                        non_bonded_lists.append(atom_id_1, i,j,component_index)


        if self._verbose:
            end_time = time()
            print '   non bonded list targets: ',len(atom_list_1),' remotes: ', len(atom_list_2),' in', "%.17g" %  (end_time-start_time), "seconds"


cdef class Fast_energy_calculator:
    cdef Constant_cache* _energy_term_cache
    cdef CDSVector[double]  *_theory_shifts
    cdef float[:] _observed_shifts
    cdef bint _verbose
    cdef int calls

    def __init__(self):
        self._energy_term_cache =  NULL
        self._theory_shifts =   NULL
        self._observed_shifts =  None
        self._verbose = False
        self.calls = 0


    def set_verbose(self,on):
        self._verbose = on

    def set_observed_shifts(self, float[:] observed_shifts):
        self._observed_shifts =  observed_shifts

    def set_calculated_shifts(self, CDSVectorFloat calculated_shifts):
        self._theory_shifts = calculated_shifts.get_data()

    def set_energy_term_cache(self, energy_term_cache ):
        self._energy_term_cache =  <Constant_cache*> <size_t> ctypes.addressof(energy_term_cache)

    cdef Constant_cache* _get_energy_terms(self, int target_atom_index):
        return &self._energy_term_cache[target_atom_index]

    cdef inline float  _get_calculated_atom_shift(self, int index):
        return self._theory_shifts[0][index]

    cdef inline float _get_observed_atom_shift(self, int index):
        return self._observed_shifts[index]

    cdef inline float  _get_shift_difference(self, int target_atom_index, int index):
        cdef float theory_shift
        cdef float observed_shift
        theory_shift = self._get_calculated_atom_shift(index)

        observed_shift = self._get_observed_atom_shift(index)

        return observed_shift - theory_shift

    cdef inline float _adjust_shift(self, float shift_diff, float flat_bottom_shift_limit):
        result  = 0.0
        if (shift_diff > 0.0):
            result = shift_diff-flat_bottom_shift_limit
        else:
            result = shift_diff + flat_bottom_shift_limit
        return result

    @cython.profile(True)
    def __call__(self,int[:] target_atom_ids, int[:] active_atom_ids=None):

        cdef double start_time = 0.0
        cdef double end_time = 0.0
        cdef float energy
        cdef int active_atom_id
        cdef int target_atom_id

        if self._verbose:
            start_time = time()

        energy = 0.0

        if active_atom_ids == None:
            for i in range(target_atom_ids.shape[0]):
                target_atom_id = target_atom_ids[i]

                energy += self._calc_one_energy(target_atom_id, i)
        else:
            for i in  range(active_atom_ids.shape[0]):
                active_atom_id = active_atom_ids[i]
                target_atom_id = target_atom_ids[active_atom_id]

                energy += self._calc_one_energy(target_atom_id, active_atom_id)

        if self._verbose:
            end_time = time()
            print '   energy calculator: ',len(target_atom_ids),' in', "%.17g" %  (end_time-start_time), "seconds"


        self.calls += 1
        return energy

    cdef inline float _calc_one_energy(self, int target_atom_index, int index):
        cdef float flat_bottom_shift_limit
        cdef float adjusted_shift_diff
        cdef float end_harmonic
        cdef float scale_harmonic
        cdef float tanh_amplitude
        cdef float tanh_elongation
        cdef float tanh_y_offset
        cdef float tanh_argument
        cdef float energy
        cdef float shift_diff


        cdef float shift_diffs
        cdef Constant_cache* energy_terms
        cdef float atom_energy

        shift_diff = self._get_shift_difference(target_atom_index, index)
        energy_terms = self._get_energy_terms(index)

        flat_bottom_shift_limit = energy_terms[0].flat_bottom_shift_limit

        energy = 0.0
        if fabs(shift_diff) > flat_bottom_shift_limit:
            adjusted_shift_diff = self._adjust_shift(shift_diff, flat_bottom_shift_limit)

            end_harmonic = energy_terms[0].end_harmonic
            scale_harmonic = energy_terms[0].scale_harmonic


            if adjusted_shift_diff < end_harmonic:
                atom_energy = (adjusted_shift_diff/scale_harmonic)**2
            else:
                tanh_amplitude = energy_terms[0].tanh_amplitude
                tanh_elongation = energy_terms[0].tanh_elongation
                tanh_y_offset = energy_terms[0].tanh_y_offset

                tanh_argument = tanh_elongation * (adjusted_shift_diff - end_harmonic)
                atom_energy = tanh_amplitude * tanh(tanh_argument) + tanh_y_offset

            energy += atom_energy * energy_terms[0].weight
        return energy

cdef class Fast_force_factor_calculator(Fast_energy_calculator):

    @cython.profile(True)
    def __call__(self, int[:] target_atom_ids, float[:] result, int[:] active_atom_ids):
        cdef double start_time =0.0
        cdef double end_time =0.0

        cdef int i

        if self._verbose:
            start_time = time()


       #TODO: shouldn't be allocated each time
        cdef int target_atom_id
        cdef int active_atom_id

        if active_atom_ids == None:
            for i in range(target_atom_ids.shape[0]):
                target_atom_id = target_atom_ids[i]
                result[i] = self._calc_one_force_factor(target_atom_id, i)
        else:
            for i in  range(active_atom_ids.shape[0]):
                active_atom_id = active_atom_ids[i]
                target_atom_id = target_atom_ids[active_atom_id]

                result[i] = self._calc_one_force_factor(target_atom_id, active_atom_id)

        if self._verbose:
            end_time = time()
            print '   force factors : ',len(target_atom_ids),' in', "%.17g" %  (end_time-start_time), "seconds"


    cdef inline float _calc_one_force_factor(self, int target_atom_id, int i):

        cdef float factor
        cdef float shift_diff
        cdef float flat_bottom_shift_limit
        cdef Constant_cache* energy_terms
        cdef float adjusted_shift_diff
        cdef float end_harmonic
        cdef float scale_harmonic
        cdef float sqr_scale_harmonic
        cdef float weight
        cdef float tanh_amplitude
        cdef float tanh_elongation

        factor = 0.0


        shift_diff = self._get_shift_difference(target_atom_id, i)
        energy_terms = self._get_energy_terms(i)


        flat_bottom_shift_limit = energy_terms[0].flat_bottom_shift_limit

        if abs(shift_diff) > flat_bottom_shift_limit:
            adjusted_shift_diff = self._adjust_shift(shift_diff, flat_bottom_shift_limit)
            end_harmonic = energy_terms[0].end_harmonic
            scale_harmonic = energy_terms[0].scale_harmonic
            sqr_scale_harmonic = scale_harmonic**2

            weight = energy_terms[0].weight

            tanh_amplitude = energy_terms[0].tanh_amplitude
            tanh_elongation = energy_terms[0].tanh_elongation

            # TODO: add factor and lambda to give fact
            fact =1.0
            if adjusted_shift_diff < end_harmonic:
                factor = 2.0 * weight * adjusted_shift_diff * fact / sqr_scale_harmonic;
            else:
                factor = weight * tanh_amplitude * tanh_elongation / (cosh(tanh_elongation * (adjusted_shift_diff - end_harmonic)))**2.0 * fact;
        return factor


cdef class Base_force_calculator:

    cdef bint _verbose
    cdef Simulation* _simulation
    cdef str _name
    cdef int[:] _active_components

    def __init__(self,potential=None, name= "not set"):
        self._verbose = False


        self._name = name


    def set_verbose(self,on):
        self._verbose = on



    def _set_components(self,components):
        self._bytes_to_components(components['ATOM'])
        self._simulation = <Simulation*><size_t>int(components['SIMU'].this)


#    TODO should most probably be a fixed array
    @cython.profile(True)
    def __call__(self, object components, int[:] component_to_result, float[:] force_factors, Out_array forces, int[:] active_components=None):

        cdef double start_time =0.0
        cdef double end_time =0.0
        if self._verbose:
            start_time = time()

        self._set_components(components)
#        TODO rename component to result to something better
        self._active_components = active_components
        self._do_calc_components(component_to_result, force_factors, forces)


        if self._verbose:
            end_time = time()
            print '   force calculator: ', self._name ,' ',len(components), ' in', "%.17g" %  (end_time-start_time), "seconds"

    cdef void _do_calc_components(self, int[:] component_to_result,  float[:] force_factors, Out_array force):
        raise Exception("this should not be called!")




cdef class Fast_distance_based_potential_force_calculator(Base_force_calculator):



    cdef int  _target_atom_index
    cdef int  _distance_atom_index_1
    cdef int  _distance_atom_index_2
    cdef int  _exponent_index
    cdef int  _coefficient_index
    cdef bint  _smoothed
    cdef float _smoothing_factor
    cdef float _cutoff
    cdef Distance_component* _compiled_components
    cdef int _num_components
    cdef object raw_data

    def __cinit__(self):
        self._compiled_components  = NULL
        self._num_components =  0

    def __init__(self, bint smoothed, name="Not set"):
        super(Fast_distance_based_potential_force_calculator, self).__init__(name=name)
        self._smoothed =  smoothed
        self._smoothing_factor =  DEFAULT_SMOOTHING_FACTOR
        self._cutoff =  DEFAULT_CUTOFF

    def set_cutoff(self, cutoff):
        self._cutoff =  cutoff

    def set_smoothing_factor(self,smoothing_factor):
        self._smoothing_factor = smoothing_factor




    def _bytes_to_components(self, data):

        self.raw_data =  data
        self._compiled_components =  <Distance_component*> <size_t> ctypes.addressof(data)
        self._num_components =  len(data)/ sizeof(Distance_component)


#   TODO: generalise this based on compiled components
    def _calc_single_force_set(self, int index, float factor, Out_array forces):
        #TODO tidy this up a hack for the test suite
        cdef Distance_component* saved_component_list = self._compiled_components
        cdef int saved_num_components = self._num_components

        self._num_components  = 1
        self._compiled_components = &saved_component_list[index]

        self._distance_calc_single_force_set(0,factor,forces)

        self._compiled_components = saved_component_list
        self._num_components = saved_num_components


    @cython.profile(False)
    cdef inline Vec3 _xyz_distances(self, int target_atom, int distance_atom):
        cdef Vec3 target_pos, distant_pos


        target_pos = self._simulation[0].atomPos(target_atom)
        distant_pos =   self._simulation[0].atomPos(distance_atom)

        return target_pos - distant_pos

    @cython.profile(False)
    cdef inline float _sum_xyz_distances_2(self, int target_atom, int distance_atom):
        cdef Vec3 target_pos, distant_pos, distance
        cdef float result =0.0
        cdef int i

        target_pos =   self._simulation[0].atomPos(target_atom)
        distant_pos =  self._simulation[0].atomPos(distance_atom)

        distance = Vec3(target_pos - distant_pos)

        for i in range(3):
            result += distance[i] * distance[i]

        return result

    def _calc_single_force_factor(self, int index, float factor):
        return self._cython_calc_single_force_factor(index, factor)

    cdef void _do_calc_components(self, int[:] component_to_result, float[:] force_factors, Out_array force):
        cdef int factor_index
        cdef int component_index

        for factor_index in range(len(self._active_components)):
            component_index = self._active_components[factor_index]
            self._distance_calc_single_force_set(component_index,force_factors[component_to_result[factor_index]],force)


    @cython.profile(False)
    cdef inline float _cython_calc_single_force_factor(self, int index, float factor):
        cdef target_distant_atom atom_ids
        cdef float exponent

        cdef float full_factor, force_factor
        cdef float ratio, pre_exponent, reduced_exponent, sum_xyz_distances_2


        sum_xyz_distances_2 = self._sum_xyz_distances_2(self._compiled_components[index].remote_atom_1, self._compiled_components[index].remote_atom_2)
        full_factor= factor *  self._compiled_components[index].coefficient

        exponent =  self._compiled_components[index].exponent
        if self._smoothed:
            ratio = sum_xyz_distances_2 / (self._cutoff**2)
            ratio =  ratio**4
            pre_exponent = exponent - (exponent + 8.0) * ratio
        else:
            pre_exponent = exponent

        reduced_exponent = (exponent - 2.0) / 2

        force_factor = full_factor *  pre_exponent * sum_xyz_distances_2 ** reduced_exponent

        return force_factor


    @cython.profile(False)
    cdef inline void _distance_calc_single_force_set(self,int index, float factor, Out_array forces):

        cdef Vec3 xyz_distances
        cdef float force_factor, distance
        cdef Vec3 result

        xyz_distances  = self._xyz_distances(self._compiled_components[index].remote_atom_1,self._compiled_components[index].remote_atom_2)

        force_factor  = self._cython_calc_single_force_factor(index, factor)

        cdef Vec3 target_forces = xyz_distances
        operator_times(target_forces, -force_factor)

        cdef Vec3 distant_forces = xyz_distances
        operator_times(distant_forces, force_factor)
#
        forces.add(self._compiled_components[index].remote_atom_1,target_forces)
        forces.add(self._compiled_components[index].remote_atom_2,distant_forces)


cdef class Fast_non_bonded_force_calculator(Fast_distance_based_potential_force_calculator):

    cdef float _nb_cutoff

    cdef Non_bonded_interaction_list _non_bonded_list


    #TODO: can we use memory views here...
    cdef object _raw_target_data
    cdef Non_bonded_target_component* _compiled_target_components
    cdef int _num_target_components


    #TODO: can we use memory views here...
    cdef object _raw_remote_data
    cdef Non_bonded_remote_atom_component* _compiled_remote_components
    cdef int _num_remote_components

    cdef object _raw_coefficient_data
    cdef Nonbonded_coefficient_component* _compiled_coefficient_components
    cdef int _num_coefficient_components

    cdef int _component_offset

    def __cinit__(self):
        self._non_bonded_list = None
        self._compiled_target_components =  NULL
        self._compiled_remote_components = NULL
        self._compiled_coefficient_components = NULL
        self._component_offset = 0

    def __init__(self, bint smoothed, str name):
        super(Fast_non_bonded_force_calculator, self).__init__(smoothed,name)
        global DEFAULT_NB_CUTOFF
        self._nb_cutoff = DEFAULT_NB_CUTOFF
        self._compiled_components =  <Distance_component*>malloc(sizeof(Distance_component))

    def set_verbose(self,on):
        self._verbose = on

    #TODO: use global version
    cdef _bytes_to_target_components(self,data):
        self._raw_target_data =  data

        self._compiled_target_components =  <Non_bonded_target_component*> <size_t> ctypes.addressof(data)
        self._num_target_components =  len(data)/ sizeof(Non_bonded_target_component)

    #TODO: use global version
    cdef _bytes_to_remote_components(self,data):
        self._raw_remote_data =  data

        self._compiled_remote_components =  <Non_bonded_remote_atom_component*> <size_t> ctypes.addressof(data)
        self._num_remote_components =  len(data)/ sizeof(Non_bonded_remote_atom_component)

    cdef _bytes_to_nonbonded_coefficient_components(self,data):
        self._raw_coefficient_data =  data

        self._compiled_coefficient_components =  <Nonbonded_coefficient_component*> <size_t> ctypes.addressof(data)
        self._num_coefficient_components =  len(data)/ sizeof(Nonbonded_coefficient_component)

    #TODO: call super?
    def _set_components(self, components):
        self._simulation = <Simulation*><size_t>int(components['SIMU'].this)
        self._non_bonded_list =  components['NBLT']
        self._bytes_to_target_components(components['ATOM'])
        self._bytes_to_remote_components(components['NBRM'])
        self._bytes_to_nonbonded_coefficient_components(components['COEF'])
        if 'ACTI' in components:
            self._active_components = components['ACTI']
        self._component_offset = components['OFFS']

    cdef inline void _set_the_component(self, int target_atom_index,int remote_atom_index, float coefficient,float exponent):
        self._compiled_components[0].target_atom =target_atom_index
        self._compiled_components[0].remote_atom_1 =target_atom_index
        self._compiled_components[0].remote_atom_2 =remote_atom_index
        self._compiled_components[0].coefficient =coefficient
        self._compiled_components[0].exponent =exponent



    cdef void _do_calc_components(self, int[:] component_to_result, float[:] force_factors, Out_array force):


        cdef int factor_index

        if self._active_components ==  None:
            for factor_index  in range(len(self._non_bonded_list)):
                non_bonded_index  = factor_index
                self._calc_one_component(factor_index, non_bonded_index, component_to_result, force_factors, force)
        else:
            for factor_index  in range(self._active_components.shape[0]):
                non_bonded_index = self._active_components[factor_index]
                self._calc_one_component(factor_index, non_bonded_index, component_to_result, force_factors, force)

    cdef void _calc_one_component(self, int factor_index, int non_bonded_index, int[:] component_to_result, float[:] force_factors, Out_array force):
        cdef double start_time = 0.0
        cdef double end_time = 0.0

        cdef Component_index_pair* non_bonded_pair

        cdef int target_component_index
        cdef int remote_component_index

        cdef int target_index
        cdef int remote_index

        cdef float distance
        non_bonded_pair  =  self._non_bonded_list.get(non_bonded_index)

        target_component_index = non_bonded_pair[0].target_index
        remote_component_index = non_bonded_pair[0].remote_index
        component_offset = self._component_offset + non_bonded_pair[0].component_index

        target_index  = self._compiled_target_components[target_component_index].target_atom_id
        remote_index = self._compiled_remote_components[remote_component_index].remote_atom_id
        distance = calc_distance_simulation(self._simulation, target_index, remote_index)
        if distance < self._nb_cutoff:
            for i in range(2):
                    self._cython_build_component(non_bonded_index, i)
                    self._distance_calc_single_force_set(0, force_factors[component_offset] , force)

    def  _build_component(self, non_bonded_index, i):
        self._cython_build_component(non_bonded_index, i)

    cdef void  _cython_build_component(self, int non_bonded_index, int i):
        cdef Component_index_pair* non_bonded_pair

        cdef int target_component_index
        cdef int remote_component_index

        cdef int target_index
        cdef int remote_index

        cdef float distance

        cdef int chem_type_id

        cdef float default_smoothing_factor = self._smoothing_factor
        cdef float smoothing_factor

        cdef float ratio


        cdef int atom_1_coefficent_offset

        cdef float coefficient,exponent

        cdef Nonbonded_coefficient_component* coefficent_component
        cdef bint result

        non_bonded_pair  =  self._non_bonded_list.get(non_bonded_index)

        target_component_index = non_bonded_pair[0].target_index
        remote_component_index = non_bonded_pair[0].remote_index

        target_index  = self._compiled_target_components[target_component_index].target_atom_id
        remote_index = self._compiled_remote_components[remote_component_index].remote_atom_id
        atom_1_coefficent_offset = self._compiled_target_components[target_component_index].atom_type_id

        chem_type_id = self._compiled_remote_components[remote_component_index].chem_type[i]

        coefficient_component = &self._compiled_coefficient_components[chem_type_id]

        exponent = coefficient_component[0].exponent

        coefficient  = coefficient_component[0].coefficients[atom_1_coefficent_offset]

        self._set_the_component(target_index,remote_index,coefficient,exponent)



cdef class Fast_dihedral_force_calculator(Base_force_calculator):

    cdef Dihedral_component* _compiled_components
    cdef int _num_components
    cdef object raw_data

    def _bytes_to_components(self, data):
        self.raw_data =  data

        self._compiled_components =  <Dihedral_component*> <size_t> ctypes.addressof(data)
        self._num_components =  len(data)/ sizeof(Dihedral_component)


    def __cinit__(self):
        self._compiled_components = NULL
        self._num_components = 0

    def __init__(self,name="not set"):
        Base_force_calculator.__init__(self,name=name)
        self.raw_data =  None



#     TODO: remove this is no longer needed
    @cython.profile(False)
    cdef inline dihedral_ids _get_dihedral_atom_ids(self, int index):
        cdef dihedral_ids result

        result.atom_id_1 = self._compiled_components[index].dihedral_atoms[0]
        result.atom_id_2 = self._compiled_components[index].dihedral_atoms[1]
        result.atom_id_3 = self._compiled_components[index].dihedral_atoms[2]
        result.atom_id_4 = self._compiled_components[index].dihedral_atoms[3]

        return result

    #TODO: remove this an correct order of parameters!
    @cython.profile(False)
    cdef inline dihedral_parameters _get_parameters(self, int index):

        cdef dihedral_parameters result

        result.param_0 = self._compiled_components[index].parameters[0]
        result.param_1 = self._compiled_components[index].parameters[3]
        result.param_2 = self._compiled_components[index].parameters[1]
        result.param_3 = self._compiled_components[index].parameters[4]
        result.param_4 = self._compiled_components[index].parameters[2]




        return result

    @cython.profile(False)
    cdef inline float _get_coefficient(self, int index):
        return self._compiled_components[index].coefficient

#    TODO make this consistent with the distance forces factor
    def _calc_single_force_factor(self, int index):
        return self._cython_calc_single_force_factor(index)

    cdef void _do_calc_components(self, int[:] component_to_result, float[:] force_factors, Out_array force):
        cdef int factor_index
        cdef int component_index

        for factor_index in range(len(self._active_components)):
            component_index = self._active_components[factor_index]
            self._dihedral_calc_single_force_set(component_index,force_factors[component_to_result[factor_index]],force)



    cdef inline float _cython_calc_single_force_factor(self, int index):

        cdef dihedral_ids dihedral_atom_ids
        cdef dihedral_parameters params

        cdef float angle,result


        dihedral_atom_ids  = dihedral_atom_ids= self._get_dihedral_atom_ids(index)

        params = self._get_parameters(index)

        angle = calc_dihedral_angle_simulation(self._simulation, dihedral_atom_ids)

        result = -3.0 * params.param_0 * sin(3.0 * angle + params.param_1) - \
                        params.param_2 * sin(angle + params.param_3)
        return result



    def _calc_single_force_set(self, int index, float factor, Out_array forces):
        #TODO tidy this up a hack for the test suite
        cdef Dihedral_component* saved_component_list = self._compiled_components
        cdef int saved_num_components = self._num_components

        self._num_components  = 1
        self._compiled_components = &saved_component_list[index]

        self._dihedral_calc_single_force_set(0,factor,forces)

        self._compiled_components = saved_component_list
        self._num_components = saved_num_components


    cdef inline void _dihedral_calc_single_force_set(self, int index, float factor, Out_array forces):
        cdef Vec3 r1, r2, r4, temp
        cdef Vec3 n1, n2
        cdef float weight
        cdef float factor_1, factor_2, factor_3
        cdef Vec3 F2, F3,
        cdef Vec3 T3, T4

        cdef float dihedral_factor = self._cython_calc_single_force_factor(index)
        cdef dihedral_ids atom_ids = self._get_dihedral_atom_ids(index)

        cdef float r2_length, r2_length_2

        v1 = self._simulation[0].atomPos(atom_ids.atom_id_1)
        v2 = self._simulation[0].atomPos(atom_ids.atom_id_2)
        v3 = self._simulation[0].atomPos(atom_ids.atom_id_3)
        v4 = self._simulation[0].atomPos(atom_ids.atom_id_4)

        r1 = v1 - v2
        r2 = v3 - v2

        cdef Vec3 r3 = r2
        operator_times(r3, -1)

        r4 = v4 - v3

        # compute normal vector to plane containing v1, v2, and v3
        n1 = cross(r1, r2)
        # compute normal vector to plane containing v2, v3, and v4
        n2 = cross(r3, r4)

        r2_length = calc_distance_simulation(self._simulation, atom_ids.atom_id_2,atom_ids.atom_id_3)
        r2_length_2 = r2_length*r2_length


        weight = factor * self._get_coefficient(index)



#        // force calculation according to Bekker, Berendsen and van Gunsteren (1995),
#        // Journal of Computational Chemistry 16, pp. 527-533:
##        // Force and virial of torsional-angle-dependent potentials.
        factor_1 = dihedral_factor * r2_length

        cdef Vec3 F1 = n1
        operator_times(F1, (-factor_1 / norm(n1)**2))

        cdef Vec3 F4 = n2
        operator_times(F4, ( factor_1 / norm(n2)**2))

        factor_2 = dot(r1, r2) / r2_length_2
        factor_3 = dot(r3, r4) / r2_length_2

        cdef Vec3 T1 = F1
        operator_times(T1,  (factor_2-1))

        cdef Vec3 T2 = F4
        operator_times(T2,  -factor_3)

        F2 = T1 + T2

        T1 = F1
        operator_times(T1,  -factor_2)

        T2 = F4
        operator_times(T2, (factor_3 -1))

        F3 = T1 + T2

        operator_times(F1, weight)
        operator_times(F2, weight)
        operator_times(F3, weight)
        operator_times(F4, weight)

        forces.add(atom_ids.atom_id_1,F1)
        forces.add(atom_ids.atom_id_2,F2)
        forces.add(atom_ids.atom_id_3,F3)
        forces.add(atom_ids.atom_id_4,F4)


cdef struct Ring_target_component:
      int target_atom_id
      int atom_type_id

cdef struct Ring_component:
    int ring_id
    int num_atoms
    int[6] atom_ids

#cdef int RING_ATOM_IDS = 1
cdef class Fast_ring_force_calculator(Base_force_calculator):


    cdef Vec3_list _centre_cache
    cdef Vec3_list _normal_cache

    cdef Ring_target_component* _compiled_components
    cdef int _num_components

    cdef Ring_component* _compiled_ring_components
    cdef int _num_ring_components
    cdef object _raw_ring_component_data

    cdef object raw_data

    cdef Coef_components _compiled_coef_components


    def __cinit__(self):
        self.raw_data =  None
        self._compiled_components = NULL
        self._num_components = 0

        self._raw_ring_component_data = None
        self._compiled_ring_components = NULL
        self._num_ring_components = 0

        #note this is not a raw array of structs it's a compiled python class
        self._compiled_coef_components = None

        self._centre_cache = None
        self._normal_cache = None

    def __init__(self,name="not set"):
        super(Fast_ring_force_calculator, self).__init__(name=name)

    def _bytes_to_components(self, data):
        self.raw_data =  data
        self._compiled_components =  <Ring_target_component*> <size_t> ctypes.addressof(data)
        self._num_components =  len(data)/ sizeof(Ring_target_component)

    cdef void _bytes_to_ring_components(self, data):
        self._raw_ring_component_data =  data
        self._compiled_ring_components =  <Ring_component*> <size_t> ctypes.addressof(data)
        self._num_ring_components =  len(data)/ sizeof(Ring_component)


    def _set_coef_components(self,coef_components, components):
        self._compiled_coef_components = Coef_components(coef_components, components)

    def _set_ring_components(self,ring_components):
        self. _bytes_to_ring_components(ring_components)



    def _set_normal_cache(self,normals):
        self._normal_cache = <Vec3_list> normals

    def _set_centre_cache(self,centres):
        self._centre_cache = <Vec3_list> centres


    @cython.profile(False)
    cdef inline Vec3* _get_ring_normal(self, int ring_id):
        return  self._normal_cache.get(ring_id)

    @cython.profile(False)
    cdef inline Vec3* _get_ring_centre(self, int ring_id):
        return  self._centre_cache.get(ring_id)




    cdef inline void _do_calc_components(self, int[:] component_to_result, float[:] force_factors, Out_array force):
        cdef int factor_index
        cdef int component_index

        if self._num_ring_components > 0:
            for factor_index in range(len(self._active_components)):
                component_index = self._active_components[factor_index]
                self._ring_calc_single_force_set(component_index,force_factors[component_to_result[factor_index]],force)

    cdef inline void _ring_calc_single_force_set(self,  int index, float force_factor, Out_array forces):
        cdef int atom_type_id  = self._compiled_components[index].atom_type_id
        cdef Component_Offsets* coeff_offset =  self._compiled_coef_components.get_id_offsets(atom_type_id)
        cdef Coef_component* coef_component

        for coef_offset in range(coeff_offset[0].offset,coeff_offset[0].offset+coeff_offset[0].length):
#           TODO: remove magic numbers or add structs
            coef_component = self._compiled_coef_components.get_component(coef_offset)
            self._calculate_single_ring_forces(self._compiled_components[index].target_atom_id, self._compiled_components[index].atom_type_id, coef_component, force_factor, forces)

    def _calculate_ring_forces(self, int atom_type_id, int ring_id, float force_factor, Python_ring_force_sub_terms python_sub_terms,Out_array forces):
        cdef Ring_force_sub_terms  terms  =  python_sub_terms.get_terms()
        self._cython_calculate_ring_forces(atom_type_id, ring_id, force_factor, terms, forces)

    cdef void _calculate_single_ring_forces(self, int target_atom_id, int atom_type_id, Coef_component* coef_component, float force_factor,Out_array forces):

        cdef Ring_force_sub_terms force_terms
        self._build_cython_force_terms(target_atom_id, coef_component.ring_id, force_terms)

        self._cython_calc_target_atom_forces(target_atom_id, force_factor * coef_component.coefficient, force_terms, forces)

#        #            #TODO: this is not how camshift does it, it uses the sum of the two ring normals
        self._cython_calculate_ring_forces(atom_type_id, coef_component.ring_id, force_factor * coef_component.coefficient, force_terms, forces)

    def _build_force_terms(self, int target_atom_id, int ring_id):
        cdef Ring_force_sub_terms terms
        self._build_cython_force_terms(target_atom_id, ring_id, terms)
        result  = Python_ring_force_sub_terms()
        result.setup(terms)
        return result

    cdef inline void _build_cython_force_terms(self, int target_atom_id, int ring_id, Ring_force_sub_terms& result):
        cdef Vec3 target_atom_pos  = self._simulation[0].atomPos(target_atom_id)

        cdef Vec3 ring_centre = self._get_ring_centre(ring_id)[0]
        cdef Vec3* ring_normal = self._get_ring_normal(ring_id)

        # distance vector between atom of interest and ring center
        cdef Vec3 d = target_atom_pos - ring_centre
        cdef float dL = norm(d)

        # squared distance of atom of interest from ring center
        cdef float dL2 = dot(d, d) #            if (dL2 < 0.5) cout << "CAMSHIFT WARNING: Distance between atom and center of ring alarmingly small at " << sqrt(dL2) << " Angstrom!" << endl;

        # calculate terms resulting from differentiating energy function with respect to query and ring atom coordinates
        cdef float dL4 = dL ** 4
        cdef float dL3 = dL ** 3
        cdef float dL6 = dL3 ** 2

        cdef float nL = norm(ring_normal[0])
        cdef float nL2 = nL ** 2
        cdef float dLnL = dL * nL
        cdef float dL3nL3 = dL3 * nL2 * nL

        cdef float dn = dot(d, ring_normal[0])
        cdef float dn2 = dn ** 2

        cdef float u = 1.0 - 3.0 * dn2 / (dL2 * nL2)

        cdef float gradUQ_factor = -6.0 * dn / (dL4 * nL2)

        #TODO: remove temporarys and operator_mu;t#
        cdef Vec3 scaled_d =  d
        operator_times(scaled_d,dn)
        cdef Vec3 scaled_normal = ring_normal[0]
        operator_times(scaled_normal,dL2)

        cdef Vec3 gradUQ = scaled_normal - scaled_d
        operator_times(gradUQ,gradUQ_factor)


        cdef float gradVQ_factor = 3.0 * dL
        cdef Vec3 gradVQ = d
        operator_times(gradVQ,gradVQ_factor)

        result.dL3nL3 = dL3nL3
        result.dLnL =dLnL
        result.dn = dn
        result.dL = dL
        result.nL = nL
        result.d = d
        result.gradUQ = gradUQ
        result.dL3 =dL3
        result.u  = u
        result.gradVQ = gradVQ
        result.dL6 = dL6
        result.ring_normal =  ring_normal[0]

##    ---
#
    def _calc_target_atom_forces(self, int target_atom_id, float force_factor, Python_ring_force_sub_terms python_sub_terms, Out_array forces):
        cdef Ring_force_sub_terms  terms  =  python_sub_terms.get_terms()
        self._cython_calc_target_atom_forces(target_atom_id, force_factor, terms, forces)

    cdef inline void _cython_calc_target_atom_forces(self, int target_atom_id, float force_factor, Ring_force_sub_terms& sub_terms, Out_array forces):
        cdef object target_force_triplet
        cdef int axis
#        target_force_triplet = self._get_or_make_target_force_triplet(forces, target_atom_id)
#
#        # update forces on query atom
        cdef float term1
        cdef Vec3 result
        for axis in range(3):
            #TODO: report this as a bug it produces an un initialized reference
            term1= (sub_terms.gradUQ[axis] * sub_terms.dL3 - sub_terms.u * sub_terms.gradVQ[axis])
            result[axis] = -force_factor * term1 / sub_terms.dL6

        forces.add(target_atom_id, result)
#             *  / sub_terms.dL6
#
#    #TODO: calculation of GradU and gradV are not consistent with force_terms for target atom correct
#    #TODO: reduce number of parameters to method

    cdef inline void _cython_calculate_ring_forces(self, int atom_type_id, int ring_id, float force_factor, Ring_force_sub_terms force_terms, Out_array forces):
        cdef Vec3 temp_normal =  Vec3(force_terms.ring_normal)
        cdef Vec3 nSum = temp_normal
        operator_times(nSum,2.0)

        #// 2 for a 5-membered ring, 3 for a 6-membered ring
        cdef int num_ring_atoms = self._compiled_ring_components[ring_id].num_atoms
        cdef int limit = num_ring_atoms - 3

#        for atom_type_id, ring_id, coefficient in coef_components:
        cdef Vec3 g = Vec3()
        cdef int ring_atom_id
        cdef int ring_atom_index
        cdef int index_1
        cdef int index_2
        cdef float pos_1
        cdef float pos_2
        cdef float pos_3
        cdef float pos_4
        cdef int offset

        cdef int[3] indices_1
        cdef int[3] indices_2

        cdef float[3] ab
        cdef float[3] c

        cdef float factor
        cdef float factor2
        cdef float one_over_num_ring_atoms
        cdef float factor3

        cdef Vec3 gradU
        cdef Vec3 gradV
        cdef Vec3 result

        cdef int axis
        cdef float sub_term
        cdef float sub_force

        for ring_atom_index in range(num_ring_atoms):
            ring_atom_id = self._compiled_ring_components[ring_id].atom_ids[ring_atom_index]
            if ring_atom_index < limit:
                for axis in range(3):
                    index_1 = (ring_atom_index + 1) % 3
                    index_2 = (ring_atom_index + 2) % 3
                    pos_1 = self._simulation[0].atomPos(self._compiled_ring_components[ring_id].atom_ids[index_1])[axis]
                    pos_2 = self._simulation[0].atomPos(self._compiled_ring_components[ring_id].atom_ids[index_2])[axis]
                    g[axis] =  pos_1 - pos_2 # atoms 3,4 (5 member) or 3,4,5 (6 member)

            else:
                if  ring_atom_index >= num_ring_atoms - limit:
                    offset = num_ring_atoms - 3 #2 for a 5-membered ring, 3 for a 6-membered ring
                    for axis in range(3):
                        index_1 = (ring_atom_index + 1 - offset) % 3 + offset
                        index_2 = (ring_atom_index + 2 - offset) % 3 + offset
                        pos_1 = self._simulation[0].atomPos(self._compiled_ring_components[ring_id].atom_ids[index_1])[axis]
                        pos_2 = self._simulation[0].atomPos(self._compiled_ring_components[ring_id].atom_ids[index_2])[axis]
                        g[axis] = pos_1 - pos_2
                else:

                    for axis in range(3):
                        pos_1 = self._simulation[0].atomPos(self._compiled_ring_components[ring_id].atom_ids[0])[axis]
                        pos_2 = self._simulation[0].atomPos(self._compiled_ring_components[ring_id].atom_ids[1])[axis]
                        pos_3 = self._simulation[0].atomPos(self._compiled_ring_components[ring_id].atom_ids[3])[axis]
                        pos_4 = self._simulation[0].atomPos(self._compiled_ring_components[ring_id].atom_ids[4])[axis]
                        g[axis] = pos_1 - pos_2 + pos_3 - pos_4

            # 0 1 2 2 1   (0+1) %3 (0+2) %3
            # 1 2 0 0 2
            # 2 0 1 10
            #atom 2 (5-membered rings)
            for axis in range(3):
                indices_1[axis] = (axis + 1) % 3
                indices_2[axis] = (axis + 2) % 3

            for axis in range(3):
                index_1 = indices_1[axis]
                index_2 = indices_2[axis]
                ab[axis] = force_terms.d[index_1] * g[index_2] - force_terms.d[index_2] * g[index_1]

            for axis in range(3):
                index_1 = indices_1[axis]
                index_2 = indices_2[axis]
                c[axis] = nSum[index_1] * g[index_2] - nSum[index_2] * g[index_1]

            factor = -6.0 * force_terms.dn / force_terms.dL3nL3
            factor2 = 0.25 * force_terms.dL / force_terms.nL
            one_over_num_ring_atoms = 1.0 / float(num_ring_atoms)
            factor3 = force_terms.nL / force_terms.dL * one_over_num_ring_atoms

            gradU = Vec3()
            for axis in range(3):
                gradU[axis] = factor * ((0.5 * ab[axis] - force_terms.ring_normal[axis] * one_over_num_ring_atoms) * force_terms.dLnL - force_terms.dn * (factor2 * c[axis] - factor3 * force_terms.d[axis]))

            gradV = Vec3()
            factor = -3 * force_terms.dL * one_over_num_ring_atoms
            for axis in range(3):
                gradV[axis] = factor * force_terms.d[axis]

#TODO: rename force terms  sub_terms or vice versa

#            ring_target_force_triplet = self._get_or_make_target_force_triplet(forces, ring_atom_id)

            for axis in range(3):
                sub_term = (gradU[axis] * force_terms.dL3 - force_terms.u * gradV[axis])
                sub_force = -force_factor * sub_term / force_terms.dL6
                result[axis] = sub_force
            forces.add(ring_atom_id,result)


cdef class Fast_non_bonded_shift_calculator(Fast_distance_shift_calculator):

    cdef float _nb_cutoff

    cdef Non_bonded_interaction_list _non_bonded_list


    #TODO: can we use memory views here...
    cdef object _raw_target_data
    cdef Non_bonded_target_component* _compiled_target_components
    cdef int _num_target_components


    #TODO: can we use memory views here...
    cdef object _raw_remote_data
    cdef Non_bonded_remote_atom_component* _compiled_remote_components
    cdef int _num_remote_components

    cdef object _raw_coefficient_data
    cdef Nonbonded_coefficient_component* _compiled_coefficient_components
    cdef int _num_coefficient_components

    cdef int _component_offset



    def __cinit__(self, simulation, bint smoothed, str name):

        self._non_bonded_list = None
        self._compiled_target_components =  NULL
        self._compiled_remote_components = NULL
        self._compiled_coefficient_components = NULL
        self. _component_offset = 0

    def __init__(self, simulation, bint smoothed, str name):
        super(Fast_non_bonded_shift_calculator, self).__init__(simulation, smoothed,name)
        self._simulation = simulationAsNative(simulation)
        global DEFAULT_NB_CUTOFF
        self._nb_cutoff = DEFAULT_NB_CUTOFF

    def set_verbose(self,on):
        self._verbose = on

    #TODO: use global version
    cdef _bytes_to_target_components(self,data):
        self._raw_target_data =  data

        self._compiled_target_components =  <Non_bonded_target_component*> <size_t> ctypes.addressof(data)
        self._num_target_components =  len(data)/ sizeof(Non_bonded_target_component)

    #TODO: use global version
    cdef _bytes_to_remote_components(self,data):
        self._raw_remote_data =  data

        self._compiled_remote_components =  <Non_bonded_remote_atom_component*> <size_t> ctypes.addressof(data)
        self._num_remote_components =  len(data)/ sizeof(Non_bonded_remote_atom_component)

    cdef _bytes_to_nonbonded_coefficient_components(self,data):
        self._raw_coefficient_data =  data

        self._compiled_coefficient_components =  <Nonbonded_coefficient_component*> <size_t> ctypes.addressof(data)
        self._num_coefficient_components =  len(data)/ sizeof(Nonbonded_coefficient_component)

    def _set_components(self, components):
        self._simulation = <EnsembleSimulation*><size_t>int(components['SIMU'].this)
        self._non_bonded_list =  components['NBLT']
        self._bytes_to_target_components(components['ATOM'])
        self._bytes_to_remote_components(components['NBRM'])
        self._bytes_to_nonbonded_coefficient_components(components['COEF'])
        self._component_offset = components['OFFS']



    @cython.profile(True)
    def __call__(self, object components, CDSSharedVectorFloat shift_cache, int[:] component_to_target, int[:] active_components):
        self._set_components(components)

        if active_components == None:
            for non_bonded_index in range(len(self._non_bonded_list)):
#                 print non_bonded_index
                self._calc_single_component_shift(non_bonded_index, non_bonded_index, shift_cache, component_to_target)
        else:
            for factor_index  in range(active_components.shape[0]):
                non_bonded_index = active_components[factor_index]
#                 print factor_index, non_bonded_index
                self._calc_single_component_shift(factor_index, non_bonded_index, shift_cache, component_to_target)


    cdef inline void  _calc_single_component_shift(self,int factor_index, int non_bonded_index, CDSSharedVectorFloat  results, int[:] component_to_target):
        cdef Component_index_pair* non_bonded_pair

        cdef int target_component_index
        cdef int remote_component_index

        cdef int target_index
        cdef int remote_index

        cdef float distance

        cdef int chem_type_id

        cdef float default_smoothing_factor = self._smoothing_factor
        cdef float smoothing_factor

        cdef float ratio

        cdef int atom_1_coefficent_offset

        cdef float coefficient,exponent

        cdef Nonbonded_coefficient_component* coefficent_component
        cdef int component_offset
        cdef int offset

        if self._verbose:
            start_time = time()


        non_bonded_pair  =  self._non_bonded_list.get(non_bonded_index)

        target_component_index = non_bonded_pair[0].target_index
        remote_component_index = non_bonded_pair[0].remote_index

        target_index  = self._compiled_target_components[target_component_index].target_atom_id
        remote_index = self._compiled_remote_components[remote_component_index].remote_atom_id

        distance = calc_distance_simulation(self._simulation, target_index, remote_index)

        component_offset = self._component_offset + non_bonded_pair[0].component_index

        if distance < self._nb_cutoff:

            atom_1_coefficent_offset = self._compiled_target_components[target_component_index].atom_type_id

            for i in range(2):
                chem_type_id = self._compiled_remote_components[remote_component_index].chem_type[i]

                coefficient_component = &self._compiled_coefficient_components[chem_type_id]

                exponent = coefficient_component[0].exponent

                coefficient  = coefficient_component[0].coefficients[atom_1_coefficent_offset]

                smoothing_factor = default_smoothing_factor
                if self._smoothed:
                    ratio = distance / self._cutoff
                    smoothing_factor = 1.0 - ratio ** 8

                offset = self.ensemble_array_offset(component_to_target[component_offset])
                results.addToResult(offset, smoothing_factor * pow(distance,  exponent) * coefficient)


        if self._verbose:
            end_time = time()
            print '   distance shift components ' ,self._name,len(self._non_bonded_list), 'in', "%.17g" % (end_time-start_time), "seconds"


