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
'''
Created on 31 Jul 2012

@author: garyt

'''

cimport cython
from vec3 import Vec3 as python_vec3
from common_constants import TARGET_ATOM_IDS_CHANGED, STRUCTURE_CHANGED
from  xplor_access cimport norm,Vec3,currentSimulation, Dihedral, Atom,  dot,  cross,  Simulation, CDSVector
from libc.math cimport cos,sin,  fabs, tanh, pow, cosh
from libc.stdlib cimport malloc, free
from libc.string cimport strcmp
from time import time
from utils import Atom_utils
from cpython cimport array

cpdef array.array allocate_array(int len, type='d'):
    result = array.array(type,[0])
    array.resize(result, len)
    array.zero(result)
    return result
    
cpdef zero_array(array.array array):
     array.zero(array)
     
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
      float[5] parameters
 
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
        
    def __dealloc__(self):
        
        if self.data != NULL:
            del self.data
            
    def __iter__(self):
        raise Exception()
    
    def __iter__(self): 
        cdef Vec3 vec3
        
        
        
        for i in range(self.data.size()):
            vec3 = self.get(i)[0]
            yield (i,python_vec3(vec3[0],vec3[1],vec3[2]))
    
    
cdef class Out_array:
    cdef long _length
    cdef double[60000] _data
    cdef int[60000] _mask
    
    def __init__ (self, length):
        self._length = length

        
    def __cinit__(self,  length):
        self._length =  60000
        
        if length > 60000: 
            raise MemoryError("couldn't allocate array beyond length %i" % 60000)
#        self._data = <double *>malloc(length * sizeof(double))
#        if not self._data:
#            raise MemoryError("couldn't allocate array of length %i" % length)
#        
#        self._mask = <long *>malloc(length * sizeof(long))
#        if not self._mask:
#            raise MemoryError("couldn't allocate array of length %i" % length)
#        
#        if  self._data  and self._mask:
#            self._clear(length)
#        print 'allocated'
            
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
#            free(self._data)
#            self._data = NULL
#            self._data = <double *>malloc(3 *(length+1) * sizeof(double))
#            if not self._data:
#                raise MemoryError("couldn't allocate array of length %i" % length)
#            if  self._data:
#                self._length =length
            self._clear(length)
        else:
            raise MemoryError("couldn't beyond of length %i" % 60000)
#
#            free(self._mask)
#            self._mask = NULL                        
#            self._mask = <int *>malloc((length+1) * sizeof(int))
#            if not self._mask:
#                raise MemoryError("couldn't allocate array of length %i" % length)
#            if  self._mask:
#                self._clear_mask(length)
            
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
#        self._check_offset(offset)
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
        if result ==  None:
            result = [None] * self._length
        
        cdef double x = 0.0
        cdef double y = 0.0
        cdef double z = 0.0
        cdef long id3
        for i in range(self._length):
            if self._mask[i] != 0:
                if result[i] ==  None:
                    result[i] = [0.0]*3
                id3 = i * 3
                x  = self._data[id3] * weight
                y  = self._data[id3+1] *weight 
                z  = self._data[id3+2] * weight  
                result[i][0] += x
                result[i][1] += y
                result[i][2] += z

                
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


#cdef float* out_array = null
#cdef long out_array_length = -1
#cdef bint allocate_out_array_for_forces(long length) except -1:
#    global out_array
#    global out_array_length
#    cdef long required_length = length*3
#    cdef bint allocate =  False
#    if out_array ==  null:
#        out_array = <float*> malloc*(required_length *sizeof(float))
#        if out_array  == null:
#            Raise
#        allocate=True
#        out_array_length = 
#        elif  len(out_array) < required_length
#    return allocate
#    
#def _clear_out_array_for_forces(self, length):
#    
#    
#    is_new  = self.allocate_out_array_for_forces(length)
#    
#    if not is_new:
#        length = length*3
#        for i in range(length):
#            self._out_array[i] =0.0
#    return self._out_array
#    
#def _copy_back_forces(self, data_dict, data_array, target_atom_ids):
#    for id in target_atom_ids:
#        result =  dict.setdefault(id,[0.0,0.0,0.0])
#        id3 = id*3
#        result[0] = data_array[id3]
#        result[1] = data_array[id3+1]
#        result[2] = data_array[id3+2]


 
cdef object vec3_as_tuple(Vec3& vec_3):
    return vec_3.x(), vec_3.y(), vec_3.z()
    
cdef class Vec3_container:
    
    cdef float[3] floats
    
    cdef set_vec3(self,Vec3& vec3):
        self.floats[0] =  vec3.x()
        self.floats[1] =  vec3.y()
        self.floats[2] =  vec3.z()
        
    @cython.profile(False)    
    cdef Vec3  get_vec3(self):
        return Vec3(self.floats[0],self.floats[1],self.floats[2])
    
    def __getitem__(self, i):
        if i > 2:
            msg = "tried to access object at %i length is %i"
            values = (i,3)
            raise IndexError(msg % values)

        return self.floats[i]
    
    
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
    
#    def __cinit__(self):
#        self._terms = NULL
        
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
    cdef int num_components
    cdef Coef_component* _components
    cdef int num_ids
    cdef Component_Offsets* _component_offsets
    
        
    def __cinit__(self, object components):
        self.num_components = len(components)
        self._components = <Coef_component *>malloc( self.num_components * sizeof(Coef_component))
        if not self._components:
            raise MemoryError()
        
        self. num_ids = len(components.get_component_atom_ids())
        self._component_offsets =  <Component_Offsets*>malloc(self.num_ids * sizeof(Coef_component))
        if not self._component_offsets:
            raise MemoryError()
        
        self.num_components = len(components)
        for i,component in enumerate(components):
            self._components[i].atom_type_id = component[0] 
            self._components[i].ring_id = component[1]
            self._components[i].coefficient = component[2]
            
        
        ids =  components.get_component_atom_ids()
        for i,id in enumerate(ids):
            start,end = components.get_component_range(id)
            self._component_offsets[i].id = id
            self._component_offsets[i].offset = start
            self._component_offsets[i].length = end-start

        
    def __init__ (self, object components):
        pass
        
    
    cdef inline Component_Offsets* get_id_offsets(self,id):
        return &self._component_offsets[id]
    
    def __dealloc__(self):
 
        
        
        if self._components != NULL: 
            free(self._components)
            self._components = NULL
            self.num_components = 0
        
        
        if self._component_offsets != NULL:
            free(self._component_offsets)
            self._component_offsets =  NULL
            self.num_ids = 0 
        
    
#    cdef int _check_offset(self, int offset) except -1:
#    
#        if not offset < self.num_components:
#            raise IndexError("tried to access object at %i length is %i" % (offset,self._num_components))
#        if self._components is NULL:
#            raise AttributeError("trying to access deallocated memory!")
#        return 0

     

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
    
    
#cdef float sum(Vec3 vec3):
#    return vec3.x()+vec3.y()+vec3.z()


@cython.profile(False)
cdef inline void operator_times (Vec3& vec3, float scale):
     vec3[0] =  vec3[0] * scale
     vec3[1] =  vec3[1] * scale
     vec3[2] =  vec3[2] * scale
     
@cython.profile(False)
cdef inline float calc_distance_simulation(Simulation* sim, int atom_index_1, atom_index_2):

    cdef Vec3 vec1 = sim[0].atomPos(atom_index_1)
    cdef Vec3 vec2 = sim[0].atomPos(atom_index_2)
    cdef Vec3 result =  vec2 - vec1
    return norm(result)



#
#
#cdef inline bint distance_less_than(int atom_index_1, atom_index_2, float target):
#    cdef float total = 0.0
#    cdef float diff
#    cdef float target2  = target*target
#    
#    cdef Vec3 vec1 = currentSimulation().atomPosArr().data(atom_index_1)
#    cdef Vec3 vec2 = currentSimulation().atomPosArr().data(atom_index_2)
#    cdef bint result  =  True
#    for i in range(3):
#        diff = vec1[i] - vec2[i]
#        total += diff * diff
#        
#        if total > target2:
#            result = False
#            break
#    return result 
        
        
    
    

#TODO how does currentSimulation().atomByID work inside
@cython.profile(False)
cdef inline float calc_dihedral_angle_simulation(Simulation* simulation, dihedral_ids dihedral_atom_ids):
    
    

    atom_1  = simulation[0].atomByID(dihedral_atom_ids.atom_id_1)
    atom_2  = simulation[0].atomByID(dihedral_atom_ids.atom_id_2)
    atom_3  = simulation[0].atomByID(dihedral_atom_ids.atom_id_3)
    atom_4  = simulation[0].atomByID(dihedral_atom_ids.atom_id_4)
    
    return Dihedral(atom_1,atom_2,atom_3,atom_4).value()



cdef float DEFAULT_CUTOFF = 5.0
cdef float DEFAULT_SMOOTHING_FACTOR = 1.0
cdef float DEFAULT_NB_CUTOFF = 5.0
    
cdef class Base_shift_calculator:
    cdef bint _verbose 
    cdef str _name
    cdef Simulation* _simulation

    def __init__(self, str name):
            self._verbose = False
            self._name = name
            self._simulation =  currentSimulation()
            
    
    def set_verbose(self,bint state):
        self._verbose =  state
        
    cdef set_simulation(self):
        self._simulation = currentSimulation()
    
    def _prepare(self,change,data):
        pass
    
        
        
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
        
    def __cinit__(self):
        self._compiled_components  = NULL
        self._num_components =  0
        
    def __dealloc__(self):
        
        self._free_compiled_components()
        
    cdef _free_compiled_components(self):
        if self._compiled_components != NULL:
            free (self._compiled_components)
            self._compiled_components = NULL   
                        
    def _prepare(self, change, data):
        if change == TARGET_ATOM_IDS_CHANGED or change == STRUCTURE_CHANGED:
             self._free_compiled_components()
                         
    def __init__(self, indices, smoothed, str name = "not set"):
        Base_shift_calculator.__init__(self, name)
        self._target_atom_index = indices.target_atom_index
        self._distance_atom_index_1 =  indices.distance_atom_index_1
        self._distance_atom_index_2 =  indices.distance_atom_index_2
        self._exponent_index  = indices.exponent_index
        self._coefficient_index  = indices.coefficient_index
        
        self._smoothed =  smoothed
        self._smoothing_factor =  DEFAULT_SMOOTHING_FACTOR
#        
        self._cutoff =  DEFAULT_CUTOFF
#    
#    def set_cutoff(self, cutoff):
#        self._cutoff =  cutoff
#    
#    def set_smoothing_factor(self,smoothing_factor):
#        self._smoothing_factor = smoothing_factor
#    
#    TODO: this needs to be removed

            
    def _set_components(self,components):
        if  self._compiled_components ==  NULL:
            self._compile_components(components)  

        
    cdef _compile_components(self,components):
        self._compiled_components = <Distance_component*>malloc(len(components) * sizeof(Distance_component))
        self._num_components = len(components)
        for i,component in enumerate(components):
            self._compiled_components[i].target_atom = components[i][self._target_atom_index]
            if len(component) == 4:
                self._compiled_components[i].remote_atom_1 = components[i][self._distance_atom_index_1]
                self._compiled_components[i].remote_atom_2 = components[i][self._distance_atom_index_2]
                self._compiled_components[i].coefficient   = components[i][self._coefficient_index]
                self._compiled_components[i].exponent      = components[i][self._exponent_index]
            elif len(component) == 5:
                self._compiled_components[i].remote_atom_1 = components[i][self._distance_atom_index_1]
                self._compiled_components[i].remote_atom_2 = components[i][self._distance_atom_index_2]
                self._compiled_components[i].coefficient   = components[i][self._coefficient_index]
                self._compiled_components[i].exponent      = components[i][self._exponent_index]
            else:
                raise Exception("bad distance component length %i should be either 4 or 5 " % len(component))
        
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
#    
    def __call__(self, object components, object results, object component_to_target):
        self.set_simulation()
        cdef double start_time = 0.0
        cdef double end_time = 0.0
        
        if self._verbose:
            start_time=time()
        
        self._set_components(components)
        cdef float smoothing_factor = self._smoothing_factor
        cdef float ratio
        cdef float result
        cdef target_distant_atom atom_indices
        cdef coefficient_exponent coef_exp
        cdef float coefficent
        cdef float exponent
        cdef object component
        cdef int target_atom_id
        cdef int distant_atom_id

        

        if self._verbose:
            start_time = time()
            
        for index in range(len(components)):
            
#             component = components[index]

            
            target_atom_id = self._compiled_components[index].remote_atom_1
            distant_atom_id  = self._compiled_components[index].remote_atom_2
            
            
            coef_exp = self._get_coefficient_and_exponent(index)
            distance =calc_distance_simulation(self._simulation, target_atom_id, distant_atom_id)
    #        Atom_utils._calculate_distance(target_atom_index, sidechain_atom_index)

            if self._smoothed:
                ratio = distance / self._cutoff
                smoothing_factor = 1.0 - ratio ** 8
            results[component_to_target[index]]  += smoothing_factor * pow(distance,  coef_exp.exponent) * coef_exp.coefficient

        if self._verbose:
            end_time = time()
            print '   distance shift components ' ,self._name,len(components), 'in', "%.17g" % (end_time-start_time), "seconds"

    

    
cdef class Fast_dihedral_shift_calculator(Base_shift_calculator):
    
    cdef Dihedral_component *_compiled_components
    cdef int _num_components
     
    def __cinit__(self):
        self._compiled_components = NULL
        self._num_components = 0
        
        
    def __init__(self, str name = "not set"):
        Base_shift_calculator.__init__(self,name)
    
    cdef _set_components(self, object components):
        if self._compiled_components ==  NULL:
            self._compile_components(components)
            
    def _compile_components(self, components):
        self._compiled_components = <Dihedral_component*>malloc(len(components) * sizeof(Dihedral_component))
        self._num_components = len(components) 
        for i,component in enumerate(components):
            self._compiled_components[i].target_atom = component[0]
            for j in range(4):
                self._compiled_components[i].dihedral_atoms[j] = component[1+j]
            self._compiled_components[i].coefficient = component[5]
            for j in range(5):
                self._compiled_components[i].parameters[j] = component[6+j]
    
    def _free_compiled_components(self):
        if self._compiled_components != NULL:
            free(self._compiled_components)
            self._compiled_components = NULL        
                   
    def _prepare(self, change, data):
         if change == TARGET_ATOM_IDS_CHANGED or change == STRUCTURE_CHANGED:
             self._free_compiled_components()  
                    
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

    def __call__(self, object components, object results, object component_to_target):
        self.set_simulation()
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
            
        self._set_components(components)
        for index in range(len(components)):
            dihedral_atom_ids = self._get_dihedral_atom_ids(index)
            
            coefficient = self._get_coefficient(index)
            
            parameters = self._get_parameters(index)
            
            angle = calc_dihedral_angle_simulation(self._simulation, dihedral_atom_ids)
    
            angle_term = parameters.param_0 * cos(3.0 * angle + parameters.param_1) + \
                         parameters.param_2 * cos(angle +  parameters.param_3) +      \
                         parameters.param_4
            shift = coefficient * angle_term
    
            results[component_to_target[index]] += shift
            
        if self._verbose:
            end_time = time()
            print '   dihedral shift components ',len(components), 'in', "%.17g" %  (end_time-start_time), "seconds"

cdef class Fast_ring_shift_calculator(Base_shift_calculator):
    cdef Vec3_list _centre_cache
    cdef Vec3_list _normal_cache
    
    cdef Ring_target_component* _compiled_components 
    cdef int _num_components
    
    cdef Ring_component* _compiled_ring_components
    cdef int _num_ring_components
    
    cdef Coef_components _compiled_coef_components
    
    def __cinit__(self):
        self._compiled_components = NULL
        self._num_components = 0
        
        self._compiled_ring_components = NULL
        self._num_ring_components = 0
        
        #note this is not a raw array of structs it's a compiled python class
        self._compiled_coef_components = None
        
        self._centre_cache = None
        self._normal_cache = None
    
    def __init__(self, str name = "not set"):
        Base_shift_calculator.__init__(self,name)

        
    def _set_components(self,components):
        if  self._compiled_components ==  NULL:
            self._compile_components(components)
            
    def _compile_components(self,components): 
        self._compiled_components = <Ring_target_component*>malloc(len(components) * sizeof(Ring_target_component))
        self._num_components = len(components)
        for i,component in enumerate(components):
            self._compiled_components[i].target_atom_id = components[i][0]
            self._compiled_components[i].atom_type_id = components[i][1]
    
    def _free_compiled_components(self):   
        if self._compiled_components != NULL:
            free(self._compiled_components)
            self._compiled_components = NULL
            
    def _free_compiled_ring_components(self):   
        if self._compiled_ring_components != NULL:
            free(self._compiled_ring_components)
            self._compiled_ring_components = NULL

    def _free_coef_components(self):
        self._compiled_coef_components = None


    def _prepare(self, change, data):
        if change == TARGET_ATOM_IDS_CHANGED or change == STRUCTURE_CHANGED:
            self._free_compiled_components() 
            self._free_compiled_ring_components()
            self._free_coef_components()

            
    def _set_normal_cache(self,normals):
        self._normal_cache = <Vec3_list> normals
        
    def _set_centre_cache(self,centres):
        self._centre_cache = <Vec3_list> centres
        
    def _set_coef_components(self,coef_components):
        self._compiled_coef_components = Coef_components(coef_components)
        
        
    #TODO: not needed ?        
    def _set_ring_components(self,ring_components):
        if  self._compiled_ring_components ==  NULL:
            self._compile_ring_components(ring_components)

    def _compile_ring_components(self,ring_components): 
        self._compiled_ring_components = <Ring_component*>malloc(len(ring_components) * sizeof(Ring_component))
        self._num_ring_components = len(ring_components)
        for i,ring_component in enumerate(ring_components):
            self._compiled_ring_components[i].ring_id  = ring_components[i][0]
            self._compiled_ring_components[i].num_atoms  = len(ring_component[1])
            for j in range(len(ring_component[1])):
                self._compiled_ring_components[i].atom_ids[j] =ring_component[1][j]
    
    def _free_compiled_ring_components(self):   
        if self._compiled_ring_components != NULL:
            free(self._compiled_ring_components)
            self._compiled_ring_components = NULL


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

    
    def __call__(self, object components, object results, object component_to_target):
        self.set_simulation()
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
        
        for index in range(self._num_components):
            target_atom_id = self._compiled_components[index].target_atom_id
            atom_type_id = self._compiled_components[index].atom_type_id

            shift = 0.0
        
            coeff_offset =  self._compiled_coef_components.get_id_offsets(atom_type_id)

            for coef_offset in range(coeff_offset[0].offset,coeff_offset[0].offset+coeff_offset[0].length):
#                TODO: remove magic numbers or add structs
                coef_component = self._compiled_coef_components.get_component(coef_offset)
                
                ring_id = coef_component[0].ring_id
                coefficient = coef_component[0].coefficient
                shift += self._calc_sub_component_shift(target_atom_id,  ring_id, coefficient)
            
            results[component_to_target[index]] += shift
        
        if self._verbose:
            end_time = time()
            print '   ring shift components ' ,self._name,len(components), 'in', "%.17g" %  (end_time-start_time), "seconds"

#   
cdef int RING_ATOM_IDS = 1
    
cdef class Fast_ring_data_calculator:
    cdef bint _verbose
    cdef Simulation* _simulation
         
    def __init__(self): 
        self._verbose = False
        self._simulation =  currentSimulation()
    
    def set_simulation(self):
        self._simulation  = currentSimulation()
        
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
        
#       
        
    
    
#    def _check_ring_size_ok(self, atom_ids):
#        num_atom_ids = len(atom_ids)
#        if num_atom_ids < 5 or num_atom_ids > 6:
#            template = "ring normals function is only implemented for 5 or six member rings i got %d atoms"
#            msg = template % num_atom_ids
#            raise Exception(msg)
#

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
        
#         cdef Vec3_container result = Vec3_container()
        self._average_2_vec_3(normal_1, normal_2, result)
#         result.set_vec3(average)
#         return result
       
    cdef inline void  _build_atom_triplet(self,atom_ids, int[3]& result):
        cdef int i
        for i in range(3):
            result[i] = atom_ids[i]
        
        
    cdef inline void _average_2_vec_3(self, Vec3& vector_1, Vec3& vector_2, Vec3* result):
        result[0] =  vector_1 + vector_2
            
        operator_times(result[0], 1.0/2.0)

    
    @cython.profile(True)        
    def __call__(self, rings, Vec3_list normals, Vec3_list centres):
        self.set_simulation()
#        cdef Vec3_container centre 
#         cdef Vec3_container normal
        cdef double start_time = 0.0 
        cdef double end_time = 0.0
        self.set_simulation()
        
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

cdef class Fast_non_bonded_calculator:
    cdef int _min_residue_seperation
    cdef float _cutoff_distance
    cdef float _jitter
    cdef float _full_cutoff_distance
    cdef bint _verbose 
    cdef Simulation* _simulation
    def __init__(self,min_residue_seperation,cutoff_distance=5.0,jitter=0.2):
        self._min_residue_seperation =  min_residue_seperation
        self._cutoff_distance =  cutoff_distance
        self._jitter = jitter
        self._full_cutoff_distance =  self._cutoff_distance + self._jitter
        self._verbose =  False
        self._simulation =  currentSimulation()
    
    def set_simulation(self):
        self._simulation = currentSimulation()
        
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
        cdef char *seg_1, *seg_2
        cdef int residue_1, residue_2
        cdef bint is_non_bonded
        cdef float distance
        
        atom_1 = currentSimulation().atomByID(atom_id_1)
        atom_2 = currentSimulation().atomByID(atom_id_2)
        
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
    
    @cython.profile(True)
    def __call__(self, atom_list_1, atom_list_2):
        
        if self._verbose:
            print '***** BUILD NON BONDED ******'
            
        cdef double start_time = 0.0 
        cdef double end_time = 0.0
        if self._verbose:
            start_time = time()


        self.set_simulation()
        non_bonded_lists = []
        cdef int i, atom_id_1, atom_id_2
        for atom_id_1 in atom_list_1:
            non_bonded_list = []
            non_bonded_lists.append(non_bonded_list)
            for i, atom_id_2 in enumerate(atom_list_2):
                if self._is_non_bonded(atom_id_1, atom_id_2):
                    non_bonded_list.append(i)
        if self._verbose:
            end_time = time()
            print '   non bonded list targets: ',len(atom_list_1),' remotes: ', len(atom_list_2),' in', "%.17g" %  (end_time-start_time), "seconds"
        return  non_bonded_lists


cdef class Fast_energy_calculator:
    cdef object _energy_term_cache 
    cdef object _theory_shifts
    cdef object _observed_shifts
    cdef bint _verbose 
    cdef Simulation* _simulation
    cdef int calls

    def __init__(self):
        self._energy_term_cache =  None
        self._theory_shifts =   None
        self._observed_shifts =  None
        self._verbose = False
        self._simulation = currentSimulation()
        self.calls = 0
    
    def set_simulation(self):
        self._simulation = currentSimulation()
        
    def set_verbose(self,on):
        self._verbose = on
        
    def set_observed_shifts(self, observed_shifts):
        self._observed_shifts =  observed_shifts
        
    def set_calculated_shifts(self, calculated_shifts):
        self._theory_shifts =  calculated_shifts
    
    def set_energy_term_cache(self, energy_term_cache ):
        self._energy_term_cache =  energy_term_cache
        
    cdef _get_energy_terms(self, int target_atom_index):
        return self._energy_term_cache[target_atom_index]
    
    cdef inline float  _get_calculated_atom_shift(self, int target_atom_index):
        return self._theory_shifts[target_atom_index]
    
    cdef inline float _get_observed_atom_shift(self, int target_atom_index):
        return self._observed_shifts.get_chemical_shift(target_atom_index)
    
    cdef inline float  _get_shift_difference(self, int target_atom_index):
        cdef float theory_shift
        cdef float observed_shift
        theory_shift = self._get_calculated_atom_shift(target_atom_index)
        
        observed_shift = self._get_observed_atom_shift(target_atom_index)
        
        return observed_shift - theory_shift

    cdef inline float _adjust_shift(self, float shift_diff, float flat_bottom_shift_limit):
        result  = 0.0
        if (shift_diff > 0.0):
            result = shift_diff-flat_bottom_shift_limit
        else:
            result = shift_diff + flat_bottom_shift_limit
        return result

    @cython.profile(True)    
    def __call__(self,int[:] target_atom_ids):
        self.set_simulation()
        cdef float energy
        cdef float flat_bottom_shift_limit
        cdef float adjusted_shift_diff
        cdef float end_harmonic
        cdef float scale_harmonic
        cdef float energy_component
        cdef float tanh_amplitude
        cdef float tanh_elongation
        cdef float tanh_y_offset
        cdef float tanh_argument
        
        cdef double start_time = 0.0
        cdef double end_time = 0.0 
        if self._verbose:
            start_time = time()
            
        energy = 0.0
        
        cdef int target_atom_index
        cdef float shift_diffs
         
        for target_atom_index in target_atom_ids:
            shift_diff = self._get_shift_difference(target_atom_index)
            energy_terms = self._get_energy_terms(target_atom_index)
            

            flat_bottom_shift_limit = energy_terms.flat_bottom_shift_limit
            
            
            if abs(shift_diff) > flat_bottom_shift_limit:
                adjusted_shift_diff = self._adjust_shift(shift_diff, flat_bottom_shift_limit)
                
                end_harmonic = energy_terms.end_harmonic
                scale_harmonic = energy_terms.scale_harmonic
                
                
                energy_component = 0.0
                if adjusted_shift_diff < end_harmonic:
                    energy_component = (adjusted_shift_diff/scale_harmonic)**2
                else:
                    tanh_amplitude = energy_terms.tanh_amplitude
                    tanh_elongation = energy_terms.tanh_elongation
                    tanh_y_offset = energy_terms.tanh_y_offset
                    
                    tanh_argument = tanh_elongation * (adjusted_shift_diff - end_harmonic)
                    energy_component = tanh_amplitude * tanh(tanh_argument) + tanh_y_offset;

                energy += energy_component

        if self._verbose:
            end_time = time()
            print '   energy calculator: ',len(target_atom_ids),' in', "%.17g" %  (end_time-start_time), "seconds"
        
        self.calls += 1    
        return energy

cdef class Fast_force_factor_calculator(Fast_energy_calculator):

    @cython.profile(True)
    def  __call__(self, int[:] target_atom_ids):
        cdef double start_time =0.0
        cdef double end_time =0.0
        if self._verbose:
            start_time = time()

        self.set_simulation()
        #TODO: shouldn't be allocated each time
        cdef float[:] result  = allocate_array(len(target_atom_ids),'f')
        cdef int i
        cdef float factor, shift_diff
        
        for i,target_atom_id in enumerate(target_atom_ids):
            factor  = 0.0
                
            
            shift_diff = self._get_shift_difference(target_atom_id)
            energy_terms = self._get_energy_terms(target_atom_id)
            

            
            flat_bottom_shift_limit = energy_terms.flat_bottom_shift_limit
            
            if abs(shift_diff) > flat_bottom_shift_limit:
                adjusted_shift_diff = self._adjust_shift(shift_diff, flat_bottom_shift_limit)
                end_harmonic = energy_terms.end_harmonic
                scale_harmonic = energy_terms.scale_harmonic
                sqr_scale_harmonic = scale_harmonic**2
                
                weight = energy_terms.weight
                
                tanh_amplitude = energy_terms.tanh_amplitude
                tanh_elongation = energy_terms.tanh_elongation
                
                # TODO: add factor and lambda to give fact
                fact =1.0
                if adjusted_shift_diff < end_harmonic:
                    factor = 2.0 * weight * adjusted_shift_diff * fact / sqr_scale_harmonic;
                else:
                    factor = weight * tanh_amplitude * tanh_elongation / (cosh(tanh_elongation * (adjusted_shift_diff - end_harmonic)))**2.0 * fact;

            result[i] = factor
            
        if self._verbose:
            end_time = time()
            print '   force factors : ',len(target_atom_ids),' in', "%.17g" %  (end_time-start_time), "seconds"

        
        return  result
cdef class Base_force_calculator:
    
    cdef bint _verbose 
    cdef Simulation* _simulation
    cdef str _name
    
    def __init__(self,potential=None, name= "not set"):
        self._verbose = False
        self._simulation =  currentSimulation()
        
        self._name = name
        
    def _prepare(self, change, data):
        pass
    
    def set_verbose(self,on):
        self._verbose = on
    
    def set_simulation(self):
        self._simulation =  currentSimulation()
    

    
        
        
        
#    TODO should most probably be a fixed array
    @cython.profile(True)
    def __call__(self, object components, int[:] component_to_result, float[:] force_factors, Out_array forces):

        cdef double start_time =0.0
        cdef double end_time =0.0
        if self._verbose:
            start_time = time()

        self._set_components(components)
        self.set_simulation()
#        TODO rename component to result to something better
        self._do_calc_components(component_to_result, force_factors, forces)

#        component_target_atom_ids = components.get_component_atom_ids()
#        for i,target_atom_id in enumerate(target_atom_ids):
#            if target_atom_id in component_target_atom_ids:
#                try:
#                    index_range = components.get_component_range(target_atom_id)
#                except ValueError:
#                    print "warning atom is not in list", target_atom_id, Atom_utils._get_atom_info_from_index(target_atom_id),self._name
#                    index_range = []
#                for index in range(*index_range):
                    
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
        
    def __cinit__(self):
        self._compiled_components  = NULL
        self._num_components =  0
        
    def __init__(self, object indices, bint smoothed, name="Not set"):
        super(Fast_distance_based_potential_force_calculator, self).__init__(name=name)
        self._target_atom_index = indices.target_atom_index
        self._distance_atom_index_1 =  indices.distance_atom_index_1
        self._distance_atom_index_2 =  indices .distance_atom_index_2
        self._exponent_index  = indices.exponent_index
        self._coefficient_index  = indices.coefficient_index
        self._smoothed =  smoothed
        self._smoothing_factor =  DEFAULT_SMOOTHING_FACTOR
        self._cutoff =  DEFAULT_CUTOFF

    def set_cutoff(self, cutoff):
        self._cutoff =  cutoff
    
    def set_smoothing_factor(self,smoothing_factor):
        self._smoothing_factor = smoothing_factor
        
    def _set_components(self,components):
        if  self._compiled_components ==  NULL:
            self._compile_components(components) 
    
    def _free_compiled_components(self):
        if self._compiled_components != NULL:
            free (self._compiled_components)
            self._compiled_components = NULL        
            
    def _prepare(self, change, data):
         if change == TARGET_ATOM_IDS_CHANGED or change == STRUCTURE_CHANGED:
             self._free_compiled_components()   

        
    cdef _compile_components(self,components):
        self._compiled_components = <Distance_component*>malloc(len(components) * sizeof(Distance_component))
        self._num_components = len(components)
        for i,component in enumerate(components):
            self._compiled_components[i].target_atom = components[i][self._target_atom_index]
            if len(component) == 4:
                self._compiled_components[i].remote_atom_1 = components[i][self._distance_atom_index_1]
                self._compiled_components[i].remote_atom_2 = components[i][self._distance_atom_index_2]
                self._compiled_components[i].coefficient   = components[i][self._coefficient_index]
                self._compiled_components[i].exponent      = components[i][self._exponent_index]
            elif len(component) == 5:
                self._compiled_components[i].remote_atom_1 = components[i][self._distance_atom_index_1]
                self._compiled_components[i].remote_atom_2 = components[i][self._distance_atom_index_2]
                self._compiled_components[i].coefficient   = components[i][self._coefficient_index]
                self._compiled_components[i].exponent      = components[i][self._exponent_index]
            else:
                raise Exception("bad distance component length %i should be either 4 or 5 " % len(component))
   
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
        for i in range(self._num_components):
            self._distance_calc_single_force_set(i,force_factors[component_to_result[i]],force)
            
    
    @cython.profile(False)
    cdef inline float _cython_calc_single_force_factor(self, int index, float factor):
        cdef target_distant_atom atom_ids
#        cdef coefficient_exponent coef_exp
        cdef float exponent 
        
        cdef float full_factor, force_factor
        cdef float ratio, pre_exponent, reduced_exponent, sum_xyz_distances_2
        
#        atom_ids = self._get_target_and_distant_atom_ids(index)
        
        sum_xyz_distances_2 = self._sum_xyz_distances_2(self._compiled_components[index].remote_atom_1, self._compiled_components[index].remote_atom_2)
#
        full_factor= factor *  self._compiled_components[index].coefficient
        
        exponent =  self._compiled_components[index].exponent
        if self._smoothed:
            ratio = sum_xyz_distances_2 / (self._cutoff**2)
            ratio =  ratio**4
            pre_exponent = exponent - (exponent + 8.0) * ratio
        else:
            pre_exponent = exponent
            
        reduced_exponent = (exponent - 2.0) / 2.0
        
        force_factor = full_factor *  pre_exponent * sum_xyz_distances_2 ** reduced_exponent

        return force_factor


 

    
#    def _calc_single_force_set(self, int index, float factor, object forces):
#        self._cython_calc_single_force_set(index, factor, forces)
        
    @cython.profile(False)
    cdef inline void _distance_calc_single_force_set(self,int index, float factor, Out_array forces):
        
#        cdef target_distant_atom atom_ids
        cdef Vec3 xyz_distances
        cdef float force_factor, distance
        cdef Vec3 result
 
#        atom_ids =  self._get_target_and_distant_atom_ids(index)
#        
##        print atom_ids.target_atom_id,atom_ids.distant_atom_id
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
    cdef float _non_bonded_cutoff
    
    def __init__(self, object indices, bint smoothed, name = "not set"):
        global DEFAULT_NB_CUTOFF
        super(Fast_non_bonded_force_calculator, self).__init__(indices,smoothed,name=name)
        self._non_bonded_cutoff = DEFAULT_NB_CUTOFF
        
#    def _calc_single_force_set(self, int index, float factor, object forces):
#        self._cython_calc_single_force_set(index, factor, forces)
    cdef void _do_calc_components(self, int[:] component_to_result, float[:] force_factors, Out_array force):
        for i in range(self._num_components):
            self._non_bonded_calc_single_force_set(i,force_factors[component_to_result[i]],force)
            
    cdef inline void  _non_bonded_calc_single_force_set(self, int index, float factor, Out_array forces):
        cdef float distance  = calc_distance_simulation(self._simulation, self._compiled_components[index].target_atom,self._compiled_components[index].remote_atom_2)
#        TODO: this should be the non bonded distance cutoff
#TODO class variable of self are not being looked up!
        if distance < 5.0:
            self._distance_calc_single_force_set(index, factor, forces)



    
cdef class Fast_dihedral_force_calculator(Base_force_calculator):
    
    cdef Dihedral_component *_compiled_components
    cdef int _num_components
    
    def __cinit__(self):
        self._compiled_components = NULL
        self._num_components = 0
        
    def __init__(self,name="not set"):
        Base_force_calculator.__init__(self,name=name)
        
        
    def _set_components(self,components):
        if self._compiled_components ==  NULL:
            self._compile_components(components)
            
    def _compile_components(self, components):
        self._compiled_components = <Dihedral_component*>malloc(len(components) * sizeof(Dihedral_component))
        self._num_components = len(components) 
        for i,component in enumerate(components):
            self._compiled_components[i].target_atom = component[0]
            for j in range(4):
                self._compiled_components[i].dihedral_atoms[j] = component[1+j]
            self._compiled_components[i].coefficient = component[5]
            for j in range(5):
                self._compiled_components[i].parameters[j] = component[6+j]
    
    def _free_compiled_components(self):
        if self._compiled_components != NULL:
            free(self._compiled_components)
            self._compiled_components = NULL        
                   
    def _prepare(self, change, data):
         if change == TARGET_ATOM_IDS_CHANGED or change == STRUCTURE_CHANGED:
             self._free_compiled_components() 
                     
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
        cdef int i
        for i in range(self._num_components):
            self._dihedral_calc_single_force_set(i,force_factors[component_to_result[i]],force)
            

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

#        print vec3_as_tuple(r1), vec3_as_tuple(r2), vec3_as_tuple(r3), vec3_as_tuple(r4)
        
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

#        print vec3_as_tuple(F1), vec3_as_tuple(F2), vec3_as_tuple(F3), vec3_as_tuple(F4)
#        // assign forces
#        force_triplet_1 = self._get_or_make_target_force_triplet(forces, atom_ids.atom_id_1)
#        force_triplet_2 = self._get_or_make_target_force_triplet(forces, atom_ids.atom_id_2)
#        force_triplet_3 = self._get_or_make_target_force_triplet(forces, atom_ids.atom_id_3)
#        force_triplet_4 = self._get_or_make_target_force_triplet(forces, atom_ids.atom_id_4)
        
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

    cdef Coef_components _compiled_coef_components

    
    def __cinit__(self):
        self._compiled_components = NULL
        self._num_components = 0
        
        self._compiled_ring_components = NULL
        self._num_ring_components = 0

        #note this is not a raw array of structs it's a compiled python class
        self._compiled_coef_components = None        
        
        self._centre_cache = None
        self._normal_cache = None

    def __init__(self,name="not set"):
        super(Fast_ring_force_calculator, self).__init__(name=name)
        

         
    def _set_components(self,components):
        if  self._compiled_components ==  NULL:
            self._compile_components(components)
            
    def _compile_components(self,components): 
        self._compiled_components = <Ring_target_component*>malloc(len(components) * sizeof(Ring_target_component))
        self._num_components = len(components)
        for i,component in enumerate(components):
            self._compiled_components[i].target_atom_id = components[i][0]
            self._compiled_components[i].atom_type_id = components[i][1]
                
        
    def _set_coef_components(self,coef_components):
        self._compiled_coef_components = Coef_components(coef_components)
            
    def _set_ring_components(self,ring_components):
        if  self._compiled_ring_components ==  NULL:
            self._compile_ring_components(ring_components)
            
    def _compile_ring_components(self,ring_components): 
        self._compiled_ring_components = <Ring_component*>malloc(len(ring_components) * sizeof(Ring_component))
        self._num_ring_components = len(ring_components)
        for i,ring_component in enumerate(ring_components):
            self._compiled_ring_components[i].ring_id  = ring_components[i][0]
            self._compiled_ring_components[i].num_atoms  = len(ring_component[1])
            for j in range(len(ring_component[1])):
                self._compiled_ring_components[i].atom_ids[j] =ring_component[1][j]
    
    def _free_compiled_ring_components(self):   
        if self._compiled_ring_components != NULL:
            free(self._compiled_ring_components)
            self._compiled_ring_components = NULL

    def _free_coef_components(self):
        self._compiled_coef_components = None

    cdef _free_compiled_components(self):
        if self._compiled_components != NULL:
            free (self._compiled_components)
            self._compiled_components = NULL   
            
    def _prepare(self, change, data):
        if change == TARGET_ATOM_IDS_CHANGED or change == STRUCTURE_CHANGED:
            self._free_compiled_components() 
            self._free_compiled_ring_components()
            self._free_coef_components()
                        
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
        

#    @cython.profile(False)
#    cdef inline target_type _get_target_and_type(self,int index):
#        component = self._get_component(index)
##        TODO: this needs to contain better names a union?
#        cdef target_type result
#        result.target_atom_id = component[0]
#        result.atom_type_id = component[1]
#        
#        return result#
    
    
    cdef inline void _do_calc_components(self, int[:] component_to_result, float[:] force_factors, Out_array force):
        for i in range(self._num_components):
            self._ring_calc_single_force_set(i,force_factors[component_to_result[i]],force)
            
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
        
#        return  dL3, u, dL6, ring_normal,  atom_type_id, coefficient, d, factor, dn, dL3nL3, dL, nL, dLnL
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
        operator_times(nSum,2.0)  #            float_type g [3], ab [3], c [3]
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
#                print AXIS_NAMES[axis],sub_force,-force_factor, gradU[axis], force_terms.dL3, force_terms.u, gradV[axis],force_terms.dL6
#            print
##        print 

cdef class Fast_non_bonded_shift_calculator(Fast_distance_shift_calculator):
    
    cdef float _nb_cutoff
    
    def __init__(self, object indices, bint smoothed, str name):
        super(Fast_non_bonded_shift_calculator, self).__init__(indices,smoothed,name)
        global DEFAULT_NB_CUTOFF
        self._nb_cutoff = DEFAULT_NB_CUTOFF
        
    def set_verbose(self,on):
        self._verbose = on
    
    @cython.profile(True)
    def __call__(self, object components, double[:] results, int[:] component_to_target):
        self._set_components(components)
        self.set_simulation()
        
        cdef float default_smoothing_factor = self._smoothing_factor
        cdef float smoothing_factor
        cdef float ratio
        cdef float result
        cdef float distance
        cdef target_distant_atom atom_indices
        cdef coefficient_exponent coef_exp
        cdef int target_atom_id
        cdef int distant_atom_id
        cdef int result_index
        cdef double value 
        
        cdef double start_time =0.0
        cdef double end_time =0.0
        if self._verbose:
            start_time = time()
            
        for index in range(self._num_components):
            target_atom_id = self._compiled_components[index].remote_atom_1
            distant_atom_id  = self._compiled_components[index].remote_atom_2
            
            distance = calc_distance_simulation(self._simulation, target_atom_id, distant_atom_id)
            if distance < self._nb_cutoff:
        
                coef_exp = self._get_coefficient_and_exponent(index)
                smoothing_factor = default_smoothing_factor
                if self._smoothed:
                    ratio = distance / self._cutoff
                    smoothing_factor = 1.0 - ratio ** 8
                results[component_to_target[index]]  += smoothing_factor * pow(distance,  coef_exp.exponent) * coef_exp.coefficient
        
        if self._verbose:
            end_time = time()
            print '   distance shift components ' ,self._name,len(components), 'in', "%.17g" % (end_time-start_time), "seconds"


