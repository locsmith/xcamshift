from  xplor_access cimport  CDSVector, EnsembleSimulation, Vec3, Simulation
from libcpp.map cimport map as cmap
from libc.stdint cimport uintptr_t

cdef int ATOM = 0
cdef int NATOM = 1
cdef int SIMU = 2
cdef int NBRM = 3
cdef int NNBRM = 4
cdef int COEF =  5
cdef int NCOEF =  6
cdef int NBLT = 7
cdef int NNBLT = 8
cdef int OFFS = 9

cdef class Out_array:
    cdef long _length
    cdef double[60000] _data
    cdef int[60000] _mask
    cdef EnsembleSimulation* _simulation
    
    cpdef get_length(Out_array self)
    cpdef add_all(Out_array self, Out_array target_array)
    cdef inline bint _check_clear_length(Out_array self, long length) nogil
    cdef inline bint _check_offset(Out_array self, long offset)
    cdef bint realloc(Out_array self,long length) nogil
    cdef inline bint _check_state(Out_array self)
    cdef void  _clear(Out_array self, int length) nogil
    cdef inline void  add(Out_array self,long offset, Vec3& value)

    

cdef struct Constant_cache:      
    int     target_atom_id      
    float   flat_bottom_shift_limit
    float   end_harmonic
    float   scale_harmonic
    float   weight
    float   tanh_amplitude
    float   tanh_elongation
    float   tanh_y_offset
    
cdef class Fast_energy_calculator_base:
    cdef Constant_cache* _energy_term_cache 
    cdef CDSVector[double]  *_theory_shifts
    cdef CDSVector[float] _observed_shifts
    cdef bint _verbose 
    cdef int calls
    
        
    cdef Constant_cache* _get_energy_terms(Fast_energy_calculator_base self, int target_atom_index) nogil
    cdef inline float  _get_calculated_atom_shift(Fast_energy_calculator_base self, int index) nogil
    cdef inline float _get_observed_atom_shift(Fast_energy_calculator_base self, int index) nogil
    cdef inline float  _get_shift_difference(Fast_energy_calculator_base self, int target_atom_index, int index) nogil
    cdef inline float _adjust_shift(Fast_energy_calculator_base self, float shift_diff, float flat_bottom_shift_limit) nogil
    cdef void set_observed_shifts(Fast_energy_calculator_base, CDSVector[float] observed_shifts) nogil
    cdef void set_calculated_shifts(self, CDSSharedVectorFloat calculated_shifts) nogil

cdef class Fast_energy_calculator(Fast_energy_calculator_base):
    cdef        float calcEnergy(Fast_energy_calculator self, CDSVector[int] target_atom_ids, CDSVector[int] *active_atom_ids=?)
    cdef inline float _calc_one_energy(Fast_energy_calculator self, int target_atom_index, int index)
    
        
cdef class Fast_force_factor_calculator(Fast_energy_calculator_base):
    cdef void calcFactors(Fast_force_factor_calculator self, CDSVector[int] target_atom_ids, float[:] result, CDSVector[int] *active_atom_ids=?) nogil
    cdef inline float _calc_one_force_factor(Fast_force_factor_calculator self, int target_atom_id, int i) nogil


cdef class CDSSharedVectorFloat:
    cdef CDSVector[double]*  data

    cdef void resize(self,int size) nogil
    cdef CDSVector[double]* get_data(self) nogil
    cdef int size(self) nogil
 
cdef class Base_shift_calculator:
    cdef bint _verbose 
    cdef str _name
    cdef EnsembleSimulation* _simulation
    
    cdef void _bytes_to_components(self, uintptr_t data, uintptr_t length) nogil
    cdef void _set_components(self, cmap[int, uintptr_t] *components) nogil
    cdef inline int ensemble_size(self) nogil
    cdef inline int ensemble_member_index(self) nogil
    cdef inline int ensemble_array_offset(self,int index) nogil
    
    cdef void calc(Base_shift_calculator self, cmap[int,uintptr_t] *components, CDSSharedVectorFloat shift_cache, int[:] component_to_target,  int[:] active_components) nogil

cdef class Base_force_calculator:
    
    cdef bint _verbose 
    cdef Simulation* _simulation
    cdef str _name
    cdef int[:] _active_components 
    
    cdef void _bytes_to_components(self, uintptr_t data, uintptr_t length) nogil
    cdef void _set_components(self, cmap[int, uintptr_t] *components) nogil
    cdef void calc(self, cmap[int,uintptr_t] *components, int[:] component_to_result, float[:] force_factors, Out_array forces, int[:] active_components=?) nogil
    cdef void _do_calc_components(self, int[:] component_to_result,  float[:] force_factors, Out_array force) nogil

cdef struct Component_index_pair:
    int target_atom_id
    int target_index
    int remote_index
    int component_index
    
cdef  class Non_bonded_interaction_list:
    cdef CDSVector[int]  *data
    cdef int length
    cdef int size_increment
    cdef int RECORD_LENGTH 
 
    cdef inline void append(self, int target_atom_id, int target_id, int remote_id, int component_index)
    cdef inline void resize(self)   
    cdef inline Component_index_pair* get(self,int offset) nogil     
    cdef int get_length(self) nogil

