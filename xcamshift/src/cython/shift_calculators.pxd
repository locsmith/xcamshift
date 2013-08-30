from  xplor_access cimport  CDSVector, EnsembleSimulation, Vec3

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
 
   