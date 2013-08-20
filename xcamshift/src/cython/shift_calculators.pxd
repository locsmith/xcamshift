from  xplor_access cimport  CDSVector

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
    cdef float[:] _observed_shifts
    cdef bint _verbose 
    cdef int calls
    
        
    cdef Constant_cache* _get_energy_terms(Fast_energy_calculator_base self, int target_atom_index) nogil
    cdef inline float  _get_calculated_atom_shift(Fast_energy_calculator_base self, int index) nogil
    cdef inline float _get_observed_atom_shift(Fast_energy_calculator_base self, int index) nogil
    cdef inline float  _get_shift_difference(Fast_energy_calculator_base self, int target_atom_index, int index) nogil
    cdef inline float _adjust_shift(Fast_energy_calculator_base self, float shift_diff, float flat_bottom_shift_limit) nogil

cdef class Fast_energy_calculator(Fast_energy_calculator_base):
    cdef        float calcEnergy(Fast_energy_calculator self, CDSVector[int] target_atom_ids, CDSVector[int] *active_atom_ids=?)
    cdef inline float _calc_one_energy(Fast_energy_calculator self, int target_atom_index, int index)
    
        
cdef class Fast_force_factor_calculator(Fast_energy_calculator_base):
    cdef void calcFactors(Fast_force_factor_calculator self, CDSVector[int] target_atom_ids, float[:] result, CDSVector[int] *active_atom_ids=?) nogil
    cdef inline float _calc_one_force_factor(Fast_force_factor_calculator self, int target_atom_id, int i) nogil


   
   