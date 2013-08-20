from xplor_access cimport  CDSVector
from libcpp.map cimport map 

cdef class Observed_shift_table(object):

    cdef map[int,float] _chemical_shifts 
    cdef CDSVector[float] _native_shifts
    cdef bint _native_shifts_set

    cdef CDSVector[float] get_native_shifts(self, CDSVector[int] target_atom_ids) nogil
