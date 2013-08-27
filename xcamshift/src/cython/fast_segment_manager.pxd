

cdef class Segment_Manager:
    cdef int _num_atoms
    cdef object _segments, _segment_info_map
    
    cdef int cython_get_number_atoms(Segment_Manager self) nogil
    