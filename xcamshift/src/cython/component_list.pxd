cdef class Component_list:
    cdef object _components
    cdef object _component_offsets
    cdef object _component_ids
    

cdef class Native_component_list(Component_list):
    
    cdef object _translator
    cdef object _preparer
    cdef object _component_struct
    
    cdef object _native_components
    cdef object _native_component_offsets
    
    cdef get_native_components(self)
