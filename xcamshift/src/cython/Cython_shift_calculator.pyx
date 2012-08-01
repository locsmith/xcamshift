'''
Created on 31 Jul 2012

@author: garyt
'''
from xplor_access cimport norm,Vec3,currentSimulation, Dihedral, Atom
from libc.math cimport cos

cdef class Fast_distance_shift_calculator:
    DEFAULT_CUTOFF  = 5.0
    DEFAULT_SMOOTHING_FACTOR = 1.0
    
    cdef int _target_atom_index
    cdef int _distance_atom_index_1
    cdef int _distance_atom_index_2
    cdef int _exponent_index
    cdef int _coefficient_index
    cdef bint _smoothed 
    cdef object _components
    cdef float _smoothing_factor 
    cdef float _cutoff
    
    def __init__(self, indices, smoothed):
        self._target_atom_index = indices.target_atom_index
        self._distance_atom_index_1 =  indices.distance_atom_index_1
        self._distance_atom_index_2 =  indices .distance_atom_index_2
        self._exponent_index  = indices.exponent_index
        self._coefficient_index  = indices.coefficient_index
        
        self._smoothed =  smoothed
        self._smoothing_factor =  self.DEFAULT_SMOOTHING_FACTOR
#        
        self._components =  None
        self._cutoff =  self.DEFAULT_CUTOFF
#    
#    def set_cutoff(self, cutoff):
#        self._cutoff =  cutoff
#    
#    def set_smoothing_factor(self,smoothing_factor):
#        self._smoothing_factor = smoothing_factor
#    
#    TODO: this needs to be removed
    cdef _set_components(self,components):
        self._components = components
    
    cdef inline _get_target_and_distant_atom_ids(self, index):
        values  = self._components.get_component(index)
        
        target_atom = values[self._distance_atom_index_1]
        distance_atom = values[self._distance_atom_index_2]
        return target_atom, distance_atom
    
    cdef inline _get_coefficient_and_exponent(self, index):
        values = self._components.get_component(index)
        
        coefficient = values[self._coefficient_index]
        exponent = values[self._exponent_index]
        
        return coefficient, exponent
#    
    def __call__(self, object components, object results):
        self._set_components(components)
        cdef float smoothing_factor = self._smoothing_factor
        cdef float ratio
        cdef float result
        for index in range(len(components)):
            target_atom_index, sidechain_atom_index = self._get_target_and_distant_atom_ids(index)
            coefficient, exponent = self._get_coefficient_and_exponent(index)
            distance =self.distance(target_atom_index, sidechain_atom_index)
    #        Atom_utils._calculate_distance(target_atom_index, sidechain_atom_index)

            if self._smoothed:
                ratio = distance / self._cutoff
                smoothing_factor = 1.0 - ratio ** 8
            results[index]  = smoothing_factor * distance ** exponent * coefficient
            

    cdef float distance(self,int atom_index_1, atom_index_2):

        vec1 = currentSimulation().atomPosArr().data(atom_index_1)
        vec2 = currentSimulation().atomPosArr().data(atom_index_2)
        cdef Vec3 result =  vec2 - vec1
        return norm(result)
    
cdef class Fast_dihedral_shift_calculator:
    cdef object _components
    
    def __init__(self):
        self._components = None
    
    cdef _set_components(self, object components):
        self._components = components
        
    cdef inline _get_component(self,int index):
        return self._components[index]
    
    cdef inline _get_dihedral_angle(self, int dihedral_1_atom_id_1, int dihedral_1_atom_id_2, 
                                          int dihedral_2_atom_id_1, int dihedral_2_atom_id_2):
        cdef Atom atom1,atom_2, atom_3, atom_4
        

        atom_1  = currentSimulation().atomByID(dihedral_1_atom_id_1)
        atom_2  = currentSimulation().atomByID(dihedral_1_atom_id_2)
        atom_3  = currentSimulation().atomByID(dihedral_2_atom_id_1)
        atom_4  = currentSimulation().atomByID(dihedral_2_atom_id_2)
        
        return Dihedral(atom_1,atom_2,atom_3,atom_4).value()


    cdef inline _get_dihedral_atom_ids(self, int index):
        component = self._get_component(index)
        
        return component[1], component[2],component[3],component[4]


    cdef inline _get_coefficient(self, int index):
        component = self._get_component(index)
        coefficient = component[5]
        return coefficient


    cdef inline _get_parameters(self, int index):
        component = self._get_component(index)
        parameters = component[6:11]
        parameter_0, parameter_1, parameter_2, parameter_3, parameter_4 = parameters
        return parameter_0, parameter_3, parameter_1, parameter_4, parameter_2

    def __call__(self, object components, object results):
        cdef float angle
        cdef float angle_term
        cdef float shift
        cdef float coefficient
        cdef float parameter_0, parameter_1, \
                   parameter_2, parameter_3, \
                   parameter_4
        cdef int dihedral_atom_id_1, dihedral_atom_id_2,dihedral_atom_id_3, dihedral_atom_id_4,  
        
        self._set_components(components)
        for index in range(len(components)):
            dihedral_atom_id_1, dihedral_atom_id_2, \
            dihedral_atom_id_3, dihedral_atom_id_4 = self._get_dihedral_atom_ids(index)
            
            coefficient = self._get_coefficient(index)
            
            parameter_0, parameter_3, parameter_1, \
            parameter_4, parameter_2               \
            = self._get_parameters(index)
            
            angle = self._get_dihedral_angle(dihedral_atom_id_1, dihedral_atom_id_2,
                                             dihedral_atom_id_3, dihedral_atom_id_4)
    
            
            angle_term = parameter_0 * cos(3.0 * angle + parameter_3) + \
                         parameter_1 * cos(angle + parameter_4) +       \
                         parameter_2
            
            shift = coefficient * angle_term
    
            results[index] = shift
