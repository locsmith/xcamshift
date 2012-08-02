'''
Created on 31 Jul 2012

@author: garyt
'''
from vec3 import Vec3 as python_vec3
from  xplor_access cimport norm,Vec3,currentSimulation, Dihedral, Atom,  dot,  cross
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
 
cdef class Fast_ring_shift_calculator:
    cdef object _components
    cdef object _coef_components
    cdef object _ring_components
    cdef object _centre_cache
    cdef object _normal_cache
    
    
    def __init__(self):
        self._components = None
        self._coef_components = None
        self._ring_components = None
        self._centre_cache = None
        self._normal_cache = None
        
    def _set_components(self,components):
        self._components = components
        
    def _set_coef_components(self,coef_components):
        self._coef_components =  coef_components
            
    def _set_ring_components(self,coef_components):
        self._ring_components =  coef_components
    
    def _set_normal_cache(self,normals):
        self._normal_cache = normals
        
    def _set_centre_cache(self,centres):
        self._centre_cache = centres

    cdef _get_coef_components(self, int atom_type_id):
        return self._coef_components.get_components_for_atom_id(atom_type_id)
            
    cdef Vec3 _get_ring_normal(self, int ring_id):
        cdef float x, y ,z 
        x,y,z =  self._normal_cache.get_component(ring_id)[1]
        return  Vec3(x,y,z)
    
    cdef Vec3 _get_ring_centre(self, int ring_id):
        cdef float x, y ,z 
        x,y,z = self._centre_cache.get_component(ring_id)[1]
        return  Vec3(x,y,z)
    
    cdef float  _calc_sub_component_shift(self, int target_atom_id, int ring_id, float coefficient):
        
        cdef Vec3 target_atom_pos
        cdef Vec3 ring_centre
        cdef Vec3 ring_normal
        cdef float lenght_normal
        cdef Vec3 direction_vector
        cdef float distance, distance3, angle, contrib
        
        target_atom_pos =  currentSimulation().atomPosArr().data(target_atom_id)
        ring_centre = self._get_ring_centre(ring_id)

        #TODO add this to a cache the same way that camshift does
        ring_normal = self._get_ring_normal(ring_id)
        length_normal = norm(ring_normal)
    
        #correct name?
        direction_vector = target_atom_pos - ring_centre
        
        distance = norm(direction_vector)
        distance3 = distance ** 3
        
        angle = dot(direction_vector, ring_normal) / (distance * length_normal)
        contrib = (1.0 - 3.0 * angle ** 2) / distance3
        
#        print Atom_utils._get_atom_info_from_index(target_atom_id), ring_id, angle, distance3, coefficient, contrib * coefficient, angle
        return contrib * coefficient

    
    def __call__(self, components, results):
        cdef int target_atom_id
        cdef int ring_id
        cdef float coefficient
        
        self._set_components(components)
        
        for index in range(len(components)):
            component = components[index]
            atom_type_id = component[1]

            shift = 0.0
        
            
            for coef_component in self._get_coef_components(atom_type_id):
#                TODO: remove magic numbers or add structs

                target_atom_id = component[0]
                ring_id = coef_component[1]
                coefficient = coef_component[2]
                shift += self._calc_sub_component_shift(target_atom_id,  ring_id, coefficient)
            
            results[index] = shift
#   
cdef int RING_ATOM_IDS = 1
    
cdef class Fast_ring_data_calculator:

    def __init__(self):
        pass
    
    cdef Vec3 _calculate_one_ring_centre(self, ring_component):

        atom_ids = ring_component[RING_ATOM_IDS]
        cdef float num_atom_ids = len(atom_ids)
        cdef Vec3 total,result

        total=Vec3(0.0,0.0,0.0)
        for atom_id in atom_ids:
            total+=currentSimulation().atomPosArr().data(atom_id)

#        TODO get operator / workin properly
        result = Vec3(total.x()/num_atom_ids,total.y()/num_atom_ids,total.z()/num_atom_ids)
        
        return result
    
    
#    def _check_ring_size_ok(self, atom_ids):
#        num_atom_ids = len(atom_ids)
#        if num_atom_ids < 5 or num_atom_ids > 6:
#            template = "ring normals function is only implemented for 5 or six member rings i got %d atoms"
#            msg = template % num_atom_ids
#            raise Exception(msg)
#

    cdef Vec3  _calculate_normal(self,int[3] atom_id_triplet):
    
        cdef Vec3 normal
        cdef Vec3 vec_1 
        cdef Vec3 vec_2

        cdef Vec3 atom_vector_1, atom_vector_2, atom_vector_3
        
        atom_vector_1 = currentSimulation().atomPosArr().data(atom_id_triplet[0])
        atom_vector_2 = currentSimulation().atomPosArr().data(atom_id_triplet[1])
        atom_vector_3 = currentSimulation().atomPosArr().data(atom_id_triplet[2])

        vec_1 =  atom_vector_1 -atom_vector_2
        vec_2 =  atom_vector_3- atom_vector_2
            
        return cross(vec_1,vec_2)
   
    #TODO could try newells method http://www.opengl.org/wiki/Calculating_a_Surface_Normal
    cdef Vec3 _calculate_one_ring_normal(self, ring_component):

        atom_ids = ring_component[RING_ATOM_IDS]
#        self._check_ring_size_ok(atom_ids)
        
        cdef int[3] atom_triplet_1
        cdef int[3] atom_triplet_2
        self._build_atom_triplet(atom_ids[:3], atom_triplet_1)
        self._build_atom_triplet(atom_ids[-3:], atom_triplet_2)
       
        
        normal_1 = self._calculate_normal(atom_triplet_1)
        normal_2 = self._calculate_normal(atom_triplet_2)
        return self. _average_2_vec_3(normal_1, normal_2)
       
    cdef inline  _build_atom_triplet(self,atom_ids, int[3]& result):
        for i in range(3):
            result[i] = atom_ids[i]
        
        
    cdef inline Vec3 _average_2_vec_3(self, Vec3& vector_1, Vec3& vector_2):
            
        cdef Vec3 sum
        
        sum =  vector_1 + vector_2
            
        return Vec3(sum.x()/2.0, sum.y()/2.0, sum.z()/2.0)

    
            
    def __call__(self, rings, normals, centres):
        for ring_component in rings:
            ring_id = ring_component[0]

            centre = self._calculate_one_ring_centre(ring_component)
            result_centre =  python_vec3(centre.x(),centre.y(),centre.z())
            centre_component = ring_id, result_centre
            centres.add_component(centre_component)
            
            normal = self._calculate_one_ring_normal(ring_component)
            result_normal =  python_vec3(normal.x(),normal.y(),normal.z())
            normal_component = ring_id, result_normal
            normals.add_component(normal_component)
            

