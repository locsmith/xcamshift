# cython: profile=True
'''
Created on 31 Jul 2012

@author: garyt
'''
from vec3 import Vec3 as python_vec3
from  xplor_access cimport norm,Vec3,currentSimulation, Dihedral, Atom,  dot,  cross
from libc.math cimport cos,sin,  fabs, tanh, pow, cosh
from libc.stdlib cimport malloc, free
from libc.string cimport strcmp

 
#cdef object vec3_as_tuple(Vec3& vec_3):
#    return vec_3.x(), vec_3.y(), vec_3.z()
    
cdef class Vec3_container:
    
    cdef float[3] floats
    
    cdef set_vec3(self,Vec3& vec3):
        self.floats[0] =  vec3.x()
        self.floats[1] =  vec3.y()
        self.floats[2] =  vec3.z()
        
        
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
        
cdef class Coef_components:
    cdef int num_components
    cdef Coef_component* _components
    
    def __init__ (self, object components):
        self.num_components = len(components)
        for i,component in enumerate(components):
            self._components[i].atom_type_id = component[0] 
            self._components[i].ring_id = component[1]
            self._components[i].coefficient = component[2]
        
    def __cinit__(self, object components):
        self._components = <Coef_component *>malloc(len(components) * sizeof(Coef_component))
        if not self._components:
            raise MemoryError()
    
    def __dealloc__(self):
        free(self._components)
        self._components = NULL
    
    cdef int _check_offset(self, int offset) except -1:
    
        if not offset < self.num_components:
            raise IndexError("tried to access object at %i length is %i" % (offset,self._num_components))
        if self._components is NULL:
            raise AttributeError("trying to access deallocated memory!")
        return 0

    
#    cdef set_component(self,int offset, coef_component& component):
#        self._check_offset(offset)
#        self._components[offset] = component

    cdef inline Coef_component* get_component(self,int offset):
        self._check_offset(offset)
        return &self._components[offset] 

cdef struct ring_atom_ids:
    int num_atoms
    int[6] atom_ids



    

cdef class Ring_atom_positions:
    cdef int num_vecs
    cdef Vec3* _vecs
    
    def __init__ (self, int num_vecs):
        self.num_vecs = num_vecs

        
    def __cinit__(self, int num_vecs):
        self._vecs = <Vec3 *>malloc(num_vecs * sizeof(Vec3))
        if not self._vecs:
            raise MemoryError()
    
    def __dealloc__(self):
        free(self._vecs)
        self._vecs = NULL
    
    cdef inline int _check_offset(self, int offset) except -1:
        if not offset < self.num_vecs:
            raise IndexError("tried to access object at %i length is %i" % (offset,self._num_vecs))
        if self._vecs is NULL:
            raise AttributeError("trying to access deallocated memory!")

    
    cdef set_vec(self,int offset, Vec3& vec):
        self._check_offset(offset)
        self._vecs[offset] = vec

    cdef Vec3* get_vec(self,int offset):
        self._check_offset(offset)
        return &self._vecs[offset] 

cdef struct dihedral_parameters:
    float param_0
    float param_1
    float param_2
    float param_3
    float param_4
    
    
#cdef float sum(Vec3 vec3):
#    return vec3.x()+vec3.y()+vec3.z()

cdef inline  Vec3 operator_times (Vec3& vec3, float scale):
    return Vec3(vec3.x() * scale, vec3.y() * scale, vec3.z() * scale)

cdef inline float calc_distance(int atom_index_1, atom_index_2):

    vec1 = currentSimulation().atomPosArr().data(atom_index_1)
    vec2 = currentSimulation().atomPosArr().data(atom_index_2)
    cdef Vec3 result =  vec2 - vec1
    return norm(result)

cdef inline float calc_dihedral_angle(dihedral_ids dihedral_atom_ids):
    
    

    atom_1  = currentSimulation().atomByID(dihedral_atom_ids.atom_id_1)
    atom_2  = currentSimulation().atomByID(dihedral_atom_ids.atom_id_2)
    atom_3  = currentSimulation().atomByID(dihedral_atom_ids.atom_id_3)
    atom_4  = currentSimulation().atomByID(dihedral_atom_ids.atom_id_4)
    
    return Dihedral(atom_1,atom_2,atom_3,atom_4).value()



cdef float DEFAULT_CUTOFF = 5.0
cdef float DEFAULT_SMOOTHING_FACTOR = 1.0
cdef float DEFAULT_NB_CUTOFF = 5.0
    
cdef class Base_shift_calculator:
    cdef bint _verbose 
    cdef str _name
    def __init__(self, str name):
            self._verbose = False
            self._name = name
            
    def is_fast(self):
        return True
    
    def set_verbose(self,bint state):
        self._verbose =  state
    
    cdef _do_verbose(self,object components, object targets):
        if self._verbose:
            msg='%s calculating %i components for %s targets'
            args =  (self._name, len(components), len(targets))
            print msg % args   
        
cdef class Fast_distance_shift_calculator(Base_shift_calculator):

    
    cdef int _target_atom_index
    cdef int _distance_atom_index_1
    cdef int _distance_atom_index_2
    cdef int _exponent_index
    cdef int _coefficient_index
    cdef bint _smoothed 
    cdef object _components
    cdef float _smoothing_factor 
    cdef float _cutoff
    
    def __init__(self, indices, smoothed, str name = "not set"):
        Base_shift_calculator.__init__(self, name)
        self._target_atom_index = indices.target_atom_index
        self._distance_atom_index_1 =  indices.distance_atom_index_1
        self._distance_atom_index_2 =  indices .distance_atom_index_2
        self._exponent_index  = indices.exponent_index
        self._coefficient_index  = indices.coefficient_index
        
        self._smoothed =  smoothed
        self._smoothing_factor =  DEFAULT_SMOOTHING_FACTOR
#        
        self._components =  None
        self._cutoff =  DEFAULT_CUTOFF
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
    
    cdef inline target_distant_atom _get_target_and_distant_atom_ids(self, int index):
        cdef object values 
        cdef target_distant_atom result 
        values  = self._components.get_component(index)
        
        result.target_atom_id = values[self._distance_atom_index_1]
        result.distant_atom_id  = values[self._distance_atom_index_2]
        return result
    
    cdef inline coefficient_exponent _get_coefficient_and_exponent(self, int index):
        cdef object values
        cdef coefficient_exponent result
        values = self._components.get_component(index)
        
        result.coefficient = values[self._coefficient_index]
        result.exponent = values[self._exponent_index]
        
        return result
#    
    def __call__(self, object components, object results, object component_to_target):
        self._do_verbose(components, results)
        
        self._set_components(components)
        cdef float smoothing_factor = self._smoothing_factor
        cdef float ratio
        cdef float result
        cdef target_distant_atom atom_indices
        cdef coefficient_exponent coef_exp
        for index in range(len(components)):
            atom_indices = self._get_target_and_distant_atom_ids(index)
            
            coef_exp = self._get_coefficient_and_exponent(index)
            distance =calc_distance(atom_indices.target_atom_id, atom_indices.distant_atom_id)
    #        Atom_utils._calculate_distance(target_atom_index, sidechain_atom_index)

            if self._smoothed:
                ratio = distance / self._cutoff
                smoothing_factor = 1.0 - ratio ** 8
            results[component_to_target[index]]  += smoothing_factor * pow(distance,  coef_exp.exponent) * coef_exp.coefficient
            
    

    
cdef class Fast_dihedral_shift_calculator(Base_shift_calculator):
    cdef object _components
    
    def __init__(self, str name = "not set"):
        Base_shift_calculator.__init__(self,name)
        self._components = None
    
    cdef _set_components(self, object components):
        self._components = components
        
    cdef inline _get_component(self,int index):
        return self._components[index]
    

    cdef inline dihedral_ids _get_dihedral_atom_ids(self, int index):
        cdef dihedral_ids result
        
        component = self._get_component(index)
        
        result.atom_id_1= component[1]
        result.atom_id_2= component[2]
        result.atom_id_3= component[3]
        result.atom_id_4= component[4]
        
        return result

    cdef inline float _get_coefficient(self, int index):
        component = self._get_component(index)
        coefficient = component[5]
        return coefficient


    cdef inline dihedral_parameters _get_parameters(self, int index):
        component = self._get_component(index)
        
        cdef dihedral_parameters result
        
        parameters = component[6:11]
        
        result.param_0 =  parameters[0]
        result.param_1 =  parameters[3]
        result.param_2 =  parameters[1]
        result.param_3 =  parameters[4]
        result.param_4 =  parameters[2]
        
            
        return result

    def __call__(self, object components, object results, object component_to_target):
        cdef float angle
        cdef float angle_term
        cdef float shift
        cdef float coefficient
        cdef dihedral_parameters parameters
        cdef dihedral_ids dihedral_atom_ids  
        
        self._do_verbose(components,results)
        
        self._set_components(components)
        for index in range(len(components)):
            dihedral_atom_ids = self._get_dihedral_atom_ids(index)
            
            coefficient = self._get_coefficient(index)
            
            parameters = self._get_parameters(index)
            
            angle = calc_dihedral_angle(dihedral_atom_ids)
    
            angle_term = parameters.param_0 * cos(3.0 * angle + parameters.param_1) + \
                         parameters.param_2 * cos(angle +  parameters.param_3) +      \
                         parameters.param_4
            shift = coefficient * angle_term
    
            results[component_to_target[index]] += shift
 
cdef class Fast_ring_shift_calculator(Base_shift_calculator):
    cdef object _components
    cdef object _coef_components
    cdef object _ring_components
    cdef object _centre_cache
    cdef object _normal_cache
    
    
    def __init__(self, str name = "not set"):
        Base_shift_calculator.__init__(self,name)
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
        cdef Vec3_container container =  self._normal_cache.get_component(ring_id)[1]
        return  container.get_vec3()
    
    cdef Vec3 _get_ring_centre(self, int ring_id):
        cdef Vec3_container container =  self._centre_cache.get_component(ring_id)[1]
        return  container.get_vec3()
    
    cpdef float  _calc_sub_component_shift(self, int target_atom_id, int ring_id, float coefficient):
        
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

    
    def __call__(self, object components, object results, object component_to_target):
        cdef int target_atom_id
        cdef int ring_id
        cdef float coefficient
        
        self._do_verbose(components,results)
        
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
            
            results[component_to_target[index]] += shift
#   
cdef int RING_ATOM_IDS = 1
    
cdef class Fast_ring_data_calculator:
    cdef bint _verbose
    
    def __init__(self): 
        self._verbose = False
    
    def set_verbose(self,on):
        self._verbose =  on
        
    cdef Vec3_container  _calculate_one_ring_centre(self, ring_component):

        atom_ids = ring_component[RING_ATOM_IDS]
        cdef float num_atom_ids = len(atom_ids)
        cdef Vec3 total,average

        total=Vec3(0.0,0.0,0.0)
        for atom_id in atom_ids:
            total+=currentSimulation().atomPosArr().data(atom_id)

#        TODO get operator / workin properly
        average = Vec3(total.x()/num_atom_ids,total.y()/num_atom_ids,total.z()/num_atom_ids)
        
        cdef Vec3_container result = Vec3_container()
        result.set_vec3(average)
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
    cdef Vec3_container _calculate_one_ring_normal(self, ring_component):

        atom_ids = ring_component[RING_ATOM_IDS]
#        self._check_ring_size_ok(atom_ids)
        
        cdef int[3] atom_triplet_1
        cdef int[3] atom_triplet_2
        self._build_atom_triplet(atom_ids[:3], atom_triplet_1)
        self._build_atom_triplet(atom_ids[-3:], atom_triplet_2)
       
        
        cdef Vec3 normal_1 = self._calculate_normal(atom_triplet_1)
        cdef Vec3 normal_2 = self._calculate_normal(atom_triplet_2)
        
        cdef Vec3_container result = Vec3_container()
        cdef Vec3 average  = self._average_2_vec_3(normal_1, normal_2)
        result.set_vec3(average)
        return result
       
    cdef inline  _build_atom_triplet(self,atom_ids, int[3]& result):
        for i in range(3):
            result[i] = atom_ids[i]
        
        
    cdef inline Vec3 _average_2_vec_3(self, Vec3& vector_1, Vec3& vector_2):
        cdef Vec3 sum
        
        sum =  vector_1 + vector_2
            
        return Vec3(sum.x()/2.0, sum.y()/2.0, sum.z()/2.0)

    
            
    def __call__(self, rings, normals, centres):
        
        cdef Vec3_container centre 
        cdef Vec3_container normal
            
        for ring_component in rings:
            ring_id = ring_component[0]

            centre = self._calculate_one_ring_centre(ring_component)
            centre_component = ring_id, centre
            centres.add_component(centre_component)
            
            normal = self._calculate_one_ring_normal(ring_component)
            normal_component = ring_id, normal
            normals.add_component(normal_component)
            

cdef class Fast_non_bonded_calculator:
    cdef int _min_residue_seperation
    cdef float _cutoff_distance
    cdef float _jitter
    cdef float _full_cutoff_distance
    cdef bint _verbose 
    def __init__(self,min_residue_seperation,cutoff_distance=5.0,jitter=0.2):
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
    
    def __call__(self, atom_list_1, atom_list_2):
        non_bonded_lists = []
        for atom_id_1 in atom_list_1:
            non_bonded_list = []
            non_bonded_lists.append(non_bonded_list)
            for i, atom_id_2 in enumerate(atom_list_2):
                if self._is_non_bonded(atom_id_1, atom_id_2):
                    non_bonded_list.append(i)
        return  non_bonded_lists


cdef class Fast_energy_calculator:
    cdef object _energy_term_cache 
    cdef object _theory_shifts
    cdef object _observed_shifts
    cdef bint _verbose 

    def __init__(self):
        self._energy_term_cache =  None
        self._theory_shifts =   None
        self._observed_shifts =  None
        self._verbose = False
    
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
        
    def __call__(self,target_atom_ids):
        
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
        
        energy = 0.0
        
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
        return energy

cdef class Fast_force_factor_calculator(Fast_energy_calculator):

    def __call__(self, target_atom_ids):
        result  = []
        for target_atom_id in target_atom_ids:
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

            result.append(factor)
        return  result
cdef class Base_force_calculator:
    
    cdef object _components
    cdef bint _verbose 
    
    def __init__(self,potential=None):
        self._components =  None
        self._verbose = False
    
    def set_verbose(self,on):
        self._verbose = on
    
    def _set_components(self,components):
        self._components = components
    
    cdef inline object _get_component(self,int index):
        return self._components.get_component(index)
    
#    TODO should most probably be a fixed array
    def __call__(self, components, target_atom_ids, force_factors, forces):
        self._set_components(components)
        component_target_atom_ids = components.get_component_atom_ids()
        for i,target_atom_id in enumerate(target_atom_ids):
            if target_atom_id in component_target_atom_ids:
                index_range = components.get_component_range(target_atom_id)
                for index in range(*index_range):
                    self._calc_single_force_set(index,force_factors[i],forces)
    
    cdef inline _get_or_make_target_force_triplet(self, object forces, object target_offset):
        target_forces = forces[target_offset]
        if target_forces == None:
            target_forces = [0.0] * 3
            forces[target_offset] = target_forces
        return target_forces

#    TODO make abstract
    def _cython_calc_single_force_set(self,int index, float force_factor, object forces):
        raise Exception("unexpected! this method should be implemented")



cdef class Fast_distance_based_potential_force_calculator(Base_force_calculator):
    


    cdef int  _target_atom_index
    cdef int  _distance_atom_index_1
    cdef int  _distance_atom_index_2
    cdef int  _exponent_index
    cdef int  _coefficient_index
    cdef bint  _smoothed 
    cdef float _smoothing_factor
    cdef float _cutoff
    
    def __init__(self, object indices, bint smoothed):
        super(Fast_distance_based_potential_force_calculator, self).__init__()
        self._target_atom_index = indices.target_atom_index
        self._distance_atom_index_1 =  indices.distance_atom_index_1
        self._distance_atom_index_2 =  indices .distance_atom_index_2
        self._exponent_index  = indices.exponent_index
        self._coefficient_index  = indices.coefficient_index
        self._components =  None
        self._smoothed =  smoothed
        self._smoothing_factor =  DEFAULT_SMOOTHING_FACTOR
        self._cutoff =  DEFAULT_CUTOFF
        
    def set_cutoff(self, cutoff):
        self._cutoff =  cutoff
    
    def set_smoothing_factor(self,smoothing_factor):
        self._smoothing_factor = smoothing_factor
        
    def _set_components(self,components):
        self._components = components
    
    cdef inline target_distant_atom _get_target_and_distant_atom_ids(self, int index):
        
        cdef object values
        cdef int target_atom_id
        cdef int distant_atom_id
        cdef target_distant_atom result
        
        values  = self._components.get_component(index)
        
        result.target_atom_id = values[self._distance_atom_index_1]
        result.distant_atom_id = values[self._distance_atom_index_2]
        
        return result
    
    cdef inline coefficient_exponent _get_coefficient_and_exponent(self, int index):
        cdef object values
        cdef float coefficient
        cdef float exponent
        cdef coefficient_exponent result
        values = self._components.get_component(index)
        


        result.coefficient = values[self._coefficient_index]
        result.exponent = values[self._exponent_index]
        
        return result
 
#    cdef inline float _distance(self, int target_atom, int distance_atom):
#        cdef Vec3 target_pos, distant_pos, distance
#        cdef float result 
#        
#        target_pos = currentSimulation().atomPosArr().data(target_atom)
#        distant_pos =  currentSimulation().atomPosArr().data(distance_atom)
#        
#        return  norm(target_pos - distant_pos)
        
    cdef inline Vec3 _xyz_distances(self, int target_atom, int distance_atom):
        cdef Vec3 target_pos, distant_pos
        
        
        target_pos = currentSimulation().atomPosArr().data(target_atom)
        distant_pos =  currentSimulation().atomPosArr().data(distance_atom)
        
        return target_pos - distant_pos
        
    cdef inline float _sum_xyz_distances_2(self, int target_atom, int distance_atom):
        cdef Vec3 target_pos, distant_pos, distance
        cdef float result =0.0
        
        target_pos = currentSimulation().atomPosArr().data(target_atom)
        distant_pos =  currentSimulation().atomPosArr().data(distance_atom)
        
        distance = Vec3(target_pos - distant_pos)
        
        for i in range(3):
            result += distance[i] * distance[i]
            
        return result 
    
    
    cpdef float _calc_single_force_factor(self, int index, float factor):
        cdef target_distant_atom atom_ids
        cdef coefficient_exponent coef_exp
        cdef float exponent 
        
        cdef float full_factor
        cdef float ratio, pre_exponent, reduced_exponent, sum_xyz_distances_2
        
        atom_ids = self._get_target_and_distant_atom_ids(index)
        
        coef_exp = self._get_coefficient_and_exponent(index)
        
        sum_xyz_distances_2 = self._sum_xyz_distances_2(atom_ids.target_atom_id, atom_ids.distant_atom_id)
#
        full_factor= factor * coef_exp.coefficient
        
        exponent = coef_exp.exponent
        if self._smoothed:
            ratio = sum_xyz_distances_2 / (self._cutoff**2)
            ratio =  ratio**4
            pre_exponent = exponent - (exponent + 8.0) * ratio
        else:
            pre_exponent = exponent
            
        reduced_exponent = (exponent - 2.0) / 2.0
        
        force_factor = full_factor *  pre_exponent * sum_xyz_distances_2 ** reduced_exponent

        return force_factor


 

    
    def _calc_single_force_set(self, int index, float factor, object forces):
        self._cython_calc_single_force_set(index, factor, forces)
        
    cdef _cython_calc_single_force_set(self, int index, float factor, object forces):
        
        cdef target_distant_atom atom_ids
        cdef Vec3 xyz_distances, target_forces, distant_forces
        cdef float force_factor, distance

        atom_ids =  self._get_target_and_distant_atom_ids(index)
        
#        print atom_ids.target_atom_id,atom_ids.distant_atom_id
        xyz_distances  = self._xyz_distances(atom_ids.target_atom_id,atom_ids.distant_atom_id)
        
        force_factor  = self._calc_single_force_factor(index, factor)
        
        
        python_target_forces = self._get_or_make_target_force_triplet(forces, atom_ids.target_atom_id)
        python_distant_forces  = self._get_or_make_target_force_triplet(forces, atom_ids.distant_atom_id)
        
        target_forces = operator_times(xyz_distances, -force_factor)
        distant_forces = operator_times(xyz_distances, force_factor)
        
        for i in range(3):
            python_target_forces[i] += target_forces[i]
            python_distant_forces [i] += distant_forces[i]

cdef class Fast_non_bonded_force_calculator(Fast_distance_based_potential_force_calculator):
    cdef float _non_bonded_cutoff
    
    def __init__(self, object indices, bint smoothed):
        global DEFAULT_NB_CUTOFF
        super(Fast_non_bonded_force_calculator, self).__init__(indices,smoothed)
        self._non_bonded_cutoff = DEFAULT_NB_CUTOFF
        
    def _calc_single_force_set(self, int index, float factor, object forces):
        self._cython_calc_single_force_set(index,  factor, forces)
        
    cdef _cython_calc_single_force_set(self, int index, float factor, object forces):
        cdef target_distant_atom atom_ids = Fast_distance_based_potential_force_calculator._get_target_and_distant_atom_ids(self,index)
        cdef float distance  = calc_distance(atom_ids.target_atom_id, atom_ids.distant_atom_id)
#        TODO: this should be the non bonded distance cutoff
#TODO class variable of self are not being looked up!
        if distance < 5.0:
            Fast_distance_based_potential_force_calculator._cython_calc_single_force_set(self,index, factor, forces)


    
cdef class Fast_dihedral_force_calculator(Base_force_calculator):
    
    def __init__(self):
        Base_force_calculator.__init__(self)

    cdef inline dihedral_ids _get_dihedral_atom_ids(self, int index):
        cdef dihedral_ids result
        
        component = self._get_component(index)
        
        result.atom_id_1= component[1]
        result.atom_id_2= component[2]
        result.atom_id_3= component[3]
        result.atom_id_4= component[4]
        
        return result
    
    cdef inline dihedral_parameters _get_parameters(self, int index):
        component = self._get_component(index)
        
        cdef dihedral_parameters result  
        
        parameters = component[6:11]
        
        result.param_0 =  parameters[0]
        result.param_1 =  parameters[3]
        result.param_2 =  parameters[1]
        result.param_3 =  parameters[4]
        result.param_4 =  parameters[2]
        
            
        return result
    
    cdef inline float _get_coefficient(self, int index):
        component = self._get_component(index)
        coefficient = component[5]
        return coefficient
    
#    TODO make this consistent with the distance forces factor
    cpdef _calc_single_force_factor(self, int index):
        
        cdef dihedral_ids dihedral_atom_ids
        cdef dihedral_parameters params
        
        
        dihedral_atom_ids  = dihedral_atom_ids= self._get_dihedral_atom_ids(index)
        
        params = self._get_parameters(index)
        
        angle = calc_dihedral_angle(dihedral_atom_ids)
        
        result = -3.0 * params.param_0 * sin(3.0 * angle + params.param_1) - \
                        params.param_2 * sin(angle + params.param_3)
        return result

    
    #TODO: is this too close?
    def _calc_single_force_set(self, int index, float factor, object forces):
        self._cython_calc_single_force_set(index, factor, forces)
        
    cdef _cython_calc_single_force_set(self, int index, float factor, object forces):
        cdef Vec3 r1, r2, r3, r4, temp
        cdef Vec3 n1, n2
        cdef float weight
        cdef float factor_1, factor_2, factor_3
        cdef Vec3 F1, F2, F3, F4
        cdef Vec3 T1, T2, T3, T4
        
        cdef float dihedral_factor = self._calc_single_force_factor(index)
        cdef dihedral_ids atom_ids = self._get_dihedral_atom_ids(index)
        
        v1 = currentSimulation().atomPosArr().data(atom_ids.atom_id_1)
        v2 = currentSimulation().atomPosArr().data(atom_ids.atom_id_2)
        v3 = currentSimulation().atomPosArr().data(atom_ids.atom_id_3)
        v4 = currentSimulation().atomPosArr().data(atom_ids.atom_id_4)
         
        r1 = v1 - v2
        r2 = v3 - v2
        temp = Vec3(r2)
        r3 = operator_times(temp, -1)
        r4 = v4 - v3

#        print vec3_as_tuple(r1), vec3_as_tuple(r2), vec3_as_tuple(r3), vec3_as_tuple(r4)
        
        # compute normal vector to plane containing v1, v2, and v3
        n1 = cross(r1, r2)
        # compute normal vector to plane containing v2, v3, and v4
        n2 = cross(r3, r4)
        
        r2_length = calc_distance(atom_ids.atom_id_2,atom_ids.atom_id_3)
        r2_length_2 = r2_length*r2_length
        

        weight = factor * self._get_coefficient(index)
        
        

#        // force calculation according to Bekker, Berendsen and van Gunsteren (1995),
#        // Journal of Computational Chemistry 16, pp. 527-533:
##        // Force and virial of torsional-angle-dependent potentials.
        factor_1 = dihedral_factor * r2_length
        temp=Vec3(n1)
        F1 = operator_times(temp, (-factor_1 / norm(n1)**2))
        temp = Vec3(n2)
        F4 = operator_times(temp, ( factor_1 / norm(n2)**2))
        
        factor_2 = dot(r1, r2) / r2_length_2
        factor_3 = dot(r3, r4) / r2_length_2

        temp = Vec3(F1)
        T1 = operator_times(temp,  (factor_2-1))
        temp = Vec3(F4)
        T2 = operator_times(temp,  -factor_3)
        
        F2 = T1 + T2
        
        temp = Vec3(F1)
        T1 = operator_times(temp,  -factor_2)
        temp =Vec3(F4)
        T2 = operator_times(temp, (factor_3 -1))

        F3 = T1 + T2

#        print vec3_as_tuple(F1), vec3_as_tuple(F2), vec3_as_tuple(F3), vec3_as_tuple(F4)
#        // assign forces
        force_triplet_1 = self._get_or_make_target_force_triplet(forces, atom_ids.atom_id_1)
        force_triplet_2 = self._get_or_make_target_force_triplet(forces, atom_ids.atom_id_2)
        force_triplet_3 = self._get_or_make_target_force_triplet(forces, atom_ids.atom_id_3)
        force_triplet_4 = self._get_or_make_target_force_triplet(forces, atom_ids.atom_id_4)
        
        self.set_force_triplet(weight,force_triplet_1,F1)
        self.set_force_triplet(weight,force_triplet_2,F2)
        self.set_force_triplet(weight,force_triplet_3,F3)
        self.set_force_triplet(weight,force_triplet_4,F4)
        
        
    
    cdef set_force_triplet(self, float weight, object force_triplet, Vec3 vector):
        for offset in range(3):
            force_component = weight * vector[offset]
            force_triplet[offset]+= force_component
            


    
    
#cdef int RING_ATOM_IDS = 1
cdef class Fast_ring_force_calculator(Base_force_calculator):

    
    cdef object _coef_components
    cdef object _ring_components
    cdef object _centre_cache
    cdef object _normal_cache
    
    
    def __init__(self):
        super(Fast_ring_force_calculator, self).__init__()
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

#    cdef _get_coef_components(self, int atom_type_id):
#        return self._coef_components.get_components_for_atom_id(atom_type_id)
            
    cdef Vec3 _get_ring_normal(self, int ring_id):
        cdef Vec3_container container =  self._normal_cache.get_component(ring_id)[1]
        return  container.get_vec3()
    
    cdef Vec3 _get_ring_centre(self, int ring_id):
        cdef Vec3_container container = self._centre_cache.get_component(ring_id)[1]
        return  container.get_vec3()       
        

    

    #TODO: use this more places
    cdef inline ring_atom_ids _get_ring_atom_ids(self, int ring_id):
        cdef object ring_component = self._ring_components.get_component(ring_id)
        cdef object python_ring_atom_ids = ring_component[RING_ATOM_IDS]
        cdef int num_atoms =  len(python_ring_atom_ids)
        cdef ring_atom_ids  ring_atoms 
        
        ring_atoms.num_atoms  = num_atoms
        for i in range(num_atoms):
            ring_atoms.atom_ids[i] = python_ring_atom_ids[i]
        return  ring_atoms
       
    

    
    cdef inline Ring_atom_positions _get_ring_atom_positions(self, int ring_id):
        cdef ring_atom_ids ring_atoms =  self._get_ring_atom_ids(ring_id)
        cdef int num_atoms = ring_atoms.num_atoms
        cdef Ring_atom_positions  positions  = Ring_atom_positions(num_atoms) 
        
        cdef Vec3* position = NULL
        for i in range(num_atoms):
            position = positions.get_vec(i)
            position[0] = currentSimulation().atomPosArr().data(ring_atoms.atom_ids[i])
            
        return positions

    cdef inline target_type _get_target_and_type(self,int index):
        component = self._get_component(index)
#        TODO: this needs to contain better names a union?
        cdef target_type result
        result.target_atom_id = component[0]
        result.atom_type_id = component[1]
        
        return result#
    
    cdef inline  Coef_components _get_coef_components(self, int atom_type_id):
        cdef object python_coef_components = self._coef_components.get_components_for_atom_id(atom_type_id)
        return Coef_components(python_coef_components)
        
    cpdef _calc_single_force_set(self, int index, float force_factor, object forces):
        cdef target_type target_atom_id_type =  self._get_target_and_type(index)
        #TODO: make array a union
        cdef Coef_components coef_components = self._get_coef_components(target_atom_id_type.atom_type_id)
        cdef Coef_component* coef_component
        for i in range(coef_components.num_components):
            coef_component = coef_components.get_component(i)
#            print coef_component.atom_type_id, coef_component.ring_id, coef_component.coefficient
    #                print 'coef_component', i, coef_component
    #                 print 'here', target_atom_id, atom_type_id,
    #                coef_component, force_factor
            self._calculate_single_ring_forces(target_atom_id_type.target_atom_id, target_atom_id_type.atom_type_id, coef_component, force_factor, forces)

    def _calculate_ring_forces(self, int atom_type_id, int ring_id, float force_factor, Python_ring_force_sub_terms python_sub_terms, object forces):
        cdef Ring_force_sub_terms  terms  =  python_sub_terms.get_terms()
        self._cython_calculate_ring_forces(atom_type_id, ring_id, force_factor, terms, forces)

    cdef _calculate_single_ring_forces(self, int target_atom_id, int atom_type_id, Coef_component* coef_component, float force_factor, object forces):
        
#    print atom_type_id, ring_id, coefficient, self._get_component_list('RING').get_components_for_atom_id(ring_id)
        cdef Ring_force_sub_terms force_terms = self._build_cython_force_terms(target_atom_id, coef_component.ring_id)
#        print coef_component
        
#        print vec3_as_tuple(force_terms.gradUQ), force_terms.dL3, force_terms.u, vec3_as_tuple(force_terms.gradVQ), force_terms.dL6
#        print Atom_utils._get_atom_info_from_index(target_atom_id)

        self._cython_calc_target_atom_forces(target_atom_id, force_factor * coef_component.coefficient, force_terms, forces)
#        #            #TODO: this is not how camshift does it, it uses the sum of the two ring normals
        self._cython_calculate_ring_forces(atom_type_id, coef_component.ring_id, force_factor * coef_component.coefficient, force_terms, forces)

    def _build_force_terms(self, int target_atom_id, int ring_id):
        cdef Ring_force_sub_terms terms  = self._build_cython_force_terms(target_atom_id, ring_id)
        result  = Python_ring_force_sub_terms()
        result.setup(terms)
        return result
        
    cdef Ring_force_sub_terms _build_cython_force_terms(self, int target_atom_id, int ring_id):
        cdef Ring_force_sub_terms result
        cdef Vec3 target_atom_pos  = currentSimulation().atomPosArr().data(target_atom_id)
        
        cdef Vec3 ring_centre = self._get_ring_centre(ring_id)
        cdef Vec3 ring_normal = self._get_ring_normal(ring_id)
        
        # distance vector between atom of interest and ring center
        cdef Vec3 d = target_atom_pos - ring_centre
        cdef float dL = norm(d)
        
        # squared distance of atom of interest from ring center
        cdef float dL2 = dot(d, d) #            if (dL2 < 0.5) cout << "CAMSHIFT WARNING: Distance between atom and center of ring alarmingly small at " << sqrt(dL2) << " Angstrom!" << endl;
        
        # calculate terms resulting from differentiating energy function with respect to query and ring atom coordinates
        cdef float dL4 = dL ** 4
        cdef float dL3 = dL ** 3
        cdef float dL6 = dL3 ** 2
        
        cdef float nL = norm(ring_normal)
        cdef float nL2 = nL ** 2
        cdef float dLnL = dL * nL
        cdef dL3nL3 = dL3 * nL2 * nL
        
        cdef float dn = dot(d, ring_normal)
        cdef float dn2 = dn ** 2
        
        cdef float u = 1.0 - 3.0 * dn2 / (dL2 * nL2)
        
        cdef float gradUQ_factor = -6.0 * dn / (dL4 * nL2)
        
        #TODO: remove temporarys and operator_mu;t#
        cdef Vec3 temp_d =  Vec3(d)
        cdef Vec3 scaled_d = operator_times(temp_d,dn)
        cdef Vec3 temp_normal = Vec3(ring_normal)
        cdef Vec3 scaled_normal = operator_times(temp_normal,dL2)
        cdef Vec3 temp_normal_distance = Vec3(scaled_normal - scaled_d)
        gradUQ = operator_times(temp_normal_distance,gradUQ_factor)
            
        
        cdef float gradVQ_factor = 3.0 * dL
        cdef Vec3 gradVQ= operator_times(d,gradVQ_factor)
        
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
        result.ring_normal =  ring_normal
        
        return result
#        return  dL3, u, dL6, ring_normal,  atom_type_id, coefficient, d, factor, dn, dL3nL3, dL, nL, dLnL
##    ---
#
    def _calc_target_atom_forces(self, int target_atom_id, float force_factor, Python_ring_force_sub_terms python_sub_terms, object forces):
        cdef Ring_force_sub_terms  terms  =  python_sub_terms.get_terms()
        self._cython_calc_target_atom_forces(target_atom_id, force_factor, terms, forces)

    cdef _cython_calc_target_atom_forces(self, int target_atom_id, float force_factor, Ring_force_sub_terms& sub_terms, object forces):
        cdef object target_force_triplet
        cdef int axis 
        target_force_triplet = self._get_or_make_target_force_triplet(forces, target_atom_id)
#    
#        # update forces on query atom
        cdef float term1
        for axis in range(3):
            #TODO: report this as a bug it produces an un initialized reference 
            term1= (sub_terms.gradUQ[axis] * sub_terms.dL3 - sub_terms.u * sub_terms.gradVQ[axis])
            target_force_triplet[axis] += -force_factor * term1 / sub_terms.dL6
#             *  / sub_terms.dL6
#
#    #TODO: calculation of GradU and gradV are not consistent with force_terms for target atom correct
#    #TODO: reduce number of parameters to method

    cdef _cython_calculate_ring_forces(self, int atom_type_id, int ring_id, float force_factor, Ring_force_sub_terms force_terms, object forces):
        cdef Vec3 temp_normal =  Vec3(force_terms.ring_normal)
        cdef Vec3 nSum = operator_times(temp_normal,2.0)  #            float_type g [3], ab [3], c [3]
        cdef ring_atom_ids ring_atoms = self._get_ring_atom_ids(ring_id)
        cdef Ring_atom_positions ring_atom_positions = self._get_ring_atom_positions(ring_id)
    #// 2 for a 5-membered ring, 3 for a 6-membered ring
        cdef int num_ring_atoms = ring_atoms.num_atoms
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
        
        cdef int axis 
        cdef float sub_term
        for ring_atom_index in range(num_ring_atoms):
            ring_atom_id = ring_atoms.atom_ids[ring_atom_index]
            if ring_atom_index < limit:
                for axis in range(3):
                    index_1 = (ring_atom_index + 1) % 3
                    index_2 = (ring_atom_index + 2) % 3
                    pos_1 = (ring_atom_positions.get_vec(index_1)[0][axis])
                    pos_2 = ring_atom_positions.get_vec(index_2)[0][axis]
                    g[axis] =  pos_1 - pos_2 # atoms 3,4 (5 member) or 3,4,5 (6 member)
            
            else:
                if  ring_atom_index >= num_ring_atoms - limit:
                    offset = num_ring_atoms - 3 #2 for a 5-membered ring, 3 for a 6-membered ring
                    for axis in range(3):
                        index_1 = (ring_atom_index + 1 - offset) % 3 + offset
                        index_2 = (ring_atom_index + 2 - offset) % 3 + offset
                        pos_1 = ring_atom_positions.get_vec(index_1)[0][axis]
                        pos_2 = ring_atom_positions.get_vec(index_2)[0][axis]
                        g[axis] = pos_1 - pos_2
                else:
                    
                    for axis in range(3):
                        pos_1 = ring_atom_positions.get_vec(0)[0][axis]
                        pos_2 = ring_atom_positions.get_vec(1)[0][axis] 
                        pos_3 = ring_atom_positions.get_vec(3)[0][axis]
                        pos_4 = ring_atom_positions.get_vec(4)[0][axis]
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

            ring_target_force_triplet = self._get_or_make_target_force_triplet(forces, ring_atom_id)
            for axis in range(3):
                sub_term = (gradU[axis] * force_terms.dL3 - force_terms.u * gradV[axis])
                sub_force = -force_factor * sub_term / force_terms.dL6
                ring_target_force_triplet[axis] += sub_force
#                print AXIS_NAMES[axis],sub_force,-force_factor, gradU[axis], force_terms.dL3, force_terms.u, gradV[axis],force_terms.dL6
#            print
##        print 
