'''
Created on 31 Jul 2012

@author: garyt
'''
from vec3 import Vec3 as python_vec3
from  xplor_access cimport norm,Vec3,currentSimulation, Dihedral, Atom,  dot,  cross
from libc.math cimport cos,fabs, tanh, pow


cdef struct target_distant_atom:
    int target_atom_id
    int distant_atom_id

cdef struct coefficient_exponent:
    float coefficient
    float exponent
    
cdef float sum(Vec3 vec3):
    return vec3.x()+vec3.y()+vec3.z()

cdef inline  Vec3 operator_times (Vec3& vec3, float scale):
    return Vec3(vec3.x() * scale, vec3.y() * scale, vec3.z() * scale)
    
cdef float DEFAULT_CUTOFF = 5.0
cdef float DEFAULT_SMOOTHING_FACTOR = 1.0
    
cdef class Fast_distance_shift_calculator:

    
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
    def __call__(self, object components, object results):
        self._set_components(components)
        cdef float smoothing_factor = self._smoothing_factor
        cdef float ratio
        cdef float result
        cdef target_distant_atom atom_indices
        cdef coefficient_exponent coef_exp
        for index in range(len(components)):
            atom_indices = self._get_target_and_distant_atom_ids(index)
            
            coef_exp = self._get_coefficient_and_exponent(index)
            distance =self.distance(atom_indices.target_atom_id, atom_indices.distant_atom_id)
    #        Atom_utils._calculate_distance(target_atom_index, sidechain_atom_index)

            if self._smoothed:
                ratio = distance / self._cutoff
                smoothing_factor = 1.0 - ratio ** 8
            results[index]  = smoothing_factor * pow(distance,  coef_exp.exponent) * coef_exp.coefficient
            
    
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
            

cdef class Fast_non_bonded_calculator:
    
    def __init__(self,min_residue_seperation,cutoff_distance=5.0,jitter=0.2):
        self._min_residue_seperation =  min_residue_seperation
        self._cutoff_distance =  cutoff_distance
        self._jitter = jitter
        self._full_cutoff_distance =  self._cutoff_distance + self._jitter
        
    cdef inline bint _filter_by_residue(self, char* seg_1, int residue_1, char* seg_2, int residue_2):
        cdef int sequence_distance
        cdef bint result
        
        result = False
        if seg_1 == seg_2:
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
    def __init__(self):
        self._energy_term_cache =  None
        self._theory_shifts =   None
        self._observed_shifts =  None
        
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

cdef class Base_force_calculator:
    
    cdef object _components
    
    def __init__(self,potential=None):
        self._components =  None
    
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
    cdef _calc_single_force_set(self,int index, float force_factor, object forces):
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
    
    cdef target_distant_atom _get_target_and_distant_atom_ids(self, int index):
        
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
 
    cdef inline float _distance(self, int target_atom, int distance_atom):
        cdef Vec3 target_pos, distant_pos, distance
        cdef float result 
        
        target_pos = currentSimulation().atomPosArr().data(target_atom)
        distant_pos =  currentSimulation().atomPosArr().data(distance_atom)
        
        return  norm(target_pos - distant_pos)
        
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
    
    
    cdef float _calc_single_force_factor(self, int index, float factor):
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


 
    cdef vec3_as_tuple(self,Vec3& vec_3):
        return vec_3.x(), vec_3.y(), vec_3.z()
    
    cdef object _calc_single_force_set(self, int index, float factor, object forces):
        
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