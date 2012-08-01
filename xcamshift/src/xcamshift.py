'''
Created on 27 Dec 2011

@author: garyt
#TODO: need to translate from_atom names
'''

#TODO: add tests to exclude atoms/distances which are not defined


from atomSel import AtomSel,intersection
from component_list import Component_list
from dihedral import Dihedral
from keys import Atom_key, Dihedral_key
from math import cos, tanh, cosh, sin
from observed_chemical_shifts import Observed_shift_table
from segment_manager import Segment_Manager
from table_manager import Table_manager
from python_utils import tupleit
from utils import Atom_utils, AXES, iter_residue_atoms,\
    iter_residues_and_segments, iter_residue_types
from vec3 import Vec3, norm, cross, dot
import abc
import sys
from common_constants import  BACK_BONE, XTRA, RANDOM_COIL, DIHEDRAL, SIDE_CHAIN, RING, NON_BONDED
import itertools
from abc import abstractmethod, ABCMeta
from cython.shift_calculators import Fast_shift_calculator

class Component_factory(object):
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def is_residue_acceptable(self, segment, residue_number, segment_manager):
        pass
    
    @abc.abstractmethod
    def create_components(self, component_list, table_source, segment,target_residue_number,selected_atoms):
        pass
        
    @abc.abstractmethod
    def get_table_name(self):
        pass
    
    def _get_table_for_residue(self, segment, target_residue_number, table_provider):
        from_residue_type = Atom_utils._get_residue_type(segment, target_residue_number)
        table = table_provider(from_residue_type)
        return table

class Residue_component_factory(Component_factory):
    
    def is_residue_acceptable(self, segment, residue_number, segment_manager):
        return True
    
    def create_components(self, component_list, table_source, segment,target_residue_number,selected_atoms):
        table = self._get_table_for_residue(segment, target_residue_number, table_source)
        self.create_residue_components(component_list, table, segment, target_residue_number)
    
    @abc.abstractmethod
    def create_residue_components(self,component_list,table, segment, residue):
        pass
    
class Atom_component_factory(Component_factory):
    __metaclass__ = abc.ABCMeta
    
    def is_residue_acceptable(self, segment, residue_number, segment_manager):
        segment_info = segment_manager.get_segment_info(segment)
        
        return residue_number > segment_info.first_residue and residue_number < segment_info.last_residue




    def create_components(self, component_list, table_provider, segment,target_residue_number,selected_atoms):
        table = self._get_table_for_residue(segment, target_residue_number, table_provider)
        self.create_atom_components(component_list, table, selected_atoms)
    
    @abc.abstractmethod
    def _build_contexts(self,atom, table):
        pass
    
    @abc.abstractmethod
    def _get_component_for_atom(self,atom, context):
        pass
    
    @abc.abstractmethod
    def  _translate_atom_name(self, atom_name,context):
        return atom_name
    
    def create_atom_components(self, component_list, table, selected_atoms):
        for atom in selected_atoms:
            contexts = self._build_contexts(atom, table)
            for context in contexts:
                if context.complete:
                    value = self._get_component_for_atom(atom, context)
                    if value != None:
                        component_list.add_component(value)

    def get_table_name(self):
        return 'ATOM'

class DihedralContext(object):
    
        def _translate_atom_name_from_table(self, residue_type, atom_name,table):
            return  table.get_translation_from_table(residue_type,atom_name)

        def _select_atom_with_translation(self, segment, residue_number_1, atom_name_1):
            
            residue_type  = Atom_utils._get_residue_type(segment, residue_number_1)
            atom_name_1 = self._translate_atom_name_from_table(residue_type, atom_name_1, self._table)
            target_atom_1 = Atom_utils.find_atom(segment, residue_number_1, atom_name_1)
        
            num_to_atom = len(target_atom_1)
            if num_to_atom > 1:
                self._get_atom_names(target_atom_1)
                raise Exception("unexpected number of to atoms selected (> 1) %d" % num_to_atom)
            return target_atom_1
    
    
        def get_atom_info(self, segment, from_atom, key_1):
            residue_number_1 = from_atom.residueNum() + key_1.offset
            atom_name_1 = key_1.atom
            target_atom_1 = self._select_atom_with_translation(segment, residue_number_1, atom_name_1)
            return target_atom_1
    
        def __init__(self, from_atom, dihedral_key ,table):
            
            self._table = table
            
            segment = from_atom.segmentName()
            
            target_atoms = []
            for key in dihedral_key.get_keys():
                target_atoms.append(self.get_atom_info(segment, from_atom, key))
            
            self.complete=True
            for target_atom in target_atoms:
                if len(target_atom) !=  1:
                    self.complete =  False
                    
            if self.complete:
                self.dihedral_indices = [target_atom[0].index()  for target_atom in target_atoms]
                self.dihedral_key = dihedral_key


class Dihedral_component_factory(Atom_component_factory):

    def get_table_name(self):
        return 'ATOM'
    
    def _translate_atom_name_to_table(self, residue_type, atom_name,table):
        
        return  table.get_translation_to_table(residue_type,atom_name)

#TODO: make this more generic and move to super class
    def _build_contexts(self, atom, table):
        contexts = []
        
        for key in table.get_dihedral_keys():
            dihedral_key = Dihedral_key(*key)
            context = DihedralContext(atom,dihedral_key,table)
            
            if context.complete:
                contexts.append(context)
                
        return contexts

    def _get_component_for_atom(self, atom, context):
        table = context._table
        
        
        from_atom_name = atom.atomName()
        from_residue_type = atom.residueName()
        from_atom_name = self._translate_atom_name_to_table(from_residue_type,from_atom_name,table)
        
        result = None
        if from_atom_name in table.get_target_atoms():
            
            dihedral_key = context.dihedral_key
            value = context._table.get_dihedral_shift(from_atom_name,dihedral_key)
            if value != None:
                from_atom_index = atom.index()
                
                dihedral_indices = context.dihedral_indices

                result = [from_atom_index]
                result.extend(dihedral_indices)
                result.append(value)
                
                for parameter_id in context._table.get_parameters():
                    parameter = context._table.get_parameter(from_atom_name,dihedral_key,parameter_id)
                    result.append(parameter)
                
                
                result.append(context._table.get_exponent())
                result = tuple(result)
                
        return result
    
    def _translate_atom_name(self, atom_name,context):
        return context._table.get_translation(atom_name)

class DistanceContext:

    def _translate_atom_name_from_table(self, residue_type, atom_name,table):
        return  table.get_translation_from_table(residue_type,atom_name)
    
    
    def __init__(self, from_atom, offset, to_atom_name ,table):
        
        self.complete = False
        self._table = table
        self.to_atom_name =  to_atom_name
        self.segment = from_atom.segmentName()
        
        self.offset = offset
        
        self.residue_number = from_atom.residueNum()
        self.residue_type = Atom_utils._get_residue_type(self.segment,self. residue_number)
        
        self.to_residue_number = self.residue_number+offset
        to_residue_type = Atom_utils._get_residue_type(self.segment, self.to_residue_number)
        
        to_atom_name = self._translate_atom_name_from_table(to_residue_type, to_atom_name, table)
        
        #TODO remove select atom with translation
        to_atom = Atom_utils._select_atom_with_translation(self.segment, self.to_residue_number, to_atom_name)
        
        if len(to_atom) == 0:
            to_atom = Atom_utils._select_atom_with_translation(self.segment, self.to_residue_number, to_atom_name)
            
        num_to_atom = len(to_atom)
        if num_to_atom > 1:
            self._get_atom_names(to_atom)
            raise Exception("unexpected number of to atoms selected (> 1) %d" % num_to_atom)
            
        if len(to_atom) == 1:
            self.to_residue_type = Atom_utils._get_residue_type(self.segment, self.to_residue_number)
            self.to_atom_index = to_atom[0].index()
            self.complete = True

class Distance_component_factory(Atom_component_factory):
    def __init__(self):
        pass
    
    #TODO: remove just kept to satisfy an abstract method declaration
    def _translate_atom_name(self, atom_name, context):
        pass

    def _translate_atom_name_to_table(self, residue_type, atom_name,table):
        
        return  table.get_translation_to_table(residue_type,atom_name)
    
    def _build_contexts(self, atom, table):
        contexts = []
        for offset in table.get_offsets():
            for to_atom_name in table.get_to_atoms():
                context = DistanceContext(atom,offset,to_atom_name,table)
                if context.complete:
                    contexts.append(context)
        return contexts
    
    def _get_component_for_atom(self, atom, context):
        table = context._table
        
        from_atom_name = atom.atomName()
        from_residue_type = atom.residueName()
        #TODO move translation into context
        from_atom_name = self._translate_atom_name_to_table(from_residue_type,from_atom_name,table)
        
        offset = context.offset
        to_atom_name = context.to_atom_name
        to_residue_type = context.to_residue_type
        #TODO move translation into context
        to_atom_name = self._translate_atom_name_to_table(to_residue_type,to_atom_name,table)
        
        
        result = None
        if from_atom_name in table.get_from_atoms():
            if to_atom_name in table.get_to_atoms():
                value = context._table.get_distance_coeeficent(from_atom_name,offset,to_atom_name)
                if value != None:
                    from_atom_index = atom.index()
                    to_atom_index = context.to_atom_index
                    if from_atom_index != to_atom_index:
                        exponent = context._table.get_exponent()
                        result = (from_atom_index,to_atom_index,value,exponent)
        return result

    
    def get_table_name(self):
        return 'ATOM'
    
#TODO: check if previous residue offsets are included
class Random_coil_context :
    def __init__(self,atom,offset,table):
        
        self._table = table
        
        segment = atom.segmentName()
        
        self.offset = offset
        
        from_residue_number = atom.residueNum()
        
        self.to_residue_number = from_residue_number + offset
        self.to_residue_type = Atom_utils._get_residue_type(segment, from_residue_number + offset)
        
        self.complete = True

class Random_coil_component_factory(Atom_component_factory):
    
    #TODO: push to base
    #TODO should random coil have a translation table
    def _translate_atom_name(self, residue_type, atom_name,table):
        
        return  table.get_translation(residue_type,atom_name)
        
        

    def _get_component_for_atom(self, atom, context):
        result = None
        
        atom_name = atom.atomName()
        residue_type = atom.residueName()
        table  =  context._table
        atom_name = self._translate_atom_name(residue_type,atom_name,table)
        
        if atom_name in context._table.get_atoms():
            value = context._table.get_random_coil_shift(context.offset, context.to_residue_type, atom_name)
            if value != None:
                result = atom.index(), value
        return result
    

    def _build_contexts(self, atom, table):
        contexts = []
        for offset in table.get_offsets():
            context = Random_coil_context(atom, offset, table)
            contexts.append(context)
        return contexts

    def get_table_name(self):
        return 'ATOM'
    
class ExtraContext(object):

    def _translate_atom_name_from_table(self, residue_type, atom_name,table):
        return  table.get_translation_from_table(residue_type,atom_name)
    
    def _select_atom_with_translation(self, segment, residue_number_1, atom_name_1):

        residue_type  = Atom_utils._get_residue_type(segment, residue_number_1)
        atom_name_1 = self._translate_atom_name_from_table(residue_type, atom_name_1, self._table)
        target_atom_1 = Atom_utils.find_atom(segment, residue_number_1, atom_name_1)

        num_to_atom = len(target_atom_1)
        if num_to_atom > 1:
            self._get_atom_names(target_atom_1)
            raise Exception("unexpected number of to atoms selected (> 1) %d" % num_to_atom)

        return target_atom_1

    def __init__(self, from_atom, key_1, key_2 ,table):
        
        self.complete = False
        self._table = table
        
        segment = from_atom.segmentName()
        
        residue_number_1 = from_atom.residueNum() + key_1.offset
        atom_name_1 = key_1.atom
        target_atom_1 = self._select_atom_with_translation(segment, residue_number_1, atom_name_1)
        
        atom_name_2 = key_2.atom
        residue_number_2 = from_atom.residueNum() + key_2.offset
        target_atom_2 =self._select_atom_with_translation(segment, residue_number_2, atom_name_2)
        
        if len(target_atom_1) == 1 and len(target_atom_2) == 1:
            self.distance_atom_index_1 = target_atom_1[0].index()
            self.distance_atom_index_2 = target_atom_2[0].index()
            self.key_1=key_1
            self.key_2=key_2
            self.complete = True


class Extra_component_factory(Atom_component_factory):

    def _translate_atom_name(self, atom_name,context):
        pass
    
    def _translate_atom_name_to_table(self, residue_type, atom_name,table):
        
        return  table.get_translation_to_table(residue_type,atom_name)
    
    def _build_contexts(self, atom, table):
        contexts = []
        for offset_1 in table.get_offsets(table.ATOM_1):
            for offset_2 in table.get_offsets(table.ATOM_2):
                for distance_atom_1 in table.get_distance_atoms(table.ATOM_1):
                    for distance_atom_2 in table.get_distance_atoms(table.ATOM_2):
                        key_1 = Atom_key(offset_1,distance_atom_1)
                        key_2 = Atom_key(offset_2,distance_atom_2)
                        context = ExtraContext(atom,key_1,key_2,table)
                        if context.complete:
                            contexts.append(context)
        return contexts
    
    def  _get_component_for_atom(self, atom, context):
        table = context._table

        from_atom_name = atom.atomName()
        from_residue_type = atom.residueName()
        #TODO move translation into context
        from_atom_name = self._translate_atom_name_to_table(from_residue_type,from_atom_name,table)
        
        result = None
        if from_atom_name in table.get_target_atoms():
            value = context._table.get_extra_shift(from_atom_name,context.key_1,context.key_2)

            if value != None:
                from_atom_index = atom.index()
                distance_index_1 = context.distance_atom_index_1
                distance_index_2 = context.distance_atom_index_2
                exponent = context._table.get_exponent()
                result = (from_atom_index,distance_index_1,distance_index_2,value,exponent)
        return result
    
    def get_table_name(self):
        return 'ATOM'
    
class Sidechain_context():
    
    
    
    def __init__(self, from_atom, residue_type, sidechain_atom_selection, table):
        
        self._table = table
        self.complete =  False
        
        from_atom_index = from_atom.index()
        
        segment = from_atom.segmentName()
        sidechain_residue_number = from_atom.residueNum()
        sidechain_atom = Atom_utils._select_atom_with_translation(segment, sidechain_residue_number, sidechain_atom_selection)
        
        self.sidechain_atom_index = sidechain_atom[0].index()
        
        if self.sidechain_atom_index != from_atom_index:
        
            self.residue_type = residue_type
            self.sidechain_atom_name = sidechain_atom[0].atomName()
            
            self.complete = True

class Sidechain_component_factory(Atom_component_factory):

    def _translate_atom_name(self, atom_name, context):
        return atom_name
    
    def _translate_atom_name_to_table(self, residue_type, atom_name,table):
        
        return  table.get_translation_to_table(residue_type,atom_name)

    
    def _build_contexts(self,atom, table):
        contexts = []
        residue_type =  atom.residueName()
        if residue_type in table.get_residue_types():
            for sidechain_atom in table.get_sidechain_atoms(residue_type):
                context = Sidechain_context(atom, residue_type, sidechain_atom,table)
                contexts.append(context)
        return contexts
    
    def  _get_component_for_atom(self, target_atom, context):
        table = context._table
        table = table

        
        target_atom_name = target_atom.atomName()
        target_atom_name = self._translate_atom_name(target_atom_name, context)
        target_residue_type = target_atom.residueName()
        #TODO move translation into context
        target_atom_name = self._translate_atom_name_to_table(target_residue_type,target_atom_name,table)
        
        result = None
        if target_atom_name in table.get_target_atoms():
            sidechain_atom_name = context.sidechain_atom_name
            residue_type = context.residue_type
            value = table.get_sidechain_coefficient(residue_type,target_atom_name,sidechain_atom_name)

            if value != None:
                target_atom_index = target_atom.index()
                sidechain_atom_index= context.sidechain_atom_index
                exponent = table.get_exponent()
                result = (target_atom_index,sidechain_atom_index,value,exponent)
        return result
    
    def get_table_name(self):
        return 'ATOM'
    
class Base_potential(object):
    
    __metaclass__ = abc.ABCMeta
    
    ALL = '(all)'
            
    def __init__(self):
        self._segment_manager = Segment_Manager()
        self._table_manager = Table_manager.get_default_table_manager()
        self._observed_shifts = Observed_shift_table()
        self._component_list_data  = {}
        self._component_factories = {}
        self._cache_list_data = {}
        
    
    def prepare(self):
        pass
    
    def get_component_table_names(self):
        result = []
        for component_factory in self._component_factories.values():
            result.append(component_factory.get_table_name())
        result.sort()
        return tuple(result)
    
    def _check_segment_length(self, segment):
        segment_info = self._segment_manager.get_segment_info(segment)
        if segment_info.segment_length < 3:
            template = "Warning: segment '%s' only contains < 3 residues (%i) and will be ignored"
            message = template % (segment, segment_info.segment_length)
            print >> sys.stderr, message
        return segment_info
    
    
    def _create_components_for_residue(self, name, segment, target_residue_number, atom_selection):
        selected_atoms = intersection(Atom_utils._select_atom_with_translation(segment, target_residue_number), atom_selection)
        
        component_factory = self._component_factories[name]
        table_source =  self._get_table_source()
        if component_factory.is_residue_acceptable(segment,target_residue_number,self._segment_manager):
            component_list =  self._component_list_data[name]
            component_factory.create_components(component_list, table_source, segment, target_residue_number,selected_atoms)
#                if component_factory.is_target_required(Atom_component_factory.ATOM):
#                    component_factory.create_atom_components(component_list, table, selected_atoms)
#                if component_factory.is_target_required(Atom_component_factory.RESIDUE):
#                    component_factory.create_residue_components(component_list, table, segment,target_residue_number)                    


    

    
    def _build_component_list(self,name,global_atom_selection):
        
        global_atom_selection = AtomSel(global_atom_selection)
        for segment in self._segment_manager.get_segments():
            
            if self._check_segment_length(segment):
                segment_info = self._segment_manager.get_segment_info(segment)
                for residue_number in range(segment_info.first_residue,segment_info.last_residue+1):
                    residue_atom_selection = Atom_utils._select_atom_with_translation(segment, residue_number)
                    target_atom_selection = intersection(residue_atom_selection,global_atom_selection)
                    self._create_components_for_residue(name,segment, residue_number, target_atom_selection)
                    
    
    def _add_component_factory(self, component_factory):
        self._component_factories[component_factory.get_table_name()] =  component_factory
    
    
    @abc.abstractmethod
    def _get_table_source(self):
        pass
        
    @abc.abstractmethod
    def get_abbreviated_name(self):
        pass

    def set_observed_shifts(self, shift_table):
        self._observed_shifts = shift_table
        
    def _get_cache_list(self,name):
        if not name in self._cache_list_data:
            self._cache_list_data[name] = Component_list()
        return self._cache_list_data[name]
    
    def clear_caches(self):
        self._cache_list_data = {}

#    TODO put 'ATOM' in a constant and rename to BB_ATOM?
    def _get_component_list(self,name=None):
        if name == None:
            name  = self._get_target_atom_list_name()
        if not name in self._component_list_data:
            self._component_list_data[name] = Component_list()
            self._build_component_list(name,"(all)")
        return self._component_list_data[name]
    
    def _get_target_atom_list_name(self):
        return 'ATOM'
    
    # TODO: make these internal
    def _get_component(self, index, name=None):

        
        components = self._get_component_list(name)
        component = components.get_component(index)
        return component
    
    def _get_components_for_id(self,list_name,lookup):
        return self._get_component_list(list_name).get_components_for_atom_id(lookup)
    
    def _get_all_components(self,name=None):
        components =  self._get_component_list(name)
        return components.get_all_components()
    
    def _get_number_components(self,name='ATOM'):
        components = self._get_component_list(name)
        
        return components.get_number_components()
    
    def _get_component_table_names(self):
        return self._component_list_data.keys()
    
    def get_target_atom_ids(self):
        components = self._get_component_list()
        return  components.get_component_atom_ids()
        
    def calc_single_atom_shift(self, target_atom_id):
        components =  self._get_component_list()

        result = 0.0
        if target_atom_id in components.get_component_atom_ids():
            index_range = components.get_component_range(target_atom_id)
            for index in range(*index_range):
                result += self._calc_component_shift(index)
        
        return result
    
    def calc_single_atom_force_set(self,target_atom_id,force_factor,forces):
        if self._have_derivative():
            components =  self._get_component_list()
        
        
            if target_atom_id in components.get_component_atom_ids():
                index_range = components.get_component_range(target_atom_id)
                for index in range(*index_range):
                    self._calc_single_force_set(index,force_factor,forces)
        return forces
        
    
    #TODO: make this just return a list in component order
    #TODO: remove
    def set_shifts(self, result):
        components =  self._get_component_list()
        
        for target_atom_id in components.get_component_atom_ids():
            result[target_atom_id] +=  self.calc_single_atom_shift(target_atom_id)
            
    def _have_derivative(self):
        return True
    
    def _get_or_make_target_force_triplet(self, forces, target_offset):
        target_forces = forces[target_offset]
        if target_forces == None:
            target_forces = [0.0] * 3
            forces[target_offset] = target_forces
        return target_forces

class Shift_calculator:
    DEFAULT_CUTOFF = 5.0
    DEFAULT_SMOOTHING_FACTOR = 1.0
    
    def __init__(self, indices, smoothed):
        self._target_atom_index = indices.target_atom_index
        self._distance_atom_index_1 =  indices.distance_atom_index_1
        self._distance_atom_index_2 =  indices .distance_atom_index_2
        self._exponent_index  = indices.exponent_index
        self._coefficient_index  = indices.coefficient_index
        self._components =  None
        self._smoothed =  smoothed
        self._smoothing_factor =  Shift_calculator.DEFAULT_SMOOTHING_FACTOR
        
        self._cutoff =  Shift_calculator.DEFAULT_CUTOFF
    
    def set_cutoff(self, cutoff):
        self._cutoff =  cutoff
    
    def set_smoothing_factor(self,smoothing_factor):
        self._smoothing_factor = smoothing_factor
        
    def set_components(self,components):
        self._components = components
    
    def _get_target_and_distant_atom_ids(self, index):
        values  = self._components.get_component(index)
        
        
        target_atom = values[self._distance_atom_index_1]
        distance_atom = values[self._distance_atom_index_2]
        return target_atom, distance_atom
    
    def _get_coefficient_and_exponent(self, index):
        values = self._components.get_component(index)
        


        
        coefficient = values[self._coefficient_index]
        exponent = values[self._exponent_index]
        
        return coefficient, exponent
    
    def __call__(self, indices, results):
        for index in indices:
            target_atom_index, sidechain_atom_index = self._get_target_and_distant_atom_ids(index)
            coefficient, exponent = self._get_coefficient_and_exponent(index)
            distance = Atom_utils._calculate_distance(target_atom_index, sidechain_atom_index)
            smoothing_factor = self._smoothing_factor
            if self._smoothed:
                ratio = distance / self._cutoff
        #            ratio2 = ratio**4
        #            for i in range(2):
        #                ratio *= ratio
        #            print ratio2,ratio
                smoothing_factor = 1.0 - ratio ** 8
            result = smoothing_factor * distance ** exponent * coefficient
        results[index] = result
    
class Distance_based_potential(Base_potential):
    
    def __init__(self,  smoothed = False, fast=False):
        super(Distance_based_potential, self).__init__()
        self._smoothed = smoothed
        #TODO: sort placement out
        self._cutoff = 5.0
        
        self.shift_calculator = Fast_shift_calculator(self._get_indices(), self._smoothed)

    class Indices(object):
        def __init__(self, target_atom_index,distance_atom_index_1,
                     distance_atom_index_2, coefficent_index, exponent_index):
            self.target_atom_index =  target_atom_index
            self.distance_atom_index_1 = distance_atom_index_1
            self.distance_atom_index_2 =  distance_atom_index_2
            self.exponent_index = exponent_index
            self.coefficient_index = coefficent_index
            
            
        def __str__(self):
            result = 'indices target= %i distance atom index 1 = %i index 2 = %i coefficent = %i exponent = %i'
            msg = result % (self.target_atom_index, self.distance_atom_index_1,
                            self.distance_atom_index_2, self.coefficient_index, self.exponent_index)
            return msg                 
            
    @abc.abstractmethod
    def _get_indices(self):
        return Distance_based_potential.Indices(target_atom_index=0,distance_atom_index_1=0,
                                                distance_atom_index_2=1,coefficent_index=2,
                                                exponent_index=3)
    
#    @abc.abstractmethod
    def _get_distance_list_name(self):
        return 'ATOM'

    def _get_target_and_distant_atom_ids(self, index):
        list_name = self._get_distance_list_name()
        values  = self._get_component(index,list_name)
        
        indices = self._get_indices()
        distance_atom_index_1 = indices.distance_atom_index_1
        distance_atom_index_2 = indices.distance_atom_index_2
        target_atom = values[distance_atom_index_1]
        distance_atom = values[distance_atom_index_2]
        return target_atom, distance_atom

    def _get_distance_components(self):
        list_name  = self._get_distance_list_name()
        return self._get_component_list(list_name)
        
        
    def _get_coefficient_and_exponent(self, index):
        list_name  = self._get_distance_list_name()
        
        values = self._get_component(index,list_name)
        
        indices = self._get_indices()
        
        coefficient_index = indices.coefficient_index
        exponent_index = indices.exponent_index
        
        coefficient = values[coefficient_index]
        exponent = values[exponent_index]
        
        return coefficient, exponent

    def _calc_single_force_factor(self,index,factor):
        
        target_atom, distance_atom = self._get_target_and_distant_atom_ids(index)
        
        coefficient, exponent = self._get_coefficient_and_exponent(index)
        
        target_pos = Atom_utils._get_atom_by_index(target_atom).pos()
        distant_pos =  Atom_utils._get_atom_by_index(distance_atom).pos()
        
        distances  = target_pos - distant_pos
        distance_2 = sum([elem**2 for elem in distances])

        factor= factor * coefficient
        
        if self._smoothed:
            ratio = distance_2 / self._cutoff**2
            ratio =  ratio**4
            pre_exponent = exponent - (exponent + 8.0) * ratio
        else:
            pre_exponent = exponent
            
        reduced_exponent = (exponent - 2.0) / 2.0
        
        force_factor = factor *  pre_exponent * distance_2 ** reduced_exponent

        return force_factor




    def _calc_single_force_set(self,index,factor, forces):
#        values  = self._get_component(index)
#        
#        indices = self._get_indices()
#        
#        distance_atom_index_1 = indices.distance_atom_index_1
#        distance_atom_index_2 =indices.distance_atom_index_2
#        
#        target_atom=values[distance_atom_index_1]
#        distance_atom = values[distance_atom_index_2]
        target_atom,distant_atom =  self._get_target_and_distant_atom_ids(index)
        
        target_pos = Atom_utils._get_atom_by_index(target_atom).pos()
        distant_pos =  Atom_utils._get_atom_by_index(distant_atom).pos()
        
        distance  = target_pos - distant_pos
        
        force_factor  = self._calc_single_force_factor(index, factor)
        

        target_offset = target_atom
        distant_offset = distant_atom
        X_OFFSET = 0
        Y_OFFSET = 1
        Z_OFFSET = 2
        
#       TODO: move to atom utils 
        OFFSETS_3 = (X_OFFSET,Y_OFFSET,Z_OFFSET)
        
        target_forces = self._get_or_make_target_force_triplet(forces, target_offset)
        distant_forces  = self._get_or_make_target_force_triplet(forces, distant_offset)
        for offset in OFFSETS_3:
            target_forces[offset] -= distance[offset] * force_factor
            distant_forces[offset] += distance[offset] * force_factor
        

    


    def _calc_component_shift(self, index):
        components = self._get_distance_components()
        self.shift_calculator.set_components(components)
        return self.shift_calculator(index)
    
class Distance_potential(Distance_based_potential):
    '''
    classdocs
    '''


    def __init__(self):
        super(Distance_potential, self).__init__()
        
        '''
        Constructor
        '''
        
        self._add_component_factory(Distance_component_factory())


    
    def _get_indices(self):
        return super(Distance_potential, self)._get_indices()
       
    def get_abbreviated_name(self):
        return BACK_BONE
    

    def _get_table_source(self):
        return self._table_manager.get_BB_Distance_Table

            
    def __str__(self):
        print self._get_number_components()
        result = []
        for from_index,to_index,value,exponent in self._get_all_components():
            from_atom = self._get_atom_name(from_index)
            to_atom = self._get_atom_name(to_index)
            
            template = '[%s] - [%s] %7.3f %7.3f'
            values = from_atom, to_atom, value, exponent
            result.append( template % values)
        return '\n'.join(result)

    


    

    
    def dump(self):
        result  = []
        for from_index,distance_index_1,value,exponent in self._get_all_components():
            sub_result  = []
            sub_result.append(Atom_utils._get_atom_info_from_index(from_index))
            
            sub_result.append(Atom_utils._get_atom_info_from_index(distance_index_1))
            
            sub_result.append(value)
            
            sub_result.append(exponent)
            
            result.append(tuple(sub_result))
        return result   


class Extra_potential(Distance_based_potential):
    def __init__(self):
        super(Extra_potential, self).__init__()
        
        self._add_component_factory(Extra_component_factory())
        
    
    def get_abbreviated_name(self):
        return XTRA
#
#    def _translate_atom_name(self, atom_name,context):
#        return context._table.get_translation(atom_name)
        
    def _get_table_source(self):
        return self._table_manager.get_extra_table
    
    


    def _get_indices(self):
        return Distance_based_potential.Indices(target_atom_index=0,distance_atom_index_1=1,
                                                distance_atom_index_2=2,coefficent_index=3,
                                                exponent_index=4)    
    
    def __str__(self):
        print len( self._components)
        result = []
        for from_index,distance_index_1,distance_index_2,value,exponent in self._components:
            from_atom = self._get_atom_name(from_index)
            distance_atom_1 = self._get_atom_name(distance_index_1)
            distance_atom_2 = self._get_atom_name(distance_index_2)
            
            template = '[%s] - [%s] - [%s] %7.3f %7.3f'
            values = from_atom, distance_atom_1, distance_atom_2, value, exponent
            result.append( template % values)
        return '\n'.join(result)

    def dump(self):
        result  = []
        for from_index,distance_index_1,distance_index_2,value,exponent in self._get_all_components():
            sub_result  = []
            sub_result.append(Atom_utils._get_atom_info_from_index(from_index))
            
            sub_result.append(Atom_utils._get_atom_info_from_index(distance_index_1))
            sub_result.append(Atom_utils._get_atom_info_from_index(distance_index_2))
            
            sub_result.append(value)
            
            result.append(tuple(sub_result))
        return result
    




    
class RandomCoilShifts(Base_potential):
    

    def __init__(self):
        super(RandomCoilShifts, self).__init__()
        
        self._add_component_factory(Random_coil_component_factory())

    
    def get_abbreviated_name(self):
        return RANDOM_COIL
    
    def _have_derivative(self):
        False
         
    
    def _get_table_source(self):
        return self._table_manager.get_random_coil_table

    
    def _calc_component_shift(self,index):
        components = self._get_component_list()
        return components.get_component(index)[1]


    def __str__(self): 
        result = []
        for i,elem in enumerate(self._get_component_list()):
            template = "%-5i  %15s % 8.3f"
            atom_index,shift = elem
            atom_name =  Atom_utils._get_atom_name(atom_index)
            string = template % (i, atom_name,shift)
            result.append(string)
        return '\n'.join(result)

class Dihedral_potential(Base_potential):

    
    def __init__(self):
        Base_potential.__init__(self)
        self._add_component_factory(Dihedral_component_factory())
    


    def get_abbreviated_name(self):
        return DIHEDRAL

    def _get_table_source(self):
        return self._table_manager.get_dihedral_table

    
                
        



    
    def dump(self):
        result  = []
        for from_index,atom_1,atom_2, atom_3, atom_4,value,param_0,param_1,param_2,param_3,param_4,exponent in self._get_all_components():
            sub_result  = []
            key = []
            sub_result.append(key)
            key.append(tuple(Atom_utils._get_atom_info_from_index(from_index)[1:]))
            
            dihedral_atom_info = []
            for atom_pair in (atom_1,atom_2), (atom_3, atom_4):
                
                info_pair  = []
                #TODO flatten this array
                info_pair.append(Atom_utils._get_atom_info_from_index(atom_pair[0])[1:])
                info_pair.append(Atom_utils._get_atom_info_from_index(atom_pair[1])[1:])
            
                dihedral_atom_info.append(tuple(info_pair))
            key.append(tuple(dihedral_atom_info))
            
            sub_result.append(value)
            sub_result.extend((param_0,param_1,param_2,param_3,param_4,))
            sub_result.append(exponent)
            
            result.append(tupleit(sub_result))
            
            
        return tuple(result)
    
    def __str__(self):
        print len( self._components)
        result = []
        for from_index, \
            atom_1,atom_2, atom_3, atom_4,                           \
            value,                                                   \
            param_0,param_1,param_2,param_3,param_4,                 \
            exponent in self._components:
                                             
            from_atom = self._get_atom_name(from_index)
            distance_atom_1 = self._get_atom_name(atom_1)
            distance_atom_2 = self._get_atom_name(atom_2)
            distance_atom_3 = self._get_atom_name(atom_3)
            distance_atom_4 = self._get_atom_name(atom_4)
                        
            template = '[%s] - [[%s] - [%s] - [%s] - [%s]] (%7.3f %7.3f %7.3f %7.3f %7.3f) %7.3f %7.3f'
            values = from_atom, distance_atom_1, distance_atom_2,         \
                                distance_atom_3, distance_atom_4,         \
                                param_0,param_1,param_2,param_3,param_4,  \
                                value, exponent
            result.append( template % values)
        return '\n'.join(result)

    def _get_dihedral_angle(self, dihedral_1_atom_id_1, dihedral_1_atom_id_2, 
                                 dihedral_2_atom_id_1, dihedral_2_atom_id_2):
        
        atom_1  = Atom_utils._get_atom_by_index(dihedral_1_atom_id_1)
        atom_2  = Atom_utils._get_atom_by_index(dihedral_1_atom_id_2)
        atom_3  = Atom_utils._get_atom_by_index(dihedral_2_atom_id_1)
        atom_4  = Atom_utils._get_atom_by_index(dihedral_2_atom_id_2)
        
        return Dihedral(atom_1,atom_2,atom_3,atom_4).value();
    
    




    def _get_dihedral_atom_ids(self, index):
        component = self._get_component(index)
        
        dihedrals = component[1:5]
        return dihedrals


    def _get_coefficient(self, index):
        component = self._get_component(index)
        coefficient = component[5]
        return coefficient


    def _get_parameters(self, index):
        component = self._get_component(index)
        parameters = component[6:11]
        parameter_0, parameter_1, parameter_2, parameter_3, parameter_4 = parameters
        return parameter_0, parameter_3, parameter_1, parameter_4, parameter_2
    
    def _get_force_parameters(self,index):
        return self._get_parameters(index)[:-1]

    def _calc_component_shift(self, index):
        
       
        dihedral_atom_ids= self._get_dihedral_atom_ids(index)
        
        coefficient = self._get_coefficient(index)
        
        parameter_0, parameter_3, parameter_1, \
        parameter_4, parameter_2               \
        = self._get_parameters(index)
        
        angle = self._get_dihedral_angle(*dihedral_atom_ids)

        
        angle_term = parameter_0 * cos(3.0 * angle + parameter_3) + \
                     parameter_1 * cos(angle + parameter_4) +       \
                     parameter_2
        
        shift = coefficient * angle_term

        return shift
    
    
#    TODO make this consistent with the distance forces factor
    def _calc_single_force_factor(self,index):
        
        dihedral_atom_ids= self._get_dihedral_atom_ids(index)
        
        
        parameter_0, parameter_3, parameter_1, \
        parameter_4 = self._get_force_parameters(index)
        
        angle = self._get_dihedral_angle(*dihedral_atom_ids)
        
        result = -3.0 * parameter_0 * sin(3.0 * angle + parameter_3) - \
                        parameter_1 * sin(angle + parameter_4)
        return result

    
    #TODO: is this too close?
    def _calc_single_force_set(self,index,factor,forces):
        dihedral_factor = self._calc_single_force_factor(index)
        dihedral_atom_ids= self._get_dihedral_atom_ids(index)
        
#        ATOM_ID_1 = 0
        ATOM_ID_2 = 1
        ATOM_ID_3 = 2
#        ATOM_ID_4 = 3
        
        positions = []
        for atom_id in dihedral_atom_ids:
            positions.append(Atom_utils._get_atom_pos(atom_id))
            
        v1,v2,v3,v4 = positions
        
        r1 = v1 - v2
        r2 = v3 - v2
        r3 = r2 * -1
        r4 = v4 - v3
        
        # compute normal vector to plane containing v1, v2, and v3
        n1 = cross(r1, r2)
        # compute normal vector to plane containing v2, v3, and v4
        n2 = cross(r3, r4)
                
        r2_length = Atom_utils._calculate_distance(dihedral_atom_ids[ATOM_ID_2], dihedral_atom_ids[ATOM_ID_3])
        r2_length_2 = r2_length**2.0
        

        weight = factor * self._get_coefficient(index)
        

#        // force calculation according to Bekker, Berendsen and van Gunsteren (1995),
#        // Journal of Computational Chemistry 16, pp. 527-533:
#        // Force and virial of torsional-angle-dependent potentials.
        factor_1 = dihedral_factor * r2_length;
        F1 = n1 *  (-factor_1 / norm(n1)**2)
        F4 = n2 *  ( factor_1 / norm(n2)**2)
        
        factor_2 = dot(r1, r2) / r2_length_2;
        factor_3 = dot(r3, r4) / r2_length_2;

        T1 = F1 * (factor_2-1);
        T2 = F4 * -factor_3
        
        F2 = T1 + T2
        
        T1 = F1 * -factor_2
        T2 = F4 * (factor_3 -1)

        F3 = T1 + T2

#        // assign forces
        X_OFFSET = 0
        Y_OFFSET = 1
        Z_OFFSET = 2
        
        OFFSETS_3 = (X_OFFSET,Y_OFFSET,Z_OFFSET)
        
        for atom_id,base_force in zip(dihedral_atom_ids,[F1,F2,F3,F4]):
            force_triplet = self._get_or_make_target_force_triplet(forces, atom_id)
            for offset in OFFSETS_3:
                force_component = weight * base_force[offset]
                force_triplet[offset]+= force_component
                

    
class Sidechain_potential(Distance_based_potential):
    
    def __init__(self):
        super(Sidechain_potential, self).__init__()
        self._add_component_factory(Sidechain_component_factory())
        
        
    def  _get_table_source(self):
        return self._table_manager.get_sidechain_table
    
    def _get_indices(self):
        return super(Sidechain_potential, self)._get_indices()
    
        

    
    def  get_abbreviated_name(self):
        return SIDE_CHAIN
    


    
    def dump(self):
        result  = []
        for from_index,sidechain_index,value,exponent in self._get_all_components():
            sub_result  = []
            
            sub_result.append(Atom_utils._get_atom_info_from_index(from_index))
            sub_result.append(Atom_utils._get_atom_info_from_index(sidechain_index))
            
            sub_result.append(value)
            sub_result.append(exponent)
            
            result.append(tuple(sub_result))
        return result
    
    def __str__(self):
        result = []
        for target_index,sidechain_index,value,exponent in self._components:
            target_atom = self._get_atom_name(target_index)
            sidechain_atom = self._get_atom_name(sidechain_index)
            
            template = '[%s] - [%s] %7.3f %7.3f'
            values = target_atom, sidechain_atom, value, exponent
            result.append( template % values)
        return '\n'.join(result)


class Backbone_atom_indexer:
    #TODO add table base to allow for varying numbers of atoms in table//

    def _get_target_atoms(self, table):
        target_atoms = table.get_target_atoms()
        return target_atoms


    def _get_table_index(self, table):
        return table._table['index']


    def _get_table_offset(self, table):
        target_atoms = self._get_target_atoms(table)
        num_target_atoms = len(target_atoms)
        offset = num_target_atoms * self._get_table_index(table)
        return offset

    def _get_atom_id(self,atom_name,table):
        target_atoms = self._get_target_atoms(table)
        offset= self._get_table_offset(table)
        return offset + target_atoms.index(atom_name)

    def _get_atom_name_offsets(self,table):
        result = []
        
        offset = self._get_table_offset(table)
        for i,atom_name in enumerate(self._get_target_atoms(table)):
            result.append((table._table['residue_type'], atom_name, offset+i,))
            
        return tuple(result)
    
    def _get_offset_names(self,table, add_to=None):
        result = add_to
        if result == None:
            result = {}
        
        for residue_type, atom_name,atom_index in self._get_atom_name_offsets(table):
            result[atom_index] = (residue_type, atom_name)
        return result
        

class Ring_backbone_context(object,Backbone_atom_indexer):
    
    def _translate_atom_name_to_table(self, residue_type, atom_name,table):
        return  table.get_translation_to_table(residue_type,atom_name)
    
    def __init__(self, atom, table):
        self.complete =  False
        
        self.target_atom_id  =  atom.index()
        
        segment = atom.segmentName()
        
        atom_name  =  Atom_utils._get_atom_name_from_index(self.target_atom_id)
        residue_number = atom.residueNum()
        residue_type = Atom_utils._get_residue_type(segment, residue_number)
        
        atom_name = self._translate_atom_name_to_table(residue_type, atom_name, table)
        
        if atom_name in table.get_target_atoms():
            self.atom_name_index =self._get_atom_id(atom_name, table)
            if self.atom_name_index >= 0:
                self.complete = True

class Ring_factory_base():
    #TODO: add a general test that the component is active for the structure

    def _residue_is_ok(self, segment, residue_number, table):
        residue_type = Atom_utils._get_residue_type(segment, residue_number)
        residue_ok = residue_type in table.get_residue_types()
        return residue_ok

    def _have_targets(self,table):
        result = False
        segment_manager =  Segment_Manager()
        for segment in segment_manager.get_segments():
            if result == True:
                break
            else:
                segment_info = segment_manager.get_segment_info(segment)
            
                for residue_number in range (segment_info.first_residue, segment_info.last_residue+1):
                    if self._residue_is_ok(segment, residue_number, table):
                        result = True
                        break
        return result
    
class Ring_backbone_component_factory (Atom_component_factory, Ring_factory_base):
    
    
    def get_table_name(self):
        return Atom_component_factory.get_table_name(self)
        
# TODO: impelement
    def _translate_atom_name(self, atom_name, context):
        return atom_name

    def _get_component_for_atom(self, atom, context):
        if context.complete:
            return context.target_atom_id, context.atom_name_index
    
    def _build_contexts(self, atom, table):
        contexts = []
        if  self._have_targets(table):
            context = Ring_backbone_context(atom,table)
            if context.complete:
                contexts.append(context)
        return contexts

class Ring_sidechain_component_factory(Residue_component_factory,Ring_factory_base):
#    TODO: make ring index a service not a mixin same with atom ids
    def __init__(self):
        self.ring_index_count = 0
        self.ring_index = {}
        self.ring_id_ring_type_map = {}
    
    def _get_ring_type_from_id(self, ring_id):
        return self.ring_id_ring_type_map[ring_id]
    
    def _get_or_make_ring_ids(self, segment,residue_number,table):
        residue_key = segment,residue_number
        
        residue_type = Atom_utils._get_residue_type(segment, residue_number)
        
        ring_ids = []
        if residue_type in table.get_residue_types():
            for ring_type in table.get_ring_types(residue_type):
                residue_key = segment, residue_number, ring_type
                
                if residue_key in self.ring_index:
                    ring_id = self.ring_index[residue_key]
                else:
                    ring_id = self.ring_index_count
                    self.ring_index[residue_key] = self.ring_index_count
                    self.ring_index_count += 1
                ring_ids.append(ring_id)
                self.ring_id_ring_type_map[ring_id] = ring_type
        return ring_ids

    def create_residue_components(self, component_list, table, segment, residue_number):
        if  self._have_targets(table):
            residue_type = Atom_utils._get_residue_type(segment,residue_number)
            
            if  residue_type in table.get_residue_types():
                
                ring_ids = self._get_or_make_ring_ids(segment,residue_number,table)
                for ring_id in ring_ids:
                    self.create_ring_components(component_list, table, segment, residue_number, ring_id)
            
    @abc.abstractmethod
    def create_ring_components(self, component_list,table, segment, residue, ring_id):
        pass
    
class Ring_sidechain_atom_factory(Ring_sidechain_component_factory):
    
    def __init__(self):
        super(Ring_sidechain_atom_factory, self).__init__()
        
        
    def get_table_name(self):
        return 'RING'

    #TODO: hang everything off ring_id?
    def create_ring_components(self, component_list,table, segment, residue_number, ring_id):
        if  self._have_targets(table):
            residue_atoms =  Atom_utils.find_atom(segment, residue_number)
            residue_type = Atom_utils._get_residue_type(segment, residue_number)
            
            
            
            atom_name_id_map = dict([(atom.atomName(),atom.index()) for atom in residue_atoms])
            
            self._get_ring_type_from_id(ring_id)
            ring_atoms = []
            ring_type = self._get_ring_type_from_id(ring_id)
            for atom_name in table.get_ring_atoms(residue_type,ring_type):
                ring_atoms.append(atom_name_id_map[atom_name]) 
            sub_result  = ring_id,tuple(ring_atoms)
            component_list.add_component(sub_result)
                    
            

# BACK_BONE-ID -> BACK_BONE-TYPE [iterate this second for find all aromatic shifts for an atom]
# BACK_BONE-TYPE ->  AROMATIC-ID [aromatic ring number] COEFF
# AROMATIC-ID -> ATOMS [iterate this first to build ring normals and centres


class Ring_coefficient_component_factory(Ring_sidechain_component_factory,Backbone_atom_indexer):
    def __init__(self):
        super(Ring_coefficient_component_factory, self).__init__()
        self._seen_tables  = set()
        
    def get_table_name(self):
        return 'COEF'

    def create_ring_components(self, component_list, table, segment, residue_number, ring_id):
#        print 'table type', table._table['residue_type']
        residue_type = Atom_utils._get_residue_type(segment, residue_number)
        ring_type = self._get_ring_type_from_id(ring_id)
        
#        print residue_number, table._table['table_type'], table._table['residue_type'],table.get_target_atoms()
        for target_atom_name in table.get_target_atoms():
            atom_id = self._get_atom_id(target_atom_name, table)
            if atom_id > -1:
                coef = table.get_ring_coefficient(target_atom_name,residue_type,ring_type)
                if coef != None:
                    coef_component = atom_id,ring_id,coef
                    component_list.add_component(coef_component)
    

    def _is_new_table_type(self, table):
        table_residue_type = table._table['residue_type']
        
        result = not table_residue_type in self._seen_tables
        self._seen_tables.add(table_residue_type)
        
        return result

    def create_residue_components(self, component_list, table, segment, residue_number):
        class residue_type_in_table(object):
            def __init__(self,accepted_residue_types=None):
                self._accepted_residue_types =  accepted_residue_types
                
            def __call__(self, segment, residue_number):
                
                result = False
                residue_type = Atom_utils._get_residue_type(segment, residue_number)
                
                if residue_type in self._accepted_residue_types:
                    result = True
                
                return result
        
        if  self._have_targets(table):
            if self._is_new_table_type(table):
#                residue_type = Atom_utils._get_residue_type(segment,residue_number)
                
#                print 'here2, residue_number type:', residue_type, table.get_residue_types(), segment, residue_number
                
                for segment,residue_number in iter_residues_and_segments(residue_type_in_table(table.get_residue_types())):
#                    print segment,residue_number
                    ring_ids = self._get_or_make_ring_ids(segment,residue_number,table)
                    for ring_id in ring_ids:
                        self.create_ring_components(component_list, table, segment, residue_number, ring_id)

    
class Ring_Potential(Base_potential):
    def __init__(self):
        super(Ring_Potential, self).__init__()
        
        self._add_component_factory(Ring_backbone_component_factory())
        self._add_component_factory(Ring_coefficient_component_factory())
        self._add_component_factory(Ring_sidechain_atom_factory())
        
        
    def get_abbreviated_name(self):
        return RING
    
    
    def _get_table_source(self):
        return self._table_manager.get_ring_table
    

    

    RING_ATOM_IDS = 1

    def _average_vec3(self, positions):
        result = Vec3(0.0,0.0,0.0)
        for position in positions:
            result += position
        
        result /= len(positions)
        return result

    def _calculate_one_ring_centre(self, ring_component):
        
        atom_ids = ring_component[self.RING_ATOM_IDS]
        positions = []
        for atom_id in atom_ids:
            positions.append(Atom_utils._get_atom_pos(atom_id))
            
        result = self._average_vec3(positions)
        
        return result
            
    
    
    def _calculate_ring_centres(self):
        ring_components = self._get_component_list('RING')
        
        result= []
        
        for ring_component in ring_components:
            result.append(self._calculate_one_ring_centre(ring_component))
            
        return result
    


    def _check_ring_size_ok(self, atom_ids):
        num_atom_ids = len(atom_ids)
        if num_atom_ids < 5 or num_atom_ids > 6:
            template = "ring normals function is only implemented for 5 or six member rings i got %d atoms"
            msg = template % num_atom_ids
            raise Exception(msg)

    #TODO could try newells method http://www.opengl.org/wiki/Calculating_a_Surface_Normal
    def _calculate_one_ring_normal(self, ring_component):
        atom_ids = ring_component[self.RING_ATOM_IDS]
        self._check_ring_size_ok(atom_ids)
        
        atom_triplets = atom_ids[:3],atom_ids[-3:]
        
        normals  = []
        for atom_triplet in atom_triplets:
            atom_vectors = []
            for atom_id in atom_triplet:
                atom_vectors.append(Atom_utils._get_atom_pos(atom_id))
            vec_1 = atom_vectors[0] -atom_vectors[1]
            vec_2 =  atom_vectors[2] - atom_vectors[1]
                
            normals.append(cross(vec_1,vec_2))
        
        result = self._average_vec3(normals)
       
        return result
    
    
    def _calculate_ring_normals(self):
        ring_components = self._get_component_list('RING')
        
        result= []
        
        for ring_component in ring_components:
            result.append(self._calculate_one_ring_normal(ring_component))
            
        return result
    
    
    def _get_ring_centre(self, ring_id):
        ring_component = self._get_component_list('RING').get_component(ring_id)
        return self._calculate_one_ring_centre(ring_component)
    
    def _get_ring_normal(self, ring_id):
        ring_component = self._get_component_list('RING').get_component(ring_id)
        return self._calculate_one_ring_normal(ring_component)
    

    def _calc_sub_component_shift(self, target_atom_id,  sub_component_id):
        
        target_atom_id, atom_type_id = self._get_component_list('ATOM')[target_atom_id]
        atom_type_id, ring_id, coefficient = self._get_component_list('COEF').get_components_for_atom_id(atom_type_id)[sub_component_id]
        
        target_atom_pos = Atom_utils._get_atom_pos(target_atom_id)
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

    def _calc_component_shift(self,component_id):

        atom_type_id = self._get_component_list('ATOM')[component_id][1]

        result = 0.0
        
        num_coef_components = len(self._get_component_list('COEF').get_components_for_atom_id(atom_type_id))
        for coef_component_id in range(num_coef_components):
            result += self._calc_sub_component_shift(component_id,  coef_component_id)
            
        return result
    

    #TODO: use this more places
    def _get_ring_atom_ids(self, ring_id):
        ring_component = self._get_component_list('RING').get_component(ring_id)
        return ring_component[self.RING_ATOM_IDS]
       
    
    
    def _get_ring_atom_positions(self, ring_id):
        ring_atom_ids =  self._get_ring_atom_ids(ring_id)
        
        result  = []
        for ring_atom_id in ring_atom_ids:
            result.append(Atom_utils._get_atom_pos(ring_atom_id))
        return result
    
    
    class Force_sub_terms(object):
        
        def _get_ring_centre(self, ring_id):
            ring_id,ring_centre = self._get_cache_list('CENT').get_component(ring_id)
            return ring_centre
        
        def _get_ring_normal(self, ring_id):
            ring_id,ring_normal = self._get_cache_list('NORM').get_component(ring_id)
            return ring_normal
        
        
    
#TODO: remove doubling of lists
        def _get_component_list(self, name):
            return self._component_list_dict[name]
        
        
        def _get_cache_list(self, name):
            return self._cache_list_data[name]
        
        def __init__(self, target_atom_id, ring_id,component_list_dict, cache_list_data):
            self.target_atom_id= target_atom_id
            self._component_list_dict = component_list_dict
            self._cache_list_data =  cache_list_data
            target_atom_pos = Atom_utils._get_atom_pos(target_atom_id)
            
            self.ring_centre = self._get_ring_centre(ring_id)
            self.ring_normal = self._get_ring_normal(ring_id)
            
            # distance vector between atom of interest and ring center
            self.d = target_atom_pos - self.ring_centre
            self.dL = norm(self.d)
            
            # squared distance of atom of interest from ring center
            dL2 = dot(self.d, self.d) #            if (dL2 < 0.5) cout << "CAMSHIFT WARNING: Distance between atom and center of ring alarmingly small at " << sqrt(dL2) << " Angstrom!" << endl;
            
            # calculate terms resulting from differentiating energy function with respect to query and ring atom coordinates
            self.dL4 = self.dL ** 4
            self.dL3 = self.dL ** 3
            self.dL6 = self.dL3 ** 2
            
            self.nL = norm(self.ring_normal)
            nL2 = self.nL ** 2
            self.dLnL = self.dL * self.nL
            self.dL3nL3 = self.dL3 * nL2 * self.nL
            
            self.dn = dot(self.d, self.ring_normal)
            dn2 = self.dn ** 2
            
            self.u = 1.0 - 3.0 * dn2 / (dL2 * nL2)
            
            factor = -6.0 * self.dn / (self.dL4 * nL2)
            self.gradUQ = [0.0] *3
            for axis in AXES:
                self.gradUQ[axis] = factor * (dL2 * self.ring_normal[axis] - self.dn * self.d[axis])
                
            
            factor = 3 *self.dL
            self.gradVQ =  [0.0] * 3
            for axis in AXES:
                self.gradVQ[axis]= factor * self.d[axis]
            
#            return  dL3, u, dL6, ring_normal,  atom_type_id, coefficient, d, factor, dn, dL3nL3, dL, nL, dLnL

    #To close? currently this is a direct port

    def _calc_target_atom_forces(self, target_atom_id, ring_id, force_factor, sub_terms, forces):
        
        X = 0
        Y = 1
        Z = 2
        AXES = X, Y, Z
        target_force_triplet = self._get_or_make_target_force_triplet(forces, target_atom_id)
    # update forces on query atom
        for axis in AXES:
            #               f.coor[pos1  ] += -fact * (gradUQx * v - u * gradVQx) / v2;
            #print "terms",-force_factor, sub_terms.gradUQ[0], sub_terms.dL3, sub_terms.u, sub_terms.gradVQ[axis], sub_terms.dL6
            target_force_triplet[axis] += -force_factor * (sub_terms.gradUQ[axis] * sub_terms.dL3 - sub_terms.u * sub_terms.gradVQ[axis]) / sub_terms.dL6
        

    #TODO: calculation of GradU and gradV are not consistent with force_terms for target atom correct
    #TODO: reduce number of parameters to method
    def _calculate_ring_forces(self, atom_type_id, ring_id, force_factor, force_terms, forces):
        nSum = force_terms.ring_normal * 2.0 #            float_type g [3], ab [3], c [3]
        ring_atoms = self._get_ring_atom_ids(ring_id)
        ring_atom_positions = self._get_ring_atom_positions(ring_id)
    #// 2 for a 5-membered ring, 3 for a 6-membered ring
        num_ring_atoms = len(ring_atoms)
        limit = num_ring_atoms - 3

#        for atom_type_id, ring_id, coefficient in coef_components:
        g = [0.0] * 3
        for ring_atom_index in range(num_ring_atoms):
            ring_atom_id = ring_atoms[ring_atom_index]
            if ring_atom_index < limit:
                for axis in AXES:
                    index_1 = (ring_atom_index + 1) % 3
                    index_2 = (ring_atom_index + 2) % 3
                    g[axis] = ring_atom_positions[index_1][axis] - ring_atom_positions[index_2][axis] # atoms 3,4 (5 member) or 3,4,5 (6 member)
            
            else:
                if  ring_atom_index >= num_ring_atoms - limit:
                    offset = num_ring_atoms - 3 #2 for a 5-membered ring, 3 for a 6-membered ring
                    for axis in AXES:
                        index_1 = (ring_atom_index + 1 - offset) % 3 + offset
                        index_2 = (ring_atom_index + 2 - offset) % 3 + offset
                        g[axis] = ring_atom_positions[index_1][axis] - ring_atom_positions[index_2][axis]
                else:
                    
                    for axis in AXES:
                        g[axis] = ring_atom_positions[0][axis] - ring_atom_positions[1][axis] + ring_atom_positions[3][axis] - ring_atom_positions[4][axis]
        
            # 0 1 2 2 1   (0+1) %3 (0+2) %3
            # 1 2 0 0 2
            # 2 0 1 10
            #atom 2 (5-membered rings)
            indices_1 = [0] * 3
            indices_2 = [0] * 3
            for axis in AXES:
                indices_1[axis] = (axis + 1) % 3
                indices_2[axis] = (axis + 2) % 3
            
            ab = [0.0] * 3
            for axis,index_1, index_2 in zip(AXES,indices_1, indices_2):
                ab[axis] = force_terms.d[index_1] * g[index_2] - force_terms.d[index_2] * g[index_1]

            c = [0.0] * 3
            for axis, index_1, index_2 in zip(AXES,indices_1, indices_2):
                c[axis] = nSum[index_1] * g[index_2] - nSum[index_2] * g[index_1]
            
            factor = -6.0 * force_terms.dn / force_terms.dL3nL3
            factor2 = 0.25 * force_terms.dL / force_terms.nL
            one_over_num_ring_atoms = 1.0 / float(num_ring_atoms)
            factor3 = force_terms.nL / force_terms.dL * one_over_num_ring_atoms
            
            gradU = [0.0] * 3
            for axis in AXES:
                gradU[axis] = factor * ((0.5 * ab[axis] - force_terms.ring_normal[axis] * one_over_num_ring_atoms) * force_terms.dLnL - force_terms.dn * (factor2 * c[axis] - factor3 * force_terms.d[axis]))
                
            gradV = [0.0] * 3
            factor = -3 * force_terms.dL * one_over_num_ring_atoms
            for axis in AXES:
                gradV[axis] = factor * force_terms.d[axis]
            
#TODO: rename force terms  sub_terms or vice versa
            ring_target_force_triplet = self._get_or_make_target_force_triplet(forces, ring_atom_id)
            for axis in AXES:
                sub_force = -force_factor * (gradU[axis] * force_terms.dL3 - force_terms.u * gradV[axis]) / force_terms.dL6
                ring_target_force_triplet[axis] += sub_force
#                print AXIS_NAMES[axis],sub_force,-force_factor, gradU[axis], force_terms.dL3, force_terms.u, gradV[axis],force_terms.dL6
#            print
#        print 

    
    def _build_ring_data_cache(self):
        #TODO: remove double normal calculation
        normals = self._get_cache_list('NORM')
        normals.clear()
        centres =  self._get_cache_list('CENT')
        centres.clear()
        
        for ring_component in self._get_component_list('RING'):
            ring_id =  ring_component[0]
            
            normal = self._calculate_one_ring_normal(ring_component)
            normal_component = ring_id,normal
            normals.add_component(normal_component)
            
            
            centre = self._calculate_one_ring_centre(ring_component)
            centre_component = ring_id,centre
            centres.add_component(centre_component)
    
    

    def _build_force_terms(self, target_atom_id, ring_id):
        force_terms = self.Force_sub_terms(target_atom_id, ring_id, self._component_list_data, self._cache_list_data)
        return force_terms


    def calculate_single_ring_forces(self, target_atom_id, atom_type_id, component, force_factor, forces):
        atom_type_id, ring_id, coefficient = component
#        print atom_type_id, ring_id, coefficient, self._get_component_list('RING').get_components_for_atom_id(ring_id)
        force_terms = self._build_force_terms(target_atom_id, ring_id)
        
#        print Atom_utils._get_atom_info_from_index(target_atom_id)
        self._calc_target_atom_forces(target_atom_id, ring_id, force_factor * coefficient, force_terms, forces)
        #            #TODO: this is not how camshift does it, it uses the sum of the two ring normals
        self._calculate_ring_forces(atom_type_id, ring_id, force_factor * coefficient, force_terms, forces)

    def calc_single_atom_force_set(self, target_atom_id, force_factor, forces):
        target_atom_components = self._get_component_list('ATOM').get_components_for_atom_id(target_atom_id)
        
        if len(target_atom_components) > 0:
            target_atom_id, atom_type_id = target_atom_components[0]
            
            #TODO: this prompts the coef_component list for ring to be created but shouldn't be needed
            # we need to populate the lists automatilly and make it part of the lists implementation
            self._get_component_list('RING')
            

            
            coef_components = self._get_component_list('COEF').get_components_for_atom_id(atom_type_id)
            for i,coef_component  in enumerate(coef_components):
#                print 'coef_component', i, coef_component
#                 print 'here', target_atom_id, atom_type_id,
#                coef_component, force_factor
                 self.calculate_single_ring_forces(target_atom_id, atom_type_id, coef_component, force_factor, forces)

        
#    def _get_component_for_atom(self, atom, context):
#        return []

#class Non_bonded_atom_context():
#    def __init__(self,atom_id,target_atom_name,table):
#        atom_name = Atom_utils._get_atom_name_from_index(atom_id)
#        
#        if atom_name in table.get_target_atoms():
#            self.complete = True
#            
#            self.atom_type_id  = self._get_atom_id(self,atom_name,table)
#            self.atom_id = atom_id
#                
#class Non_bonded_backbone_component_factory(Atom_component_factory, Backbone_atom_indexer):
#    
#    def get_table_name(self):
#        return 'ATOM'
#    
#    def _build_contexts(self, atom, table):
#        contexts = []
#        
#        for target_atom in table.get_target_atoms():
#            
#            context = Non_bonded_atom_context(atom,target_atom,table)
#            
#            if context.complete:
#                contexts.append(context)
#                
#        return contexts
#
#    def _get_component_for_atom(self, atom, context):
#        table = context._table
#        
#        from_atom_name = atom.atomName()
#        from_atom_name = self._translate_atom_name(from_atom_name, context)
#        
#        result = None
#        if from_atom_name in table.get_target_atoms():
#            
#            dihedral_key = context.dihedral_key
#            value = context._table.get_dihedral_shift(from_atom_name,dihedral_key)
#            if value != None:
#                from_atom_index = atom.index()
#                
#                dihedral_indices = context.dihedral_indices
#
#                result = [from_atom_index]
#                result.extend(dihedral_indices)
#                result.append(value)
#                
#                for parameter_id in context._table.get_parameters():
#                    parameter = context._table.get_parameter(from_atom_name,dihedral_key,parameter_id)
#                    result.append(parameter)
#                
#                
#                result.append(context._table.get_exponent())
#                result = tuple(result)
#                
#        return result
#    
#    def _translate_atom_name(self, atom_name,context):
#        return context._table.get_translation(atom_name)

# TODO: it should be the default that all residues are accepted
class Non_bonded_backbone_component_factory(Ring_backbone_component_factory):
    def _residue_is_ok(self, segment, residue_number, table):
        return True

class Non_bonded_remote_component_factory(Atom_component_factory):
    

    def __init__(self):
            self._table_manager = Table_manager.get_default_table_manager()
            self._residue_types = self._table_manager.get_residue_types_for_table(NON_BONDED)

            self._all_spheres = self._get_all_spheres(self._table_manager, self._residue_types)
            self._chem_type_indexer = Non_bonded_potential.Chem_type_indexer(Table_manager.get_default_table_manager())

    def is_residue_acceptable(self, segment, residue_number, segment_manager):
        return True

    class Remote_non_bonded_context(Backbone_atom_indexer):

                
#        def _get_target_atom_id_coefficient_map(self, chem_type, sphere, table):
#            
#            non_bonded_type = table.get_chem_type_translation(chem_type)
##            print 'non bonded type', non_bonded_type
#            target_atom_coefficents = {}
#            for atom_name in table.get_target_atoms():
#                atom_index = self._get_atom_id(atom_name, table)
##                print non_bonded_type
#                coefficient = table.get_non_bonded_coefficient(atom_name, sphere, *non_bonded_type)
#                target_atom_coefficents[atom_index] = coefficient
#                
#            max_key = max(target_atom_coefficents.keys())
#            
#            result = [0.0]* (max_key+1)
#            
#            target_atom_ids = target_atom_coefficents.keys()
#            target_atom_ids.sort()
#            for target_atom_id in target_atom_ids:
#                result[target_atom_id] = target_atom_coefficents[target_atom_id]
#            return result
#
#        def _get_sphere_id(self,sphere,table):
#            sphere_id = table.get_spheres().index(sphere)
#            
#            return sphere_id
        
        def __init__(self, atom, spheres, table):
            
            chem_type_indexer = Non_bonded_potential.Chem_type_indexer(Table_manager.get_default_table_manager())
        #            self.exponent  =  table.get_exponent(sphere)
            
            raw_chem_type = Atom_utils._get_chem_type(atom)
            # TODO this is not efficient
            segid,residue_type, atom_name = Atom_utils._get_atom_info_from_index(atom.index())
            residue_type = Atom_utils._get_residue_type_from_atom_id(atom.index())
            raw_chem_type = table.get_chem_type_conversion(residue_type,atom_name,raw_chem_type)            
            chem_type = table.get_chem_type_translation(raw_chem_type)
            self.atom_index = atom.index()

#            key = bb_index, sphere, tuple(chem_type)
            
            self.chem_type_ids = []
            for sphere in spheres:
                chem_type_key = chem_type,sphere
                chem_type_id  = chem_type_indexer.get_index_for_key(chem_type_key)
                self.chem_type_ids.append(chem_type_id)
#            if raw_chem_type != chem_type:
#                print raw_chem_type, chem_type
            
#            self.target_atom_coefficents = self._get_target_atom_id_coefficient_map(raw_chem_type, sphere, table)
            
            
#            HEADER_SIZE = 3  # remote_atom_id  exponent
#            OFFSET_1 = 1
            
#            REMOTE_ATOM_ID = 0
#            SPHERE_ID = 1
#            EXPONENT = 2
            
#            component = [None] * (max_key+HEADER_SIZE+OFFSET_1)
#            component[REMOTE_ATOM_ID] = atom_id.index()
#            self.sphere_id = self._get_sphere_id(sphere,table)
            
            self.complete = True
#            component[SPHERE_ID] = 
#            component[EXPONENT] = exponent
            
#            results.append(component)
    
    def _get_all_spheres(self, table_manager, residue_types):
        seen_residue_types = set()
        spheres = []
        for residue_type in residue_types:
            table = table_manager.get_non_bonded_table(residue_type)
            real_residue_type = table._table['residue_type']
            
            
            if not real_residue_type in seen_residue_types:
                seen_residue_types.add(real_residue_type)
                spheres.extend(table.get_spheres())
        return tuple(set(spheres))
        
    def _build_contexts(self, atom, table):
        results = []

#        print atom.residueNum(),atom.atomName(),chem_type,table.is_non_bonded_chem_type(chem_type)

#        print self._all_spheres
#        for residue_type in self._residue_types:
#            bb_table = self._table_manager.get_non_bonded_table(residue_type)
#            bb_table_index  = bb_table._table['index']
#            for sphere in self._all_spheres:
#                if table.is_non_bonded_chem_type(chem_type):
        results.append(self.Remote_non_bonded_context(atom, self._all_spheres, table))
                
#                remote_atom_types = table.get_remote_atom_types(sphere)
        return results
    
    def _get_component_for_atom(self, atom, context):
        result = [atom.index()]
        result.extend(context.chem_type_ids)
#        result.extend(context.target_atom_coefficents)
        return tuple(result)
    #TODO is this still needed 
    def _translate_atom_name(self, atom_name, context):
        return atom_name
    
    def get_table_name(self):
        return 'NBRM'
    
    
class Non_bonded_list(object):
    def __init__(self,cutoff_distance= 5.0,jitter=0.2,update_frequency=5, min_residue_separation = 2):
        self._cutoff_distance = cutoff_distance
        self._jitter = jitter
        self._update_frequency =update_frequency
        
        self._box_update_count = update_frequency
        self._non_bonded = []
        self._non_bonded_calculation_count = 0
        self._non_bonded_call_count = 0

        self._min_residue_seperation = min_residue_separation
        self._pos_cache = {}

    def _get_cutoff_distance_2(self):
        return (self._cutoff_distance+self._jitter)**2
    

    def _filter_by_residue(self, seg_1, residue_1, seg_2, residue_2):
        result = False
        
        if seg_1 == seg_2:
            distance =abs(residue_1-residue_2)
            if distance < self._min_residue_seperation:
                result =True
        return result
    
    def update(self):
        self._box_update_count +=1
        
    def get_boxes(self,component_list_1, component_list_2,target_component_list,coefficient_list):
#        print self._box_update_count, self._update_frequency,self._box_update_count >= self._update_frequency
        self._non_bonded_call_count += 1
        
        if self._box_update_count >= self._update_frequency:
            target_component_list.clear()
            self._box_update_count = -1
            self._non_bonded_calculation_count += 1
            self._build_boxes(component_list_1, component_list_2, target_component_list,coefficient_list)
            
        return target_component_list
        
    
    def _get_cached_pos(self,atom_id):
        if atom_id in  self._pos_cache:
            result = self._pos_cache[atom_id]
        else:
            result = Atom_utils._get_atom_pos(atom_id)
            self._pos_cache = result
        return result
    
        
    def _is_non_bonded(self, atom_id_1, atom_id_2):
        
        seg_1, residue_1 = Atom_utils._get_atom_info_from_index(atom_id_1)[:2]
        seg_2, residue_2 = Atom_utils._get_atom_info_from_index(atom_id_2)[:2]


        cutoff_distance = self._cutoff_distance + self._jitter
        is_non_bonded = True
        if self._filter_by_residue(seg_1, residue_1, seg_2, residue_2):
            is_non_bonded = False
        else:
            is_non_bonded =  Atom_utils._calculate_distance(atom_id_1, atom_id_2) < cutoff_distance

        return is_non_bonded
    
    class NonBoolean:
        def __nonzero__(self):
            raise Exception("internal error this object should never be used in a boolean comparison!!")
    
    def _build_boxes(self, component_list_1, component_list_2, target_component_list, coefficient_list):
        
#        print
#        for elem in component_list_2:
#            print elem
#        print
        self._pos_cache = {}
        for component_1 in component_list_1:
            for component_2 in component_list_2: 
                
                atom_id_1, atom_1_coefficent_offset = component_1
                atom_id_2 = component_2[0]
                
                is_non_bonded = Non_bonded_list.NonBoolean()
                for i, chem_type_id in enumerate(component_2[1:]):
                    #TODO replace with a more direct lookup, we shouldn't need multiple values here
                    #so only list chem type ids by know chem_types in the sturtcure not all chem types
#                    print chem_type_id
#                    print coefficient_list
#                    print chem_type_id, sphere_id, exponent 
                    
                    
                    if i == 0:
                        is_non_bonded = self._is_non_bonded(atom_id_1, atom_id_2)
                        
                    if is_non_bonded:
                        coefficient_data  = coefficient_list.get_components_for_atom_id(chem_type_id)
    #                    print coefficient_data
                        coefficient_data = coefficient_data[0]
                        coefficients = coefficient_data[3:]
                        
                        chem_type_id, sphere_id, exponent = coefficient_data[:3]
#                    print chem_type_id, sphere_id, exponent 
                    
                        try:
                            result_component = atom_id_1,atom_id_2,coefficients[atom_1_coefficent_offset],exponent 
                        except:
                            print 'offset', atom_1_coefficent_offset
                        target_component_list.add_component(result_component)


class Null_component_factory(object):
   
    def __init__(self,name):
        self._name = name
    
    def is_residue_acceptable(self, segment, residue_number, segment_manager):
        False
    
    def create_components(self, component_list, table_source, segment,target_residue_number,selected_atoms):
        return []
        
    def get_table_name(self):
        return self._name



class Non_bonded_coefficient_factory(Atom_component_factory):
    
    def __init__(self):
        self._seen_spheres_and_chem_types =  set()
        self._chem_type_indexer = Non_bonded_potential.Chem_type_indexer(Table_manager.get_default_table_manager())
        self._table_manager = Table_manager.get_default_table_manager() 
        
    
    class Non_bonded_coefficient_context(object):
        
        def __init__(self, chem_type, sphere, non_bonded_tables):
            table_manager = Table_manager.get_default_table_manager()
            chem_type_indexer  = Non_bonded_potential.Chem_type_indexer(table_manager)
            sphere_indexer =  Non_bonded_potential.Sphere_indexer(table_manager)
            
            self.sphere_id =  sphere_indexer.get_index_for_key(sphere)
            
            chem_type_key = chem_type,sphere
            self.chem_type_id = chem_type_indexer.get_index_for_key(chem_type_key)

            self.exponent  =  non_bonded_tables[0].get_exponent(sphere)
            
            self.coefficients = []
            
            for table in non_bonded_tables:
                for target_atom in table.get_target_atoms():
                    self.coefficients.append(table.get_non_bonded_coefficient(target_atom,sphere,*chem_type))
            
#            print self.chem_type_id, self.exponent, self.coefficients
        
            self.complete =  True
            
            
    def is_residue_acceptable(self, segment, residue_number, segment_manager):
        return True
    

    def get_translated_chem_type(self, atom, table):
        raw_chem_type = Atom_utils._get_chem_type(atom)
        # TODO this is not efficient
        segid,residue_type, atom_name = Atom_utils._get_atom_info_from_index(atom.index())
        residue_type = Atom_utils._get_residue_type_from_atom_id(atom.index())
        raw_chem_type = table.get_chem_type_conversion(residue_type,atom_name,raw_chem_type)
        chem_type = table.get_chem_type_translation(raw_chem_type)
        return tuple(chem_type)

    def _build_contexts(self, atom, table):
        results = []
        residue_types = self._table_manager.get_all_known_residue_types()
        non_bonded_tables  = [self._table_manager.get_non_bonded_table(residue_type) for residue_type in residue_types]
        chem_type = self.get_translated_chem_type(atom, table)
        for sphere in table.get_spheres():
            
            key =  sphere, chem_type
            if not key in self._seen_spheres_and_chem_types:
                context = Non_bonded_coefficient_factory.Non_bonded_coefficient_context(chem_type, sphere, non_bonded_tables)
                results.append(context)
                self._seen_spheres_and_chem_types.add(key)
#            self.exponent  =  table.get_exponent(sphere)
            
        return results
        

    def _get_component_for_atom(self, atom, context):
        result = [context.chem_type_id, context.sphere_id, context.exponent]
        result.extend(context.coefficients)

        return tuple(result)
    
    def _translate_atom_name(self, atom_name, context):
        return atom_name
    
    def get_table_name(self):
        return 'COEF'





# target_atom_id, target_atom_type_id
# remote_atom_id  remote_atom_type_id 
# remote_atom_type_id exponent coefficient_by target_atom_id
class Non_bonded_potential(Distance_based_potential):

    def __init__(self,smoothed=True):
        super(Non_bonded_potential, self).__init__(smoothed=smoothed)
        
        self._add_component_factory(Non_bonded_backbone_component_factory())
        self._add_component_factory(Non_bonded_remote_component_factory())
        self._add_component_factory(Non_bonded_coefficient_factory())
        self._add_component_factory(Null_component_factory('NBLT'))
        
        self._non_bonded_list = Non_bonded_list()
    
    def _get_table_source(self):
        return Table_manager.get_default_table_manager().get_non_bonded_table
        
    def get_abbreviated_name(self):
        return NON_BONDED
    
    def _get_indices(self):
        return Distance_based_potential.Indices(target_atom_index=0, distance_atom_index_1=0, distance_atom_index_2=1, coefficent_index=2, exponent_index=3)
    
    def _get_target_atom_list_name(self):
        return 'NBLT'
    
    def _get_distance_list_name(self):
        return 'NBLT'

    def _get_non_bonded_list(self):
        non_bonded_list = self._get_component_list('NBLT')
        
        target_atom_list = self._get_component_list('ATOM')
        remote_atom_list = self._get_component_list('NBRM')
        coefficient_list  = self._get_component_list('COEF')
        self._non_bonded_list.get_boxes(target_atom_list, remote_atom_list, non_bonded_list, coefficient_list)
        
        return non_bonded_list
     
    def get_target_atom_ids(self):
        self.update_non_bonded_list()
        return super(Non_bonded_potential, self).get_target_atom_ids()

    def update_non_bonded_list(self):
        self._get_non_bonded_list()
        self._non_bonded_list.update()
        

    def calc_single_atom_shift(self, target_atom_id):
#        self.update_non_bonded_list()
        return Distance_based_potential.calc_single_atom_shift(self, target_atom_id)
    
    def calc_single_atom_force_set(self, target_atom_id, force_factor, forces):
#        self.update_non_bonded_list()
        return Distance_based_potential.calc_single_atom_force_set(self, target_atom_id, force_factor, forces)
    
    def _calc_single_force_set(self,index,factor, forces):
        target_atom_index,distant_atom_index = self._get_target_and_distant_atom_ids(index)
        
        result  = 0.0
        distance  = Atom_utils._calculate_distance(target_atom_index, distant_atom_index)
#        TODO: this should be the non bonded distance cutoff
        if distance < 5.0:
            result = Distance_based_potential._calc_single_force_set(self,index, factor, forces)

    
    def _calc_component_shift(self, index):
        target_atom_index,distant_atom_index = self._get_target_and_distant_atom_ids(index)
        
        result  = 0.0
        distance  = Atom_utils._calculate_distance(target_atom_index, distant_atom_index)
        if distance < 5.0:
            result = super(Non_bonded_potential, self)._calc_component_shift(index)
        return result
    
    class Base_indexer(object):
        
        __metaclass__ = ABCMeta
        
        def __init__(self,table_manager):
            self._index = {}
            self._inverted_index = {}
            self._max_index = 0
            
            self._build_index(table_manager)
            
        @abstractmethod
        def iter_keys(self,table):
            pass
        
        def _get_tables_by_index(self, table_manager):
            table_index_map = {}
            
            for residue_type in table_manager.get_residue_types_for_table(NON_BONDED):
                table = table_manager.get_non_bonded_table(residue_type)
                index  = table._table['index']
                table_index_map[index]=table
            
            keys = table_index_map.keys()
            keys.sort()
            result = []
            
            for key in keys:
                result.append(table_index_map[key])
            
            return tuple(result) 
                   
        def _build_index(self, table_manager):
            table_manager.get_non_bonded_table('base')
            tables_by_index = self._get_tables_by_index(table_manager)
            
            i = 0
            for table in tables_by_index:
                for key in self.iter_keys(table):
                    
#                    for chem_type in table.get_remote_atom_types(sphere):
                    if not key in self._index:
                        self._index[key] = i
                        self._inverted_index[i] = key
                        self._max_index = i
                        i += 1
                        

        def get_index_for_key(self,key):
            return self._index[key]
        
        def get_key_for_index(self,index):
            return self._inverted_index[index]
        
        def get_max_index(self):
            return self._max_index
            

        def _flatten(self, lst):
            for el in lst:
                if hasattr(el, '__iter__') and not isinstance(el, basestring):
                    for x in self._flatten(el):
                        yield x
                else:
                    yield el


        def get_key_str(self, index):
            key = self.get_key_for_index(index)
            flat_key = tuple([elem for elem in self._flatten(key)])
            format_string = '[%s, %-4s, %8s]'
            key_str = format_string % flat_key
            return key_str

        def __str__(self):
            result  = []
            
            result.append('non bonded chem_type index (%i entries)' % (self.get_max_index() +1))
            result.append('')
            keys = self._inverted_index.keys()
            keys.sort()
            for index in keys:
                key_str = self.get_key_str(index)
                index_str = '%3i.' %( index+1) 
                result.append(index_str + ' '  + key_str)
            return '\n'.join(result)

    
    class Sphere_indexer(Base_indexer):
        def __init__(self,table_manager):
            super(Non_bonded_potential.Sphere_indexer, self).__init__(table_manager)
            
        def iter_keys(self,table):
            for sphere in table.get_spheres():
                yield sphere           
        
        
    class Chem_type_indexer(Base_indexer):

        def __init__(self,table_manager):
            super(Non_bonded_potential.Chem_type_indexer, self).__init__(table_manager)
            
            
          
        def iter_keys(self, table):
            i = 0
            for sphere in table.get_spheres():
                for chem_type in table.get_remote_atom_types(sphere):  
                    i+= 1
                    yield chem_type,sphere

                         
            
                
        
    #TODO centralise table manager (each sub potential has its own table manager at the moment)
    #TODO complete
    def __str__(self):
        non_bonded_indexer = Non_bonded_potential.Chem_type_indexer(Table_manager.get_default_table_manager())
        atom_list = self._get_component_list('ATOM')
        
        indexer  = Backbone_atom_indexer()
        table_manager = self._table_manager
        residue_types = table_manager.get_all_known_residue_types()   
        
        offset_names = {}
        
        for residue_type in residue_types:
            table =  table_manager.get_non_bonded_table(residue_type)
            indexer._get_offset_names(table, offset_names)                             

        result = []
        result.append('Non bonded table')
        result.append('')
        result.append('residue types: %s' %  ', '.join(residue_types))
        result.append('')
        result.append(non_bonded_indexer.__str__())
        result.append('')
        result.append('type table:')
        result.append('')

        offset_items = offset_names.items()
        offset_items.sort()
        
        result.append('\n'.join(['%2i. %-4s %-2s' % ((i+1,) + elem) for i,elem in offset_items]))
        result.append('')
        
        result.append('')
        result.append('atom list (%i entries)' % len(atom_list) )
        result.append('')
        
        for i, atom_elem in enumerate(atom_list):
            atom_id, atom_type = atom_elem
            
            atom_sel = Atom_utils._get_atom_name(atom_id)
            i_elem =  (`i+1`+'.',) + atom_elem + (atom_sel,)
            
            all_elem=i_elem + offset_names[atom_type]
            result.append('%-5s [%-4i %2i] : %s - [%-4s,%-2s]' % all_elem)
        
        non_bonded_remote_list = self._get_component_list('NBRM')
        result.append('')
        result.append('non bonded remote list (%i entries)' % len(non_bonded_remote_list) )
        for remote_elem in non_bonded_remote_list:
            atom_id = remote_elem[0]
            atom_sel = Atom_utils._get_atom_name(atom_id)
            sub_result_1 = '%i (%s) ' % (atom_id,atom_sel)
            sub_result_2 = []
            sub_result_3 = []
            for chem_type_index in remote_elem[1:]:
                sub_result_2.append(non_bonded_indexer.get_key_str(chem_type_index))
                sub_result_3.append('%2i'  % chem_type_index)
                
            result.append('%s %s %s'  % (sub_result_1,' '.join(sub_result_3),  '.'.join(sub_result_2)))
            

        non_bonded_coefficient_list  = self._get_component_list('COEF')
        result.append('')
        result.append('coefficient table (%i entries)' % len(non_bonded_coefficient_list))
        result.append('')
    
        for i, elem in enumerate(non_bonded_coefficient_list):
            chem_type_id,sphere_id,exponent = elem[:3]
            coefficients  = list(elem[3:])
            for j,coefficient in enumerate(coefficients):
                if coefficient == None:
                    coefficients[j] = float('nan')
            num_atom_types = len(offset_names)
#           
            pre_string = '%-3i %-3i - %- 8s %+2.1f - [%i]' % (i, chem_type_id, non_bonded_indexer.get_key_str(chem_type_id), exponent, sphere_id)
            clear_prestring = ' ' * len(pre_string)
            
            if i == 0:
                headers_format = '%12s ' * (num_atom_types /  len(residue_types))
                n_headers = num_atom_types / len(residue_types)
                headers_data = []
                for i in range(n_headers):
                    headers_data.append(offset_names[i][1])
                header = headers_format % tuple(headers_data)
                result.append(clear_prestring + '      ' + header)
            for i,coefficient_offset in enumerate(range(0,num_atom_types, (num_atom_types /  len(residue_types)))):
                coeff_format = '% +12.7g ' * (num_atom_types / len(residue_types))
                coefficients_set =  tuple(coefficients)[coefficient_offset:coefficient_offset+(num_atom_types / len(residue_types))]
                coefficient_string =  coeff_format % coefficients_set
                residue_types_offset = coefficient_offset / (num_atom_types / len(residue_types))
#                print residue_types_offset, coefficient_offset
                if i == 0:
                    line = pre_string + (' %4s ' % residue_types[residue_types_offset]) + coefficient_string
                else:
                    line = clear_prestring + (' %4s ' % residue_types[residue_types_offset]) + coefficient_string
                line = line.replace('+',' ')
                line = line.replace('nan','   ')
                result.append(line)
        return '\n'.join(result)

class Xcamshift():
    def __init__(self):
        self.potential = [
                          RandomCoilShifts(),
                          Distance_potential(),
                          Extra_potential(),
                          Dihedral_potential(),
                          Sidechain_potential(),
                          Ring_Potential(),
                          Non_bonded_potential()
                          ]
        self._shift_table = Observed_shift_table()
        self._shift_cache = {}
    
    def clear_shift_cache(self):
        self._shift_cache = {}
        
    def get_sub_potential_names(self):
        return [potential.get_abbreviated_name() for potential in self.potential]
    
    def get_named_sub_potential(self,name):
        result =  None
        for potential in self.potential:
            if potential.get_abbreviated_name() ==  name:
                result = potential
                break
        if potential ==  None:
            template = "cannot find potential with the name %s in %s"
            all_potential_names  =  ",".join(self.get_sub_potential_names())
            msg = template % (name,all_potential_names)
            raise Exception(msg)
        return result
#    
    def print_shifts(self):
        result  = [0] * Segment_Manager().get_number_atoms()
        
        result_elements = {}
#        keys=[]
#        for potential in self.potential:
#            num_atoms = len(result)
#            sub_result  = [0.0] * num_atoms
#            potential.set_shifts(sub_result)
#            key = potential.get_abbreviated_name()
#            keys.append(key)
#            result_elements[key] = sub_result
#        
        
        keys = []
        for potential in self.potential:
            num_atoms = len(result)
            sub_result  = [0.0] * num_atoms
            potential.set_shifts(sub_result)
            key = potential.get_abbreviated_name()
            keys.append(key)
            result_elements[key] = sub_result
        
        
        total = [sum(elems) for elems in zip(*result_elements.values())]
        keys.append('TOTL')
        
        result_elements['TOTL'] = total
        residues  = []
        atoms = []
        result_elements['ATOM'] =  atoms
        result_elements['RESD'] =  residues
        
        keys.insert(0, 'RESD')
        keys.insert(0,'ATOM')
        
        for i in range(num_atoms):
            segments,residue,atom = Atom_utils._get_atom_info_from_index(i)
            residues.append(('%5i' %residue))
            atoms.append(atom)
        
        for key in keys:
            values = []
            print key.ljust(5),
            for total,value in zip(result_elements['TOTL'],result_elements[key]):
                if total > 0:
                    if isinstance(value, float):
                        string  = '%- 7.4f' % value
                    else:
                        string = value
                    string = string.rjust(9)
                    
                    values.append(string)
            print ' '.join(values)
            
        
    def _get_target_atom_ids(self):
        result  = set()
        for potential in self.potential:
            for target_atom_id in potential.get_target_atom_ids():
                result.add(target_atom_id)
        result = list(result)
        result.sort()
        return result
    
    
    def set_shifts(self, result):
        target_atom_ids =  self._get_target_atom_ids()
        for target_atom_id in target_atom_ids:
            shift = self.calc_single_atom_shift(target_atom_id)
            result[target_atom_id] = shift
        
        return result
    
    def set_observed_shifts(self, shift_table):
        self._shift_table  =  shift_table
        
    def calc_single_atom_shift(self,atom_index):
        result  = 0.0
        if atom_index in self._shift_cache:
            result = self._shift_cache[atom_index]
        else:
            for potential in self.potential:
                result += potential.calc_single_atom_shift(atom_index)
        return result
    
    
    
    def _get_flat_bottom_shift_limit(self, residue_type, atom_name):
        table_manager = Table_manager.get_default_table_manager()
        constants_table = table_manager.get_constants_table(residue_type)
        
        atom_name = constants_table.get_translation_to_table(residue_type,atom_name)
        flat_bottom_limit = constants_table.get_flat_bottom_limit(atom_name)
        #TODO:  move to class body this is a global scaling for the 
        # whole force field
        flat_bottom_constant  = constants_table.get_flat_bottom_constant()
        return flat_bottom_limit * flat_bottom_constant
    

    def _adjust_shift(self, shift_diff, flat_bottom_shift_limit):
        result  = 0.0
        if (shift_diff > 0.0):
            result = shift_diff-flat_bottom_shift_limit
        else:
            result = shift_diff + flat_bottom_shift_limit
        return result
    
#    TODO: investigate need for atom name and residue_type for these lookups do they need caching?
#          how does cacmshift do it 
    def _get_end_harmonic(self, residue_type, atom_name):
        table_manager = Table_manager.get_default_table_manager()
        constants_table = table_manager.get_constants_table(residue_type)
        
        atom_name = constants_table.get_translation_to_table(residue_type,atom_name)
        
        return constants_table.get_end_harmonic(atom_name)
    
    
    def _get_scale_harmonic(self, residue_type, atom_name):
        table_manager = Table_manager.get_default_table_manager()
        constants_table = table_manager.get_constants_table(residue_type)
        
        atom_name = constants_table.get_translation_to_table(residue_type,atom_name)
                
        return constants_table.get_scale_harmonic(atom_name)
    
    
    def _get_tanh_y_offset(self, residue_type, atom_name):
        table_manager = Table_manager.get_default_table_manager()
        constants_table = table_manager.get_constants_table(residue_type)
        
        atom_name = constants_table.get_translation_to_table(residue_type,atom_name)
                
        return constants_table.get_tanh_y_offset(atom_name)
    

        

    def get_shift_difference(self, target_atom_index):
        
        theory_shift = self.calc_single_atom_shift(target_atom_index)
        
        observed_shift = self._shift_table.get_chemical_shift(target_atom_index)
        
        return observed_shift - theory_shift
        

    def _calc_single_atom_energy(self, target_atom_index):
        
        energy = 0.0
        if target_atom_index in self._shift_table.get_atom_indices():
            
            shift_diff = self.get_shift_difference(target_atom_index)
            
            residue_type = Atom_utils._get_residue_type_from_atom_id(target_atom_index)
            atom_name = Atom_utils._get_atom_name_from_index(target_atom_index)
            
            flat_bottom_shift_limit = self._get_flat_bottom_shift_limit(residue_type, atom_name)
            
            if abs(shift_diff) > flat_bottom_shift_limit:
                adjusted_shift_diff = self._adjust_shift(shift_diff, flat_bottom_shift_limit)
                
                end_harmonic = self._get_end_harmonic(residue_type, atom_name)
                scale_harmonic = self._get_scale_harmonic(residue_type, atom_name)
                
                
                energy_component = 0.0
                if adjusted_shift_diff < end_harmonic:
                    energy_component = (adjusted_shift_diff/scale_harmonic)**2
                else:
                    tanh_amplitude = self._get_tanh_amplitude(residue_type,atom_name)
                    tanh_elongation = self._get_tanh_elongation(residue_type, atom_name)
                    tanh_y_offset = self._get_tanh_y_offset(residue_type, atom_name)
                    
                    tanh_argument = tanh_elongation * (adjusted_shift_diff - end_harmonic)
                    energy_component = tanh_amplitude * tanh(tanh_argument) + tanh_y_offset;

                energy += energy_component
        return energy
    
        
    def _calc_single_atom_force_set_with_potentials(self, target_atom_id, forces, potentials_list):
#        print 'forces _calc_single_atom_force_set_with_potentials', forces 
#        self._calc_single_force_factor(target_atom_id, forces)
        for potential in potentials_list:
            factor  = self._calc_single_factor(target_atom_id)
            potential.calc_single_atom_force_set(target_atom_id,factor,forces)
        return forces

    def _calc_single_atom_force_set(self,target_atom_id,forces,potentials=None):
#        print 'forces xcamshift._calc_single_atom_force_set',forces
        if potentials ==  None:
            potentials =  self.potential
        self._calc_single_atom_force_set_with_potentials(target_atom_id, forces, potentials)
        
        return forces

    def _get_weight(self, residue_type, atom_name):
        table_manager = Table_manager.get_default_table_manager()
        constants_table = table_manager.get_constants_table(residue_type)
        
        atom_name = constants_table.get_translation_to_table(residue_type,atom_name)
        
        return constants_table.get_weight(atom_name)
    
    #TODO: maybe call this a scaling??

    
    def _get_tanh_amplitude(self, residue_type, atom_name):
        table_manager = Table_manager.get_default_table_manager()
        constants_table = table_manager.get_constants_table(residue_type)
        
        atom_name = constants_table.get_translation_to_table(residue_type,atom_name)
        
        return constants_table.get_tanh_amplitude(atom_name)
    
    
    def _get_tanh_elongation(self, residue_type, atom_name):
        table_manager = Table_manager.get_default_table_manager()
        constants_table = table_manager.get_constants_table(residue_type)
        
        atom_name = constants_table.get_translation_to_table(residue_type,atom_name)
        
        return constants_table.get_tanh_elongation(atom_name)
    
    
    def _calc_single_factor(self, target_atom_id):
        
        factor = 0.0
        if target_atom_id in self._shift_table.get_atom_indices():
            
            theory_shift = self.calc_single_atom_shift(target_atom_id)
            observed_shift = self._shift_table.get_chemical_shift(target_atom_id)
            
            shift_diff = observed_shift - theory_shift
            
            residue_type = Atom_utils._get_residue_type_from_atom_id(target_atom_id)
            atom_name = Atom_utils._get_atom_name_from_index(target_atom_id)
            
            flat_bottom_shift_limit = self._get_flat_bottom_shift_limit(residue_type, atom_name)
            
            if abs(shift_diff) > flat_bottom_shift_limit:
                adjusted_shift_diff = self._adjust_shift(shift_diff, flat_bottom_shift_limit)
                end_harmonic = self._get_end_harmonic(residue_type, atom_name)
                scale_harmonic = self._get_scale_harmonic(residue_type, atom_name)
                sqr_scale_harmonic = scale_harmonic**2
                
                weight = self._get_weight(residue_type,atom_name)
                
                tanh_amplitude = self._get_tanh_amplitude(residue_type,atom_name)
                tanh_elongation = self._get_tanh_elongation(residue_type,atom_name)
                
                # TODO: add factor and lambda to give fact
                fact =1.0
                if adjusted_shift_diff < end_harmonic:
                    factor = 2.0 * weight * adjusted_shift_diff * fact / sqr_scale_harmonic;
                else:
                    factor = weight * tanh_amplitude * tanh_elongation / (cosh(tanh_elongation * (adjusted_shift_diff - end_harmonic)))**2.0 * fact;
                
        else:
            msg = "requested factor for target [%s] which is not in shift table"
            target_atom_info = Atom_utils._get_atom_info_from_index(target_atom_id)
            raise Exception(msg % target_atom_info)
        return factor

    def _calc_single_force_factor(self,target_atom_index,forces):
        factor = 1.0
        for potential in self.potential:
            if potential._have_derivative():
                potential.set_observed_shifts(self._shift_table)
                
                if target_atom_index in self._shift_table.get_atom_indices():
                    factor = self._calc_single_factor(target_atom_index)
#                    potential._calc_single_factor(target_atom_index,factor,forces)
        return factor
    

    def _get_active_target_atom_ids(self):
        target_atom_ids = set(self._get_target_atom_ids())
        observed_shift_atom_ids = self._shift_table.get_atom_indices()
        active_target_atom_ids = target_atom_ids.intersection(observed_shift_atom_ids)
        active_target_atom_ids = list(active_target_atom_ids)
        return active_target_atom_ids
    
    def prepare(self):
        for potential in self.potential:
            potential.prepare()
            
    def calcEnergy(self, prepare =  True):
        if prepare:
            self.prepare()
            
        active_target_atom_ids = self._get_active_target_atom_ids()
        
        energy = 0.0
        for target_atom_id in active_target_atom_ids:
            energy+= self._calc_single_atom_energy(target_atom_id)
        
        return energy
    
    

    def _calc_derivs(self, derivs,potentials=None):
        
        active_target_atom_ids = self._get_active_target_atom_ids()
        
        for target_atom_id in active_target_atom_ids:
            self._calc_single_atom_force_set(target_atom_id, derivs, potentials)

    def calcEnergyAndDerivs(self,derivs):
        self.prepare()
        energy = self.calcEnergy(prepare=False)
        
        self._calc_derivs(derivs)
        
        return energy
        

