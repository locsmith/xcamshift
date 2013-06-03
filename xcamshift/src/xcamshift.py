#-------------------------------------------------------------------------------
# Copyright (c) 2013 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 27 Dec 2011

@author: garyt
#TODO: need to translate from_atom names
'''

#TODO: add tests to exclude atoms/distances which are not defined


from atomSel import AtomSel,intersection
from component_list import Component_list, Native_component_list
from dihedral import Dihedral
from keys import Atom_key, Dihedral_key
from observed_chemical_shifts import Observed_shift_table
from cython.fast_segment_manager import Segment_Manager
from table_manager import Table_manager
from python_utils import tupleit
from utils import Atom_utils, iter_residues_and_segments
import sys
from common_constants import  BACK_BONE, XTRA, RANDOM_COIL, DIHEDRAL, SIDE_CHAIN, RING, NON_BONDED, DISULPHIDE
from common_constants import  TARGET_ATOM_IDS_CHANGED, ROUND_CHANGED, STRUCTURE_CHANGED, SHIFT_DATA_CHANGED
from abc import abstractmethod, ABCMeta
from cython.shift_calculators import Fast_distance_shift_calculator, Fast_dihedral_shift_calculator, \
                                     Fast_ring_shift_calculator, Fast_ring_data_calculator,          \
                                     Fast_non_bonded_calculator, Fast_energy_calculator,             \
                                     Fast_distance_based_potential_force_calculator,                 \
                                     Fast_dihedral_force_calculator,                                 \
                                     Fast_ring_force_calculator,                                     \
                                     Fast_force_factor_calculator,                                   \
                                     Out_array, Vec3_list, allocate_array, zero_array,               \
                                     Fast_non_bonded_shift_calculator,                               \
                                     Fast_non_bonded_force_calculator,                               \
                                     Non_bonded_interaction_list,                                    \
                                     Fast_random_coil_shift_calculator
from time import time
import array#

#TODO: REMOVE!
class Index_translator(object):
    def __init__(self, indices):
        self._indices=indices
        
    def __call__(self, component):
        result = []
        result.append(component[self._indices.target_atom_index])
        if len(component) == 4:
            result.append(component[self._indices.distance_atom_index_1])
            result.append(component[self._indices.distance_atom_index_2])
            result.append(component[self._indices.coefficient_index])
            result.append(component[self._indices.exponent_index])
        elif len(component) == 5:
            result.append(component[self._indices.distance_atom_index_1])
            result.append(component[self._indices.distance_atom_index_2])
            result.append(component[self._indices.coefficient_index])
            result.append(component[self._indices.exponent_index])
        else:
            raise Exception("bad distance component length %i should be either 4 or 5 " % len(component))
        
        return result
    
class Component_factory(object):
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def is_residue_acceptable(self, segment, residue_number, segment_manager):
        pass
    
    @abstractmethod
    def create_components(self, component_list, table_source, segment,target_residue_number,selected_atoms):
        pass
        
    @abstractmethod
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
    
    @abstractmethod
    def create_residue_components(self,component_list,table, segment, residue):
        pass
    
class Atom_component_factory(Component_factory):
    __metaclass__ = ABCMeta
    
    def is_residue_acceptable(self, segment, residue_number, segment_manager):
        segment_info = segment_manager.get_segment_info(segment)
        
        return residue_number > segment_info.first_residue and residue_number < segment_info.last_residue




    def create_components(self, component_list, table_provider, segment,target_residue_number,selected_atoms):
        table = self._get_table_for_residue(segment, target_residue_number, table_provider)
        self.create_atom_components(component_list, table, selected_atoms)
    
    @abstractmethod
    def _build_contexts(self,atom, table):
        pass
    
    @abstractmethod
    def _get_component_for_atom(self,atom, context):
        pass
    
    @abstractmethod
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

class Disulphide_context :
    def __init__(self,atom,table):
        self.complete = False
        
        self._table = table
        
        segment = atom.segmentName()
        
        from_residue_number = atom.residueNum()
        from_residue_type = Atom_utils._get_residue_type(segment, from_residue_number )
        if from_residue_type ==  'CYS':
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
        for (distance_atom_1,offset_1), (distance_atom_2,offset_2)in table.get_offsets_and_target_atoms():
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
        
        if len(sidechain_atom) > 0:
            self.sidechain_atom_index = sidechain_atom[0].index()
            if self.sidechain_atom_index != from_atom_index:
            
                self.residue_type = residue_type
                self.sidechain_atom_name = sidechain_atom[0].atomName()
                
                self.complete = True
        else:
            message = "NOTICE: couldn't find the atom %s:%i@%s"
            data = segment, sidechain_residue_number, sidechain_atom_selection
            print >> sys.stderr, message % data

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
    
    __metaclass__ = ABCMeta
    
    ALL = '(all)'
            
    def __init__(self):
        self._segment_manager = Segment_Manager.get_segment_manager()
        self._table_manager = Table_manager.get_default_table_manager()
        self._observed_shifts = Observed_shift_table()
        self._component_list_data  = {}
        self._component_factories = {}
        self._cache_list_data = {}
        self._freeze  = False
    
        #TODO: this can go in the end... we just need some more clever logic in the get
        self._shift_calculator = None
        
        
    
        self._component_to_result = None
        self._active_components = None
        self._verbose = False
        
    def set_frozen(self,on):
        self._freeze =  (on == True)
        
    def set_verbose(self, on=True):
        self._verbose = on
        if self._shift_calculator != None:
            self._shift_calculator.set_verbose(on)
            
    def _get_force_calculator(self):
        raise Exception("Not used!")
        return Base_force_calculator(self)



    def _build_active_components_list(self, target_atom_ids, components=None, force=False):
        if force or not self._freeze  or  self._active_components == None:
            if self._verbose:
                print '   filtering components %s' % self.get_abbreviated_name(),
            
            if target_atom_ids == None:
                test = lambda x:True
            else:
                target_atom_ids = sorted(set(target_atom_ids))
                test = lambda x:x[0] in target_atom_ids
            if components == None: 
                components =  self._get_component_list() 
                 
            result  =  components.build_selection_list(test)

            if self._verbose:
                print ' %i reduced to %i' % (len(components),len(self._active_components))
        else:
            result = self._active_components
        
        return result
    
    def _get_components(self):
        return self._get_component_list().get_native_components()
    
    def _calc_component_shift(self, index):
        
        component_to_result = array.array('i', [0])
        results = array.array('d',[0.0])
        components = self._get_components()
        active_components = array.array('i', [index])
        self._shift_calculator(components,results,component_to_result, active_components=active_components)
        
        return results[0]
      
    

    def _build_component_to_result(self, active_components, target_atom_ids, components):
        result = allocate_array(len(active_components), 'i')
        
        for i, component_index in enumerate(active_components):
            out_atom_id = components[component_index][0]
            if out_atom_id in target_atom_ids:
                out_index = target_atom_ids.index(out_atom_id)
            else:
                out_index =  -1
            result[i] = out_index
        
        return result

    def _update_component_to_result(self, target_atom_ids):
        if not self._freeze:
            self._component_to_result = self._build_component_to_result(self._active_components, target_atom_ids, self._get_component_list())
         

    def _prepare(self,change, data):
        if change == STRUCTURE_CHANGED:
            self.clear_caches()
            
        if change == TARGET_ATOM_IDS_CHANGED:
            self._active_components = self._build_active_components_list(data)
            self._update_component_to_result(data)
        

    
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
    
    
    @abstractmethod
    def _get_table_source(self):
        pass
        
    @abstractmethod
    def get_abbreviated_name(self):
        pass

    def set_observed_shifts(self, shift_table):
        self._observed_shifts = shift_table
        
    def _get_cache_list(self,name):
        if not name in self._cache_list_data:
            if name in ('NORM', 'CENT'):
                self._cache_list_data[name] = Vec3_list()
            else:
                self._cache_list_data[name] = Component_list()
        return self._cache_list_data[name]
    
    def clear_caches(self):
        self._cache_list_data = {}
        self._active_components = None

#    TODO put 'ATOM' in a constant and rename to BB_ATOM?
    def _get_component_list(self,name=None):
        if name == None:
            name  = self._get_target_atom_list_name()
        if not name in self._component_list_data:
            self._component_list_data[name] = self._create_component_list(name)
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
        
    def _calc_single_atom_shift(self, target_atom_id):
        components =  self._get_component_list()

        result = 0.0
        if target_atom_id in components.get_component_atom_ids():
            index_range = components.get_component_range(target_atom_id)
            for index in range(*index_range):
                result += self._calc_component_shift(index)
        
        return result
    
    def calc_force_set(self,target_atom_ids,force_factors,forces):
        if self._have_derivative():
            components = self._get_components()
            self._force_calculator(components, self._component_to_result, force_factors, forces, active_components=self._get_active_components())
            
    
    def calc_single_atom_force_set(self,target_atom_id,force_factor,forces):
        target_atom_ids = array.array('i',[target_atom_id])
        force_factors = [force_factor]
        self.calc_force_set(target_atom_ids,force_factor,forces)
    
    #TODO: make this just return a list in component order
    #TODO: remove
    def set_shifts(self, result):
        components =  self._get_component_list()
        
        for target_atom_id in components.get_component_atom_ids():
            result[target_atom_id] +=  self._calc_single_atom_shift(target_atom_id)
            
    def _have_derivative(self):
        return True
    
#    TODO remove or make batched now in force calculaor
    def _get_or_make_target_force_triplet(self, forces, target_offset):
        target_forces = forces[target_offset]
        if target_forces == None:
            target_forces = [0.0] * 3
            forces[target_offset] = target_forces
        return target_forces

    def _create_component_list(self,name):
        return Component_list()
    
    def _get_active_components(self):
        return self._active_components
    
    #TODO: unify with ring random coil and disuphide shift calculators
    def calc_shifts(self, target_atom_ids, results):
       
        components = self._get_components()
        self._shift_calculator(components,results,self._component_to_result, active_components=self._get_active_components())
             

            

            


class Distance_based_potential(Base_potential):
    
    def __init__(self,  smoothed = False, fast=False):
        super(Distance_based_potential, self).__init__()
        self._smoothed = smoothed
        #TODO: sort placement out, move to non bonded
        #
        self._cutoff = 5.0
        
        #TODO: move to base potential
        self._shift_calculator = self._get_shift_calculator()
        self._force_calculator = self._get_force_calculator()
        
        
    #TODO: move to base potential
    def _get_shift_calculator(self):
#        if self._fast:
        result  = Fast_distance_shift_calculator(self._smoothed, name=self.get_abbreviated_name())
        result.set_verbose(self._verbose)
#       else:
#        result = self._shift_calculator = Distance_shift_calculator(self._get_indices(), self._smoothed)
        return result
    
    def _get_force_calculator(self):
#        if self._fast:
        result = Fast_distance_based_potential_force_calculator(self._smoothed, name=self.get_abbreviated_name())
#        else:
#        result = Distance_based_potential_force_calculator(self._get_indices(), self._smoothed)
        result.set_verbose(self._verbose)
        return result 
                
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
            
    @abstractmethod
    def _get_indices(self):
        return Distance_based_potential.Indices(target_atom_index=0,distance_atom_index_1=0,
                                                distance_atom_index_2=1,coefficent_index=2,
                                                exponent_index=3)
    
#    @abstractmethod
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

    #TODO move this to su classes when required   
    def _create_component_list(self, name):

            
        if name == 'ATOM':
            translator = Index_translator(self._get_indices())
            return Native_component_list(translator=translator,format='iiiff')
        else:
            print 'note non native component list for %s %s' % (self.get_abbreviated_name(), name)
            return Component_list()
    
 
    
class Distance_potential(Distance_based_potential):
    '''
    classdocs
    '''


    def __init__(self, fast=False):
        super(Distance_potential, self).__init__(smoothed = False, fast=fast)
        
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
        self._shift_calculator = self._get_shift_calculator()

    
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
    
    def _create_component_list(self, name):
        if name == 'ATOM':
            return Native_component_list(format='if')
    
    def _get_shift_calculator(self):
        result = Fast_random_coil_shift_calculator(name=self.get_abbreviated_name())
        result.set_verbose(self._verbose)
        
        return result
        
class Disulphide_shift_calculator(Base_potential):
    

    def __init__(self):
        super(Disulphide_shift_calculator, self).__init__()
        
        self._add_component_factory(Disulphide_shift_component_factory())

    
    def get_abbreviated_name(self):
        return DISULPHIDE
    
    def _have_derivative(self):
        False
         
    
    def _get_table_source(self):
        return self._table_manager.get_disulphide_table
    
    #TODO: this is the same as for the random coild shifts
    def calc_shifts(self, target_atom_ids, result):
        components = self._get_component_list()
        for i, component_index in enumerate(self._active_components):
            result[self._component_to_result[i]] += components[component_index][1]
            
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
    
        self._shift_calculator = self._get_shift_calculator()
        self._force_calculator = self._get_force_calculator()
    
#    TODO: remove only used for compatability during upgrades
    def calc_single_atom_force_set(self,target_atom_id,force_factor,forces):
        if self._have_derivative():
            components =  self._get_component_list()
        
        
            if target_atom_id in components.get_component_atom_ids():
                index_range = components.get_component_range(target_atom_id)
                for index in range(*index_range):
                    self._calc_single_force_set(index,force_factor,forces)
        return forces
        
        
#     def calc_shifts(self, target_atom_ids, results):cal
#         components  = self._filter_components(target_atom_ids)
#         self._shift_calculator(components,results,self._component_to_result)
        
    def get_abbreviated_name(self):
        return DIHEDRAL

    def _get_shift_calculator(self):
#        if self._fast:
        result = Fast_dihedral_shift_calculator(name=self.get_abbreviated_name())
        result.set_verbose(self._verbose)
#        else:
#            result = Dihedral_shift_calculator()
        return result
    
    def _get_table_source(self):
        return self._table_manager.get_dihedral_table

    def _get_distance_components(self):
        return self._get_component_list()
    

        



    
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
    
    def _get_force_calculator(self):
#        if self._fast:
        result = Fast_dihedral_force_calculator(name=self.get_abbreviated_name())
#        else:
#        result  = Dihedral_force_calculator()
        result.set_verbose(self._verbose)
        return result

    def _create_component_list(self, name):
        if name == 'ATOM':
            return Native_component_list(format='iiiiifffffff')
        
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

class Disulphide_shift_component_factory(Atom_component_factory, Ring_factory_base):
    
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
            value = context._table.get_disulphide_shift(atom_name)
            if value != None:
                result = atom.index(), value
        return result
    

    def _build_contexts(self, atom, table):
        contexts = [Disulphide_context(atom, table)]
        return contexts

    def get_table_name(self):
        return 'ATOM'

    
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
            
    @abstractmethod
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

    
    RING_ATOM_IDS = 1
    
    def __init__(self):
        super(Ring_Potential, self).__init__()
        
        self._add_component_factory(Ring_backbone_component_factory())
        self._add_component_factory(Ring_coefficient_component_factory())
        self._add_component_factory(Ring_sidechain_atom_factory())
        self._shift_calculator = self._get_shift_calculator()
        self._ring_data_calculator = self._get_ring_data_calculator()
        self._force_calculator = self._get_force_calculator()
      
    def set_verbose(self,on):
        super(Ring_Potential, self).set_verbose(on)
        self._ring_data_calculator.set_verbose(on)

        
    def _get_force_calculator(self):
#        if self._fast:
        result = Fast_ring_force_calculator(name=self.get_abbreviated_name())
#        else:
#            result =  Ring_force_calculator()
        result.set_verbose(self._verbose)
        return result
            
    def _prepare(self, change, target_atom_ids): 
        super(Ring_Potential, self)._prepare(change, target_atom_ids)
        
#         if change == STRUCTURE_CHANGED:
#             self._clear_ring_cache()
            
        if change  == ROUND_CHANGED or change == TARGET_ATOM_IDS_CHANGED:
            if self._verbose:
                print "update ring cache"
                
            self._build_ring_data_cache()
            self._setup_ring_calculator(self._shift_calculator)
            self._setup_ring_calculator(self._force_calculator)

         

    def _setup_ring_calculator(self,calculator):
        calculator._set_coef_components(self._get_component_list('COEF').get_native_components(), self._get_component_list('COEF').get_native_component_offsets())
        calculator._set_ring_components(self._get_component_list('RING').get_native_components())
        calculator._set_normal_cache(self._get_cache_list('NORM'))
        calculator._set_centre_cache(self._get_cache_list('CENT'))
    
    def _get_ring_data_calculator(self):
#        if self._fast:
        result = Fast_ring_data_calculator()
#        else:
#            result = Ring_data_calculator()
        result.set_verbose(self._verbose)
        return result
    
        
        
    def get_abbreviated_name(self):
        return RING

    def _get_shift_calculator(self):
#        if self._fast:
        result = Fast_ring_shift_calculator(name=self.get_abbreviated_name())
        result.set_verbose(self._verbose)
#        else:
#            result = Ring_shift_calculator()
        return result
    
        
    def _get_table_source(self):
        return self._table_manager.get_ring_table
    
    def _get_distance_components(self):
        return self._get_component_list('ATOM')




    



    
    def _calc_component_shift(self, index):
        components = self._create_component_list('ATOM')
        component = self._get_distance_components()[index]
        components.add_component(component)
        results = array.array('d',[0.0])
        self._setup_ring_calculator(self._shift_calculator)
        self._setup_ring_calculator(self._shift_calculator)
        component_to_result  = array.array('i',[0])
        self._shift_calculator(components.get_native_components(),results,component_to_result, active_components=array.array('i',[0]))
        return results[0]
    
    

#     def _clear_ring_cache(self):
#         normals = self._get_cache_list('NORM')
# #         normals.clear()
#         centres = self._get_cache_list('CENT')
# #         centres.clear()
#         return normals, centres
#     def _clear_ring_cache(self):
#         normals = self._get_cache_list('NORM')
# #         normals.clear()
#         centres = self._get_cache_list('CENT')
# #         centres.clear()
#         return normals, centres

    def _build_ring_data_cache(self):
        #TODO: remove double normal calculation
#         self._clear_ring_cache()
        
        normals = self._get_cache_list('NORM')
        centres = self._get_cache_list('CENT')
        rings = self._get_component_list('RING')

        self._ring_data_calculator(rings, normals, centres)
    
#    def _calc_single_force_set(self, index, force_factor, forces):
#        pass 


#    def _calc_single_force_set(self, index, force_factor, forces ):
#        component = self._get_component_list().get_component(index)
#        target_atom_id, atom_type_id = component
#        coef_components = self._get_component_list('COEF').get_components_for_atom_id(atom_type_id)
#        for coef_component in coef_components:
#    #                print 'coef_component', i, coef_component
#    #                 print 'here', target_atom_id, atom_type_id,
#    #                coef_component, force_factor
#            self.calculate_single_ring_forces(target_atom_id, atom_type_id, coef_component, force_factor, forces)
#                    
#    def calc_single_atom_force_set(self, target_atom_id, force_factor, forces):
#        target_atom_ids = [target_atom_id]
#        force_factors = [force_factor]
#        self.calc_force_set(target_atom_ids, force_factors, forces)

        
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

    def _create_component_list(self, name):
        if name == 'ATOM':
            return Native_component_list(format='ii')
        elif name == 'RING':
            def ring_translator(component):
                result = list()
                result.append(component[0])
                result.append(len(component[1]))
                result.extend(component[1])
                for i in range(8-len(result)):
                    result.append(0)
                return result
            
            return Native_component_list(format='iiiiiiii', translator=ring_translator)
        elif name == 'COEF':
            return Native_component_list(format='iif')
        else:
            raise Exception('unexpected component %s' % name)

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
        self._min_residue_seperation = min_residue_separation
        
        self._reset()
        

        self._verbose = False 
        self._non_bonded_list_calculator = self._get_non_bonded_calculator()

        
    def _reset(self):
        self._box_update_count = self._update_frequency
        self._non_bonded = []
        self._non_bonded_calculation_count = 0
        self._non_bonded_call_count = 0
        self._interaction_count = 0
        self._pos_cache = {}
        self._num_target_atoms = None
        
    def get_cutoff_distance(self):
        return self._cutoff_distance
    
    def _get_non_bonded_calculator(self):
#        if self._fast:
        result = Fast_non_bonded_calculator(self._min_residue_seperation, self._cutoff_distance, self._jitter)
#        else:
#        result  = Non_bonded_calculator(self._min_residue_seperation, self._cutoff_distance, self._jitter)
        return result
    
    def _get_cutoff_distance_2(self):
        return (self._cutoff_distance+self._jitter)**2
    

    def set_verbose(self,on):
        self._verbose = on
        self._non_bonded_list_calculator.set_verbose(on)

    
    def update(self):
        self._box_update_count +=1
        if self._verbose:
             print "  BOX COUNT UPDATED TO: ", self._box_update_count
        
    def get_boxes(self,component_list_1, component_list_2,target_component_list,coefficient_list):
#        print self._box_update_count, se        if self._verbose:
        updated = False
        if self._verbose:
            print '  update, box_count= ',self._box_update_count

        if self._box_update_count < self._update_frequency:
            if self._verbose:
                print '  skip boxes call count = ', self._box_update_count
        else:
            if self._verbose:
                print '  update boxes call count = ', self._box_update_count

            target_component_list.clear()
            self._box_update_count = -1
            self._non_bonded_calculation_count += 1
            
            native_target_atom_list = component_list_1.get_native_components()
            native_remote_atom_list = component_list_2.get_native_components()
            
            self._num_target_atoms = len(component_list_1.get_component_atom_ids())
            self._non_bonded_list_calculator(native_target_atom_list, native_remote_atom_list, target_component_list)
            
        
            updated = True

        return updated
    
    def get_num_target_atoms(self):
        return self._num_target_atoms
    
    def _get_cached_pos(self,atom_id):
        if atom_id in  self._pos_cache:
            result = self._pos_cache[atom_id]
        else:
            result = Atom_utils._get_atom_pos(atom_id)
            self._pos_cache = result
        return result
    
        

    class NonBoolean:
        def __nonzero__(self):
            raise Exception("internal error this object should never be used in a boolean comparison!!")
    



            
    def __str__(self):
        result  = []
        
        result.append('non_bonded_list (%i elements)' % len(self._non_bonded))
        result.append('')
        result.append('cutoff: %4.2f jitter: %4.2f')
        result.append('')
        result.append('update frequency: %i ' % self._update_frequency)
        result.append('current update: %i ' % self._box_update_count)
        result.append('calculation count: %7.4f call count: %7.4f' % (self._non_bonded_calculation_count, self._non_bonded_call_count))
            
        return '\n'.join(result)
        
        
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
            
            
            for table in self._sort_tables_by_index(non_bonded_tables):
                for target_atom in table.get_target_atoms():
                    self.coefficients.append(table.get_non_bonded_coefficient(target_atom,sphere,*chem_type))
            
#            print self.chem_type_id, self.exponent, self.coefficients
        
            self.complete =  True
            
        def _sort_tables_by_index(self,tables):
            table_map = {}
            
            for table in tables:
                table_map[table.get_table_index()] =  table
                
            result = []
            for key in sorted(table_map.keys()):
                result.append(table_map[key])
            return result
    
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
        self._component_set = None
        self._selected_components =  None
    
    def _get_shift_calculator(self):
        result  = Fast_non_bonded_shift_calculator(smoothed=self._smoothed, name = self.get_abbreviated_name())
        result.set_verbose(self._verbose)
        return result
    
    def _get_force_calculator(self):
#        if self._fast:
        result = Fast_non_bonded_force_calculator( smoothed=self._smoothed, name = self.get_abbreviated_name())
#        else:
#        result = Non_bonded_force_calculator(self._get_indices(), smoothed=self._smoothed)
        result.set_verbose(self._verbose)
        return result
    
    def set_verbose(self, on=True):
        Distance_based_potential.set_verbose(self,on)
        self._non_bonded_list.set_verbose(on)
            
    def _get_table_source(self):
        return Table_manager.get_default_table_manager().get_non_bonded_table
        
    def get_abbreviated_name(self):
        return NON_BONDED
    
    def _get_indices(self):
        return Distance_based_potential.Indices(target_atom_index=0, distance_atom_index_1=0, distance_atom_index_2=1, coefficent_index=2, exponent_index=3)
    
    def _get_target_atom_list_name(self):
        return 'ATOM'
    
    def _get_distance_list_name(self):
        return 'NBLT'

    def _get_non_bonded_list(self):
        non_bonded_list = self._get_component_list('NBLT')
        
        target_atom_list = self._get_component_list('ATOM')
        
        remote_atom_list = self._get_component_list('NBRM')
        
        coefficient_list  = self._get_component_list('COEF')
        updated = self._non_bonded_list.get_boxes(target_atom_list, remote_atom_list, non_bonded_list, coefficient_list)
        
        return non_bonded_list
    
    def _prepare(self, change, target_atom_ids):
        super(Non_bonded_potential, self)._prepare(change, target_atom_ids)
        
        if change == STRUCTURE_CHANGED:
            self._reset_non_bonded_list()
        
        if change == ROUND_CHANGED:
            self.update_non_bonded_list()
        
        if change == TARGET_ATOM_IDS_CHANGED:
            self._selected_components = self._build_selected_components(target_atom_ids)

        
    def get_target_atom_ids(self):
        #TODO: this  call is required check why....
        self.update_non_bonded_list(increment=False)
        components = self._get_component_list('ATOM')
        return components.get_component_atom_ids()
    
    def _reset_non_bonded_list(self, increment=True):
        self._non_bonded_list._reset()
            
    def update_non_bonded_list(self, increment=True):
        if increment:
            self._non_bonded_list.update()
        self._get_non_bonded_list()
        

    
    def calc_single_atom_force_set(self, target_atom_id, force_factor, forces):
#        self.update_non_bonded_list()
        return Distance_based_potential.calc_single_atom_force_set(self, target_atom_id, force_factor, forces)
    
#    def _calc_single_force_set(self,index,factor, forces):
#        print 'here'
#        target_atom_index,distant_atom_index = self._get_target_and_distant_atom_ids(index)
#        
#        distance  = Atom_utils._calculate_distance(target_atom_index, distant_atom_index)
##        TODO: this should be the non bonded distance cutoff
#        if distance < 5.0:
#            self._force_calculator._set_components(self._get_component_list())
#            self._force_calculator._calc_single_force_set(index, factor, forces)

    

    
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
        
        @abstractmethod
        def get_name(self):
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
            
            name = self.get_name()
            result.append('%s index (%i entries)' % (name,self.get_max_index()))
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
                
        def get_name(self):
            return 'sphere'
        
        
    class Chem_type_indexer(Base_indexer):

        def __init__(self,table_manager):
            super(Non_bonded_potential.Chem_type_indexer, self).__init__(table_manager)
            
        def get_name(self):
            return 'chem type'
          
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
        residue_types = table_manager.get_all_known_table_residue_types()   
        
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
                
        if self._verbose:
            result.append(self._non_bonded_list.__str__())
            result.append('')
            
            non_bonded_list = self._get_component_list('NBLT')
            result.append('distances (%i)' %  (len(non_bonded_list)/2))
            result.append('')
            
            seen_distances  = set()
            for i,elem in enumerate(non_bonded_list):
                
                target_atom_id = elem[0]
                remote_atom_id = elem[1]
                
                distance_key =  target_atom_id,remote_atom_id
                if distance_key not in seen_distances:
                    distance = Atom_utils._calculate_distance(target_atom_id, remote_atom_id)
                    active = ' ' 
                    #TODO add a getter for private variable
                    if distance < self._non_bonded_list.get_cutoff_distance():
                        active = 'A'
                    data = i, Atom_utils._get_atom_info_from_index(target_atom_id), Atom_utils._get_atom_info_from_index(remote_atom_id),distance, active
                    result.append('%6i %-17s - %-17s %- 7.3f  %s' % data)
                    seen_distances.add(distance_key)
            
        return '\n'.join(result)

    def _create_component_list(self, name):
        if name  == 'NBLT':
            return Non_bonded_interaction_list()
        elif name == "ATOM":
            return Native_component_list(format='ii')
        elif name == 'NBRM':
            return Native_component_list(format='iii')
        elif name == 'COEF':
            def coef_translator(component):
                result = list()
                result.extend(component)
                COEF_COMP_LENGTH = 3 + (3 * 7)
                extension_length =  COEF_COMP_LENGTH - len(component)
                result.extend([0.0]*extension_length)
                for i,elem in enumerate(result):
                    if elem ==  None:
                        result[i]=0.0
                return result
            
            def coef_expander(components):
                keyed_components = {}
                for component in components:
                    keyed_components[component[0]]= component
                
                result = []
                last_component_index  =  components[-1][0]
                for i in range(0,last_component_index+1):
                    if i in keyed_components:
                        result.append(keyed_components[i])
                    else:
                        num_floats = len(components[0])-2
                        result.append([0,0] + ([0.0] * num_floats))                
                return result
            
            return Native_component_list('iif'+ ('f'*3*7), translator=coef_translator, preparer=coef_expander)
        else:
            return Component_list()


            

    def _get_components(self):
        if self._component_set ==  None:
            target_component_list = self._get_component_list('ATOM')
            native_target_atom_list = target_component_list.get_native_components()
            
            remote_component_list = self._get_component_list('NBRM')
            native_remote_atom_list = remote_component_list.get_native_components()
            
            coefficient_list = self._get_component_list('COEF')
            native_coefficient_list = coefficient_list.get_native_components()
    
            non_bonded_list = self._get_component_list('NBLT')
            
            self._component_set = {'NBLT':non_bonded_list, 'ATOM':native_target_atom_list, 
                          'NBRM':native_remote_atom_list, 'COEF':native_coefficient_list,
                          'OFFS' : 0}
        
        return self._component_set
    
    def _build_selected_components(self, target_atom_ids):
                
        target_component_list = self._get_component_list('ATOM')
        non_bonded_list =  self._get_component_list('NBLT')

        
        num_target_atoms = self._non_bonded_list.get_num_target_atoms()
        
        self._get_components()['OFFS'] = 0
        
        if num_target_atoms != len(target_atom_ids):
            active_components  = self._build_active_components_list(target_atom_ids, non_bonded_list)
            if len(active_components) > 0:
                self._get_components()['OFFS'] = -non_bonded_list[active_components[0]][3]
        else:
            active_components = None
            
        return  active_components

    #TODO: this is no longer as canonical as the other versions 
    # NB calculator doesn't have a 1:1 correspondence between nb lists elememnts and distance components....
    def _calc_single_atom_shift(self, target_atom_id):
        components =  self._get_component_list('NBLT')
        
        result = 0.0
        for i,component in enumerate(components):
            if component[0] ==  target_atom_id:
                result += self._calc_component_shift(i)
        
        return result    
    
    def _calc_component_shift(self, index):
         
        calc = self._get_shift_calculator()
            
        components = self._get_components()
         
        target_component_list = self._get_component_list('ATOM')
        non_bonded_list =  components['NBLT']
         
        target_atom_ids = [non_bonded_list[index][0],]
        
        active_components = array.array('i',[index,])
        components['OFFS'] = - non_bonded_list[index][3]
        
        component_to_result = array.array('i',[0,])
        
        results = array.array('d',[0.0])

        calc(components,results,component_to_result, active_components=active_components)
 
        return results[0]
    
    def _get_active_components(self):
        return self._selected_components
    
#     def calc_shifts(self, target_atom_ids, results):
#          
#         components = self._get_components()
#         
#         self._shift_calculator(components,results,self._component_to_result, active_components=self._selected_components)
# 
# 
# 
#     def calc_force_set(self,target_atom_ids,force_factors,forces):
# 
#            
#         components = self._get_components()
#            
        self._force_calculator(components, self._component_to_result, force_factors, forces, active_components=self._selected_components)
         


from pyPot import PyPot

class Xcamshift(PyPot):

    
    
    


    def _get_force_factor_calculator(self):
#        if self._fast:
        result  =  Fast_force_factor_calculator()
#        else:
#            result = Force_factor_calculator()
        return result
    
            
            
    def __init__(self, verbose=False):
        super(Xcamshift, self).__init__("xcamshift")
        self.potential = [
                          RandomCoilShifts(),
                          Distance_potential(),
                          Extra_potential(),
                          Dihedral_potential(),
                          Sidechain_potential(),
                          Ring_Potential(),
                          Non_bonded_potential(),
                          Disulphide_shift_calculator()
                          ]
        self._verbose=verbose
        self._shift_table = Observed_shift_table()
        self._shift_cache = None
        self._out_array =  None
        self._energy_term_cache = None
        self._energy_calculator = self._get_energy_calculator()
        self._force_factor_calculator =  self._get_force_factor_calculator()

        self.set_verbose(verbose)
        self._freeze = False
        self._active_target_atom_ids = None
        self._factors = None

    
    
    

    def set_verbose(self,on=True):
        for potential in self.potential:
            potential.set_verbose(on)
        self._force_factor_calculator.set_verbose(on)
        self._energy_calculator.set_verbose(on)
        self._verbose=on
    
    def _get_energy_term_cache(self):
        if self._energy_term_cache == None:
            self._energy_term_cache = self._create_energy_term_cache()
        return self._energy_term_cache
    
    def _update_calculator(self, calculator):
        calculator.set_calculated_shifts(self._shift_cache)
        calculator.set_observed_shifts(self._shift_table.get_native_shifts(self._get_active_target_atom_ids()))
        calculator.set_energy_term_cache(self._get_energy_term_cache().get_native_components())
    
    def update_energy_calculator(self):
        self._update_calculator(self._energy_calculator)

    def update_force_factor_calculator(self):
        self._update_calculator(self._force_factor_calculator)
    
    #TODO:  might be better to have a prepapre shift cache function    
    def clear_shift_cache(self):
        if self._shift_cache != None:
            zero_array(self._shift_cache)
        
    def get_sub_potential_names(self):
        return [potential.get_abbreviated_name() for potential in self.potential]
    
    def _get_energy_calculator(self):
#        if self._fast:
        result  = Fast_energy_calculator()
#        else:
#            result  = Energy_calculator()
        result.set_verbose(self._verbose)
        
        return result
    
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
            key = potential.L()
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
            
        
    def _get_all_component_target_atom_ids(self):
        result  = set()
        for potential in self.potential:
            result.update(potential.get_target_atom_ids())
        result = list(result)
        result.sort()
        return result
    
    def _calc_shift_cache(self,target_atom_ids):
        if self._verbose:
            start_time =  time()
            
        if self._shift_cache == None:
            self._shift_cache =  allocate_array(len(target_atom_ids),'d')
        elif len(self._shift_cache) != len(target_atom_ids):
            resize_array(self._shift_cache, len(target_atom_ids))
        
        self.clear_shift_cache()
        self._calc_shifts(target_atom_ids, self._shift_cache)
        
        if self._verbose:
            end_time = time()
            
            print "shifts completed in ", "%.17g " %  (end_time-start_time),"seconds"

    def calc_shifts(self, target_atom_ids=None, result=None):
        def target_atom_ids_as_selection_strings(target_atom_ids):
            result  = []
            for target_atom_id in target_atom_ids:
                result.append(Atom_utils._get_atom_info_from_index(target_atom_id))
            return tuple(result)

        if target_atom_ids  == None:
            target_atom_ids = self._get_all_component_target_atom_ids() 
            
  
        if result == None or len(result) < len(target_atom_ids):
            result = array.array('d',[0.0] *len(target_atom_ids))


        
        self._calc_shifts(target_atom_ids, result)
        
        #TODO review whole funtion and resturn types
        return (target_atom_ids_as_selection_strings(target_atom_ids), tuple(result))
        
    #TODO: add standalone mode flag or wrap in a external wrapper that call round changed etc    
    def _calc_shifts(self, target_atom_ids=None, result=None):
                
        if self._verbose:
            start_time =  time()
            print 'start calc shifts'
            
            

        #TODO: currently needed to generate a component to result remove by making component to result lazy?
        if not self._freeze:
            self._prepare(TARGET_ATOM_IDS_CHANGED, target_atom_ids)
            
        for potential in self.potential:
            potential.calc_shifts(target_atom_ids, result)

       
        if self._verbose:
            end_time =  time()
            print 'end calc shifts (calculated shifts in  %.17g seconds) \n' % (end_time-start_time)
        
                           
    #TODO: deprecated remove use _calc_shifts
    def set_shifts(self, result):
        print 'deprecated remove use _calc_shifts'
        target_atom_ids =  self._get_all_component_target_atom_ids()
        self._prepare(TARGET_ATOM_IDS_CHANGED, target_atom_ids)
        result_shifts  = allocate_array(len(target_atom_ids))
        self._calc_shifts(target_atom_ids, result_shifts)
        for target_atom_id, shift in zip(target_atom_ids,result_shifts):
            result[target_atom_id] =  shift
        
        return result
    
    def set_observed_shifts(self, shift_table):
        self._shift_table  =  shift_table
        #TODO: could be better
        self._shift_cache =  None
        self._energy_term_cache =  None
        self._prepare(SHIFT_DATA_CHANGED,None)
        

    def _calc_single_atom_shift(self,target_atom_id):
        #TODO: could be better
        self._shift_cache = allocate_array(1)
        self._calc_shift_cache([target_atom_id,])
        return  self._shift_cache [0]
    
    
    
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
        
        theory_shift = self._calc_single_atom_shift(target_atom_index)
        
        observed_shift = self._shift_table.get_chemical_shift(target_atom_index)
        
        return observed_shift - theory_shift
        
#    def _calc_single_atom_shift(self,target_atom_index):
#        calculated_shift  = self._calc_single_atom_shift(target_atom_index)
#        self._shift_cache[target_atom_index] =  calculated_shift
        
    def _calc_single_atom_energy(self, target_atom_index):
        
        target_atom_ids = array.array('i',[target_atom_index])
        self._prepare(TARGET_ATOM_IDS_CHANGED,target_atom_ids)
        
        #TODO: setting a single atom shift with self._calc_single_atom_shift(target_atom_index) doesn't work as the 
        # and then using a single index doesn't work as the energy terms are indexed as well
        self._calc_shift_cache(self._get_active_target_atom_ids())
#         _single_atom_shift(target_atom_index)
        self.update_energy_calculator()

        #TODO: this is a hack remove!        
        active_target_atom_ids =  self._get_active_target_atom_ids()
        if active_target_atom_ids == None:
            
            active_target_indices = array.array('i',[0])
            active_target_atom_ids =  array.array('i',[target_atom_index])
        else:
            active_target_indices =  array.array('i',[active_target_atom_ids.index(target_atom_index)])
        return self._energy_calculator(active_target_atom_ids,active_target_indices)
        
    
    def _calc_force_set_with_potentials(self, target_atom_ids, forces, potentials_list):
        if self._factors == None or len(self._factors) < len(target_atom_ids):
            self._factors = allocate_array(len(target_atom_ids),'f')
        elif len(self._factors) > len(target_atom_ids):
            self._factors =  resize_array(self._factors,len(target_atom_ids))
        
        self._calc_factors(target_atom_ids, self._factors)
        for potential in potentials_list:
            potential.calc_force_set(target_atom_ids,self._factors,forces)
        
             
    def _calc_single_atom_force_set_with_potentials(self, target_atom_id, forces, potentials_list):
#        print 'forces _calc_single_atom_force_set_with_potentials', forces 
#        self._calc_single_force_factor(target_atom_id, forces)
        target_atom_ids = array.array('i',[target_atom_id])
        self._calc_force_set_with_potentials(target_atom_ids, forces, potentials_list)

    def _calc_force_set(self,target_atom_ids,forces,potentials=None):
        if potentials ==  None:
            potentials =  self.potential
        self._calc_force_set_with_potentials(target_atom_ids, forces, potentials)
         
    def _calc_single_atom_force_set(self,target_atom_id,forces,potentials=None):
#        print 'forces xcamshift._calc_single_atom_force_set',forces
        target_atom_ids = [target_atom_id]
        self._calc_force_set(target_atom_ids, forces, potentials)
        

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
        target_atom_ids = array.array('i',[target_atom_id])
        
        if not target_atom_id in self._shift_table.get_atom_indices():
            msg = "requested factor for target [%s] which is not in shift table"
            target_atom_info = Atom_utils._get_atom_info_from_index(target_atom_id)
            raise Exception(msg % target_atom_info)
        
        result = allocate_array(1,'f')
        self. _calc_factors(target_atom_ids, result)
        return  result[0]
    
    def _calc_factors(self, target_atom_ids, factors):
        #TODO move to prepare or function called by prepare
        active_components = None
        active_target_atom_ids = self._get_active_target_atom_ids()
        if len(target_atom_ids) != len(active_target_atom_ids):
            active_components =  []
            for target_atom_id in target_atom_ids:
                active_components.append(active_target_atom_ids.index(target_atom_id))
            active_components = array.array('i',active_components)
            
        self._force_factor_calculator(active_target_atom_ids, factors, active_components)
        return factors
    

    def _calc_single_force_factor(self,target_atom_index,forces):
        factor = 1.0
        for potential in self.potential:
            if potential._have_derivative():
                potential.set_observed_shifts(self._shift_table)
                
                if target_atom_index in self._shift_table.get_atom_indices():
                    factor = self._calc_single_factor(target_atom_index)
#                    potential._calc_single_factor(target_atom_index,factor,forces)
        return factor
    


    def get_selected_atom_ids(self):
       
        
        return   self._shift_table.get_atom_indices()

    def _get_active_target_atom_ids(self):
        
        if self._active_target_atom_ids == None:
            target_atom_ids = set(self._get_all_component_target_atom_ids())
            observed_shift_atom_ids = self.get_selected_atom_ids()
            active_target_atom_ids = target_atom_ids.intersection(observed_shift_atom_ids)
            active_target_atom_ids = list(active_target_atom_ids)
            
            self._active_target_atom_ids =  array.array('i',sorted(active_target_atom_ids))
        
         
        return self._active_target_atom_ids
    


    def _prepare_potentials(self, change, data):
        
        for potential in self.potential:
            potential._prepare(change, data)

    def _prepare(self, change, data):
        if self._verbose:
            print 'do prepare %s' % change
            
        if change == STRUCTURE_CHANGED:
            self.clear_shift_cache()
            self._active_target_atom_ids = None
            
        if change == TARGET_ATOM_IDS_CHANGED:
            self._active_target_atom_ids = None
            if data ==  None:
                data  = self._get_active_target_atom_ids()
        
        if change ==  SHIFT_DATA_CHANGED:
            self._active_target_atom_ids = None
            self.update_energy_calculator()
        
        self._prepare_potentials(change, data)
         
        
    
    def _get_constants(self, residue_type, atom_name):
        return                                                             \
            self._get_flat_bottom_shift_limit(residue_type, atom_name),    \
            self._get_end_harmonic(residue_type, atom_name),               \
            self._get_scale_harmonic(residue_type, atom_name),             \
            self._get_weight(residue_type, atom_name),                     \
            self._get_tanh_amplitude(residue_type,atom_name),              \
            self._get_tanh_elongation(residue_type, atom_name),            \
            self._get_tanh_y_offset(residue_type, atom_name)
        
    
    def _create_energy_term_cache(self):
        cache = Native_component_list(format = 'i' + ('f'*7))
        seen_types = {}
        table_manager =  Table_manager.get_default_table_manager()
        for target_atom_index in self._get_active_target_atom_ids():
            residue_type = Atom_utils._get_residue_type_from_atom_id(target_atom_index)
            atom_name = Atom_utils._get_atom_name_from_index(target_atom_index)
            

            table = table_manager.get_constants_table(residue_type)
            table_key = table.get_table_residue_type(), atom_name
            if not table_key in seen_types:
                seen_types[table_key] =  self._get_constants(residue_type, atom_name)
            cache_data = [target_atom_index]
            cache_data.extend(seen_types[table_key])
            cache_data = tuple(cache_data)
            cache.add_component(cache_data)

        return cache
    

    def _get_out_array(self):
        num_atoms = Segment_Manager.get_segment_manager().get_number_atoms()
        result = self._out_array
        if result == None:
            result = Out_array(num_atoms)
        else:
            result.realloc(num_atoms)
        return result
        


        
    def calcEnergy(self, prepare =  True, active_target_atom_ids = None):
        if self._verbose:
            start_time = time()

        if active_target_atom_ids == None:
            active_target_atom_ids = self._get_active_target_atom_ids()
            
        if prepare:
            self._prepare(TARGET_ATOM_IDS_CHANGED,active_target_atom_ids)
            self._prepare(ROUND_CHANGED,None)
            self._calc_shift_cache(active_target_atom_ids)
            

        self.update_energy_calculator()
        
        if self._verbose:
            end_time =  time()
            print 'energy calculation completed in %.17g seconds ' % (end_time-start_time), "seconds"
            
        return self._energy_calculator(active_target_atom_ids)
    
    

    def _calc_derivs(self, derivs, active_target_atom_ids ,potentials=None):
        
        out_array = self._get_out_array()
        self.update_force_factor_calculator()
        self._calc_force_set(active_target_atom_ids, out_array, potentials)
        out_array.add_forces_to_result(derivs)

    def calcEnergyAndDerivs(self,derivs):
        if self._verbose:
            print "start energy and derivatives"
            start_time = time()
        
        if self._verbose:
            print 'start_prepare'
            prepare_start_time =  time()
            
        target_atom_ids = self._get_active_target_atom_ids()
        self._prepare(ROUND_CHANGED, None)
        self._calc_shift_cache(target_atom_ids)
        if self._verbose:
            prepare_end_time = time()
            print 'prepare completed in %.17g seconds' %(prepare_end_time-prepare_start_time)

        energy = self.calcEnergy( active_target_atom_ids=target_atom_ids, prepare=False)

        self._calc_derivs(derivs, target_atom_ids)
        
        
        if self._verbose:
            end_time = time()
        
            print "energy and derivatives completed in %.17g " %  (end_time-start_time)," second\ns"
        return energy
    
    def _set_frozen(self,on=True):
        for potential in self.potential:
            potential.set_frozen(on)
        self._freeze = True
            
    def setup(self):
        #TODO: do we need a STRUCTURE_CHANGED here
        self._prepare(TARGET_ATOM_IDS_CHANGED, self._get_active_target_atom_ids())
        self._set_frozen()
    
    def reset(self):
        self._prepare(STRUCTURE_CHANGED, None)
        self._prepare(TARGET_ATOM_IDS_CHANGED, None)
        #TODO: do we need a  set froze here

    


