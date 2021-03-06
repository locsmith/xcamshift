#-------------------------------------------------------------------------------
# Copyright (c) 2013-2015 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 11 Feb 2012

@author: garyt
'''
from bisect import insort
from collections import defaultdict
from struct import Struct
import ctypes
from copy import copy
import array

class Component_list(object):
    def __init__(self):
        self._components = []
        self._component_offsets = None
        self._component_ids =  set()

        

    def _clear_component_offsets(self):
        self._component_offsets = None


    def _basic_add_component(self,component):
        insort(self._components, component)
        self._component_ids.add(component[0])

    def add_component(self, component):
        self._basic_add_component(component)
        self._clear_component_offsets()
        
    def add_components(self,components):
        for component in components:
            self._basic_add_component(component)
        self._clear_component_offsets()
    
    def replace_component(self,index,component):
        atom_id = self._components[index][0]
        range = self.get_component_range(atom_id)
        if range[1]-range[0] == 1:
            self._component_ids.remove(atom_id)
        del self._components[index]
        self.add_component(component)
        
        
        
    def get_all_components(self):
        return tuple(self._components)
    
#    TODO rename to get_target_atom_ids
    def get_component_atom_ids(self):
        result = sorted(list(self._component_ids))
        return tuple(result)
    

    def report_bad_atom_id(self, atom_id):
        component_ranges = self._get_component_ranges()
        
        template = "There is no component with the atom id %i, the list of known atom ids is: [%s]"
        keys = component_ranges.keys()
        keys.sort()
        
        message = template % (atom_id, ','.join(keys))
        
        raise Exception(message)

    #TODO: make private
    def get_component_range(self, atom_id):
        component_ranges = self._get_component_ranges()
        start, end = 0,0
        if atom_id in component_ranges:
            start, end = component_ranges[atom_id]
        return start, end

    def get_components_for_atom_id(self,atom_id):
        start, end = self.get_component_range(atom_id)
        result = Component_list()
        result.add_components(self._components[start:end])
        return result
    
    def get_component(self,i):
        return self._components[i]
        
    def _index_components(self, components):
#         print components
        component_id_list =[]
        for component in components:
            component_id_list.append(component[0])
            
        reversed_component_id_list =  list(component_id_list)
        reversed_component_id_list.reverse()
        
        result = {}
#         print self._component_ids
        for component_id in self._component_ids:
#             print component_id_list, component_id
            start = component_id_list.index(component_id)
            end = len(component_id_list) - reversed_component_id_list.index(component_id)
            result[component_id] = (start,end)
        
        return result

    def _get_component_ranges(self):
        if self._component_offsets  == None:
            self._component_offsets = self._index_components(self._components)
        return self._component_offsets
    
    def get_number_components(self):
        return len(self._components)
    
    def get_number_target_atoms(self):
        return len(self._component_ids)
    
    def __iter__(self):
        return self.next()
    
    def next(self):
        for component in self._components:
            yield component
            
    def __len__(self):
        return len(self._components)
    
    def __getitem__(self,index):
        return self._components[index]
    
    def __str__(self):
        return self._components.__str__()
    
    def clear(self):
        self._components =[]
        self._component_ids = set()
    
    def build_selection_list(self, accept = lambda x: True):
        result = []
        for i,elem in enumerate(self._components):
            if accept(elem):
                result.append(i)
        return array.array('i',result)
                            

class Native_component_list(Component_list):
    def __init__(self, format, translator = lambda x : x, preparer = lambda x : x):
        super(Native_component_list, self).__init__()
        #TODO remove and simplifiy!
        self.set_translator(translator)
        self.set_format(format)    
        self._native_components =  None
        self._native_component_offsets = None
        self.set_preparer(preparer)

   
    def set_translator(self, translator):
        self._translator =  translator
        
    def set_preparer(self,preparer):
        self._preparer =  preparer
        
    def set_format(self,format):
        self._component_struct =  Struct(format)
    
    def _reset(self):
        self._native_components =  None
        self._native_component_offsets = None
                   
    def _basic_add_component(self,component):
        super(Native_component_list, self)._basic_add_component(component)
        self._reset()
        
    def clear(self):
        super(Native_component_list, self).clear()
        self._reset()
            
    def _translate_to_native_component(self, component):
        return self._translator(component)
    
    def _prepare_component_list(self, components):
        return self._preparer(components)

    def _build_native_components(self):
        
        target_components =  self._prepare_component_list(self._components)
        num_components = len(target_components)
        if num_components > 0 :
            format_length = len(self._component_struct.format)
            data_length =  len(self._translate_to_native_component(target_components[0]))
            
            if format_length != data_length:
                msg = 'data and format must have the same length, got len(fomat) = %i and len(data) = %i'
                raise Exception(msg % (format_length, data_length))
        
        struct_size = self._component_struct.size
        bytes = ctypes.create_string_buffer(struct_size * num_components)
        for i in range(num_components):
            native_component = self._translate_to_native_component(target_components[i])
            self._component_struct.pack_into(bytes, struct_size * i, *native_component)
        
        return bytes

    def get_native_components(self):
        if self._native_components == None:
            self._native_components = self._build_native_components()
        bytes = self._native_components 
        
        return bytes
       
    def _build_native_component_offsets(self):
        struct =  Struct('iii')
        
        ids =  self.get_component_atom_ids()
        num_ids = 0
        if len(ids) > 0:
            num_ids = max(ids)+1
        struct_size = struct.size
        bytes =  ctypes.create_string_buffer(struct_size * num_ids)
         
        for id in range(num_ids):
            if id in ids:
                start,end = self.get_component_range(id)
                native_component = id,start,end-start
            else:
                native_component = id,-1,-1
            struct.pack_into(bytes, struct_size * id, *native_component)
        
        return bytes

    def get_native_component_offsets(self):
        if self._native_component_offsets == None:
            self._native_component_offsets = self._build_native_component_offsets()
        bytes = self._native_component_offsets
        return bytes
    
    
#TODO: not a useful string
#    def __str__(self):
#        result = []
#        for target_atom_id in self._component_ids:
#            result.append("%i %s" %  (target_atom_id, `self.get_component_range(target_atom_id)`))
#        return ",".join(result)
    
