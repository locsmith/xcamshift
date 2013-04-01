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
Created on 11 Feb 2012

@author: garyt
'''
from bisect import insort

class Component_list():
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
        return self._components[start:end]
    
    def get_component(self,i):
        return self._components[i]
        
    def _index_components(self, components):
        component_id_list =[]
        for component in components:
            component_id_list.append(component[0])
            
        reversed_component_id_list =  list(component_id_list)
        reversed_component_id_list.reverse()
        
        result = {}
        for component_id in self._component_ids:
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
#TODO: not a useful string
#    def __str__(self):
#        result = []
#        for target_atom_id in self._component_ids:
#            result.append("%i %s" %  (target_atom_id, `self.get_component_range(target_atom_id)`))
#        return ",".join(result)
    
