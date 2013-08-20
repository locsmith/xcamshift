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
Created on 24 Jan 2012

@author: garyt
'''

cdef extern from "instantiate.hh":
    pass

from python_utils import tupleit
from utils import Atom_utils
from cython.shift_calculators import allocate_array
from xplor_access cimport  CDSVector
from libcpp.map cimport map 

cdef class Observed_shift_table(object):
    '''
    classdocs
    '''

    def __cinit__(self):
        self._native_shifts_set = False
    
    def __init__(self,shift_data={}):
        self._chemical_shifts = self.process_observed_shifts(shift_data)
        
    def process_observed_shifts(self,shift_data):
        cdef map[int,float]  result   
        for key in shift_data:
            if len(key) == 2:
                search_key = '*',key[0],key[1]
            else:
                search_key = tuple(key)
            if len(search_key) ==  3:
                #TODO move Base_potential components to utilities
                atoms = Atom_utils.find_atom(*search_key)
                for atom in atoms:
                    result[atom.index()] = shift_data[key]
            else:
                template = """key with unexpected length should be either 'segid,residue_no,atom_name' 
                              or 'residue_no,atom_name' but i got '%s'"""
                msg = template % `key`
                raise Exception(msg)
                    
        return result
    
    def get_chemical_shift(self,atom_index):
        return self._chemical_shifts[atom_index]
    
    def get_atom_indices(self):
        return self._chemical_shifts.keys()
    
    def get_indices_for_atom_id(self):
        result ={}
        for i, atom_id in enumerate(self._chemical_shifts.keys()):
            result[atom_id] = i
        return result
        
    #TODO add more generic dump capabilities
    def dump_observed_shifts(self):
        results = []
        for atom_index, value in self._chemical_shifts:
            sub_result  = []
            results.append(sub_result)
            sub_result.append(Atom_utils._get_atom_info_from_index(atom_index))
            sub_result.append(value)
        return tupleit(results)
    
    cdef CDSVector[float] get_native_shifts(self, CDSVector[int] target_atom_ids):
        if not self._native_shifts_set:
            self._native_shifts.resize(target_atom_ids.size())
            for i in range(target_atom_ids.size()):
                target_atom_id = target_atom_ids[i]
                self._native_shifts[i] = self._chemical_shifts[target_atom_id]
                
        return self._native_shifts
            
    def py_get_native_shifts(self, target_atom_ids):
        cdef CDSVector[int] native_target_atom_ids
        cdef CDSVector[float] native_result
        
        native_target_atom_ids.resize(len(target_atom_ids))
                      
        for i in range(len(target_atom_ids)):
           native_target_atom_ids[i] = target_atom_ids[i]
        
        native_result = self.get_native_shifts(native_target_atom_ids)
        
        result = []
        for i in range(len(target_atom_ids)):
            result.append(native_result[i])
        
        return result
        
        
    def __str__(self): 
        
        result  = ["shift table", "-----------",""]
        for (segid,resid,atom_name),shift in self.dump_observed_shifts():
            result.append("[%4s]:%i@%s  %7.4f" % (segid,resid,atom_name,shift))
        return "\n".join(result)
            