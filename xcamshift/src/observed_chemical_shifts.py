#-------------------------------------------------------------------------------
# Copyright (c) 2013 gary thompson.
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
from python_utils import tupleit
from utils import Atom_utils

class Observed_shift_table(object):
    '''
    classdocs
    '''

    def __init__(self,shift_data={}):
        self._chemical_shifts = self.process_observed_shifts(shift_data)
        
    def process_observed_shifts(self,shift_data):
        result  = {}
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
        for atom_index in self._chemical_shifts:
            sub_result  = []
            results.append(sub_result)
            sub_result.append(Atom_utils._get_atom_info_from_index(atom_index))
            sub_result.append(self._chemical_shifts[atom_index])
        return tupleit(results)
    
    def __str__(self): 
        
        result  = ["shift table", "-----------",""]
        for (segid,resid,atom_name),shift in self.dump_observed_shifts():
            result.append("[%4s]:%i@%s  %7.4f" % (segid,resid,atom_name,shift))
        return "\n".join(result)
            
