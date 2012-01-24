'''
Created on 24 Jan 2012

@author: garyt
'''
from xcamshift import Base_potential
import sys
from utils import tupleit

class Observed_shift_table(object):
    '''
    classdocs
    '''

    def __init__(self,shift_data):
        self._chemical_shifts = self.process_observed_shifts(shift_data)
        
    def process_observed_shifts(self,shift_data):
        result  = {}
        for key in shift_data:
            search_key = key
            if len(key) == 2:
                search_key = list(key)
                search_key.insert(0,"*")
                search_key =  tuple(search_key)
                #TODO move Base_potential components to utilities
                atoms = Base_potential.find_atom(*search_key)
                for atom in atoms:
                    result[atom.index()] = shift_data[key]
        return result
    
    #TODO add more generic dump capabilities
    def dump_observed_shifts(self):
        results = []
        for atom_index in self._chemical_shifts:
            sub_result  = []
            results.append(sub_result)
            sub_result.append(Base_potential._get_atom_info_from_index(atom_index))
            sub_result.append(self._chemical_shifts[atom_index])
        return tupleit(results)
        