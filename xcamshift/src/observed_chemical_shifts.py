'''
Created on 24 Jan 2012

@author: garyt
'''
from utils import tupleit, Atom_utils

class Observed_shift_table(object):
    '''
    classdocs
    '''

    def __init__(self,shift_data={}):
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
                atoms = Atom_utils.find_atom(*search_key)
                for atom in atoms:
                    result[atom.index()] = shift_data[key]
        return result
    
    def get_chemical_shift(self,atom_index):
        return self._chemical_shifts[atom_index]
    
    def get_atom_indices(self):
        return self._chemical_shifts.get_atom_indices()
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
        result  = []
        for (segid,resid,atom_name),shift in self.dump_observed_shifts():
            result.append("[%4s]:%i@%s  %7.4f" % (segid,resid,atom_name,shift))
        return "\n".join(result)
            