'''
Created on 30 Dec 2011

@author: garyt
'''

class Random_coil_table(object):
    
    DATA = 'data'
    def __init__(self, table):
        self.table = table
        

    def __check_residue_keys(self, residue, atom):
        data = self.table[self.DATA]
        
        if residue not in data:
            message = "residue %s not in random coil shift table" % residue
            raise KeyError(message)
        
    
    
    def get_random_coil_shift(self,residue,atom):
        self.__check_residue_keys(residue,atom)
        
        result = None
        
        atom_random_coil_table  = self.table[self.DATA][residue]
        if atom in atom_random_coil_table:
            result =  atom_random_coil_table[atom]
            
        return result