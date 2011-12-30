'''
Created on 30 Dec 2011

@author: garyt
'''

class Random_coil_table(object):
    
    DATA = 'data'
    def __init__(self, table):
        self.table = table
        

    def __check_random_coil_keys(self, residue, atom):
        data = self.table[self.DATA]
        
        if residue not in data:
            message = "residue %s not in random coil shift table" % residue
            raise KeyError(message)
        
        if atom not in data[residue]:
            template = "atom %s not found in random coil table for residue %s"
            message = template % (atom,residue)
            raise KeyError(message)
            
    
    
    def get_random_coil_shift(self,residue,atom):
        self.__check_random_coil_keys(residue,atom)
        
        return self.table[self.DATA][residue][atom]