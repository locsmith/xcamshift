'''
Created on 30 Dec 2011

@author: garyt
'''

TRANSLATIONS = 'translations'
        
class Random_coil_table(object):
    
    DATA = 'data'
    def __init__(self, table):
        self._table = table
        
        self.ATOMS = 'atoms'
        self.OFFSETS='offsets'
        
    def get_offsets(self):
        return self._table[self.OFFSETS]
    
    def get_atoms(self):
        return self._table[self.ATOMS]
    
    def __check_residue_keys(self, offset, residue, atom):
        
        if not offset in self.get_offsets():
            offset_strings = [`offset` for offset in self.get_offsets()]
            message_values = offset, ', '.join(offset_strings)
            message = "offset %s is not in %s" % message_values
            raise KeyError(message)
        
        if atom not in self.get_atoms():
            message = "atom %s not in random coil shift _table atoms %s" % (atom,', '.join(self.get_atoms()))
            raise KeyError(message)
        
    def get_random_coil_shift(self,offset, residue, atom):
        self.__check_residue_keys(offset,residue,atom)
        
        result = None
        
        atom_random_coil_table  = self._table[self.DATA][offset][residue]
        if atom in atom_random_coil_table:
            result =  atom_random_coil_table[atom]
            
        return result
    
    def get_translation(self,residue,atom):
        result =  atom
        
        key = residue, atom
        if key in self._table[TRANSLATIONS]:
            result =  self._table[TRANSLATIONS][key]
        return result