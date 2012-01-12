'''
Created on 30 Dec 2011

@author: garyt
'''

class Extra_table(object):
    
    DATA = 'data'

    EXPONENT = 'exponent'
    ATOM_1=0
    ATOM_2=1
    
    EXTRA_ATOMS='extra_atoms_1','extra_atoms_2'
    OFFSETS='extra_offsets_1','extra_offsets_2'
    

    TARGET_ATOMS='target_atoms'
    
    def __init__(self, table):
        self.table = table
    
    def get_translation(self,atom):
        return atom
    
    def get_distance_atoms(self,atom_index):
        return self.table[self.EXTRA_ATOMS[atom_index]]
    
    def get_offsets(self,atom_index):
        return self.table[self.OFFSETS[atom_index]]
    
    def get_target_atoms(self):
        return self.table[self.TARGET_ATOMS]
    
    def get_exponent(self):
        return self.table[self.EXPONENT]

    def _check_offset(self, atom_index, offset):
        
        offsets = self.get_offsets(atom_index)
        
        if not offset in offsets:
            offset_strings = [`offset` for offset in offsets]
            message_values = atom_index, offset, ', '.join(offset_strings)
            message = "atom_%i offset %s is not in %s" % message_values
            raise KeyError(message)
    
    def _check_extra_atom(self, atom_index, atom_name):
        
        atoms = self.get_distance_atoms(atom_index)
        if atom_name not in atoms:
            message = "atom_name %s not in extra shift table atoms_%i (%s)" % (atom_name, atom_index, ', '.join(atoms))
            raise KeyError(message)

    def _check_distance_key(self, atom_index, offset, atom_name):
        
        self._check_offset(atom_index, offset)
        
        self._check_extra_atom(atom_index, atom_name)
    
    def _check_target_atom_key(self, atom):
        
        atoms = self._get_target_atoms()
        if not atom in atoms:
            message = "atom %s not in target atoms %s" % (atom,', '.join(atoms))
            raise KeyError(message)
    
    
    def get_extra_shift(self,target_atom,key_values_1,key_values_2):
        self._check_distance_key(self.ATOM_1,key_values_1.offset,key_values_1.atom)
        self._check_distance_key(self.ATOM_2, key_values_2.offset,key_values_2.atom)
        
        result = None
        
        key_1 = key_values_1.get_key()
        key_2 = key_values_2.get_key()
        if key_1 in self.table[self.DATA]:
            data_1 = self.table[self.DATA][key_values_1.get_key()]
            if key_2 in data_1:
                target_atoms = data_1[key_2]
                if target_atom in target_atoms:
                    result  = target_atoms[target_atom]
            
        return result