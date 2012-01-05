'''
Created on 30 Dec 2011

@author: garyt
'''

class Extra_table(object):
    
    DATA = 'data'

    ATOM_1=0
    ATOM_2=1
    
    EXTRA_ATOMS='extra_atoms_1','extra_atoms_2'
    OFFSETS='extra_offsets_1','extra_offsets_2'
    

    TARGET_ATOMS='target_atoms'
    
    def __init__(self, table):
        self.table = table
    
    def _get_atoms(self,atom_index):
        return self.table[self.EXTRA_ATOMS[atom_index]]
    
    def _get_offsets(self,atom_index):
        return self.table[self.OFFSETS[atom_index]]
    
    def _get_target_atoms(self):
        return self.table[self.TARGET_ATOMS]
    

    def _check_offset(self, atom_index, offset):
        
        offsets = self._get_offsets(atom_index)
        
        if not offset in offsets:
            offset_strings = [`offset` for offset in offsets]
            message_values = atom_index, offset, ', '.join(offset_strings)
            message = "atom_%i offset %s is not in %s" % message_values
            raise KeyError(message)
    
    def _check_extra_atom(self, atom_index, atom_name):
        
        atoms = self._get_atoms(atom_index)
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

    def get_extra_shift(self,offset_1,atom_1,offset_2,atom_2,target_atom):
        self._check_distance_key(self.ATOM_1,offset_1,atom_1)
        self._check_distance_key(self.ATOM_2, offset_2,atom_2)
        self._check_target_atom_key(target_atom)
        
        result = None
        key_1 = atom_1,offset_1
        key_2 = atom_2,offset_2
        
        extra_table  = self.table[self.DATA]
        
        if key_1 in extra_table:
            data_1 = extra_table[key_1]
            if key_2 in data_1:
                target_atoms = data_1[key_2]
                if target_atom in target_atoms:
                    result  = target_atoms[target_atom]
            
        return result