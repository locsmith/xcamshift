'''
Created on 30 Dec 2011

@author: garyt
'''

TRANSLATIONS = 'translations'

class Extra_table(object):
    
    DATA = 'data'

    EXPONENT = 'exponent'
    ATOM_1=0
    ATOM_2=1
    
    EXTRA_ATOMS='extra_atoms_1','extra_atoms_2'
    OFFSETS='extra_offsets_1','extra_offsets_2'
    

    TARGET_ATOMS='target_atoms'
    
    def __init__(self, table):
        self._table = table
        self._translation_to_table = self._build_from_translation_table(table)
        self._translation_from_table = self._build_from_translation_table(table)
        
    def _build_to_translation_table(self,data_table):
        result  = {}
        if TRANSLATIONS in data_table:
            result  =  data_table[TRANSLATIONS]
        return result
    
    def _build_from_translation_table(self,data_table):
        result  = {}
        if TRANSLATIONS in data_table:
            for key,to_atom in data_table[TRANSLATIONS].items():
                residue, from_atom = key
                new_key = (residue,to_atom)
                result[new_key]= from_atom
        return result
    def get_translation(self,atom):
        return atom
    
    def get_distance_atoms(self,atom_index):
        return self._table[self.EXTRA_ATOMS[atom_index]]
    
    def get_offsets(self,atom_index):
        return self._table[self.OFFSETS[atom_index]]
    
    def get_target_atoms(self):
        return self._table[self.TARGET_ATOMS]
    
    def get_exponent(self):
        return self._table[self.EXPONENT]

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
            message = "atom_name %s not in extra shift _table atoms_%i (%s)" % (atom_name, atom_index, ', '.join(atoms))
            raise KeyError(message)

    def _check_distance_key(self, atom_index, offset, atom_name):
        
        self._check_offset(atom_index, offset)
        
        self._check_extra_atom(atom_index, atom_name)
    
    def _check_target_atom(self, atom):
        
        atoms = self._get_target_atoms()
        if not atom in atoms:
            message = "atom %s not in target atoms %s" % (atom,', '.join(atoms))
            raise KeyError(message)
    
    
    def get_extra_shift(self,target_atom,key_values_1,key_values_2):
        self._check_distance_key(self.ATOM_1,key_values_1.offset,key_values_1.atom)
        self._check_distance_key(self.ATOM_2, key_values_2.offset,key_values_2.atom)
        
        result = None
        
        key_1 = key_values_1.get_atom_key()
        key_2 = key_values_2.get_atom_key()
        if key_1 in self._table[self.DATA]:
            data_1 = self._table[self.DATA][key_values_1.get_atom_key()]
            if key_2 in data_1:
                target_atoms = data_1[key_2]
                if target_atom in target_atoms:
                    result  = target_atoms[target_atom]
            
        return result
        #TODO: move to base
    #TODO: make distance base
    def get_translation_to_table(self,residue,atom):
        result =  atom
        
        key = residue, atom
        if key in self._table[TRANSLATIONS]:
            result =  self._table[TRANSLATIONS][key]
        return result

    def get_translation_from_table(self, residue, atom):
        result = atom
        
        key = residue,atom
        if key in self._translation_from_table:
            result = self._translation_from_table[key]
        return result
