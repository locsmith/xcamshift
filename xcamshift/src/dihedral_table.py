'''
Created on 30 Dec 2011

@author: garyt
'''

class Dihedral_key(object):

    MAX_INDEX = 3
    def __init__(self, atom_key_1, atom_key_2,atom_key_3,atom_key_4):
        self.atom_keys=[atom_key_1, atom_key_2,atom_key_3,atom_key_4]
    
    def get_key(self,index):
        if index > self.MAX_INDEX:
            raise IndexError("dihedrals only contain 4 atoms you requested atom %d" % index +1)
        
    def get_keys(self):
        return tuple(self.atom_keys)
        
class Dihedral_table(object):
    
    DATA = 'data'

    EXPONENT = 'exponent'
    OFFSETS='offsets'
    TARGET_ATOMS='target_atoms'
    DIHEDRAL_ATOMS='dihedral_atoms'
    TRANSLATIONS='translations'
    
    def __init__(self, table):
        self.table = table
    
    def get_translation(self,atom):
        return atom
#    
    def get_dihedral_atoms(self):
        return self.table[self.DIHEDRAL_ATOMS]
    
    def get_offsets(self):
        return self.table[self.OFFSETS]
    
    def get_target_atoms(self):
        return self.table[self.TARGET_ATOMS]

    def get_exponent(self):
        return self.table[self.EXPONENT]

    def _check_offset(self, offset):
        
        offsets = self.get_offsets()
        
        if not offset in offsets:
            message = "offset %s is not in %s" % (offset,", ".join(offsets))
            raise KeyError(message)
    
    def _check_dihedral_atom(self, atom_name):
        
        atoms = self.get_dihedral_atoms()
        if atom_name not in atoms:
            message = "atom_name %s not in dihedral shift table atoms (%s)" % (atom_name, ', '.join(atoms))
            raise KeyError(message)

    def _check_dihedral_key(self,dihedral_key):
        
        for key in dihedral_key.get_keys():
            self._check_dihedral_atom(key.atom)
            self._check_offset(key.offset)
            
    def _check_atom_key(self, atom_index, offset, atom_name):
        
        self._check_offset(atom_index, offset)
        
        self._check_extra_atom(atom_index, atom_name)
    
    def _check_target_atom(self, atom):
        
        atoms = self.get_target_atoms()
        if not atom in atoms:
            message = "atom %s not in target atoms %s" % (atom,', '.join(atoms))
            raise KeyError(message)
    
    
    def get_dihedral_shift(self,target_atom,dihedral_key):
        self._check_target_atom(target_atom)
        self._check_dihedral_key(dihedral_key)
        
        result = None
        
        raw_dihedral_key = dihedral_key.get_keys()
        if raw_dihedral_key in self.table[self.DATA]:
            data = self.table[self.DATA][raw_dihedral_key]
            if target_atom in data:
                result  = data[target_atom]
#            
        return result