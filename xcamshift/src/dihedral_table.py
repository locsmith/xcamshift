'''
Created on 30 Dec 2011

@author: garyt
'''
from keys import Atom_key

class Dihedral_table_base(object):
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
    
    def get_dihedral_keys(self):
        keys = self.table[self.DATA].keys()
        
        result = []
        for index in range(len(keys)):
            atom_keys = [Atom_key(offset,atom) for atom,offset in keys[index]]
            result.append(tuple(atom_keys))
            
        return result
        
class Dihedral_table(Dihedral_table_base):
    

    def __init__(self,table):   
        super(Dihedral_table, self).__init__(table)    
            
    
    def get_dihedral_shift(self,target_atom,dihedral_key):
        self._check_target_atom(target_atom)
        self._check_dihedral_key(dihedral_key)
        
        result = None
        
        raw_dihedral_key = dihedral_key.get_dihedral_key()

        if raw_dihedral_key in self.table[self.DATA]:
            data = self.table[self.DATA][raw_dihedral_key]
            if target_atom in data:
                result  = data[target_atom]
#            
        return result

class Dihedral_parameter_table(Dihedral_table_base):
    
    PARAMETERS='parameters'
    

    def __init__(self,table):   
        super(Dihedral_parameter_table, self).__init__(table)


    def get_parameters(self):
        return self.table[self.PARAMETERS]
    
    
    def check_parameter_key(self, parameter):
        parameters = self.get_parameters()
        if not parameter in parameters:
            message = "parameter %s not in parameters %s" % (parameter,', '.join(parameters))
            raise KeyError(message)
    
    
    def get_dihedral_parameter(self,target_atom,dihedral_key,parameter):
        self._check_target_atom(target_atom)
        self._check_dihedral_key(dihedral_key)
        self.check_parameter_key(parameter)
        
        result = None
        
        raw_dihedral_key = dihedral_key.get_dihedral_key()

        if raw_dihedral_key in self.table[self.DATA]:
            data = self.table[self.DATA][raw_dihedral_key]
            if parameter in data:
                parameter = data[parameter]
                if target_atom in parameter:
                    result  = parameter[target_atom]
#            
        return result    