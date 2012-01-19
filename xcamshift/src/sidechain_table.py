'''
Created on 19 Jan 2012

@author: garyt
'''


class Sidechain_table(object):


    TARGET_ATOMS = "target_atoms"
    EXPONENT = "exponent"
    SIDECHAIN_ATOMS = "sidechain_atoms"
    DATA = 'data'
    
    def __init__(self, table):
        self.table = table
    
    def get_exponent(self):
        return self.table[self.EXPONENT]
    
    def get_target_atoms(self):
        return self.table[self.TARGET_ATOMS]
    
    def get_residue_types(self):
        return self.table[self.SIDECHAIN_ATOMS].keys()
    
    def get_sidechain_atoms(self,residue_type):
        self._check_residue_type(residue_type)
        
        return self.table[self.DATA][residue_type].keys()
    

    def _check_target_atom(self, atom_name):
        atoms = self.get_target_atoms()
        if atom_name not in atoms:
            template = "atom_name %s is not in sidechain shift table target atoms (%s)"
            message = template % (atom_name, ', '.join(atoms))
            raise KeyError(message)

    
    
    def _check_residue_type(self, residue_type):
        residue_types = self.get_residue_types()
        if residue_type not in residue_types    :
            template = "residue type %s not in sidechain shift table residue types (%s)"
            message = template % (residue_type, ', '.join(residue_types))
            raise KeyError(message)
    
    
    def _check_sidechain_atom(self, residue_type, atom_type):
        sidechain_atoms = self.get_sidechain_atoms(residue_type)
        if atom_type not in sidechain_atoms:
            template = "sidechain atom type %s for residue type %s is not in sidechain shift table atoms (%s)"
            message = template % (atom_type, residue_type, ', '.join(sidechain_atoms))
            raise KeyError(message)        
    
    
    def get_sidechain_coefficient(self,residue_type,target_atom,sidechain_atom):
        self._check_target_atom(target_atom)
        self._check_sidechain_atom(residue_type, sidechain_atom)
        
        return self.table[self.DATA][residue_type][sidechain_atom][target_atom]
    
#        self.ATOMS = 'atoms'
#        self.OFFSETS='offsets'