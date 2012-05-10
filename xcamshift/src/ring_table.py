'''
Created on 19 Jan 2012

@author: garyt
'''

TRANSLATIONS = 'translations'

class Ring_table(object):


    TARGET_ATOMS = "target_atoms"
    RINGS = 'rings'
    DATA = 'data'
    
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

        
    def get_target_atoms(self):
        return self._table[self.TARGET_ATOMS]
    
    def get_residue_types(self):
        return self._table[self.RINGS].keys()
    
    def get_rings(self,residue_type,ring_type):
        self._check_residue_type(residue_type)
        self._check_ring_type(residue_type, ring_type)
        
        return self._table[self.DATA][residue_type][ring_type].keys()
    
    def get_ring_atoms(self,residue_type,ring_type):
        self._check_residue_type(residue_type)
        self._check_ring_type(residue_type, ring_type)
        
        return self._table[self.RINGS][residue_type][ring_type]
    
    def get_ring_types(self,residue_type):
        self._check_residue_type(residue_type)
        
        return self._table[self.RINGS][residue_type].keys()
        
    def _check_ring_type(self,residue_type,ring_type):
        self._check_residue_type(residue_type)
        
        ring_types = self.get_ring_types(residue_type)
        
        if not ring_type in ring_types:
            template = "ring type %s is not in ring types %s for residue (%s)"
            message = template % (ring_type, ', '.join(ring_types),residue_type)
            raise KeyError(message)
        
#        TODO: much in here adn elsewhere can be pushed to a super class
    def _check_target_atom(self, atom_name):
        atoms = self.get_target_atoms()
        if atom_name not in atoms:
            template = "atom_name %s is not in ring shift _table target atoms (%s)"
            message = template % (atom_name, ', '.join(atoms))
            raise KeyError(message)

    
    
    def _check_residue_type(self, residue_type):
        residue_types = self.get_residue_types()
        if residue_type not in residue_types    :
            template = "residue type %s not in sidechain shift _table residue types (%s)"
            message = template % (residue_type, ', '.join(residue_types))
            raise KeyError(message)
    
    
    def _check_ring_atom(self, residue_type, atom_type):
        sidechain_atoms = self.get_rings_atoms(residue_type)
        if atom_type not in sidechain_atoms:
            template = "sidechain atom type %s for residue type %s is not in rung shift _table atoms (%s)"
            message = template % (atom_type, residue_type, ', '.join(sidechain_atoms))
            raise KeyError(message)        
    
    
    def get_ring_coefficient(self,target_atom,residue_type,ring_type):
        self._check_target_atom(target_atom)
        self._check_residue_type(residue_type)
        self._check_ring_type(residue_type, ring_type)
        
        return self._table[self.DATA][residue_type,ring_type][target_atom]*1000

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
