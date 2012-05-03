'''
Created on 30 Dec 2011

@author: garyt
'''

TRANSLATIONS = 'translations'

class Distance_table(object):
    '''
    classdocs
    '''


    def __init__(self,table):
        '''
        Constructor
        '''
        self._table = table
        self.DATA = 'data'
        self.EXPONENT = 'exponent'
        self.FROM_ATOMS = 'from_atoms'
        self.TO_ATOMS = 'to_atoms'
        self.OFFSETS='offsets'
        
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
    
    def get_offsets(self):
        return self._table[self.OFFSETS]
    
    def get_from_atoms(self):
        return self._table[self.FROM_ATOMS]
    
    def get_to_atoms(self):
        return self._table[self.TO_ATOMS]
        
    def get_exponent(self):
        return self._table[self.EXPONENT]
    
    
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
    

    def __check_distance_coefficient_keys(self, to_atom, offset, from_atom):
        if not from_atom in self.get_from_atoms():
            message_values = from_atom, ', '.join(self.get_from_atoms())
            message = "from atom key %s is not in %s" % message_values
            raise KeyError(message)
        
        if not offset in self.get_offsets():
            message_values = offset, ', '.join(self.get_offsets())
            message = "offset %s is not in %s" % message_values
            raise KeyError(message)
    
        if not to_atom in self.get_to_atoms():
            message_values = to_atom, ', '.join(self.get_to_atoms())
            message = "to atom key %s is not in %s" % message_values
            raise KeyError(message)
        
        
    def get_distance_coeeficent(self,from_atom,offset,to_atom):
        '''
            get the distance coefficient from atom_i to the the atom_j 
            in the residue at offset offset
        '''
        
        self.__check_distance_coefficient_keys(to_atom,offset,from_atom)
        
        return self._table[self.DATA][offset][to_atom][from_atom]

