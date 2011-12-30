'''
Created on 30 Dec 2011

@author: garyt
'''

class Distance_table(object):
    '''
    classdocs
    '''


    def __init__(self,table):
        '''
        Constructor
        '''
        self.table = table
        self.DATA = 'data'
        self.EXPONENT = 'exponent'
        self.FROM_ATOMS = 'from_atoms'
        self.TO_ATOMS = 'to_atoms'
        self.OFFSETS='offsets'
        
    def get_offsets(self):
        return self.table[self.OFFSETS]
    
    def get_from_atoms(self):
        return self.table[self.FROM_ATOMS]
    
    def get_to_atoms(self):
        return self.table[self.TO_ATOMS]
        
    def get_exponent(self):
        return self.table[self.EXPONENT]
    

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
        
        return self.table[self.DATA][offset][to_atom][from_atom]

