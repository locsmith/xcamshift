#-------------------------------------------------------------------------------
# Copyright (c) 2013 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 30 Dec 2011

@author: garyt
'''

from table_base import Table_base
    
class Distance_table(Table_base):
    '''
    classdocs
    '''


    def __init__(self,table):
        '''
        Constructor
        '''
        
        super(Distance_table, self).__init__(table)
        
        self.DATA = 'data'
        self.EXPONENT = 'exponent'
        self.FROM_ATOMS = 'from_atoms'
        self.TO_ATOMS = 'to_atoms'
        self.OFFSETS='offsets'
        

    def get_offsets(self):
        return self._table[self.OFFSETS]
    
    def get_from_atoms(self):
        return self._table[self.FROM_ATOMS]
    
    def get_target_atoms(self):
        return self.get_from_atoms()
    
    def get_to_atoms(self):
        return self._table[self.TO_ATOMS]
        
    def get_exponent(self):
        return self._table[self.EXPONENT]
    
    
    #TODO: move to base
    #TODO: make distance base
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

