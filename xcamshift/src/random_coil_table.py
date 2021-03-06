#-------------------------------------------------------------------------------
# Copyright (c) 2013-2015 Gary Thompson & The University of Leeds.
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
from table_base import Table_base, TRANSLATIONS

        
class Random_coil_table(Table_base):
    
    DATA = 'data'
    def __init__(self, table):
        super(Random_coil_table, self).__init__(table)
        
        self.ATOMS = 'atoms'
        self.OFFSETS='offsets'
        
    def get_offsets(self):
        return self._table[self.OFFSETS]
    
    def get_atoms(self):
        return self._table[self.ATOMS]
    
    def get_target_atoms(self):
        return self.get_atoms()
    
    def __check_residue_keys(self, offset, residue, atom):
        
        if not offset in self.get_offsets():
            offset_strings = [`offset` for offset in self.get_offsets()]
            message_values = offset, ', '.join(offset_strings)
            message = "offset %s is not in %s" % message_values
            raise KeyError(message)
        
        if atom not in self.get_atoms():
            message = "atom %s not in random coil shift _table atoms %s" % (atom,', '.join(self.get_atoms()))
            raise KeyError(message)
        
    def get_random_coil_shift(self,offset, residue, atom):
        self.__check_residue_keys(offset,residue,atom)
        
        result = None
        
        atom_random_coil_table  = self._table[self.DATA][offset][residue]
        if atom in atom_random_coil_table:
            result =  atom_random_coil_table[atom]
            
        return result
    
    def get_translation(self,residue,atom):
        result =  atom
        
        key = residue, atom
        if key in self._table[TRANSLATIONS]:
            result =  self._table[TRANSLATIONS][key]
        return result
