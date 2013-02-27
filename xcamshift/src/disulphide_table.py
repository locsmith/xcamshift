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

from table_base import Table_base, TRANSLATIONS

DISU = 'DISU'
CYS = 'CYS'
        

DATA = 'data'
TARGET_ATOMS = 'target_atoms'        
RESIDUE_TYPE_OFFSET = 1
class Disulphide_table(Table_base):
    
    def __init__(self, table):
        super(Disulphide_table, self).__init__(table)

    def get_target_atoms(self):
        return self._table[TARGET_ATOMS]
    
    def get_residue_types(self):
        return [entry[RESIDUE_TYPE_OFFSET] for entry in self._table[DATA]]
#        
    def get_disulphide_shift(self, atom):
        self.__check_keys(atom)
        
        result = None
        

        key = (DISU, CYS)
        atom_disulphide_table  = self._table[DATA][key]
        if atom in atom_disulphide_table:
            result =  atom_disulphide_table[atom]
            
        return result
    
    def __check_keys(self, atom):
        
        if atom not in self.get_target_atoms():
            message = "atom %s not in disulphide_shift_table atoms %s" % (atom,', '.join(self.get_target_atoms()))
            raise KeyError(message)
    

