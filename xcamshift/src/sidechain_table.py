#-------------------------------------------------------------------------------
# Copyright (c) 2013 gary thompson.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v2.1
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 19 Jan 2012

@author: garyt
'''

DATA = 'data'
TRANSLATIONS = 'translations'

class Sidechain_table(object):


    TARGET_ATOMS = "target_atoms"
    EXPONENT = "exponent"
    SIDECHAIN_ATOMS = "sidechain_atoms"
    

    def __init__(self, table):
        self._table = table
        self._sidechain_atoms = self._build_sidechain_atoms(table)

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

    def _build_sidechain_atoms(self, table):
        result = {}
        for residue_key in self._table[DATA]:
            result[residue_key] = self._table[DATA][residue_key].keys()
            
        return result
    
    def get_exponent(self):
        return self._table[self.EXPONENT]
    
    def get_target_atoms(self):
        return self._table[self.TARGET_ATOMS]
    
    def get_residue_types(self):
        return self._sidechain_atoms.keys()
    
    def get_sidechain_atoms(self,residue_type):
        self._check_residue_type(residue_type)
        
        return self._sidechain_atoms[residue_type]
    

    def _check_target_atom(self, atom_name):
        atoms = self.get_target_atoms()
        if atom_name not in atoms:
            template = "atom_name %s is not in sidechain shift _table target atoms (%s)"
            message = template % (atom_name, ', '.join(atoms))
            raise KeyError(message)

    
    
    def _check_residue_type(self, residue_type):
        residue_types = self.get_residue_types()
        if residue_type not in residue_types    :
            template = "residue type %s not in sidechain shift _table residue types (%s)"
            message = template % (residue_type, ', '.join(residue_types))
            raise KeyError(message)
    
    
    def _check_sidechain_atom(self, residue_type, atom_type):
        sidechain_atoms = self.get_sidechain_atoms(residue_type)
        if atom_type not in sidechain_atoms:
            template = "sidechain atom type %s for residue type %s is not in sidechain shift _table atoms (%s)"
            message = template % (atom_type, residue_type, ', '.join(sidechain_atoms))
            raise KeyError(message)        
    
    
    def get_sidechain_coefficient(self,residue_type,target_atom,sidechain_atom):
        self._check_target_atom(target_atom)
        self._check_sidechain_atom(residue_type, sidechain_atom)
        
        return self._table[DATA][residue_type][sidechain_atom][target_atom]
    
#        self.ATOMS = 'atoms'
#        self.OFFSETS='offsets'

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
