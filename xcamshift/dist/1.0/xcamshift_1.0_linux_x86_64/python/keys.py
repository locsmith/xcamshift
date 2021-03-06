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
Created on 13 Jan 2012

@author: garyt
'''

class Atom_key(object):
    def __init__(self, offset,atom):
        self.offset = offset
        self.atom = atom
        
    def __str__(self):
        message = "atom key (offset: %1i, atom: %s)"
        return message % (self.offset,self.atom)
    
    def get_atom_key(self):
        return self.atom,self.offset
    
class Dihedral_key(object):

    MAX_INDEX = 3
    def __init__(self, *atom_keys):
        self._atom_keys=[Atom_key(atom_key.offset, atom_key.atom)  for atom_key in atom_keys]
    
    def get_atom_key(self,index):
        if index > self.MAX_INDEX:
            raise IndexError("dihedrals only contain 4 atoms you requested atom %d" % index +1)
        return self._atom_keys[index]
        
    def get_keys(self):
        return tuple(self._atom_keys)
    
    def get_dihedral_key(self):
        return tuple([atom_key.get_atom_key() for atom_key in self._atom_keys])
