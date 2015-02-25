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
Created on 24 Jan 2012

@author: garyt
'''
from python_utils import tupleit
from utils import Atom_utils
from  cython.shift_calculators import allocate_array
from UserDict import DictMixin

class Observed_shift_table(DictMixin):
    '''
    classdocs
    '''

    def __init__(self,shift_data={},format='map'):
        if format ==  'map':
            self._chemical_shifts = self._process_observed_shift_map(shift_data)
        elif format == 'xplor':
            self._chemical_shifts =  self._process_xplor_shifts(shift_data)
        else:
            raise Exception("unexpected shift format %s" % format)

        self._native_shifts = None

    def _process_xplor_shifts(self,shift_data):
        result = {}
        for datum in shift_data:
            result[datum.atom_id] = datum.shift
        return result

    def _process_observed_shift_map(self,shift_data):
        result  = {}
        for key in shift_data:
            if len(key) == 2:
                search_key = '*',key[0],key[1]
            else:
                search_key = tuple(key)
            if len(search_key) ==  3:
                #TODO move Base_potential components to utilities
                atoms = Atom_utils.find_atom(*search_key)
                for atom in atoms:
                    result[atom.index()] = shift_data[key]
            else:
                template = """key with unexpected length should be either 'segid,residue_no,atom_name'
                              or 'residue_no,atom_name' but i got '%s'"""
                msg = template % `key`
                raise Exception(msg)

        return result

    def add_shifts(self,shift_table):
        for atom_id in shift_table:
            self._chemical_shifts[atom_id] =  shift_table[atom_id]
        self._native_shifts = None


    def get_chemical_shift(self,atom_index):
        return self._chemical_shifts[atom_index]

    def get_atom_indices(self):
        return self._chemical_shifts.keys()

    def get_indices_for_atom_id(self):
        result ={}
        for i, atom_id in enumerate(self._chemical_shifts.keys()):
            result[atom_id] = i
        return result

    #TODO add more generic dump capabilities
    def dump_observed_shifts(self):
        results = []
        for atom_index in self._chemical_shifts:
            sub_result  = []
            results.append(sub_result)
            sub_result.append(Atom_utils._get_atom_info_from_index(atom_index))
            sub_result.append(self._chemical_shifts[atom_index])
        return tupleit(results)

    def get_native_shifts(self, target_atom_ids):
        if self._native_shifts == None:
            self._native_shifts = allocate_array(len(target_atom_ids),'f')
            for i,target_atom_id in enumerate(target_atom_ids):
                self._native_shifts[i] = self.get_chemical_shift(target_atom_id)

        return self._native_shifts

    def __len__(self):
        return len(self._chemical_shifts)

    def __delitem__(self,index):
        del self._chemical_shifts[index]
        self._native_shifts = None

    def __contains__(self, atom_id):
        return atom_id in self._chemical_shifts.keys()

    def keys(self):
        return self._chemical_shifts.keys()

    def __getitem__(self,index):
        result = None
        if index in self._chemical_shifts:
            result  = self._chemical_shifts[index]
        else:
            raise KeyError('shift with atom id %i not found' % index)
        return result

    def __setitem__(self,index,value):
        self._chemical_shifts[index]=value

    def __str__(self):

        result  = ["shift table", "-----------",""]
        for (segid,resid,atom_name),shift in self.dump_observed_shifts():
            result.append("[%4s]:%i@%s  %7.4f" % (segid,resid,atom_name,shift))
        return "\n".join(result)

