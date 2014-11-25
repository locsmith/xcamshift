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
Created on 17 Jan 2012

@author: garyt
'''

from atom import Atom
from atomSel import AtomSel
from simulation import currentSimulation
from vec3 import norm
from cython.fast_segment_manager import Segment_Manager
from cache_manager import Cache_manager

X = 0
Y = 1
Z = 2
AXES = X,Y,Z

#TODO: look at replacing caching with cython code
#TODO: add SegmentManager to cache infractrusture

RESIDUE_TYPE =  'RESIDUE_TYPE'
def _build_residue_type_cache():
    result = {}
    for atom_id in iter_atom_ids():
        residue_type = Atom_utils._get_residue_type_from_atom_id(atom_id)
        segid,residue_number,atom_name =  Atom_utils._get_atom_info_from_index(atom_id)
        result[segid,residue_number] = residue_type
    return result

Cache_manager.get_cache_manager().add_cache_builder(RESIDUE_TYPE,_build_residue_type_cache)

def _get_residue_type_cache():
    return  Cache_manager.get_cache_manager().get_cache(RESIDUE_TYPE)


ATOM_SELECTION =  'ATOM_SELECTION'
ATOM_ID = 'ATOM_ID'
def _build_empty_cache():
    return {}

Cache_manager.get_cache_manager().add_cache_builder(ATOM_SELECTION,_build_empty_cache)

def _get_atom_selection_cache():
    return  Cache_manager.get_cache_manager().get_cache(ATOM_SELECTION)

Cache_manager.get_cache_manager().add_cache_builder(ATOM_ID,_build_empty_cache)

def _get_atom_id_cache():
    return  Cache_manager.get_cache_manager().get_cache(ATOM_ID)



class Atom_utils(object):
    @staticmethod
    def _calculate_distance(distance_atom_id_1, distance_atom_id_2):

        from_atom_pos = Atom_utils._get_atom_pos(distance_atom_id_1)
        to_atom_pos = Atom_utils._get_atom_pos(distance_atom_id_2)

        xyz_distance = from_atom_pos - to_atom_pos

        distance = norm(xyz_distance)

        return distance

    @staticmethod
    def clear_cache():
        Cache_manager.get_cache_manager().clear_cache()

    @staticmethod
    def find_all_atoms():
        cache = _get_atom_selection_cache()
        result = None
        key = "(all)"
        if key in cache:
            result = cache[key]
        else:
            result = AtomSel(key)
            cache[key]=result
        return result

    @staticmethod
    def find_atom(segment='*', residue_number='#', atom='*'):
        cache = _get_atom_selection_cache()
        result = None
        key = segment, residue_number, atom
        if key in cache:
            result = cache[key]
        else:
            selection = '(segid "%s" and resid %s and name %s)' % key
            result = AtomSel(selection)
            cache[key]=result
        return result

    @staticmethod
    def find_atom_ids(segment='*', residue_number='#', atom='*'):
        atoms  = Atom_utils.find_atom(segment,residue_number,atom)
        return [atom.index() for atom in atoms]

    #TODO what is the difference between this and the called method
    @staticmethod
    def _select_atom_with_translation(segment='*', residue_number='#',atom='*'):
        return Atom_utils.find_atom(segment, residue_number, atom)

    @staticmethod
    def _get_residue_type(segment, residue_number):
        _residue_type_cache = _get_residue_type_cache()

        key = segment, residue_number
        if key in _residue_type_cache:
            result = _residue_type_cache[key]
        else:
            residue_atoms = Atom_utils._select_atom_with_translation(segment, residue_number)
            result = residue_atoms[0].residueName()
            _residue_type_cache[key] = result

        return result

    @staticmethod
    def _get_residue_type_from_atom_id(atom_index):
        atom = Atom_utils._get_atom_by_index(atom_index)
        return atom.residueName()

    @staticmethod
    def _get_atom_by_index(atom_index):
        atom_by_index_cache = _get_atom_id_cache()
        if atom_index in atom_by_index_cache:
            result = atom_by_index_cache[atom_index]
        else:
            result = AtomSel("(id %i)" % (atom_index + 1))[0]
            atom_by_index_cache[atom_index] =  result
        return result

    @staticmethod
    def _get_atom_info_from_index(atom_index):
        atom = Atom_utils._get_atom_by_index(atom_index)
        return atom.segmentName(), atom.residueNum(), atom.atomName()

    @staticmethod
    def _get_atom_name_from_index(atom_index):
        atom = Atom_utils._get_atom_by_index(atom_index)
        return atom.atomName()

    @staticmethod
    def _get_atom_name(atom_index, template="%-5i '%4s' %3i [%3s] %-4s"):
        atom = Atom_utils._get_atom_by_index(atom_index)

        segid = atom.segmentName()
        residue_number = atom.residueNum()
        residue_type = atom.residueName()
        atom_name = atom.atomName()

        return template % (atom_index, segid, residue_number, residue_type, atom_name)

    @staticmethod
    def _format_pretty_atom_name(segid,residue_num,atom_name):
        return "%s:%i@%s" % (segid,residue_num,atom_name)

    @staticmethod
    def _get_pretty_atom_name_from_index(atom_index):
        data = Atom_utils._get_atom_info_from_index(atom_index)
        return Atom_utils._format_pretty_atom_name(*data)

    @staticmethod
    def _get_atom_names(atoms,template  = "%-5i '%4s' %i [%3s] %-4s", joiner=", "):
        atom_names = []
        for atom in atoms:
            atom_names.append(Atom_utils._get_atom_name(atom))
        return joiner.join(atom_names)

    @staticmethod
    def _get_atom_pos(atom_id):
        atom = Atom(currentSimulation(), atom_id)
        atom_pos = atom.pos()
        return atom_pos

    @staticmethod
    def _get_chem_type(atom):
        return atom.chemType()

    @staticmethod
    def _get_bonded_atom_ids(atom_index):
        return currentSimulation().select("bondedto index %i" % atom_index)

def return_true(*args, **kwargs):
    return True


def iter_atoms(predicate=return_true):
    atoms = Atom_utils.find_all_atoms()
    for atom in atoms:
        if predicate(atom):
            yield atom


def iter_atom_ids(predicate=return_true):
    for atom in iter_atoms():
        atom_id = atom.index()
        if predicate(atom_id):
            yield atom_id






def iter_residue_types(predicate=return_true):
    _residue_type_cache = _get_residue_type_cache()
    residue_types = list(set(_residue_type_cache.values()))
    residue_types.sort()
    for residue_type in residue_types:
        if predicate(residue_type):
            yield residue_type


def iter_residue_atoms(residue_predicate = return_true):
    current_residue_atoms = None
    current_residue_segment = None,None

    for atom in iter_atoms():
        residue_number = atom.residueNum()
        segid = atom.segmentName()

        key = segid, residue_number

        if not key == current_residue_segment:

            if current_residue_atoms != None  and residue_predicate(current_residue_atoms):
                yield tuple(current_residue_atoms)

            current_residue_segment = key
            current_residue_atoms = []

        current_residue_atoms.append(atom)

    if current_residue_atoms != None  and residue_predicate(current_residue_atoms):
        yield tuple(current_residue_atoms)



def iter_residue_atom_ids(predicate =  return_true):

    for atoms in iter_residue_atoms():
        atom_ids = tuple([atom.index() for atom in atoms])
        if predicate(atom_ids):
            yield atom_ids

def iter_residues_and_segments(predicate =  return_true):
    segment_manager =  Segment_Manager()

    for segment in segment_manager.get_segments():
        segment_info = segment_manager.get_segment_info(segment)
        for residue_num in range (segment_info.first_residue, segment_info.last_residue+1):
            if predicate(segment,residue_num):
                yield segment,residue_num

