'''
Created on 12 Aug 2012

@author: garyt
'''
from utils import Atom_utils
from unittest2.case import TestCase

translations = {('GLY','HA')  :'HA1',
                ('ILE','CD')  : 'CD1',
                ('ILE','HD1') : 'HD11',
                ('ILE','HD2') : 'HD12',
                ('ILE','HD3') : 'HD13'}
inverted_translations  =  None

def invert_translations(translations):
    result  ={}
    
    for key in translations:
        translation = translations[key]
        new_key = key[0],translation
        result[new_key] =  key[1]
    return result

def get_key_for_atom_index( atom_id):
    global inverted_translations
    if inverted_translations == None:
        inverted_translations = invert_translations(translations)
        

    atom_key = Atom_utils._get_atom_info_from_index(atom_id)
    res_type = Atom_utils._get_residue_type_from_atom_id(atom_id)
    key = res_type,atom_key[2]
    if key in inverted_translations:
        atom_key = list(atom_key)
        atom_key[2] = inverted_translations[key]
        atom_key =  tuple(atom_key)
    return atom_key


def get_atom_index(key):
    residue_type = Atom_utils._get_residue_type(key[0], key[1])

    atom_name = key[2]
    
    residue_key = list(key)
    residue_key[2] =  translations.setdefault((residue_type,atom_name),atom_name)
    residue_key = tuple(residue_key)
    
    target_atom_index = Atom_utils.find_atom(*residue_key)[0].index()
    return target_atom_index

def _shift_cache_as_result(shift_cache,target_atom_ids):
    result = []
    for atom_id in target_atom_ids:
        result.append(shift_cache[atom_id])
    return result

DEFAULT_DECIMAL_PLACES =  5
def _check_shift_results(active_target_atoms_ids, result, expected_shifts):
    global DEFAULT_DECIMAL_PLACES
    for i, target_atom_index in enumerate(active_target_atoms_ids):
        key = get_key_for_atom_index(target_atom_index)
        expected_shift_diff = expected_shifts[key]
        if abs(result[i] - expected_shift_diff) >  1.0**-(DEFAULT_DECIMAL_PLACES-2):
            raise AssertionError("diff too big %f - %f " % (result[i], expected_shift_diff))

class Empty_loader(object):
    def __init__(self):
        pass
    def length(self):
        return 0
    
    def get_item(self, key):
        raise IndexError("index out of range; this sequence is always empty")
    

class Virtual_list(object):
    def __init__(self, loader):
        self._loader = loader
        
    def __len__(self):
        return self._loader.length()
    
    def __getitem__(self, key):
        return self._loader.get_item(key)