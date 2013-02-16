'''
Created on 12 Aug 2012

@author: garyt
'''
import os
from utils import Atom_utils
from unittest2.case import TestCase
from table_builders.yaml_patches import add_access_to_yaml_list_based_keys

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


def translate_atom_key(atom_key):
    global inverted_translations
    if inverted_translations == None:
        inverted_translations = invert_translations(translations)
    res_type = Atom_utils._get_residue_type(*atom_key[:2])
    key = res_type, atom_key[2]
    if key in inverted_translations:
        atom_key = list(atom_key)
        atom_key[2] = inverted_translations[key]
        atom_key = tuple(atom_key)
    return atom_key

def get_key_for_atom_index(atom_id):

        

    atom_key = Atom_utils._get_atom_info_from_index(atom_id)
    atom_key = translate_atom_key(atom_key)
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

class Yaml_loader(object):
    def __init__(self, files, root=None):
        self._files = files
        self._root = root
    
    def length(self):
        return len(self._files)
    

    def get_item(self, key):
        path = self._add_root_path_to_key(key)
       
        return self._load(path)
    
    #TODO: isolate an test
    def _add_root_path_to_key(self, key):
        path = self._files[key]
        if self._root != None:
            path = os.path.join(self._root, path)
        return path

    def _open_stream(self, file):
        return  open(file, "r")


    def  _load(self, file):
        from yaml import load, dump
        try:
           from yaml import CLoader as Loader, CDumper as Dumper
        except ImportError:
           print 'warning: using slow native python yaml loader'
           from yaml import Loader, Dumper
        add_access_to_yaml_list_based_keys(Loader)

        stream = self._open_stream(file)
        result = load(stream, Loader=Loader)
        stream.close()

        return result 
    
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