'''
Created on 17 Jan 2012

@author: garyt
'''

from atom import Atom
from atomSel import AtomSel
from simulation import currentSimulation
from vec3 import norm
from cython.fast_segment_manager import Segment_Manager
from table_manager import Table_manager


X = 0
Y = 1
Z = 2
AXES = X,Y,Z


#TODO: add more general caching mechanism IMPORTANT and add as part of force field
#TODO: 
cache = {}
atom_by_index_cache = {}
_residue_type_cache = None
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
        global cache
        global _residue_type_cache
        _residue_type_cache =  None
        cache = {}
        atom_by_index_cache = {}
        Table_manager.get_default_table_manager().reset_default_table_manager()
        
    @staticmethod  
    def find_all_atoms():
        global cache
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
        global cache
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
        
        
        if residue_number in _residue_type_cache:
            result = _residue_type_cache[residue_number]
        else:
            residue_atoms = Atom_utils._select_atom_with_translation(segment, residue_number)
            result = residue_atoms[0].residueName()
            _residue_type_cache[residue_number] = result
            
        return result
    
    @staticmethod
    def _get_residue_type_from_atom_id(atom_index):
        atom = Atom_utils._get_atom_by_index(atom_index)
        return atom.residueName()

    @staticmethod
    def _get_atom_by_index(atom_index):
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
        




def _get_residue_type_cache():
    global _residue_type_cache
    
    if _residue_type_cache == None:
        _residue_type_cache = {}
        for atom_id in iter_atom_ids():
            residue_type = Atom_utils._get_residue_type_from_atom_id(atom_id)
            residue_number = Atom_utils._get_atom_info_from_index(atom_id)[1]
            _residue_type_cache[residue_number] = residue_type
    return _residue_type_cache


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
    
