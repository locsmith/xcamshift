'''
Created on 17 Jan 2012

@author: garyt
'''
from simulation import currentSimulation
from atom import Atom
from vec3 import  norm
from atomSel import AtomSel

def tupleit(t):
    return tuple(map(tupleit, t)) if isinstance(t, (list, tuple)) else t

class Atom_utils(object):
    @staticmethod
    def _calculate_distance(distance_atom_id_1, distance_atom_id_2):
        
        from_atom_pos = Atom_utils._get_atom_pos(distance_atom_id_1)
        to_atom_pos = Atom_utils._get_atom_pos(distance_atom_id_2)
        
        xyz_distance = from_atom_pos - to_atom_pos
        
        distance = norm(xyz_distance)
        
        return distance
    
    @staticmethod
    def find_atom(segment='*', residue_number='#', atom='*'):
        selection = '(segid "%s" and resid %i and name %s)' % (segment, int(residue_number), atom)
        residue_atoms = AtomSel(selection)
        return residue_atoms

    @staticmethod
    def find_atom_id(segment='*', residue_number='#', atom='*'):
        atoms  = Atom_utils.find_atom(segment,residue_number,atom)
        return [atom.index() for atom in atoms]
    
    #TODO what is the difference between this and the called method
    @staticmethod
    def _select_atom_with_translation(segment='*', residue_number='#',atom='*'):
        residue_atoms = Atom_utils.find_atom(segment, residue_number, atom)
        return residue_atoms

    @staticmethod
    def _get_residue_type(segment, residue_number):
        residue_atoms = Atom_utils._select_atom_with_translation(segment, residue_number)
        return residue_atoms[0].residueName()
    
    @staticmethod
    def _get_residue_type_from_atom_id(atom_index):
        atom = Atom_utils._get_atom_by_index(atom_index)
        return atom.residueName()

    @staticmethod
    def _get_atom_by_index(atom_index):
        return AtomSel("(id %i)" % (atom_index + 1))[0]
    
    @staticmethod
    def _get_atom_info_from_index(atom_index):
        atom = Atom_utils._get_atom_by_index(atom_index)
        return atom.segmentName(), atom.residueNum(), atom.atomName()
    
    @staticmethod
    def _get_atom_name_from_index(atom_index):
        atom = Atom_utils._get_atom_by_index(atom_index)
        return atom.atomName()
    
    @staticmethod
    def _get_atom_name(atom_index, template="%-5i '%4s' %i [%3s] %-4s"):
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
