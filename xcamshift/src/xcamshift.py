'''
Created on 27 Dec 2011

@author: garyt
'''
from atomSel import AtomSel

from protocol import initStruct
from pdbTool import PDBTool
from textwrap import dedent
from segment_manager import Segment_Manager
from Table_manager import Table_manager
import sys
from atom import Atom

        

class Distance_potential(object):
    '''
    classdocs
    '''
    distance_set = {
        (-1,'CA')  : {'HA': 0.1049538},
        ( 0,'CA')  : 2,
        ( 1,'CA')  : 3  
    }

    def __init__(self):
        '''
        Constructor
        '''
        
        self.distance_index = {}
        self.create_distances("(all)")
    
    def print_atom(self,atom):
        return "%i. [%s]:%i[%s]@%s" % (atom.index(),atom.segmentName(), atom.residueNum(),atom.residueName(), atom.atomName())
        
    def create_distances(self,selection):
        seg_residue_atom = {}
        
        atom_selection = AtomSel(selection)
        
        for atom in atom_selection:
            key = (atom.segmentName(),atom.residueNum(),atom.atomName())
            seg_residue_atom[key]=atom
        
        for atom in atom_selection:
            
            segment = atom.segmentName()
            residue_num = atom.residueNum()
            for distance in self.distance_set:
                offset, atom_name = distance
                
                key = (segment,residue_num+offset,atom_name)
                if key in seg_residue_atom:
                    print self.print_atom(atom),self.print_atom(seg_residue_atom[key])




class RandomCoilShifts(object):
    
    def __init__(self):
        
        self.__segment_manager = Segment_Manager()
        self.__table_manager = Table_manager.get_default_table_manager()
        self.__shifts_list = self.__create_shift_list('(all)')
        

    def __select_atoms(self, segment='*', residue_number='#',atom='*'):
        selection = '(segid "%s" and resid %i and name %s)' % (segment, residue_number, atom)
        residue_atoms = AtomSel(selection)
        return residue_atoms


    def __get_residue_type(self, segment, residue_number):
        residue_atoms = self.__select_atoms(segment, residue_number)
        return residue_atoms[0].residueName()

    def __create_shift_list(self,atom_selection):
        
        result  = []
        
        random_coil_table = self.__table_manager.get_random_coil_table()
        
        atom_selection = AtomSel(atom_selection)
        for segment in self.__segment_manager.get_segments():
            segment_info = self.__segment_manager.get_segment_info(segment)
            
            if segment_info.segment_length < 3:
                template = "Warning: segment '%s' only contains < 3 residues (%i) and will be ignored"
                message = template % (segment,segment_info.segment_length)
                print >> sys.stderr, message
                continue
            
            for residue_number in range(segment_info.first_residue+1,segment_info.last_residue):
        
                residue_type = self.__get_residue_type(segment, residue_number)
                for atom in self.__select_atoms(segment, residue_number):
                    atom_name = atom.atomName()
                    atom_index = atom.index()
                    value  = random_coil_table.get_random_coil_shift(residue_type,atom_name)
                    
                    if value != None:
                        result.append((atom_index,value))
        return result
        
    def set_shifts(self,shift_list):
        for atom_index,shift in self.__shifts_list:
            shift_list[atom_index] = shift
            

    def __get_atom_name(self, atom_index, template="%-5i '%4s' %i [%3s] %-4s"):
        atom = AtomSel("(id %i)" % (atom_index+1))[0]
        
        segid = atom.segmentName()
        residue_number = atom.residueNum()
        residue_type = atom.residueName()
        atom_name = atom.atomName()
        
        return template % (atom_index, segid, residue_number, residue_type, atom_name)
    
    
    def __str__(self): 
        result = []
        for atom_index,shift in self.__shifts_list:
            template = "%s %7.3f"
            atom_name =  self.__get_atom_name(atom_index)
            string = template % (atom_name,shift)
            result.append(string)
        return '\n'.join(result)
    
if __name__ == "__main__":
    initStruct("test_data/3_ala/3ala.psf")
    PDBTool("test_data/3_ala/3ala.pdb").read()
    RandomCoilShifts()
#    Distance_potential()
    
        