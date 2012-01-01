'''
Created on 27 Dec 2011

@author: garyt
'''
from atomSel import AtomSel

from protocol import initStruct
from pdbTool import PDBTool
from segment_manager import Segment_Manager
from Table_manager import Table_manager
from simulation import currentSimulation
from atom import Atom
from vec3 import  norm
import sys

        
class Base_potential(object):
    def __init__(self):
        self._segment_manager = Segment_Manager()
        self._table_manager = Table_manager.get_default_table_manager()
        
    def _select_atoms(self, segment='*', residue_number='#',atom='*'):
        selection = '(segid "%s" and resid %i and name %s)' % (segment, residue_number, atom)
        residue_atoms = AtomSel(selection)
        return residue_atoms


    def _get_residue_type(self, segment, residue_number):
        residue_atoms = self._select_atoms(segment, residue_number)
        return residue_atoms[0].residueName()
    
    
    def _get_atom_name(self, atom_index, template="%-5i '%4s' %i [%3s] %-4s"):
        atom = AtomSel("(id %i)" % (atom_index+1))[0]
        
        segid = atom.segmentName()
        residue_number = atom.residueNum()
        residue_type = atom.residueName()
        atom_name = atom.atomName()
        
        return template % (atom_index, segid, residue_number, residue_type, atom_name)
    
#    def _print_atom(self,atom):
#        return "%i. [%s]:%i[%s]@%s" % (atom.index(),atom.segmentName(), atom.residueNum(),atom.residueName(), atom.atomName())
        
class Distance_potential(Base_potential):
    '''
    classdocs
    '''
    distance_set = {
        (-1,'CA')  : {'HA': 0.1049538},
        ( 0,'CA')  : 2,
        ( 1,'CA')  : 3  
    }

    def __init__(self):
        super(Distance_potential, self).__init__()
        
        '''
        Constructor
        '''
        
        self._distances = self.create_distances("(all)")
    

    def create_distances(self,atom_selection):
        
        result  = []
        
        atom_selection = AtomSel(atom_selection)
        for segment in self._segment_manager.get_segments():
            segment_info = self._segment_manager.get_segment_info(segment)
            
            if segment_info.segment_length < 3:
                template = "Warning: segment '%s' only contains < 3 residues (%i) and will be ignored"
                message = template % (segment,segment_info.segment_length)
                print >> sys.stderr, message
                continue
            
            for from_residue_number in range(segment_info.first_residue+1,segment_info.last_residue):
                residue_type = self._get_residue_type(segment, from_residue_number)
        
                distance_table = self._table_manager.get_BB_Distance_Table(residue_type)
                
                exponent = distance_table.get_exponent()
                from_atoms  =  distance_table.get_from_atoms()
                
                for offset in distance_table.get_offsets():
                    for from_atom in self._select_atoms(segment, from_residue_number):
                        from_atom_name = from_atom.atomName()
                        from_atom_index = from_atom.index()
                        if from_atom_name in from_atoms:
                            for to_atom_name in distance_table.get_to_atoms():
                                to_residue_number =  from_residue_number+offset
                                to_atoms =  self._select_atoms(segment,to_residue_number, to_atom_name)
                                for to_atom in to_atoms:
                                    value = None
                                    to_atom_name = to_atom.atomName()
                                    to_atom_index =  to_atom.index()
                                    value = distance_table.get_distance_coeeficent(from_atom_name,offset,to_atom_name)
    #                            value  = distance_table.get_random_coil_shift(residue_type,atom_name)
                                
                                    if value != None:
                                        result.append((from_atom_index,to_atom_index,value,exponent))
                return result
            
    def __str__(self):
        print len( self._distances)
        result = []
        for from_index,to_index,value,exponent in self._distances:
            from_atom = self._get_atom_name(from_index)
            to_atom = self._get_atom_name(to_index)
            
            template = '[%s] - [%s] %7.3f %7.3f'
            values = from_atom, to_atom, value, exponent
            result.append( template % values)
        return '\n'.join(result)

    

    def _get_atom_pos(self, atom_id):
        atom = Atom(currentSimulation(), atom_id)
        atom_pos = atom.pos()
        return atom_pos

    def set_shifts(self, result):
        
        for from_atom_id,to_atom_id,coefficent,exponent in self._distances:
            
            from_atom_pos = self._get_atom_pos(from_atom_id)
            to_atom_pos = self._get_atom_pos(to_atom_id)
            
            xyz_distance = from_atom_pos - to_atom_pos
            distance  = norm(xyz_distance)
            
            shift = distance * coefficent ** exponent
            print >> sys.stderr, self._get_atom_name(from_atom_id), self._get_atom_name(to_atom_id)
            print >> sys.stderr, shift, coefficent, exponent,distance,shift
            result[from_atom_id] += shift
            
        print >> sys.stderr
        return shift
#    
#        seg_residue_atom = {}
#        
#        atom_selection = AtomSel(atom_selection)
#        
#        for atom in atom_selection:
#            key = (atom.segmentName(),atom.residueNum(),atom.atomName())
#            seg_residue_atom[key]=atom
#        
#        for atom in atom_selection:
#            
#            segment = atom.segmentName()
#            residue_num = atom.residueNum()
#            for distance in self.distance_set:
#                offset, atom_name = distance
#                
#                key = (segment,residue_num+offset,atom_name)
#                if key in seg_residue_atom:
#                    print self.print_atom(atom),self.print_atom(seg_residue_atom[key])




class RandomCoilShifts(Base_potential):
    
    def __init__(self):
        super(RandomCoilShifts, self).__init__()

        self._shifts_list = self._create_shift_list('(all)')
        

    def _create_shift_list(self,atom_selection):
        
        result  = []
        
        random_coil_table = self._table_manager.get_random_coil_table()
        
        atom_selection = AtomSel(atom_selection)
        for segment in self._segment_manager.get_segments():
            segment_info = self._segment_manager.get_segment_info(segment)
            
            if segment_info.segment_length < 3:
                template = "Warning: segment '%s' only contains < 3 residues (%i) and will be ignored"
                message = template % (segment,segment_info.segment_length)
                print >> sys.stderr, message
                continue
            
            for residue_number in range(segment_info.first_residue+1,segment_info.last_residue):
        
                residue_type = self._get_residue_type(segment, residue_number)
                for atom in self._select_atoms(segment, residue_number):
                    atom_name = atom.atomName()
                    atom_index = atom.index()
                    value  = random_coil_table.get_random_coil_shift(residue_type,atom_name)
                    
                    if value != None:
                        result.append((atom_index,value))
        return result
        
    def set_shifts(self,shift_list):
        for atom_index,shift in self._shifts_list:
            shift_list[atom_index] = shift
            

    def __str__(self): 
        result = []
        for atom_index,shift in self._shifts_list:
            template = "%s %7.3f"
            atom_name =  self._get_atom_name(atom_index)
            string = template % (atom_name,shift)
            result.append(string)
        return '\n'.join(result)
    
if __name__ == "__main__":
    initStruct("test_data/3_ala/3ala.psf")
    PDBTool("test_data/3_ala/3ala.pdb").read()
#    RandomCoilShifts()
    Distance_potential()
    
        