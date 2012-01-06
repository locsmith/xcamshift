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
from atomSel import intersection
from vec3 import  norm
import sys
from bzrlib.osutils import stat

        
class Base_potential(object):
    ALL = '(all)'
            
    def __init__(self):
        self._segment_manager = Segment_Manager()
        self._table_manager = Table_manager.get_default_table_manager()
    
    @staticmethod
    def _select_atoms(segment='*', residue_number='#',atom='*'):
        selection = '(segid "%s" and resid %i and name %s)' % (segment, residue_number, atom)
        residue_atoms = AtomSel(selection)
        return residue_atoms

    @staticmethod
    def _get_residue_type(segment, residue_number):
        residue_atoms = Base_potential._select_atoms(segment, residue_number)
        return residue_atoms[0].residueName()
    
    @staticmethod
    def _get_atom_name(atom_index, template="%-5i '%4s' %i [%3s] %-4s"):
        atom = AtomSel("(id %i)" % (atom_index+1))[0]
        
        segid = atom.segmentName()
        residue_number = atom.residueNum()
        residue_type = atom.residueName()
        atom_name = atom.atomName()
        
        return template % (atom_index, segid, residue_number, residue_type, atom_name)
    
    @staticmethod
    def _get_atom_pos(atom_id):
        atom = Atom(currentSimulation(), atom_id)
        atom_pos = atom.pos()
        return atom_pos
#    def _print_atom(self,atom):
#        return "%i. [%s]:%i[%s]@%s" % (atom.index(),atom.segmentName(), atom.residueNum(),atom.residueName(), atom.atomName())

    def _check_segment_length(self, segment):
        segment_info = self._segment_manager.get_segment_info(segment)
        if segment_info.segment_length < 3:
            template = "Warning: segment '%s' only contains < 3 residues (%i) and will be ignored"
            message = template % (segment, segment_info.segment_length)
            print >> sys.stderr, message
        return segment_info
        
class Distance_potential(Base_potential):
    '''
    classdocs
    '''


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

    



    def set_shifts(self, result):
        
        for from_atom_id,to_atom_id,coefficent,exponent in self._distances:
            
            from_atom_pos = self._get_atom_pos(from_atom_id)
            to_atom_pos = self._get_atom_pos(to_atom_id)
            
            xyz_distance = from_atom_pos - to_atom_pos
            distance  = norm(xyz_distance)
            
            shift = distance ** exponent * coefficent 
            result[from_atom_id] += shift
            
        return shift

class ExtraDistanceShifts(Base_potential):
    def __init__(self):
        Base_potential.__init__(self)
        
        self._shifts_list =  self._create_shift_list()

    def _create_shift_list(self, all):
        pass
    
    

class RandomCoilShifts(Base_potential):
    

    def __init__(self):
        super(RandomCoilShifts, self).__init__()

        #TODO bad name
        self._shifts_list = self._create_shift_list(self.ALL)


    def _get_table(self, from_residue_type):
        return self._table_manager.get_random_coil_table(from_residue_type)


    def _get_coefficents(self, atom, context):
        value = None
        atom_name = atom.atomName()
        if atom_name in context.table.get_atoms():
            value = context.table.get_random_coil_shift(context.offset, context.to_residue_type, atom_name)
        return value
    
    class ResidueOffsetContext :
        def __init__(self,atom,offset,table):
            
            self.table = table
            
            self.segment = atom.segmentName()
            
            self.offset = offset
            
            self.from_residue_number = atom.residueNum()
            self.from_residue_type = Base_potential._get_residue_type(self.segment, self.from_residue_number)
            
            self.to_residue_number = self.from_residue_number+1
            self.to_residue_type = Base_potential._get_residue_type(self.segment, self.from_residue_number + offset)


    def _build_contexts(self, atom, table):
        contexts = []
        for offset in table.get_offsets():
            context = RandomCoilShifts.ResidueOffsetContext(atom, offset, table)
            contexts.append(context)
        return contexts

    def add_components_for_residue(self, segment, target_residue_number, atom_selection):
        result  = []
        
        from_residue_type = self._get_residue_type(segment, target_residue_number)
        random_coil_table = self._get_table(from_residue_type)
        selected_atoms = intersection(self._select_atoms(segment, target_residue_number), atom_selection)
        
        for atom in selected_atoms:
            
            contexts = self._build_contexts(atom, random_coil_table)
            
            for context in contexts:
                
                value = None
                value = self._get_coefficents(atom, context)
                    
                if value != None:
                    atom_index = atom.index()
                    result.append((atom_index, value))
                    
        return result




    def _create_shift_list(self,global_atom_selection):
        
        result  = []
        
        
        global_atom_selection = AtomSel(global_atom_selection)
        for segment in self._segment_manager.get_segments():
            
            if self._check_segment_length(segment):
                segment_info = self._segment_manager.get_segment_info(segment)
                for residue_number in range(segment_info.first_residue+1,segment_info.last_residue):
                    residue_atom_selection = self._select_atoms(segment, residue_number)
                    target_atom_selection = intersection(residue_atom_selection,global_atom_selection)
                    result.extend(self.add_components_for_residue(segment, residue_number, target_atom_selection))
                
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
    
        