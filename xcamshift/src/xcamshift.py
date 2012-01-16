'''
Created on 27 Dec 2011

@author: garyt
TODO need to translate from_atom names
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
import abc 
from keys import Atom_key, Dihedral_key
        
class Base_potential(object):
    
    __metaclass__ = abc.ABCMeta
    
    ALL = '(all)'
            
    def __init__(self):
        self._segment_manager = Segment_Manager()
        self._table_manager = Table_manager.get_default_table_manager()
        
    @staticmethod
    def _calculate_distance(distance_atom_id_1, distance_atom_id_2):
        
        from_atom_pos = Base_potential._get_atom_pos(distance_atom_id_1)
        to_atom_pos = Base_potential._get_atom_pos(distance_atom_id_2)
        
        xyz_distance = from_atom_pos - to_atom_pos
        
        distance = norm(xyz_distance)
        
        return distance
    
    @staticmethod
    def _select_atom_with_translation(segment='*', residue_number='#',atom='*'):
        selection = '(segid "%s" and resid %i and name %s)' % (segment, residue_number, atom)
        residue_atoms = AtomSel(selection)
        return residue_atoms

    @staticmethod
    def _get_residue_type(segment, residue_number):
        residue_atoms = Base_potential._select_atom_with_translation(segment, residue_number)
        return residue_atoms[0].residueName()
    

    @staticmethod
    def _get_atom_by_index(atom_index):
        return AtomSel("(id %i)" % (atom_index + 1))[0]
    
    @staticmethod
    def _get_atom_info_from_index(atom_index):
        atom = Base_potential._get_atom_by_index(atom_index)
        return atom.segmentName(), atom.residueNum(), atom.atomName()
            
    @staticmethod
    def _get_atom_name(atom_index, template="%-5i '%4s' %i [%3s] %-4s"):
        atom = Base_potential._get_atom_by_index(atom_index)
        
        segid = atom.segmentName()
        residue_number = atom.residueNum()
        residue_type = atom.residueName()
        atom_name = atom.atomName()
        
        return template % (atom_index, segid, residue_number, residue_type, atom_name)
    
    @staticmethod
    def _get_atom_names(atoms,template  = "%-5i '%4s' %i [%3s] %-4s", joiner=", "):
        atom_names = []
        for atom in atoms:
            atom_names.append(Base_potential._get_atom_name(atom))
        return joiner.join(atom_names)
            
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
    
    def add_components_for_residue(self, segment, target_residue_number, atom_selection):
        result  = []
        
        from_residue_type = self._get_residue_type(segment, target_residue_number)
        random_coil_table = self._get_table(from_residue_type)
        selected_atoms = intersection(self._select_atom_with_translation(segment, target_residue_number), atom_selection)
        
        for atom in selected_atoms:
            
            contexts = self._build_contexts(atom, random_coil_table)
            
            for context in contexts:
                
                value = None
                if context.complete:
                    value = self._get_component_for_atom(atom, context)
                        
                    if value != None:
                        result.append(value)
                    
        return result
    
    def _create_component_list(self,global_atom_selection):
        
        result  = []
        
        
        global_atom_selection = AtomSel(global_atom_selection)
        for segment in self._segment_manager.get_segments():
            
            if self._check_segment_length(segment):
                segment_info = self._segment_manager.get_segment_info(segment)
                for residue_number in range(segment_info.first_residue+1,segment_info.last_residue):
                    residue_atom_selection = self._select_atom_with_translation(segment, residue_number)
                    target_atom_selection = intersection(residue_atom_selection,global_atom_selection)
                    result.extend(self.add_components_for_residue(segment, residue_number, target_atom_selection))
                
        return result
    
    @abc.abstractmethod
    def _get_component_for_atom(self, atom, context):
        pass
    
    @abc.abstractmethod
    def _build_contexts(self, atom, table):
        contexts = []
        for offset in table.get_offsets():
            context = RandomCoilShifts.ResidueOffsetContext(atom, offset, table)
            contexts.append(context)
        return contexts

    @abc.abstractmethod
    def set_shifts(self,shift_list):
        pass
    
    @abc.abstractmethod
    def _get_table(self, from_residue_type):
        pass


    @abc.abstractmethod
    def _translate_atom_name(self,atom_name):
        return atom_name
        
class Distance_potential(Base_potential):
    '''
    classdocs
    '''


    def __init__(self):
        super(Distance_potential, self).__init__()
        
        '''
        Constructor
        '''
        
        self._distances = self._create_component_list("(all)")
        
    class ResidueAtomOffsetContext:

        
        
        def __init__(self, from_atom, offset, to_atom_name ,table):
            
            self.complete = False
            self.table = table
            
            self.segment = from_atom.segmentName()
            
            self.offset = offset
            
#            self.from_atom_name = from_atom.atomName()
#            self.from_atom_name = self.table.get_translation(self.from_atom_name)
#            self.from_atom_index = from_atom.index()
            from_residue_number = from_atom.residueNum()
#            self.from_residue_type = Base_potential._get_residue_type(self.segment, self.from_residue_number)
            
            self.to_atom_name = to_atom_name
            self.to_residue_number = from_residue_number+offset
            to_atom = Base_potential._select_atom_with_translation(self.segment, self.to_residue_number, self.to_atom_name)
            
            if len(to_atom) == 0:
                self.to_atom_name =  self.table.get_translation(self.to_atom_name)
                to_atom = Base_potential._select_atom_with_translation(self.segment, self.to_residue_number, self.to_atom_name)
                
            num_to_atom = len(to_atom)
            if num_to_atom > 1:
                self._get_atom_names(to_atom)
                raise Exception("unexpected number of to atoms selected (> 1) %d" % num_to_atom)
                
            if len(to_atom) == 1:
                self.to_residue_type = Base_potential._get_residue_type(self.segment, self.to_residue_number)
                self.to_atom_index = to_atom[0].index()
                self.complete = True
            
            
    def _translate_atom_name(self, atom_name, context):
        return context.table.get_translation(atom_name)
    
    def _build_contexts(self, atom, table):
        contexts = []
        for offset in table.get_offsets():
            for to_atom_name in table.get_to_atoms():
                context = Distance_potential.ResidueAtomOffsetContext(atom,offset,to_atom_name,table)
                if context.complete:
                    contexts.append(context)
        return contexts
                
    def _get_component_for_atom(self, atom, context):
        table = context.table
        
        from_atom_name = atom.atomName()
        from_atom_name = self._translate_atom_name(from_atom_name, context)
        offset = context.offset
        to_atom_name = context.to_atom_name
        to_atom_name = self._translate_atom_name(to_atom_name,context)
        
        
        result = None
        if from_atom_name in table.get_from_atoms():
            if to_atom_name in table.get_to_atoms():
                value = context.table.get_distance_coeeficent(from_atom_name,offset,to_atom_name)
                if value != None:
                    from_atom_index = atom.index()
                    to_atom_index = context.to_atom_index
                    exponent = context.table.get_exponent()
                    result = (from_atom_index,to_atom_index,value,exponent)
        return result

    def _get_table(self, from_residue_type):
        return self._table_manager.get_BB_Distance_Table(from_residue_type)

#    def _create_component_list(self,atom_selection):
#        
#        result  = []
#        
#        atom_selection = AtomSel(atom_selection)
#        for segment in self._segment_manager.get_segments():
#            segment_info = self._segment_manager.get_segment_info(segment)
#            
#            if segment_info.segment_length < 3:
#                template = "Warning: segment '%s' only contains < 3 residues (%i) and will be ignored"
#                message = template % (segment,segment_info.segment_length)
#                print >> sys.stderr, message
#                continue
#            
#            for from_residue_number in range(segment_info.first_residue+1,segment_info.last_residue):
#                residue_type = self._get_residue_type(segment, from_residue_number)
#        
#                distance_table = self._table_manager.get_BB_Distance_Table(residue_type)
#                
#                exponent = distance_table.get_exponent()
#                from_atoms  =  distance_table.get_from_atoms()
#                
#                for offset in distance_table.get_offsets():
#                    for from_atom in self._select_atom_with_translation(segment, from_residue_number):
#                        from_atom_name = from_atom.atomName()
#                        from_atom_index = from_atom.index()
#                        if from_atom_name in from_atoms:
#                            for to_atom_name in distance_table.get_to_atoms():
#                                to_residue_number =  from_residue_number+offset
#                                to_atoms =  self._select_atom_with_translation(segment,to_residue_number, to_atom_name)
#                                for to_atom in to_atoms:
#                                    value = None
#                                    to_atom_name = to_atom.atomName()
#                                    to_atom_index =  to_atom.index()
#                                    value = distance_table.get_distance_coeeficent(from_atom_name,offset,to_atom_name)
#    #                            value  = distance_table.get_random_coil_shift(residue_type,atom_name)
#                                
#                                    if value != None:
#                                        result.append((from_atom_index,to_atom_index,value,exponent))
#                return result
            
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

class Extra_potential(Base_potential):
    def __init__(self):
        Base_potential.__init__(self)
        
        self._distances =  self._create_component_list(self.ALL)
    
    def _translate_atom_name(self, atom_name,context):
        return context.table.get_translation(atom_name)
        
    def _get_table(self, residue_type):
        return self._table_manager.get_extra_table(residue_type)
    
    class ExtraContext(object):

        def _select_atom_with_translation(self, segment, residue_number_1, atom_name_1):
            target_atom_1 = Base_potential._select_atom_with_translation(segment, residue_number_1, atom_name_1)
            if len(target_atom_1) == 0:
                atom_name_1 = self.table.get_translation(atom_name_1)
                target_atom_1 = Base_potential._select_atom_with_translation(segment, residue_number_1, atom_name_1)
            num_to_atom = len(target_atom_1)
            if num_to_atom > 1:
                self._get_atom_names(target_atom_1)
                raise Exception("unexpected number of to atoms selected (> 1) %d" % num_to_atom)
            return target_atom_1

        def __init__(self, from_atom, key_1, key_2 ,table):
            
            self.complete = False
            self.table = table
            
            segment = from_atom.segmentName()
            
            residue_number_1 = from_atom.residueNum() + key_1.offset
            atom_name_1 = key_1.atom
            target_atom_1 = self._select_atom_with_translation(segment, residue_number_1, atom_name_1)
            
            atom_name_2 = key_2.atom
            residue_number_2 = from_atom.residueNum() + key_2.offset
            target_atom_2 =self._select_atom_with_translation(segment, residue_number_2, atom_name_2)
            
            if len(target_atom_1) == 1 and len(target_atom_2) == 1:
                self.distance_atom_index_1 = target_atom_1[0].index()
                self.distance_atom_index_2 = target_atom_2[0].index()
                self.key_1=key_1
                self.key_2=key_2
                self.complete = True


        
    def _build_contexts(self, atom, table):
        contexts = []
        for offset_1 in table.get_offsets(table.ATOM_1):
            for offset_2 in table.get_offsets(table.ATOM_2):
                for distance_atom_1 in table.get_distance_atoms(table.ATOM_1):
                    for distance_atom_2 in table.get_distance_atoms(table.ATOM_2):
                        key_1 = Atom_key(offset_1,distance_atom_1)
                        key_2 = Atom_key(offset_2,distance_atom_2)
                        context = Extra_potential.ExtraContext(atom,key_1,key_2,table)
                        if context.complete:
                            contexts.append(context)
        return contexts
    
    def  _get_component_for_atom(self, atom, context):
        table = context.table

        
        from_atom_name = atom.atomName()
        from_atom_name = self._translate_atom_name(from_atom_name, context)
        
        result = None
        if from_atom_name in table.get_target_atoms():
            value = context.table.get_extra_shift(from_atom_name,context.key_1,context.key_2)
#            print self._get_atom_name(atom.index()),value
            if value != None:
                from_atom_index = atom.index()
                distance_index_1 = context.distance_atom_index_1
                distance_index_2 = context.distance_atom_index_2
                exponent = context.table.get_exponent()
                result = (from_atom_index,distance_index_1,distance_index_2,value,exponent)
        return result
    
    def __str__(self):
        print len( self._distances)
        result = []
        for from_index,distance_index_1,distance_index_2,value,exponent in self._distances:
            from_atom = self._get_atom_name(from_index)
            distance_atom_1 = self._get_atom_name(distance_index_1)
            distance_atom_2 = self._get_atom_name(distance_index_2)
            
            template = '[%s] - [%s] - [%s] %7.3f %7.3f'
            values = from_atom, distance_atom_1, distance_atom_2, value, exponent
            result.append( template % values)
        return '\n'.join(result)

    def dump(self):
        result  = []
        for from_index,distance_index_1,distance_index_2,value,exponent in self._distances:
            sub_result  = []
            sub_result.extend(self._get_atom_info_from_index(from_index)[1:])
            
            sub_result.extend(self._get_atom_info_from_index(distance_index_1)[1:])
            sub_result.extend(self._get_atom_info_from_index(distance_index_2)[1:])
            
            sub_result.append(value)
            
            result.append(tuple(sub_result))
        return result
    



    def _calc_single_shift(self, index):
        distance_atom_id_1, distance_atom_id_2, coefficient, exponent = self._distances[index][1:]
        
        distance = self._calculate_distance(distance_atom_id_1, distance_atom_id_2)
        
        shift = distance ** exponent * coefficient
        
        return shift

    def set_shifts(self,result):
        for index in range(len(self._distances)):
            shift = self._calc_single_shift(index)
            from_atom_id = self._distances[index][0]
            result[from_atom_id] += shift

        return result
    
class RandomCoilShifts(Base_potential):
    

    def __init__(self):
        super(RandomCoilShifts, self).__init__()

        #TODO bad name
        self._shifts_list = self._create_component_list(self.ALL)

    def _translate_atom_name(self, atom_name):
        #print >> sys.stderr, "WARNING no atom name translations in radom coil shifts yet!"
        return super(RandomCoilShifts, self)._translate_atom_name(atom_name)
    
    def _get_table(self, from_residue_type):
        return self._table_manager.get_random_coil_table(from_residue_type)


    def _get_component_for_atom(self, atom, context):
        result = None
        atom_name = atom.atomName()
        atom_name = self._translate_atom_name(atom_name)
        if atom_name in context.table.get_atoms():
            value = context.table.get_random_coil_shift(context.offset, context.to_residue_type, atom_name)
            if value != None:
                result = atom.index(), value
        return result
    
    class ResidueOffsetContext :
        def __init__(self,atom,offset,table):
            
            self.table = table
            
            segment = atom.segmentName()
            
            self.offset = offset
            
            from_residue_number = atom.residueNum()
            
            self.to_residue_number = from_residue_number + offset
            self.to_residue_type = Base_potential._get_residue_type(segment, from_residue_number + offset)
            
            self.complete = True


    def _build_contexts(self, atom, table):
        contexts = []
        for offset in table.get_offsets():
            context = RandomCoilShifts.ResidueOffsetContext(atom, offset, table)
            contexts.append(context)
        return contexts
        
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

class Dihedral_potential(Base_potential):

    def __init__(self):
        Base_potential.__init__(self)
        
        self._distances =  self._create_component_list(self.ALL)
    
    def _translate_atom_name(self, atom_name,context):
        return context.table.get_translation(atom_name)
        
    def _get_table(self, residue_type):
        return self._table_manager.get_dihedral_table(residue_type)
#    
    class DihedralContext(object):

        def _select_atom_with_translation(self, segment, residue_number_1, atom_name_1):
            target_atom_1 = Base_potential._select_atom_with_translation(segment, residue_number_1, atom_name_1)
            if len(target_atom_1) == 0:
                atom_name_1 = self.table.get_translation(atom_name_1)
                target_atom_1 = Base_potential._select_atom_with_translation(segment, residue_number_1, atom_name_1)
            num_to_atom = len(target_atom_1)
            if num_to_atom > 1:
                self._get_atom_names(target_atom_1)
                raise Exception("unexpected number of to atoms selected (> 1) %d" % num_to_atom)
            return target_atom_1


        def get_atom_info(self, segment, from_atom, key_1):
            residue_number_1 = from_atom.residueNum() + key_1.offset
            atom_name_1 = key_1.atom
            target_atom_1 = self._select_atom_with_translation(segment, residue_number_1, atom_name_1)
            return target_atom_1

        def __init__(self, from_atom, dihedral_key ,table):
            
            self.table = table
            
            segment = from_atom.segmentName()
            
            target_atoms = []
            for key in dihedral_key.get_keys():
                target_atoms.append(self.get_atom_info(segment, from_atom, key))
            
            self.complete=True
            for target_atom in target_atoms:
                if len(target_atom) !=  1:
                    self.complete =  False
                    
            if self.complete:
                self.dihedral_atom_index_pair_1 = (target_atoms[0][0].index(),
                                                   target_atoms[1][0].index())
                self.distance_atom_index_pair_2 = (target_atoms[2][0].index(),
                                                   target_atoms[3][0].index())
                
                self.dihedral_key = dihedral_key
        
    def _build_contexts(self, atom, table):
        contexts = []
        
        for key in table.get_dihedral_keys():
            dihedral_key = Dihedral_key(*key)
            context = Dihedral_potential.DihedralContext(atom,dihedral_key,table)
            
            if context.complete:
                contexts.append(context)
                
        return contexts


#    
    def  _get_component_for_atom(self, atom, context):
        pass
#        table = context.table
#
#        
#        from_atom_name = atom.atomName()
#        from_atom_name = self._translate_atom_name(from_atom_name, context)
#        
#        result = None
#        if from_atom_name in table.get_target_atoms():
#            value = context.table.get_extra_shift(from_atom_name,context.key_1,context.key_2)
##            print self._get_atom_name(atom.index()),value
#            if value != None:
#                from_atom_index = atom.index()
#                distance_index_1 = context.distance_atom_index_1
#                distance_index_2 = context.distance_atom_index_2
#                exponent = context.table.get_exponent()
#                result = (from_atom_index,distance_index_1,distance_index_2,value,exponent)
#        return result
#    
#    def __str__(self):
#        print len( self._distances)
#        result = []
#        for from_index,distance_index_1,distance_index_2,value,exponent in self._distances:
#            from_atom = self._get_atom_name(from_index)
#            distance_atom_1 = self._get_atom_name(distance_index_1)
#            distance_atom_2 = self._get_atom_name(distance_index_2)
#            
#            template = '[%s] - [%s] - [%s] %7.3f %7.3f'
#            values = from_atom, distance_atom_1, distance_atom_2, value, exponent
#            result.append( template % values)
#        return '\n'.join(result)
#
#    def dump(self):
#        result  = []
#        for from_index,distance_index_1,distance_index_2,value,exponent in self._distances:
#            sub_result  = []
#            sub_result.extend(self._get_atom_info_from_index(from_index)[1:])
#            
#            sub_result.extend(self._get_atom_info_from_index(distance_index_1)[1:])
#            sub_result.extend(self._get_atom_info_from_index(distance_index_2)[1:])
#            
#            sub_result.append(value)
#            
#            result.append(tuple(sub_result))
#        return result
#    
#
#
#
#    def _calc_single_shift(self, index):
#        distance_atom_id_1, distance_atom_id_2, coefficient, exponent = self._distances[index][1:]
#        
#        distance = self._calculate_distance(distance_atom_id_1, distance_atom_id_2)
#        
#        shift = distance ** exponent * coefficient
#        
#        return shift
#
    def set_shifts(self,result):
        pass
#        for index in range(len(self._distances)):
#            shift = self._calc_single_shift(index)
#            from_atom_id = self._distances[index][0]
#            result[from_atom_id] += shift
#
#        return result

if __name__ == "__main__":
    initStruct("test_data/3_ala/3ala.psf")
    PDBTool("test_data/3_ala/3ala.pdb").read()
#    RandomCoilShifts()
    Distance_potential()
    
        