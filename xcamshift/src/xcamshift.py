'''
Created on 27 Dec 2011

@author: garyt
TODO need to translate from_atom names
'''



from atomSel import AtomSel, intersection
from dihedral import Dihedral
from keys import Atom_key, Dihedral_key
from math import cos, tanh
from observed_chemical_shifts import Observed_shift_table
from pdbTool import PDBTool
from protocol import initStruct
from segment_manager import Segment_Manager
from table_manager import Table_manager
from utils import tupleit, Atom_utils
from vec3 import norm
import abc
import sys

        
class Base_potential(object):
    
    __metaclass__ = abc.ABCMeta
    
    ALL = '(all)'
            
    def __init__(self):
        self._segment_manager = Segment_Manager()
        self._table_manager = Table_manager.get_default_table_manager()
        self._observed_shifts = Observed_shift_table()
        
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
        
        from_residue_type = Atom_utils._get_residue_type(segment, target_residue_number)
        random_coil_table = self._get_table(from_residue_type)
        selected_atoms = intersection(Atom_utils._select_atom_with_translation(segment, target_residue_number), atom_selection)
        
        for atom in selected_atoms:
            
            contexts = self._build_contexts(atom, random_coil_table)
            
            for context in contexts:
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
                    residue_atom_selection = Atom_utils._select_atom_with_translation(segment, residue_number)
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
        
    @abc.abstractmethod
    def get_abbreviated_name(self):
        pass
    #TODO should be abc
    def set_observed_shifts(self, shift_table):
        self._observed_shifts = shift_table

class Distance_based_potential(Base_potential):
    
    class Indices(object):
        def __init__(self, target_atom_index,distance_atom_index_1,
                     distance_atom_index_2, coefficent_index, exponent_index):
            self.target_atom_index =  target_atom_index
            self.distance_atom_index_1 = distance_atom_index_1
            self.distance_atom_index_2 =  distance_atom_index_2
            self.exponent_index = exponent_index
            self.coefficient_index = coefficent_index
            
        def __str__(self):
            result = 'indices target= %i distance atom index 1 = %i index 2 = %i coefficent = %i exponent = %i'
            msg = result % (self.target_atom_index, self.distance_atom_index_1,
                            self.distance_atom_index_2, self.coefficient_index, self.exponent_index)
            return msg                 
            
#    @abc.abstractmethod
    def _get_indices(self):
        return Distance_based_potential.Indices(target_atom_index=0,distance_atom_index_1=0,
                                                distance_atom_index_2=1,coefficent_index=2,
                                                exponent_index=3)
    
#    @abc.abstractmethod
    def _get_data_table(self):
        pass
    
    def _calc_single_force_factor(self,index,factor):
        values  = self._get_data_table()[index]
        
        indices = self._get_indices()
        
        distance_atom_index_1 = indices.distance_atom_index_1
        distance_atom_index_2 =indices.distance_atom_index_2
        coefficient_index = indices.coefficient_index
        exponent_index  = indices.exponent_index
        
        
        target_atom=values[distance_atom_index_1]
        distance_atom = values[distance_atom_index_2]
        coefficient  =  values[coefficient_index]
        exponent =  values[exponent_index]
        
        target_pos = Atom_utils._get_atom_by_index(target_atom).pos()
        distant_pos =  Atom_utils._get_atom_by_index(distance_atom).pos()
        
        distance  = target_pos - distant_pos
        distance_2 = sum([elem**2 for elem in distance])

        factor= factor * coefficient
        modified_exponent = (exponent - 2.0) / 2.0
        force_factor = factor *  exponent * distance_2**modified_exponent

        return force_factor

    def _calc_single_force_set(self,index,factor, forces):
        values  = self._get_data_table()[index]
        
        indices = self._get_indices()
        
        distance_atom_index_1 = indices.distance_atom_index_1
        distance_atom_index_2 =indices.distance_atom_index_2
        
        target_atom=values[distance_atom_index_1]
        distance_atom = values[distance_atom_index_2]
        
        target_pos = Atom_utils._get_atom_by_index(target_atom).pos()
        distant_pos =  Atom_utils._get_atom_by_index(distance_atom).pos()
        
        distance  = target_pos - distant_pos
        
        force_factor  = self._calc_single_force_factor(index, factor)
        
        DIM_3 = 3
        target_offset = target_atom * DIM_3
        distant_offset = distance_atom * DIM_3
        X_OFFSET = 0
        Y_OFFSET = 1
        Z_OFFSET = 2
        
        OFFSETS_3 = (X_OFFSET,Y_OFFSET,Z_OFFSET)
        
        for offset in OFFSETS_3:
#            print distance
#            print forces
            forces[target_offset+offset] -= distance[offset] * force_factor
            forces[distant_offset+offset] += distance[offset] * force_factor
        
        return forces
    
class Distance_potential(Distance_based_potential):
    '''
    classdocs
    '''


    def __init__(self):
        super(Distance_potential, self).__init__()
        
        '''
        Constructor
        '''
        
        self._distances = self._create_component_list("(all)")

    def _get_data_table(self):
        return self._distances
    
    def _get_indices(self):
        return super(Distance_potential, self)._get_indices()
       
    def get_abbreviated_name(self):
        return "BB  "
    
    class ResidueAtomOffsetContext:

        
        
        def __init__(self, from_atom, offset, to_atom_name ,table):
            
            self.complete = False
            self._table = table
            
            self.segment = from_atom.segmentName()
            
            self.offset = offset
            
#            self.from_atom_name = from_atom.atomName()
#            self.from_atom_name = self._table.get_translation(self.from_atom_name)
#            self.from_atom_index = from_atom.index()
            from_residue_number = from_atom.residueNum()
#            self.from_residue_type = Base_potential._get_residue_type(self.segment, self.from_residue_number)
            
            self.to_atom_name = to_atom_name
            self.to_residue_number = from_residue_number+offset
            to_atom = Atom_utils._select_atom_with_translation(self.segment, self.to_residue_number, self.to_atom_name)
            
            if len(to_atom) == 0:
                self.to_atom_name =  self._table.get_translation(self.to_atom_name)
                to_atom = Atom_utils._select_atom_with_translation(self.segment, self.to_residue_number, self.to_atom_name)
                
            num_to_atom = len(to_atom)
            if num_to_atom > 1:
                self._get_atom_names(to_atom)
                raise Exception("unexpected number of to atoms selected (> 1) %d" % num_to_atom)
                
            if len(to_atom) == 1:
                self.to_residue_type = Atom_utils._get_residue_type(self.segment, self.to_residue_number)
                self.to_atom_index = to_atom[0].index()
                self.complete = True
            
            
    def _translate_atom_name(self, atom_name, context):
        return context._table.get_translation(atom_name)
    
    def _build_contexts(self, atom, table):
        contexts = []
        for offset in table.get_offsets():
            for to_atom_name in table.get_to_atoms():
                context = Distance_potential.ResidueAtomOffsetContext(atom,offset,to_atom_name,table)
                if context.complete:
                    contexts.append(context)
        return contexts
                
    def _get_component_for_atom(self, atom, context):
        table = context._table
        
        from_atom_name = atom.atomName()
        from_atom_name = self._translate_atom_name(from_atom_name, context)
        offset = context.offset
        to_atom_name = context.to_atom_name
        to_atom_name = self._translate_atom_name(to_atom_name,context)
        
        
        result = None
        if from_atom_name in table.get_from_atoms():
            if to_atom_name in table.get_to_atoms():
                value = context._table.get_distance_coeeficent(from_atom_name,offset,to_atom_name)
                if value != None:
                    from_atom_index = atom.index()
                    to_atom_index = context.to_atom_index
                    if from_atom_index != to_atom_index:
                        exponent = context._table.get_exponent()
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

    

    def _calc_single_shift(self,i):
        from_atom_id,to_atom_id,coefficent,exponent = self._distances[i]
        from_atom_pos = Atom_utils._get_atom_pos(from_atom_id)
        to_atom_pos = Atom_utils._get_atom_pos(to_atom_id)
        
        xyz_distance = from_atom_pos - to_atom_pos
        distance  = norm(xyz_distance)
        
        return  distance ** exponent * coefficent 
    

#        print target_atom_info, distant_atom_info,distance_atom,distance 
#        for elem in 
#        
#        pass
#    #TODO move to Base_potential
    def set_shifts(self, result):
        for index in range(len(self._distances)):
            shift = self._calc_single_shift(index)
            from_atom_id = self._distances[index][0]
            result[from_atom_id] += shift
#        for from_atom_id,to_atom_id,coefficent,exponent in self._distances:
#            
#            from_atom_pos = Atom_utils._get_atom_pos(from_atom_id)
#            to_atom_pos = Atom_utils._get_atom_pos(to_atom_id)
#            
#            xyz_distance = from_atom_pos - to_atom_pos
#            distance  = norm(xyz_distance)
#            
#            shift = distance ** exponent * coefficent 
#            result[from_atom_id] += shift
            
        return result
    

    
#    def get_shift(self):
#        pass
#    def get_energies(self, data):
#        pass
#    
#    def get_derivatives(self,data):
#        pass
#

    
    def dump(self):
        result  = []
        for from_index,distance_index_1,value,exponent in self._distances:
            sub_result  = []
            sub_result.append(Atom_utils._get_atom_info_from_index(from_index))
            
            sub_result.append(Atom_utils._get_atom_info_from_index(distance_index_1))
            
            sub_result.append(value)
            
            sub_result.append(exponent)
            
            result.append(tuple(sub_result))
        return result   

    

        
#        self.set_shifts(result)
    
    
    
    
    
    

class Extra_potential(Distance_based_potential):
    def __init__(self):
        Base_potential.__init__(self)
        
        self._distances =  self._create_component_list(self.ALL)
    
    def get_abbreviated_name(self):
        return "XTRA"

    def _translate_atom_name(self, atom_name,context):
        return context._table.get_translation(atom_name)
        
    def _get_table(self, residue_type):
        return self._table_manager.get_extra_table(residue_type)
    
    def _get_data_table(self):
        return self._distances
    
    class ExtraContext(object):

        def _select_atom_with_translation(self, segment, residue_number_1, atom_name_1):
            target_atom_1 = Atom_utils._select_atom_with_translation(segment, residue_number_1, atom_name_1)
            if len(target_atom_1) == 0:
                atom_name_1 = self._table.get_translation(atom_name_1)
                target_atom_1 = Atom_utils._select_atom_with_translation(segment, residue_number_1, atom_name_1)
            num_to_atom = len(target_atom_1)
            if num_to_atom > 1:
                self._get_atom_names(target_atom_1)
                raise Exception("unexpected number of to atoms selected (> 1) %d" % num_to_atom)
            return target_atom_1

        def __init__(self, from_atom, key_1, key_2 ,table):
            
            self.complete = False
            self._table = table
            
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


    def _get_indices(self):
        return Distance_based_potential.Indices(target_atom_index=0,distance_atom_index_1=1,
                                                distance_atom_index_2=2,coefficent_index=3,
                                                exponent_index=4)    
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
        table = context._table

        
        from_atom_name = atom.atomName()
        from_atom_name = self._translate_atom_name(from_atom_name, context)
        
        result = None
        if from_atom_name in table.get_target_atoms():
            value = context._table.get_extra_shift(from_atom_name,context.key_1,context.key_2)
#            print self._get_atom_name(atom.index()),value
            if value != None:
                from_atom_index = atom.index()
                distance_index_1 = context.distance_atom_index_1
                distance_index_2 = context.distance_atom_index_2
                exponent = context._table.get_exponent()
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
            sub_result.append(Atom_utils._get_atom_info_from_index(from_index))
            
            sub_result.append(Atom_utils._get_atom_info_from_index(distance_index_1))
            sub_result.append(Atom_utils._get_atom_info_from_index(distance_index_2))
            
            sub_result.append(value)
            
            result.append(tuple(sub_result))
        return result
    



    def _calc_single_shift(self, index):
        distance_atom_id_1, distance_atom_id_2, coefficient, exponent = self._distances[index][1:]
        
        distance = Atom_utils._calculate_distance(distance_atom_id_1, distance_atom_id_2)
        
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
    
    def get_abbreviated_name(self):
        return "RC  "
    
    def _translate_atom_name(self, atom_name):
        #print >> sys.stderr, "WARNING no atom name translations in radom coil shifts yet!"
        return super(RandomCoilShifts, self)._translate_atom_name(atom_name)
    
    def _get_table(self, from_residue_type):
        return self._table_manager.get_random_coil_table(from_residue_type)


    def _get_component_for_atom(self, atom, context):
        result = None
        atom_name = atom.atomName()
        atom_name = self._translate_atom_name(atom_name)
        if atom_name in context._table.get_atoms():
            value = context._table.get_random_coil_shift(context.offset, context.to_residue_type, atom_name)
            if value != None:
                result = atom.index(), value
        return result
    
    class ResidueOffsetContext :
        def __init__(self,atom,offset,table):
            
            self._table = table
            
            segment = atom.segmentName()
            
            self.offset = offset
            
            from_residue_number = atom.residueNum()
            
            self.to_residue_number = from_residue_number + offset
            self.to_residue_type = Atom_utils._get_residue_type(segment, from_residue_number + offset)
            
            self.complete = True


    def _build_contexts(self, atom, table):
        contexts = []
        for offset in table.get_offsets():
            context = RandomCoilShifts.ResidueOffsetContext(atom, offset, table)
            contexts.append(context)
        return contexts
        
    def set_shifts(self,shift_list):
        for atom_index,shift in self._shifts_list:
            shift_list[atom_index] += shift
            

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

    def get_abbreviated_name(self):
        return "DHA "
    
    def _translate_atom_name(self, atom_name,context):
        return context._table.get_translation(atom_name)
        
    def _get_table(self, residue_type):
        return self._table_manager.get_dihedral_table(residue_type)
#    
    class DihedralContext(object):

        def _select_atom_with_translation(self, segment, residue_number_1, atom_name_1):
            target_atom_1 = Atom_utils._select_atom_with_translation(segment, residue_number_1, atom_name_1)
            if len(target_atom_1) == 0:
                atom_name_1 = self._table.get_translation(atom_name_1)
                target_atom_1 = Atom_utils._select_atom_with_translation(segment, residue_number_1, atom_name_1)
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
            
            self._table = table
            
            segment = from_atom.segmentName()
            
            target_atoms = []
            for key in dihedral_key.get_keys():
                target_atoms.append(self.get_atom_info(segment, from_atom, key))
            
            self.complete=True
            for target_atom in target_atoms:
                if len(target_atom) !=  1:
                    self.complete =  False
                    
            if self.complete:
                self.dihedral_indices = [target_atom[0].index()  for target_atom in target_atoms]
                self.dihedral_key = dihedral_key
                
        
    def _build_contexts(self, atom, table):
        contexts = []
        
        for key in table.get_dihedral_keys():
            dihedral_key = Dihedral_key(*key)
            context = Dihedral_potential.DihedralContext(atom,dihedral_key,table)
            
            if context.complete:
                contexts.append(context)
                
        return contexts


    def  _get_component_for_atom(self, atom, context):
        table = context._table
        
        from_atom_name = atom.atomName()
        from_atom_name = self._translate_atom_name(from_atom_name, context)
        
        result = None
        if from_atom_name in table.get_target_atoms():
            
            dihedral_key = context.dihedral_key
            value = context._table.get_dihedral_shift(from_atom_name,dihedral_key)
            if value != None:
                from_atom_index = atom.index()
                
                dihedral_indices = context.dihedral_indices

                result = [from_atom_index]
                result.extend(dihedral_indices)
                result.append(value)
                
                for parameter_id in context._table.get_parameters():
                    parameter = context._table.get_parameter(from_atom_name,dihedral_key,parameter_id)
                    result.append(parameter)
                
                
                result.append(context._table.get_exponent())
                result = tuple(result)
                
        return result
    
    def dump(self):
        result  = []
        for from_index,atom_1,atom_2, atom_3, atom_4,value,param_0,param_1,param_2,param_3,param_4,exponent in self._distances:
            sub_result  = []
            key = []
            sub_result.append(key)
            key.append(tuple(Atom_utils._get_atom_info_from_index(from_index)[1:]))
            
            dihedral_atom_info = []
            for atom_pair in (atom_1,atom_2), (atom_3, atom_4):
                
                info_pair  = []
                info_pair.append(Atom_utils._get_atom_info_from_index(atom_pair[0])[1:])
                info_pair.append(Atom_utils._get_atom_info_from_index(atom_pair[1])[1:])
            
                dihedral_atom_info.append(tuple(info_pair))
            key.append(tuple(dihedral_atom_info))
            
            sub_result.append(value)
            sub_result.extend((param_0,param_1,param_2,param_3,param_4,))
            sub_result.append(exponent)
            
            result.append(tupleit(sub_result))
            
            
        return tuple(result)
    
    def __str__(self):
        print len( self._distances)
        result = []
        for from_index, \
            atom_1,atom_2, atom_3, atom_4,                           \
            value,                                                   \
            param_0,param_1,param_2,param_3,param_4,                 \
            exponent in self._distances:
                                             
            from_atom = self._get_atom_name(from_index)
            distance_atom_1 = self._get_atom_name(atom_1)
            distance_atom_2 = self._get_atom_name(atom_2)
            distance_atom_3 = self._get_atom_name(atom_3)
            distance_atom_4 = self._get_atom_name(atom_4)
                        
            template = '[%s] - [[%s] - [%s] - [%s] - [%s]] (%7.3f %7.3f %7.3f %7.3f %7.3f) %7.3f %7.3f'
            values = from_atom, distance_atom_1, distance_atom_2,         \
                                distance_atom_3, distance_atom_4,         \
                                param_0,param_1,param_2,param_3,param_4,  \
                                value, exponent
            result.append( template % values)
        return '\n'.join(result)

    def _get_dihedral_angle(self, dihedral_1_atom_id_1, dihedral_1_atom_id_2, 
                                 dihedral_2_atom_id_1, dihedral_2_atom_id_2):
        
        atom_1  = Atom_utils._get_atom_by_index(dihedral_1_atom_id_1)
        atom_2  = Atom_utils._get_atom_by_index(dihedral_1_atom_id_2)
        atom_3  = Atom_utils._get_atom_by_index(dihedral_2_atom_id_1)
        atom_4  = Atom_utils._get_atom_by_index(dihedral_2_atom_id_2)
        
        return Dihedral(atom_1,atom_2,atom_3,atom_4).value();
    
    
    def _calc_single_shift(self, index):
        
        data = self._distances[index]
        
        dihedral_1_atom_id_1, dihedral_1_atom_id_2, \
        dihedral_2_atom_id_1, dihedral_2_atom_id_2 = data[1:5]
        
        coefficient =  data[5]
        
        parameter_0, parameter_1,parameter_2,parameter_3,parameter_4 = data[6:11]
        
        angle = self._get_dihedral_angle(dihedral_1_atom_id_1, dihedral_1_atom_id_2, 
                                        dihedral_2_atom_id_1, dihedral_2_atom_id_2)

        
        angle_term = parameter_0 * cos(3.0 * angle + parameter_3) + \
                     parameter_1 * cos(angle + parameter_4) +       \
                     parameter_2
        
        shift = coefficient * angle_term

        return shift
    
    # TODO move to base_class
    def set_shifts(self,result):
        for index in range(len(self._distances)):
            shift = self._calc_single_shift(index)
            from_atom_id = self._distances[index][0]
            result[from_atom_id] += shift

        return result

    
class Sidechain_potential(Distance_based_potential):
    
    def __init__(self):
        Base_potential.__init__(self)
        
        self._distances =  self._create_component_list(self.ALL)    
        
    def  _get_table(self, residue_type):
        return self._table_manager.get_sidechain_table(residue_type)
    
    def _get_data_table(self):
        return self._distances
    
    def _get_indices(self):
        return super(Sidechain_potential, self)._get_indices()
    
    def _build_contexts(self,atom, table):
        contexts = []
        residue_type =  atom.residueName()
        if residue_type in table.get_residue_types():
            for sidechain_atom in table.get_sidechain_atoms(residue_type):
                context = Sidechain_potential.Sidechain_context(atom, residue_type, sidechain_atom,table)
                contexts.append(context)
        return contexts
    
    class Sidechain_context():
        
        
        
        def __init__(self, from_atom, residue_type, sidechain_atom, table):
            
            self._table = table
            
            segment = from_atom.segmentName()
            sidechain_residue_number = from_atom.residueNum()
            sidechain_atom = Atom_utils._select_atom_with_translation(segment, sidechain_residue_number, sidechain_atom)
            
            
            self.sidechain_atom_index = sidechain_atom[0].index()
            
            self.residue_type = residue_type
            self.sidechain_atom_name = sidechain_atom[0].atomName()
            
            self.complete = True
        

    
    def  get_abbreviated_name(self):
        return "SDCH"
    
    def  _get_component_for_atom(self, target_atom, context):
        table = context._table
        table = table

        
        target_atom_name = target_atom.atomName()
        target_atom_name = self._translate_atom_name(target_atom_name, context)
        
        result = None
        if target_atom_name in table.get_target_atoms():
            sidechain_atom_name = context.sidechain_atom_name
            residue_type = context.residue_type
            value = table.get_sidechain_coefficient(residue_type,target_atom_name,sidechain_atom_name)

            if value != None:
                target_atom_index = target_atom.index()
                sidechain_atom_index= context.sidechain_atom_index
                exponent = table.get_exponent()
                result = (target_atom_index,sidechain_atom_index,value,exponent)
        return result

    def _translate_atom_name(self, atom_name,context):
        return atom_name
    
    def  set_shifts(self,result):
        for index in range(len(self._distances)):
            shift = self._calc_single_shift(index)
            from_atom_id = self._distances[index][0]
            result[from_atom_id] += shift

        return result

    def _calc_single_shift(self, index):
        data = self._distances[index]
        
        target_atom_index,sidechain_atom_index,value,exponent =  data
        
        
        distance = Atom_utils._calculate_distance(target_atom_index, sidechain_atom_index)

        return distance ** exponent * value
    
    def dump(self):
        result  = []
        for from_index,sidechain_index,value,exponent in self._distances:
            sub_result  = []
            
            sub_result.append(Atom_utils._get_atom_info_from_index(from_index))
            sub_result.append(Atom_utils._get_atom_info_from_index(sidechain_index))
            
            sub_result.append(value)
            sub_result.append(exponent)
            
            result.append(tuple(sub_result))
        return result
    
    def __str__(self):
        result = []
        for target_index,sidechain_index,value,exponent in self._distances:
            target_atom = self._get_atom_name(target_index)
            sidechain_atom = self._get_atom_name(sidechain_index)
            
            template = '[%s] - [%s] %7.3f %7.3f'
            values = target_atom, sidechain_atom, value, exponent
            result.append( template % values)
        return '\n'.join(result)
    
class Xcamshift():
    def __init__(self):
        self.potential = [RandomCoilShifts(),
                          Distance_potential(),
                          Extra_potential(),
                          Dihedral_potential(),
                          Sidechain_potential()]
        self._shift_table = Observed_shift_table()
                
    def print_shifts(self):
        result  = [0] * Segment_Manager().get_number_atoms()
        
        result_elements = {}
        keys=[]
        for potential in self.potential:
            num_atoms = len(result)
            sub_result  = [0.0] * num_atoms
            potential.set_shifts(sub_result)
            key = potential.get_abbreviated_name()
            keys.append(key)
            result_elements[key] = sub_result
        
        
        keys = []
        for potential in self.potential:
            num_atoms = len(result)
            sub_result  = [0.0] * num_atoms
            potential.set_shifts(sub_result)
            key = potential.get_abbreviated_name()
            keys.append(key)
            result_elements[key] = sub_result
        
        
        total = [sum(elems) for elems in zip(*result_elements.values())]
        keys.append('TOTL')
        
        result_elements['TOTL'] = total
        residues  = []
        atoms = []
        result_elements['ATOM'] =  atoms
        result_elements['RESD'] =  residues
        
        keys.insert(0, 'RESD')
        keys.insert(0,'ATOM')
        
        for i in range(num_atoms):
            segments,residue,atom = Atom_utils._get_atom_info_from_index(i)
            residues.append(('%5i' %residue))
            atoms.append(atom)
        
        for key in keys:
            values = []
            print key.ljust(5),
            for total,value in zip(result_elements['TOTL'],result_elements[key]):
                if total > 0:
                    if isinstance(value, float):
                        string  = '%- 7.4f' % value
                    else:
                        string = value
                    string = string.rjust(9)
                    
                    values.append(string)
            print ' '.join(values)
            
        
    def set_shifts(self, result):
        keys = []
        result_elements = {}
        for potential in self.potential:
            num_atoms = len(result)
            sub_result  = [0.0] * num_atoms
            potential.set_shifts(sub_result)
            key = potential.get_abbreviated_name()
            keys.append(key)
            result_elements[key] = sub_result
        
        
        result = [sum(elems) for elems in zip(*result_elements.values())]
        
        return result
    
    def set_observed_shifts(self, shift_table):
        self._shift_table  =  shift_table
        
    #TODO grossly inefficient!
    def _calc_single_shift(self,atom_index):
        shifts = [0.0] * Segment_Manager().get_number_atoms()
        shifts = self.set_shifts(shifts)
        return shifts[atom_index]
    
    
    def _get_flat_bottom_shift_limit(self, residue_type, atom_name):
        table_manager = Table_manager.get_default_table_manager()
        constants_table = table_manager.get_constants_table(residue_type)
        
        flat_bottom_limit = constants_table.get_flat_bottom_limit(atom_name)
        #TODO move to class body this is a global scaling for the 
        # whole force field
        flat_bottom_constant  = constants_table.get_flat_bottom_constant()
        return flat_bottom_limit * flat_bottom_constant
    

    def _adjust_shift(self, shift_diff, flat_bottom_shift_limit):
        result  = 0.0
        if (shift_diff > 0.0):
            result = shift_diff-flat_bottom_shift_limit
        else:
            result = shift_diff + flat_bottom_shift_limit
        return result
    
    
    def _get_end_harmonic(self, residue_type, atom_name):
        table_manager = Table_manager.get_default_table_manager()
        constants_table = table_manager.get_constants_table(residue_type)
        
        return constants_table.get_end_harmonic(atom_name)
    
    
    def _get_scale_harmonic(self, residue_type, atom_name):
        table_manager = Table_manager.get_default_table_manager()
        constants_table = table_manager.get_constants_table(residue_type)
        
        return constants_table.get_scale_harmonic(atom_name)
    
    
    def _calc_single_energy(self, target_atom_index):
        
        energy = 0.0
        if target_atom_index in self._shift_table.get_atom_indices():
            
            theory_shift = self._calc_single_shift(target_atom_index)
            observed_shift = self._shift_table.get_chemical_shift(target_atom_index)
            
            shift_diff = observed_shift - theory_shift
            
            residue_type = Atom_utils._get_residue_type_from_atom_id(target_atom_index)
            atom_name = Atom_utils._get_atom_name_from_index(target_atom_index)
            
            flat_bottom_shift_limit = self._get_flat_bottom_shift_limit(residue_type, atom_name)
            
            if abs(shift_diff) > flat_bottom_shift_limit:
                adjusted_shift_diff = self._adjust_shift(shift_diff, flat_bottom_shift_limit)
                
                end_harmonic = self._get_end_harmonic(residue_type, atom_name)
                scale_harmonic = self._get_scale_harmonic(residue_type, atom_name)
                
                
                energy_component = 0.0
                if adjusted_shift_diff < end_harmonic:
                    energy_component = (adjusted_shift_diff/scale_harmonic)**2
                else:
#                    tanh_amplitude = self._get_tanh_amplitude(residue_type,atom_name)
#                    tanh_elongation = self._get_tanh_elongation(residue_type, atom_name)
#                    tan_y_offset = self._get_tan_y_offset(residue_type, atom_name)
#
#                    tanh_argument = tanh_elongation * (adjusted_shift_diff - end_harmonic)
#                    energy_component = tanh_amplitude * tanh(tanh_argument) + tan_y_offset;

                    raise Exception("not implemented")
                energy += energy_component
        return energy
    

    def _get_weight(self, residue_type, atom_name):
        table_manager = Table_manager.get_default_table_manager()
        constants_table = table_manager.get_constants_table(residue_type)
        
        return constants_table.get_weight(atom_name)
    
    #TODO maybe call this a scaling??
    def _calc_single_factor(self, target_atom_index):
        
        factor = 0.0
        if target_atom_index in self._shift_table.get_atom_indices():
            
            theory_shift = self._calc_single_shift(target_atom_index)
            observed_shift = self._shift_table.get_chemical_shift(target_atom_index)
            
            shift_diff = observed_shift - theory_shift
            
            residue_type = Atom_utils._get_residue_type_from_atom_id(target_atom_index)
            atom_name = Atom_utils._get_atom_name_from_index(target_atom_index)
            
            flat_bottom_shift_limit = self._get_flat_bottom_shift_limit(residue_type, atom_name)
            
            if abs(shift_diff) > flat_bottom_shift_limit:
                adjusted_shift_diff = self._adjust_shift(shift_diff, flat_bottom_shift_limit)
                end_harmonic = self._get_end_harmonic(residue_type, atom_name)
                scale_harmonic = self._get_scale_harmonic(residue_type, atom_name)
                sqr_scale_harmonic = scale_harmonic**2
                
                weight = self._get_weight(residue_type,atom_name)
                
                # TODO add factor and lambda to give fact
                fact =1.0
                if adjusted_shift_diff < end_harmonic:
                    factor = 2.0 * weight * adjusted_shift_diff * fact / sqr_scale_harmonic;
                else:
                    raise Exception("not implemented")
                
        else:
            msg = "requested factor for target [%s] which is not in shift table"
            target_atom_info = Atom_utils._get_atom_info_from_index(target_atom_index)
            raise Exception(msg % target_atom_info)
        return factor

    def _calc_single_force_factor(self,target_atom_index,forces):
        for potential in self.potential:
            if isinstance(potential, Distance_potential):
                print >> sys.stderr, "warning forces only set for distance potential!"
                potential.set_observed_shifts(self._shift_table)
                
                if target_atom_index in self._shift_table.get_atom_indices():
                    factor = self._calc_single_factor(target_atom_index)
                    potential._calc_single_factor(target_atom_index,factor,forces)

                          
if __name__ == "__main__":
    initStruct("test_data/3_ala/3ala.psf")
    PDBTool("test_data/3_ala/3ala.pdb").read()
#    RandomCoilShifts()
    Distance_potential()
    
        