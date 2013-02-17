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
Created on 31 Dec 2011

@author: garyt
'''

from atomSel import AtomSel
from textwrap import dedent

class Segment_Manager():

    def __init__(self):
        raise Exception("do not use")
        (self.__segments, self.__segment_info_map) =  self.__build_segments_and_residues()
    
    __default = None
    
    @staticmethod 
    def get_segment_manager():
        if Segment_Manager.__default == None:
            Segment_Manager.__default = Segment_Manager()
        return Segment_Manager.__default
    
    @staticmethod
    def reset_segment_manager():
        Segment_Manager.__default = None
        
    
    class Segment_info():
        def __init__(self, first_residue, last_residue, 
                     first_atom_index, last_atom_index):
            
            self.first_atom_index = first_atom_index
            self.last_atom_index=last_atom_index
            
            self.first_residue=first_residue
            self.last_residue=last_residue
            
            self.segment_length = self.last_residue -  first_residue + 1
            
        def __str__(self):
            template =  ''' 
                            first residue: %i
                            last residue  %i
                            
                            first atom index: %i
                            last atom index: %i'''
            template = dedent(template)
            
            data = (self.first_residue, self.last_residue,
                    self.first_atom_index, self.last_atom_index)
            
            return template % data
        
    def get_number_atoms(self):
        return self._num_atoms
    
    def __build_segments_and_residues(self):
        
        atom_selection = AtomSel('(all)')
        
        self._num_atoms = len(atom_selection)
        
        segments  =  set()
        segment_residues = {}
        
        for atom in atom_selection:
            atom_index = atom.index()
            segment = atom.segmentName()
            segments.add(segment)
            
            residue_number = atom.residueNum()
            residues  = segment_residues.setdefault(segment,{})
            residues.setdefault(residue_number,[]).append(atom_index)
            
        segment_info_map = {}
        for segment in segment_residues:
            residues = segment_residues[segment].keys()
            residues.sort()
            
            first_residue = min(residues)
            last_residue = max(residues)
            
            first_atom_index =  min(segment_residues[segment][first_residue])
            last_atom_index = max(segment_residues[segment][last_residue])
            
            segment_info = Segment_Manager.Segment_info(first_residue,last_residue,
                                                        first_atom_index,last_atom_index)
            segment_info_map[segment] =  segment_info
        return segments,segment_info_map
    
    def get_segments(self):
        
        return tuple(self.__segments)
    

    def __check_segment(self, segment):
        if segment not in self.__segments:
            template = "segment %s is not defined  should be one of %s"
            message = template % (segment, ', '.join(self.segments))
            raise KeyError(message)

    def get_segment_info(self,segment):
        
        self.__check_segment(segment)
        
        return self.__segment_info_map[segment]
    
    
         
