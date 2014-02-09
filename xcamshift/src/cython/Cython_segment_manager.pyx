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
# cython: c_string_encoding=ascii
'''
Created on 11 Aug 2012

@author: garyt
'''

'''
Created on 31 Dec 2011

@author: garyt
'''

cdef extern from "instantiate.hh":
    pass

from textwrap import dedent
from  xplor_access cimport currentSimulation, Atom

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
 
cdef object _default = None

cdef class Segment_Manager:
    cdef int _num_atoms
    cdef object _segments, _segment_info_map
    def __init__(self):
        self._num_atoms=0
        (self._segments, self._segment_info_map) =  self._build_segments_and_residues()
    
    
    
    @classmethod 
    def get_segment_manager(cls):
        global _default
        if _default == None:
            _default = Segment_Manager()
        return _default
    
    @classmethod
    def reset_segment_manager(cls):
        global _default
        _default = None
        
    
        
    cpdef get_number_atoms(self):
        return self._num_atoms
    
    def _build_segments_and_residues(self):
        
        cdef Atom atom
        cdef int atom_id
        cdef str segment
        cdef int residue_number
        self._num_atoms  =  currentSimulation().numAtoms()
        
        segments  =  set()
        segment_residues = {}
        
        for atom_id in range(self._num_atoms):
            
            atom = currentSimulation().atomByID(atom_id)
            segment = <char*> atom.segmentName()
            segments.add(segment)
            
            residue_number = atom.residueNum()
            residues  = segment_residues.setdefault(segment,{})
            residues.setdefault(residue_number,[]).append(atom_id)
            
        segment_info_map = {}
        for segment in segment_residues:
            residues = segment_residues[segment].keys()
            residues.sort()
            
            first_residue = min(residues) 
            last_residue = max(residues)
            
            first_atom_index =  min(segment_residues[segment][first_residue])
            last_atom_index = max(segment_residues[segment][last_residue])
            
            segment_info = Segment_info(first_residue,last_residue,
                                                        first_atom_index,last_atom_index)
            segment_info_map[segment] =  segment_info
        return segments,segment_info_map
    
    def get_segments(self):
        
        return tuple(self._segments)
    

    def __check_segment(self, segment):
        if segment not in self._segments:
            template = "segment %s is not defined  should be one of %s"
            message = template % (segment, ', '.join(self._segments))
            raise KeyError(message)

    def get_segment_info(self,segment):
        
        self.__check_segment(segment)
        
        return self._segment_info_map[segment]
    
    
#         
