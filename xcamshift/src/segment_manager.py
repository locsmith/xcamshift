'''
Created on 31 Dec 2011

@author: garyt
'''

from atomSel import AtomSel

class Segment_Manager():

    def __init__(self):
        (self.__segments, self.__segment_info_map) =  self.__build_segments_and_residues()
    
    class Segment_info():
        def __init__(self, first_residue, last_residue, 
                     first_atom_index, last_atom_index):
            
            self.first_atom_index = first_atom_index
            self.last_atom_index=last_atom_index
            
            self.first_residue=first_residue
            self.last_residue=last_residue
            
        def __str__(self):
            template =  ''' 
                            first residue: %i
                            last ressidue  %i
                            
                            first atom index: %i
                            last atom index: %i'''
            template = dedent(template)
            
            data = (self.first_residue, self.last_residue,
                    self.first_atom_index, self.last_atom_index)
            
            return template % data
    
    def __build_segments_and_residues(self):
        
        atom_selection = AtomSel('(all)')
       
        segments  =  set()
        segment_residues = {}
        
        for atom in atom_selection:
            atom_index = atom.index()
            segment = atom.segmentName()
            segments.add(segment)
            
            residue_num = atom.residueNum()
            residues  = segment_residues.setdefault(segment,{})
            residues.setdefault(residue_num,[]).append(atom_index)
            
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
    
    def get_segment_info(self,segment):
        
        if segment not in self.__segments:
            template = "segment %s is not defined  should be one of %s"
            message = template % (segment, ', '.join(self.segments))
            raise KeyError(message)
        
        return self.__segment_info_map[segment]
         