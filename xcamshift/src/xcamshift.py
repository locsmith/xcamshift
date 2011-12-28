'''
Created on 27 Dec 2011

@author: garyt
'''
from atomSel import AtomSel

from protocol import initStruct
from pdbTool import PDBTool

class Distance_potential(object):
    '''
    classdocs
    '''
    distance_set = {
        (-1,'CA')  : 1,
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
            

if __name__ == "__main__":
    initStruct("test_data/gb3/gb3.psf")
    PDBTool("test_data/gb3/gb3_refined_II.pdb").read()
    Distance_potential()
    
        