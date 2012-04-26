'''
Created on 24 Apr 2012

@author: garyt





'''
import sys
import re
from collection_backport import OrderedDict

class Camshift_output_parser(object):
    
    def __init__(self, source):
        self._source  = source
        
    def _get_num_residues_or_raise(self,lines):
        residues = None
        for line in lines:
            line = line.strip()
            if re.match('>> Number of fragments', line):
                elems = line.split()
                residues = int(elems[4])
                
        if residues ==  None:
            raise Exception("can't find a line of the form '>> Number of fragments:   3'")

            
        return residues
    
    def _open(self):
        return open(self._source)
    
    def _close(self,file_h):
        file_h.close()
        
    def parse(self):
        
        file_h  =  self._open()
        
        result = None
        try:
            lines = file_h.readlines()
            
            residues = self._get_num_residues_or_raise(lines)
            
            shifts = lines[-residues:]
            header = lines[-residues - 2]
            atom_types = header.strip().split()[2:]
            
            result = OrderedDict()
            for shift in shifts:
                elems = shift.strip().split()
                segid = ''
                residue_number = elems[0]
                for atom, shift_string in zip(atom_types, elems[2:]):
                    key = segid, residue_number, atom
                    result[key] = float(shift_string)
        finally:
            self._close(file_h)
            
        return result
        
if __name__ == '__main__':
    
    parser = Camshift_output_parser(sys.argv[1])
    print parser.parse()
   
            
        