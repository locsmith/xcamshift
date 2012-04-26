'''
Created on 24 Apr 2012

@author: garyt


'''
import sys
import re
from collection_backport import OrderedDict

def get_num_residues_or_exit(lines):
    residues = None
    for line in lines:
        line = line.strip()
        if re.match('>> Number of fragments', line):
            elems = line.split()
            residues = int(elems[4])
            
    if residues ==  None:
        print >> sys.stderr, "can't find a line of the form '>> Number of fragments:   3'"
        print >> sys.stderr, "exiting..."
        sys.exit(-1)
        
    return residues
            
if __name__ == '__main__':
    file_h = open(sys.argv[1])
    
    
    lines = file_h.readlines()
    
    residues  = get_num_residues_or_exit(lines)
    shifts = lines[-residues:]
    header = lines[-residues-2]
    atom_types = header.strip().split()[2:]
    print atom_types
    
    result = OrderedDict()
    
    for shift in shifts:
        elems  = shift.strip().split()
        
        segid=''
        residue_number = elems[0]
        for atom,shift_string in zip(atom_types,elems[2:]):
            key =  segid,residue_number, atom
            result[key] = float(shift_string)
    print result
            
        