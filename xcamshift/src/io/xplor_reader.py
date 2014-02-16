'''
Created on 11 Feb 2014

@author: garyt
'''
import re
from collections import namedtuple
from atomSel import AtomSel
import sys

Shift_data = namedtuple('Shift_data', ('atom_id','shift','error','weight','atom_class'))
DEFAULT_CLASS='DEFAULT'
DEFAULT_WEIGHT=1.0

class Xplor_reader:
    def __init__(self,file_name):
        self._file_name = file_name
       

    def get_file(self):
        return open(self._file_name, 'r')

    def read(self):
        msg = "WARNING reading cs data from %s with a very primitive xplor data reader"
        print >> sys.stderr, msg % self._file_name
        print >> sys.stderr, "WARNING insufficient data validation being done!!"
        
        results = []
        with self.get_file() as file:
            weight = DEFAULT_WEIGHT
            atom_class = DEFAULT_CLASS
            for line_no,line in enumerate(file):
                
                assign_keyword_re = re.compile('^\s*[Aa][Ss][Ss][Ii][Gg]?[Nn]? (\(.*\))(.+)')
                
                matches = assign_keyword_re.match(line).groups()
                
                if len(matches) > 1:
                    selection =()
                    try:
                        selection = AtomSel(matches[0])
                    except:
                        pass
                    
                    if len(selection) > 1 or len(selection) == 0:
                        msg = "atom selection must select a single atom , got %i from %s at line %i in file %s:\n\t%s"
                        raise Exception(msg(len(selection),matches[0],line_no,self._file_name,line))
                    else:
                        atom_id  =  selection[0].index()
                        
                    shift_fields  = matches[1].split()
                    
                    for i,field in enumerate(shift_fields):
                        shift_fields[i] = float(field)
                    result  = Shift_data(atom_id=atom_id,shift=shift_fields[0],error=shift_fields[1],weight=weight,atom_class=DEFAULT_CLASS)

                    results.append(result)
        return results
#                     last_bracket_index = line.find(')')
#                     print line[:last_bracket_index], line[last_bracket_index:]
            
        