'''
Created on 11 Feb 2014

@author: garyt
'''
import re
from atomSel import AtomSel
import sys

DEFAULT_WEIGHT=1.0

class Xplor_reader:
    def __init__(self,file_name):
        self._file_name = file_name
        
    def read(self):
        msg = "WARNING reading cs data from %s with a very primitive xplor data reader"
        print >> sys.stderr, msg % self._file_name
        print >> sys.stderr, "WARNING insufficient data validation being done!!"
        
        result = []
        with open(self._file_name,'r') as file:
            weight = DEFAULT_WEIGHT
            for line in file:
                assign_keyword_re = re.compile('^\s*[Aa][Ss][Ss][Ii][Gg]?[Nn]? (\(.*\))(.+)')
                
                matches = assign_keyword_re.match(line).groups()
                
                if len(matches) > 1:
                    selection =()
                    try:
                        selection = AtomSel(matches[0])
                    except:
                        pass
                    
                    if len(selection) > 1:
                        msg = "atom selection must select a single atom , got %i from %s"
                        raise Exception(msg(len(selection),matches[0]))
                    else:
                        atom_id  =  selection[0].index()
                    distance_fields  = matches[1].split()
                    
                    for i,field in enumerate(distance_fields):
                        distance_fields[i] = float(field)
                    line_data = [atom_id]
                    line_data.extend(distance_fields)
                    line_data.append(weight)
                    result.append(line_data)
        return result
#                     last_bracket_index = line.find(')')
#                     print line[:last_bracket_index], line[last_bracket_index:]
            
        