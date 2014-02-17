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
       

    def get_line_source(self):
        return open(self._file_name, 'r')


    def raise_exception(self, msg):
        line_position = 'line %i in file \'%s\', line data:\n|%s|' % (self.line_index+1, self._file_name, self.line)
        msg = '%s %s' % (msg, line_position)
        print msg
        raise Exception(msg)



    def get_selection_or_raise(self, matches):
        selection = ()
        
        try:
            selection = AtomSel(matches[0])
        except:
            pass
        
        if len(selection) > 1 or len(selection) == 0:
            msg = "atom selection must select a single atom, got %i atoms at" % len(selection)
            msg = self.raise_exception(msg)
        else:
            atom_id = selection[0].index()
            
        return atom_id

    def read(self):
        msg = "WARNING reading cs data from %s with a very primitive xplor data reader"
        print >> sys.stderr, msg % self._file_name
        print >> sys.stderr, "WARNING insufficient data validation being done!!"
        
        results = []
        with self.get_line_source() as lines:
            weight = DEFAULT_WEIGHT
            atom_class = DEFAULT_CLASS
            for self.line_index,self.line in enumerate(lines):
                
                assign_keyword_re = re.compile('^\s*[Aa][Ss][Ss][Ii][Gg]?[Nn]? (\(.*\))(.+)')
                
                matches = assign_keyword_re.match(self.line).groups()
                
                if len(matches) > 1:
                    atom_id = self.get_selection_or_raise(matches)
                        
                    shift_fields  = matches[1].split()
                    if len(shift_fields) <1 or len(shift_fields) > 2:
                        msg = "no shifts or too many shifts [%i] found at line %i in lines %s:\n\t%s"
                        raise Exception(msg % (len(shift_fields), line_index,self._file_name,line))
                                         
                    for i,field in enumerate(shift_fields):
                        try:
                            shift_fields[i] = float(field)
                        except:
                            msg = "couldn't convert field %i [%s] to a float at line %i in lines %s:\n\t%s"
                            raise Exception(msg % (i,shift_fields[i], line_index, self._file_name,line))

                    if len(shift_fields) == 1:
                        print >> sys.stderr, "warning error set to 0.1, not a sensible value"
                        error = 0.1
                    elif len(shift_fields) == 2:
                        error = shift_fields[1]

                    result  = Shift_data(atom_id=atom_id,shift=shift_fields[0],error=error,weight=weight,atom_class=DEFAULT_CLASS)
                    results.append(result)
        return results
#                     last_bracket_index = line.find(')')
#                     print line[:last_bracket_index], line[last_bracket_index:]
            
        