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
    def __init__(self):
        pass


    def raise_exception(self, msg):
        line_position = 'line %i in input, line data:\n|%s|' % (self.line_index+1, self.line)
        msg = '%s %s' % (msg, line_position)
        
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


    def check_number_shift_fields_or_raise(self, shift_fields):
        if len(shift_fields) < 1 or len(shift_fields) > 2:
            self.raise_exception("no shifts or too many shifts [found = %i, expected 1-2] found at" % len(shift_fields))


    def convert_fields_to_float_or_raise(self, shift_fields):
        for field_index, field in enumerate(shift_fields):
            try:
                shift_fields[field_index] = float(field)
            except:
                msg = "couldn't convert field %i [value=|%s|] to a float" % (field_index+1, shift_fields[field_index])
                self.raise_exception(msg)
        return shift_fields


    def get_error_or_default_error(self, shift_fields):
        if len(shift_fields) == 1:
            error = None
        elif len(shift_fields) == 2:
            error = shift_fields[1]
        return error

    def read(self,lines):
        msg = "WARNING reading cs data from %s with a very primitive xplor data reader"
        print >> sys.stderr, "WARNING insufficient data validation being done!!"
        
        results = []

        weight = DEFAULT_WEIGHT
        atom_class = DEFAULT_CLASS
        assign_keyword_re = re.compile('^\s*[Aa][Ss][Ss][Ii][Gg]?[Nn]?\s+(\(.*\))\s+(.+)')
        for self.line_index,self.line in enumerate(lines.strip().split("\n")):
            
            
            matches = assign_keyword_re.match(self.line).groups()
            
            if matches == None:
                self.raise_exception("line doesn't contain a correctly formated assign statement")
            elif len(matches) > 1 and len(matches) < 3:
                atom_id = self.get_selection_or_raise(matches)
                    
                shift_fields  = matches[1].split()
                
                self.check_number_shift_fields_or_raise(shift_fields)
                 
                shif_fields = self.convert_fields_to_float_or_raise(shift_fields)

                error = self.get_error_or_default_error(shift_fields)

                result  = Shift_data(atom_id=atom_id,shift=shift_fields[0],error=error,weight=weight,atom_class=DEFAULT_CLASS)
                
                results.append(result)
            else:
                self.raise_exception("internal error unexpected number of matches [%i] for assign statement" % len(matces))
        return results

            
        