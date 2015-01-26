'''
Created on 11 Feb 2014

@author: garyt
'''
import re
from collections import namedtuple
from atomSel import AtomSel
import sys
from __builtin__ import file

Shift_data = namedtuple('Shift_data', ('atom_id','shift','error','weight','atom_class', 'comment'))
DEFAULT_CLASS='DEFAULT'
DEFAULT_WEIGHT=1.0

class IOException (Exception):
    def __init__(self,message,line_number,line_value,file_name):
        message = 'Error reading shifts [%s] at line %i in %s\n\tline data: |%s|' % (message,line_number,file_name,line_value)
        super(IOException, self).__init__(message)
        self._file_name =  file_name
        self._line_number = line_number
        self._message = message
        self._line_value = line_value

class Xplor_reader:
    def __init__(self,file_name = '\'unknown file\''):
        self._file_name = file_name


    def raise_exception(self, msg):
        raise IOException(msg,self.line_index+1, self.line,self._file_name)



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
        result = []
        for field_index, field in enumerate(shift_fields):
            try:
                result.append(float(field))
            except:
                msg = "couldn't convert field %i [value=|%s|] to a float" % (field_index+1, shift_fields[field_index])
                self.raise_exception(msg)
        return result


    def get_error_or_default_error(self, shift_fields):
        if len(shift_fields) == 1:
            error = None
        elif len(shift_fields) == 2:
            error = shift_fields[1]
        return error


    def raise_bad_statement(self,name):
        msg = "line doesn't contain a correctly formated statement for %s" %name
        return self.raise_exception(msg)

    def read(self,lines):
        msg = "WARNING reading cs data from %s with a very primitive xplor data reader"
        print >> sys.stderr, "WARNING insufficient data validation being done!!"

        results = []

        weight = DEFAULT_WEIGHT
        atom_class = DEFAULT_CLASS


        assign_keyword_re = re.compile('^\s*[Aa][Ss][Ss][Ii][Gg]?[Nn]?\s+(\(.*\))\s+(.+)')
        weight_keyword_re = re.compile('^\s*[Ww][Ee][Ii][Gg][Hh]?[Tt]?\s+(.+)')
        class_keyword_re = re.compile('^\s*[Cc][Ll][Aa][Ss][Ss]?\s+(.+)')

        matchers = (assign_keyword_re,'ASSI',(1,3)),(weight_keyword_re,'WEIG',(0,2)),(class_keyword_re,'CLAS',(0,2))
        for self.line_index,self.line in enumerate(lines.strip().split("\n")):
            line_complete = False
            comment  = None

            self.line = self.line.strip()
            comment_start =  self.line.find('!')
            if comment_start > -1:
                comment = self.line[comment_start+1:]
                self.line =  self.line[:comment_start]

            if len(self.line.strip()) == 0:
                continue

            for matcher,name,range in matchers:

                if not self.line.upper().startswith(name):
                    continue

                match = matcher.match(self.line)
                matches = None
                if match:
                    matches = match.groups()


                if name == 'ASSI':
                    if matches != None and (len(matches) > range[0] and len(matches) < range[1]):
                        atom_id = self.get_selection_or_raise(matches)

                        shift_fields  = matches[1].split()

                        self.check_number_shift_fields_or_raise(shift_fields)

                        shift_fields = self.convert_fields_to_float_or_raise(shift_fields)

                        error = self.get_error_or_default_error(shift_fields)

                        result  = Shift_data(atom_id=atom_id,shift=shift_fields[0],error=error,weight=weight,atom_class=atom_class, comment=comment)

                        results.append(result)
                        line_complete = True
                        break

                    else:
                        self.raise_bad_statement('assign')



                elif name == 'WEIG':
                    if matches != None and (len(matches) > range[0] and len(matches) < range[1]):
                        weight = self.convert_fields_to_float_or_raise(matches)[0]

                        line_complete = True
                        break
                    else:
                        self.raise_bad_statement('weight')

                elif name == 'CLAS':
                    if matches != None and (len(matches) > range[0] and len(matches) < range[1]):
                        atom_class = matches[0].split()
                        if len(atom_class) > 1:
                            raise Exception()
                        else:
                            atom_class = atom_class[0]
                        line_complete = True
                        break
                    else:
                        self.raise_bad_statement('class')

            if line_complete != True:
                self.raise_exception("unexpected expression \'%s\'"  % self.line)
        return results


