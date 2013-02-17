#-------------------------------------------------------------------------------
# Copyright (c) 2013 gary thompson.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
#!/usr/bin/python
'''
Created on 21 Jul 2012

@author: garyt
'''

from common_constants import h_keys
import os

class Make_zero_shifts(object):
    '''
    classdocs
    '''

    template = '%sshifts.dat'

    def __init__(self, path, length, glycines=()):
        '''
        Constructor
        '''
        self.path = path
        self.range = [2,length]
        self.glycines=glycines
        
    def write(self):
        for atom_type in h_keys:
            file_name = self.template % atom_type
            file_path = os.path.join(self.path,file_name)
            
            file_h = open(file_path,'w+')
            
            for i in range(*self.range):
                if atom_type == 'CB' and i in self.glycines:
                    continue
                print >> file_h, "%-4i %9.8f" % (i, 0.000001) 
            
            file_h.close()
        
if __name__ == "__main__":
    path = 'data'
    
    shift_maker = Make_zero_shifts(path,56,(9,14,38,41))
    shift_maker.write()

    
