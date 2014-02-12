'''
Created on 11 Feb 2014

@author: garyt
'''

class Xplor_reader:
    def __init__(self,file_name):
        self._file_name = file_name
        
    def read(self):
        with open(self._file_name,'r') as file:
            for line in file:
                print line
            
        