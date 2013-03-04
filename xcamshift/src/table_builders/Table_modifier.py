'''
Created on 1 Mar 2013

@author: garyt
'''

from collection_backport import OrderedDict

class Table_modifier(object):
    '''
    classdocs
    '''


    def __init__(self, program):
        '''
        Constructor
        '''
        self._program=program
    
    
    def run(self,table):
        for expression in self._program:
            table  = self._execute_expression_or_raise(expression,table)
        
        return table
    
    def _execute_expression_or_raise(self,expression,table):
        action = '_process_' + expression[0]
        if hasattr(self, action):
            func = getattr(self, action)
            func(expression, table)
        else:
            raise Exception("action %s is not defined" % action)
        return table
            
        
    def _process_append(self,expression,table):
        name,value  = expression[1]
        
        table[name]=value
        
        return table
    
    def _process_append_to(self,expression,table):
        path= expression[1]
        name,value=expression[2]
        
        target =  self._get_path_or_raise(table, path)
        
        target[name]=value
    
        return table

    def _get_path_or_raise(self,table,path):
        current = table
        for elem in path:
            current = current[elem]
        return current        

        
     