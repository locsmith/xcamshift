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
        if self._program != None:
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
        for name,value  in expression[1:]:
            
            table[name]=value
        
        return table
    
    def _process_append_to(self,expression,table):
        path= expression[1]
        
        target =  self._get_path_or_raise(table, path)
        
        for name,value  in expression[2:]: 
            target[name]=value
    
        return table
        
    def _do_prepend_to_dict_multiple(self, target, name_values):
        temp = OrderedDict()
        for key in target:
            temp[key] = target[key]
        
        for key in target:
            del target[key]
            
        for name,value in name_values:
            target[name] = value
            
        for key in temp:
            target[key] = temp[key]

    def _do_prepend_to_dict(self, target, name, value):
        temp = OrderedDict()
        for key in target:
            temp[key] = target[key]
        
        for key in target:
            del target[key]
        
        target[name] = value
        for key in temp:
            target[key] = temp[key]

    def _process_prepend_to(self,expression,table):
        path= expression[1]
        target =  self._get_path_or_raise(table, path)
        
        self._do_prepend_to_dict_multiple(target, expression[2:])
        
        return table
    
    def _process_prepend(self,expression,table):
        
        self._do_prepend_to_dict_multiple(table, expression[1:])
        
        return table
    
    def _process_replace(self,expression,table):
        path =  expression[1][:-1]
        name = expression[1][-1]
        value =  expression[2]
        
        target =  self._get_path_or_raise(table, path)
        target[name]=value

        return table

    def _process_add(self,expression,table):
        path =  expression[1][:-1]
        name = expression[1][-1]
        value =  expression[2]
        
        target =  self._get_path_or_raise(table, path)
        target[name]+=value

        return table        
        
    def _get_path_or_raise(self,table,path):
        current = table
        for elem in path:
            current = current[elem]
        return current        

        
     