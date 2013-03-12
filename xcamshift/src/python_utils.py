#-------------------------------------------------------------------------------
# Copyright (c) 2013 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on Apr 10, 2012

@author: garyt
'''
from UserDict import DictMixin

def tupleit(t):
    """ 
    Converts lists to tuples recursively
    """
    return tuple(map(tupleit, t)) if isinstance(t, (list, tuple)) else t

def IsMappingType(obj, require_set=False):
    """
    Returns ``True`` if an object appears to
    support the ``MappingType`` protocol.

    If ``require_set`` is ``True`` then the
    object must support ``__setitem__`` as
    well as ``__getitem__`` and ``keys``.
    """
    if require_set and not hasattr(obj, '__setitem__'):
        return False
    if hasattr(obj, 'keys') and hasattr(obj, '__getitem__'):
        return True
    else:
        return False

def IsSequenceType(obj, require_set=False):
    """
    Returns ``True`` if an object appears to
    support the ``SequenceType`` protocol.

    If ``require_set`` is ``True`` then the
    object must support ``__setitem__`` as
    well as ``__getitem__``.

    The object must not have a ``keys``
    method.
    """
    if require_set and not hasattr(obj, '__setitem__'):
        return False
    if (not hasattr(obj, 'keys')) and hasattr(obj, '__getitem__'):
        return True
    else:
        return False

class Dict_walker():
    
    def __init__(self, include_dicts=False):
        self._seen_keys = None
        self._key_path = None
        self._call_back = None
        self._include_dicts = include_dicts
    
    def walk_dict(self, target_dict,call_back=None):
      
        self._setup(call_back)
        self._walk(target_dict)
    
    def _setup(self, call_back):
        self._seen_keys = []
        self._key_path = []
        self._call_back = call_back
        
    def _walk(self,target_dict):
        
        for key in target_dict.keys():
            self._key_path.append(key)
            
            is_mapping_type = IsMappingType(target_dict[key])
            if is_mapping_type:
                if not key in self._seen_keys:
                    self._seen_keys.append(key)
                    
                    self._walk(target_dict[key])
                    
            if not is_mapping_type or self._include_dicts and is_mapping_type:
                self._call_back(target_dict[key],tuple(self._key_path))
            
            
            self._key_path.pop()
            
    class Aggregator():
        def __init__(self):
            self._result = {}
            
        def __call__(self,elem,path):
            self._result[path]  = elem
            
        def get_result(self):
            return self._result
        
        
def value_from_key_path(target_dict,key_path):
    
    result  =  target_dict
    for key in key_path:
        if key in result:
            result = result[key]

    return result



        
def filter_dict(target_dict, pred, invert = False):
    if invert:
        to_delete = [k for k,v in target_dict.iteritems() if not pred(k,v)]
    else:
        to_delete = [k for k,v in target_dict.iteritems() if pred(k,v)]
    for key in to_delete:
        del target_dict[key]
        
    return target_dict




PARENT_KWARG = 'parent'
class Hierarchical_dict(object,DictMixin):
    '''
    classdocs
    '''
    
    def __init__(self,*args, **kwargs):
        '''
        Constructor
        '''
        
        if PARENT_KWARG in kwargs:
            self._parent = kwargs[PARENT_KWARG]
            del kwargs[PARENT_KWARG]
        else:
            self._parent = None
            
        self._data  = dict(*args,**kwargs)
        
    def set_parent(self,parent):
        self._parent = parent
        
    def get_parent(self):
        return self._parent
    
    def get_hierarchy(self, result =  None):
        if result == None:
            result  = []
            
        
        parent = self.get_parent()
        if parent != None:
            try:
                parent.get_hierarchy(result)
            except:
                result.insert(0,parent)
        result.insert(0,self)
            
        return result
        
    def __setitem__(self, key, value):
        self._data[key]=value
        
    def __getitem__(self,key):
        result = None
        if key in self._data:
            result  = self._data[key]
        elif self._parent != None:
            result = self._parent[key]
        else:
            raise KeyError(key)
        return result
            
    def __delitem__(self, key):
        del self._data[key]
        
    def keys(self):
        
        keys = []
        if self._parent != None:
            keys = self._parent.keys()
        keys.extend(self._data.keys())
        
        return keys
            
    def __str__(self):
        result = [self._data.__str__(), self._parent.__str__()]
        return ' -> '.join(result)
