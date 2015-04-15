'''
Created on 21 Feb 2013

@author: garyt
'''


class Cache_manager(object):
    '''
    classdocs
    '''

    _cache_manager = None
    
    def __init__(self):
        '''
        Constructor
        '''
        
        self._reset()
        self._cache_builders = {}
        
    def _reset(self):    
        self._cache =  {}

        
    @staticmethod
    def get_cache_manager():
        if Cache_manager._cache_manager == None:
            Cache_manager._cache_manager = Cache_manager()
        return Cache_manager._cache_manager

            
    def clear_cache(self):
        self._reset()
#
    def add_cache_builder(self,name,builder):
        if name not in self._cache_builders:
            self._cache_builders[name] =  builder
        else:
            raise Exception("there is already a builder for %s" % name)
        
    def get_cache(self,name):
        if name not in self._cache:
            if name in self._cache_builders:
                builder = self._cache_builders[name]
                self._cache[name] =  builder()
            else:
                raise Exception("caching is not enabled for %s" % name)
        return self._cache[name]