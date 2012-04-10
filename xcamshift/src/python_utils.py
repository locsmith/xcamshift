'''
Created on Apr 10, 2012

@author: garyt
'''

def tupleit(t):
    return tuple(map(tupleit, t)) if isinstance(t, (list, tuple)) else t