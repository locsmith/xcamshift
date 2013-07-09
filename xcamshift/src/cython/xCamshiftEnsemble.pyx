'''
Created on 17 Jun 2013

@author: garyt
'''
from cython.xPyPot import PyPot

from vec3 import norm
from cython.pyEnsemblePot import PyEnsemblePot



class XCamshiftEnsemble(PyPot):
    
    
    '''
    classdocs
    '''
    def __init__(self,name):
        '''constructor - force constant is optional.'''
        super(XCamshiftEnsemble, self).__init__(name)
        self._proxyPot = PyEnsemblePot(name,"camshift",self.simulation())
        
    def calcEnergy(self):
        return self._proxyPot.calcEnergy()
    
    def calcEnergyAndDerivList(self, derivList):
        return self._proxyPot.calcEnergyAndDerivList(derivList)

