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
        print 'calc energy'
#         result = self._proxyPot.calcEnergy()
        return 0.0
    
    def calcEnergyAndDerivList(self, derivList):
#        return self._proxyPot.calcEnergyAndDerivList(derivList)
        return 0.0

