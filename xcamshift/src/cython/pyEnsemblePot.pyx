

from cython.xPyPot import PyPot
cimport cpython.ref as cpy_ref
from  xplor_access cimport Simulation, String, DerivList, EnsembleSimulation, EnsembleMemberSimulation
from libc.stdio cimport printf
# 
cdef class PyEnsemblePotData:
    cdef PyEnsemblePotProxy* ensemblePotProxy
    cdef String* _instance_name
    cdef String* _potential_name
    cdef Simulation *_simulation 
# 
    def __init__(self, instance_name, potential_name='test'):
        self._instance_name = new String(instance_name, len(instance_name))
        self._potential_name = new String(potential_name, len(potential_name))
        self._simulation = <Simulation*><size_t>self.simulation()
        self.ensemblePotProxy = new PyEnsemblePotProxy(self._instance_name[0], self._potential_name[0],self._simulation, <cpy_ref.PyObject*>self)
     
    def calcEnergyAndDerivList(self,derivList):
        pointer = int(derivList.this)
        result = self.ensemblePotProxy[0].calcEnergyAndDerivs((<DerivList*><size_t>pointer)[0])
        return result
    
class PyEnsemblePot(PyPot,PyEnsemblePotData):
    def __init__(self,name):
        PyPot.__init__(self,name)
        PyEnsemblePotData.__init__(self, self.instanceName(),self.potName())
      
    def calcEnergyAndDerivsMaybe0(self, Py_ssize_t derivListPtr, Py_ssize_t ensembleSimulationPtr, bint calcDerivatives):
        return 0.0
  
    def calcEnergyAndDerivsMaybe1(self, Py_ssize_t derivListPtr, Py_ssize_t ensembleSimulationPtr, bint calcDerivatives):
        return 0.0
    
    def calcEnergyAndDerivsMaybe2(self, Py_ssize_t derivListPtr, Py_ssize_t ensembleSimulationPtr, bint calcDerivatives):
        return 0.0
 
    def calcEnergyAndDerivsMaybe3(self, Py_ssize_t derivListPtr, Py_ssize_t ensembleSimulationPtr, bint calcDerivatives):
        return 0.0
 
    def calcEnergyAndDerivsMaybe4(self, Py_ssize_t derivListPtr, Py_ssize_t ensembleSimulationPtr, bint calcDerivatives):
        return 0.0
# 
cdef public api double cy_call_calc_energy_and_derivs_maybe(object self, int i, DerivList* derivList, EnsembleSimulation* esim, bint calcDerivs, float* result,  int *error):
 
    if i == 0:
        result[0] = self.calcEnergyAndDerivsMaybe0(<Py_ssize_t> derivList, <Py_ssize_t> esim, calcDerivs)
    elif i == 1:
        result[0] = self.calcEnergyAndDerivsMaybe1(<Py_ssize_t> derivList, <Py_ssize_t> esim, calcDerivs)
    elif i == 2:
        result[0] = self.calcEnergyAndDerivsMaybe2(<Py_ssize_t> derivList, <Py_ssize_t> esim, calcDerivs)        
    elif i == 3:
        result[0] = self.calcEnergyAndDerivsMaybe3(<Py_ssize_t> derivList, <Py_ssize_t> esim, calcDerivs)  
    elif i == 4:
        result[0] = self.calcEnergyAndDerivsMaybe4(<Py_ssize_t> derivList, <Py_ssize_t> esim, calcDerivs)


