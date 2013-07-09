cimport pyEnsemblePot 
cimport cpython.ref as cpy_ref
from  xplor_access cimport Simulation, String, DerivList

from cython.operator cimport dereference as deref

cdef class PyEnsemblePot:  
    cdef PyEnsemblePotProxy* ensemblePotProxy
    cdef String* _instance_name
    cdef String* _potential_name
    cdef Simulation *_simulation
    
    def __init__(self, instance_name="test1", potential_name="test2", simulation=None):
        print ("[cy]  PyEnsemblePot.__init__")
        
        self._instance_name = new String(instance_name, len(instance_name))
        self._potential_name = new String(potential_name, len(potential_name))
        self._simulation = <Simulation*> <size_t> simulation
        self.ensemblePotProxy = new PyEnsemblePotProxy(self._instance_name[0], self._potential_name[0],self._simulation, <cpy_ref.PyObject*>self)
#,
#            <cpy_ref.PyObject*>self)
        
    
    def calcEnergyAndDerivList(self,derivList):
        pointer = int(derivList.this)
        result = self.ensemblePotProxy[0].calcEnergyAndDerivs((<DerivList*><size_t>pointer)[0])
        return result

        
    #------- non-virutal methods --------
    def doRun(self):
        print ("[cy]  PyEnsemblePot.doRun()")
        #self.ensemblePotProxy.run()


#------- virtual methods --------

cdef public api void cy_call_run(object self, int *error):
    print("[cy]  cy_call_run()")
    try:
        func = self.virtual_call
    except AttributeError:
        error[0] = 1
    else:
        error[0] = 0
        func()


#
#
print("[cy]  pyEnsemblePot module imported")

 

