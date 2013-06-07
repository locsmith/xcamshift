from cpython.ref cimport PyObject
from xplor_access cimport String, Simulation, currentSimulation
cdef extern from "PyEnsemblePotProxy.hh":

    

    cdef cppclass PyEnsemblePotProxy:
        #PyObject *obj
        #RunFct fct
        
        PyEnsemblePotProxy(String& potName,  String& instanceName, Simulation* simulation, PyObject *obj)
        void doRun()
        void run()

    


