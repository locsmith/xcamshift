from cpython.ref cimport PyObject
from xplor_access cimport String, Simulation, currentSimulation, float_type, DerivList, EnsembleSimulation


    
cdef extern from "PyEnsemblePotProxy.hh":

    
    cdef cppclass PyEnsemblePotProxy:
        PyEnsemblePotProxy(String& potName,  String& instanceName, Simulation* simulation, PyObject *object)
        float_type calcEnergyAndDerivs(DerivList&) nogil
        float calcEnergy() nogil
        EnsembleSimulation *ensembleSimulation() nogil
        
