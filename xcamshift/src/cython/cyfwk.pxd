from cpython.ref cimport PyObject
from xplor_access cimport String, Simulation, currentSimulation
cdef extern from "lib.h":

    

    cdef cppclass CyAlgBase:
        #PyObject *obj
        #RunFct fct
        
        CyAlgBase(String& potName,  String& instanceName, Simulation* simulation, PyObject *obj)
        void doRun()
        void run()

    


