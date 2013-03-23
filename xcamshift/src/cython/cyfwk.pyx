import xplorWrap
cimport cyfwk 
cimport cpython.ref as cpy_ref
from  xplor_access cimport Simulation, String

from cython.operator cimport dereference as deref
import inspect
#cimport cpython.cobject as cpy_cobj 

cdef class AlgBase:
    cdef CyAlgBase* alg
    cdef String* _instance_name
    cdef String* _potential_name
    cdef Simulation *_simulation
    
    def __init__(self, instance_name, potential_name, simulation):
        print ("[cy]  AlgBase.__init__")
        
        self._instance_name = new String(instance_name, len(instance_name))
        self._potential_name = new String(potential_name, len(potential_name))
        self._simulation = currentSimulation()
        self.alg = new CyAlgBase(self._instance_name[0], self._potential_name[0],self._simulation,
            <cpy_ref.PyObject*>self)
            
    #------- non-virutal methods --------
    def doRun(self):
        print ("[cy]  AlgBase.doRun()")
        self.alg.run()


#------- virutal methods --------

cdef public api void cy_call_run(object self, int *error):
    print("[cy]  cy_call_run()")
    try:
        func = self.virtual_call
    except AttributeError:
        error[0] = 1
    else:
        error[0] = 0
        func()



#------- implementation --------

cdef class MyAlg(AlgBase):
    cdef int i
    def __init__(self, i=-1):
        print("[cy]  MyAlg.__init__()")
        self.i = i
        AlgBase.__init__(self)
    def run(self):
        print("[cy]  MyAlg.run() => %i" % self.i)

#
#cdef public api IAlg* cy_create_alg():
#    cdef AlgBase cyalg = MyAlg(42)
#    return <IAlg*>(cyalg.alg)

#------- experiments --------

# this is not being used
cdef public api class Test[object TestObject, type TestType]:
    cpdef jump(self):
        print "[cy]  jump"


# the following is a "generic" lookup, but itonly works for methods that match this signature
cdef public api void cy_call_fct(void *ptr, char* method, int *error):
    """lookup and execute a method. called from c++"""
    print("[cy]  cy_call_fct()")
    cdef AlgBase alg = <AlgBase>(ptr)
    try:
        func = getattr(alg, method)
    except AttributeError:
        print "attr error"
        error[0] = 1
    else:
        error[0] = 0
        func()

print("[cy]  cyfwk module imported")



