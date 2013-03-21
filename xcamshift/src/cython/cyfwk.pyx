cimport cyfwk 
cimport cpython.ref as cpy_ref
from cython.operator cimport dereference as deref
import inspect
#cimport cpython.cobject as cpy_cobj

cdef class AlgBase:
    cdef CyAlgBase* alg
    def __init__(self):
        print ("[cy]  AlgBase.__init__")
        self.alg = new CyAlgBase(
            <cpy_ref.PyObject*>self)
            
    #------- non-virutal methods --------
    def doRun(self):
        print ("[cy]  AlgBase.doRun()")
        self.alg.doRun()

    def doJump(self, int i):
        print ("[cy]  AlgBase.doJump(%d)" % i)
        self.alg.doJump(i)

#------- virutal methods --------

cdef public api void cy_call_run(object self, int *error):
    print("[cy]  cy_call_run()")
    try:
        func = self.run
    except AttributeError:
        error[0] = 1
    else:
        error[0] = 0
        func()

cdef public api void cy_call_jump(object self, int i, int *error):
    print("[cy]  cy_call_jump(%d)" % i)
    try:
        func = self.jump
    except AttributeError:
        error[0] = 1
    else:
        error[0] = 0
        func(i)

#------- implementation --------

cdef class MyAlg(AlgBase):
    cdef int i
    def __init__(self, i=-1):
        print("[cy]  MyAlg.__init__()")
        self.i = i
        AlgBase.__init__(self)
    def run(self):
        print("[cy]  MyAlg.run() => %i" % self.i)


cdef public api IAlg* cy_create_alg():
    cdef AlgBase cyalg = MyAlg(42)
    return <IAlg*>(cyalg.alg)

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



