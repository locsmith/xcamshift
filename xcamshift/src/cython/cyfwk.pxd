from cpython.ref cimport PyObject

cdef extern from "lib.h":
    cdef cppclass IAlg:
        void doRun()
        void run()
        void doJump()
        void jump()

    cdef IAlg* create_cy_alg()

    cdef cppclass CyAlgBase:
        #PyObject *obj
        #RunFct fct
        
        CyAlgBase(PyObject *obj)
        void doRun()
        void run()
        void doJump(int i)
        void jump(int i)
    
cdef public api IAlg* cy_create_alg()


