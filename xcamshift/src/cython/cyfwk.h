#ifndef __PYX_HAVE__cyfwk
#define __PYX_HAVE__cyfwk

struct TestObject;

/* "cyfwk.pyx":64
 * 
 * # this is not being used
 * cdef public api class Test[object TestObject, type TestType]:             # <<<<<<<<<<<<<<
 *     cpdef jump(self):
 *         print "[cy]  jump"
 */
struct TestObject {
  PyObject_HEAD
  struct __pyx_vtabstruct_5cyfwk_Test *__pyx_vtab;
};

#ifndef __PYX_HAVE_API__cyfwk

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

__PYX_EXTERN_C DL_IMPORT(PyTypeObject) TestType;

__PYX_EXTERN_C DL_IMPORT(void) cy_call_run(PyObject *, int *);
__PYX_EXTERN_C DL_IMPORT(void) cy_call_fct(void *, char *, int *);

#endif /* !__PYX_HAVE_API__cyfwk */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initcyfwk(void);
#else
PyMODINIT_FUNC PyInit_cyfwk(void);
#endif

#endif /* !__PYX_HAVE__cyfwk */
