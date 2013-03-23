#ifndef __PYX_HAVE_API__cyfwk
#define __PYX_HAVE_API__cyfwk
#include "Python.h"
#include "cyfwk.h"

static PyTypeObject *__pyx_ptype_5cyfwk_Test = 0;
#define TestType (*__pyx_ptype_5cyfwk_Test)

static void (*__pyx_f_5cyfwk_cy_call_run)(PyObject *, int *) = 0;
#define cy_call_run __pyx_f_5cyfwk_cy_call_run
static void (*__pyx_f_5cyfwk_cy_call_fct)(void *, char *, int *) = 0;
#define cy_call_fct __pyx_f_5cyfwk_cy_call_fct














































































































































































































































































































































































































#if !defined(__Pyx_PyIdentifier_FromString)
#if PY_MAJOR_VERSION < 3
  #define __Pyx_PyIdentifier_FromString(s) PyString_FromString(s)
#else
  #define __Pyx_PyIdentifier_FromString(s) PyUnicode_FromString(s)
#endif
#endif

























































































































































































































































































































































































































































































#ifndef __PYX_HAVE_RT_ImportModule
#define __PYX_HAVE_RT_ImportModule
static PyObject *__Pyx_ImportModule(const char *name) {
    PyObject *py_name = 0;
    PyObject *py_module = 0;

    py_name = __Pyx_PyIdentifier_FromString(name);
    if (!py_name)
        goto bad;
    py_module = PyImport_Import(py_name);
    Py_DECREF(py_name);
    return py_module;
bad:
    Py_XDECREF(py_name);
    return 0;
}
#endif

#ifndef __PYX_HAVE_RT_ImportFunction
#define __PYX_HAVE_RT_ImportFunction
static int __Pyx_ImportFunction(PyObject *module, const char *funcname, void (**f)(void), const char *sig) {
    PyObject *d = 0;
    PyObject *cobj = 0;
    union {
        void (*fp)(void);
        void *p;
    } tmp;

    d = PyObject_GetAttrString(module, (char *)"__pyx_capi__");
    if (!d)
        goto bad;
    cobj = PyDict_GetItemString(d, funcname);
    if (!cobj) {
        PyErr_Format(PyExc_ImportError,
            "%s does not export expected C function %s",
                PyModule_GetName(module), funcname);
        goto bad;
    }
#if PY_VERSION_HEX >= 0x02070000 && !(PY_MAJOR_VERSION==3&&PY_MINOR_VERSION==0)
    if (!PyCapsule_IsValid(cobj, sig)) {
        PyErr_Format(PyExc_TypeError,
            "C function %s.%s has wrong signature (expected %s, got %s)",
             PyModule_GetName(module), funcname, sig, PyCapsule_GetName(cobj));
        goto bad;
    }
    tmp.p = PyCapsule_GetPointer(cobj, sig);
#else
    {const char *desc, *s1, *s2;
    desc = (const char *)PyCObject_GetDesc(cobj);
    if (!desc)
        goto bad;
    s1 = desc; s2 = sig;
    while (*s1 != '\0' && *s1 == *s2) { s1++; s2++; }
    if (*s1 != *s2) {
        PyErr_Format(PyExc_TypeError,
            "C function %s.%s has wrong signature (expected %s, got %s)",
             PyModule_GetName(module), funcname, sig, desc);
        goto bad;
    }
    tmp.p = PyCObject_AsVoidPtr(cobj);}
#endif
    *f = tmp.fp;
    if (!(*f))
        goto bad;
    Py_DECREF(d);
    return 0;
bad:
    Py_XDECREF(d);
    return -1;
}
#endif





























































































































































































































































































































































































































#ifndef __PYX_HAVE_RT_ImportType
#define __PYX_HAVE_RT_ImportType
static PyTypeObject *__Pyx_ImportType(const char *module_name, const char *class_name,
    size_t size, int strict)
{
    PyObject *py_module = 0;
    PyObject *result = 0;
    PyObject *py_name = 0;
    char warning[200];

    py_module = __Pyx_ImportModule(module_name);
    if (!py_module)
        goto bad;
    py_name = __Pyx_PyIdentifier_FromString(class_name);
    if (!py_name)
        goto bad;
    result = PyObject_GetAttr(py_module, py_name);
    Py_DECREF(py_name);
    py_name = 0;
    Py_DECREF(py_module);
    py_module = 0;
    if (!result)
        goto bad;
    if (!PyType_Check(result)) {
        PyErr_Format(PyExc_TypeError,
            "%s.%s is not a type object",
            module_name, class_name);
        goto bad;
    }
    if (!strict && (size_t)((PyTypeObject *)result)->tp_basicsize > size) {
        PyOS_snprintf(warning, sizeof(warning),
            "%s.%s size changed, may indicate binary incompatibility",
            module_name, class_name);
        #if PY_VERSION_HEX < 0x02050000
        if (PyErr_Warn(NULL, warning) < 0) goto bad;
        #else
        if (PyErr_WarnEx(NULL, warning, 0) < 0) goto bad;
        #endif
    }
    else if ((size_t)((PyTypeObject *)result)->tp_basicsize != size) {
        PyErr_Format(PyExc_ValueError,
            "%s.%s has the wrong size, try recompiling",
            module_name, class_name);
        goto bad;
    }
    return (PyTypeObject *)result;
bad:
    Py_XDECREF(py_module);
    Py_XDECREF(result);
    return NULL;
}
#endif

static int import_cyfwk(void) {
  PyObject *module = 0;
  module = __Pyx_ImportModule("cyfwk");
  if (!module) goto bad;
  if (__Pyx_ImportFunction(module, "cy_call_run", (void (**)(void))&__pyx_f_5cyfwk_cy_call_run, "void (PyObject *, int *)") < 0) goto bad;
  if (__Pyx_ImportFunction(module, "cy_call_fct", (void (**)(void))&__pyx_f_5cyfwk_cy_call_fct, "void (void *, char *, int *)") < 0) goto bad;
  Py_DECREF(module); module = 0;
  __pyx_ptype_5cyfwk_Test = __Pyx_ImportType("cyfwk", "Test", sizeof(struct TestObject), 1); if (!__pyx_ptype_5cyfwk_Test) goto bad;
  return 0;
  bad:
  Py_XDECREF(module);
  return -1;
}

#endif /* !__PYX_HAVE_API__cyfwk */
