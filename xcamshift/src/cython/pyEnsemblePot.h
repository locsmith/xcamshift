#ifndef __PYX_HAVE__pyEnsemblePot
#define __PYX_HAVE__pyEnsemblePot


#ifndef __PYX_HAVE_API__pyEnsemblePot

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

__PYX_EXTERN_C DL_IMPORT(double) cy_call_calc_energy_and_derivs_maybe(PyObject *, int, DerivList *, EnsembleSimulation *, int, float *, int *);

#endif /* !__PYX_HAVE_API__pyEnsemblePot */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initpyEnsemblePot(void);
#else
PyMODINIT_FUNC PyInit_pyEnsemblePot(void);
#endif

#endif /* !__PYX_HAVE__pyEnsemblePot */
