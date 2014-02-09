#ifndef __PYX_HAVE__shift_calculators
#define __PYX_HAVE__shift_calculators

struct Component_index_pair;

/* "Cython_shift_calculator.pyx":104
 *                       compiled_components[i].parameters[6]
 * 
 * cdef public struct Component_index_pair:             # <<<<<<<<<<<<<<
 *     int index_1
 *     int index_2
 */
struct Component_index_pair {
  int index_1;
  int index_2;
};

#ifndef __PYX_HAVE_API__shift_calculators

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

#endif /* !__PYX_HAVE_API__shift_calculators */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initshift_calculators(void);
#else
PyMODINIT_FUNC PyInit_shift_calculators(void);
#endif

#endif /* !__PYX_HAVE__shift_calculators */
