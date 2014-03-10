
#ifndef __pyPot_hh__
#define __pyPot_hh__

//
// provide mechanism for writing potential terms in Python
//

#include <pot.hh>
#include <cdsString.hh>

#include <Python.h>

class DerivList;
class Simulation;


class PyPot: public Pot {

  PyObject*   pobj;
  PyObject*   energyMeth;
  PyObject*   energyAndDerivsMeth;
  PyObject*   energyAndDerivListMeth;
  PyObject*   new_VecVec3;
  Simulation* sim;      //is this member needed?

public:

  PyPot(const String&     instanceName,
	      PyObject*   pyObject,
	      Simulation* defaultSimulation);
  ~PyPot();

  float_type calcEnergy();
  float_type calcEnergyAndDerivs(DerivList&);

  //don't return anything useful yet
  float_type rms()        { return -1.; }
  float_type violations() { return -1; }
  int        numRestraints() { return -1; }

  // return the python potential object
  PyObject* pythonPot() { return pobj; }

  size_t simulation() { return  (size_t) sim; }
};

void 
fromPy(PyPot*& ret, 
       PyObject*         obj);

#endif /* __pyPot_hh__ */
