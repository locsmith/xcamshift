#include "Python.h"
#include "PyEnsemblePotProxy.hh"
#include "stdio.h"

#include "ensembleSimulation.hh"
#include "pyEnsemblePot_api.h"
#include "derivList.hh"

PyEnsemblePotProxy::PyEnsemblePotProxy(const String& potName, const String& instanceName, Simulation* simulation, PyObject *obj):
  EnsemblePot(potName, instanceName, simulation),  m_obj(obj)
{
  // do we need to un register?
  registerTo(esim);


  if (import_pyEnsemblePot()) {
    printf("[c+]  error in import_pyEnsemblePot!\n");
  } else {
    Py_XINCREF(this->m_obj);
  }

}

PyEnsemblePotProxy::~PyEnsemblePotProxy()
{
  Py_XDECREF(this->m_obj);
  unRegister(esim);
}


#if OLD_ENSEMBLE_INTERFACE > 0

float_type PyEnsemblePotProxy::energyMaybeDerivs0(DerivList& derivList, bool calcDerivs) {return callCyEnergyMaybeDerivs(&derivList, calcDerivs, 0);}
float_type PyEnsemblePotProxy::energyMaybeDerivs1(DerivList& derivList, bool calcDerivs) {return callCyEnergyMaybeDerivs(&derivList, calcDerivs, 1);}
float_type PyEnsemblePotProxy::energyMaybeDerivs2(DerivList& derivList, bool calcDerivs) {return callCyEnergyMaybeDerivs(&derivList, calcDerivs, 2);}
float_type PyEnsemblePotProxy::energyMaybeDerivs3(DerivList& derivList, bool calcDerivs) {return callCyEnergyMaybeDerivs(&derivList, calcDerivs, 3);}
float_type PyEnsemblePotProxy::energyMaybeDerivs4(DerivList& derivList, bool calcDerivs) {return callCyEnergyMaybeDerivs(&derivList, calcDerivs, 4);}

#else /* OLD_ENSEMBLE_INTERFACE */

float_type PyEnsemblePotProxy::energyMaybeDerivs0(DerivList* derivList){return callCyEnergyMaybeDerivs(derivList, true, 0);}
float_type PyEnsemblePotProxy::energyMaybeDerivs1(DerivList* derivList){return callCyEnergyMaybeDerivs(derivList, true, 1);}
float_type PyEnsemblePotProxy::energyMaybeDerivs2(DerivList* derivList){return callCyEnergyMaybeDerivs(derivList, true, 2);}
float_type PyEnsemblePotProxy::energyMaybeDerivs3(DerivList* derivList){return callCyEnergyMaybeDerivs(derivList, true, 3);}
float_type PyEnsemblePotProxy::energyMaybeDerivs4(DerivList* derivList){return callCyEnergyMaybeDerivs(derivList, true, 4);}

#endif /* OLD_ENSEMBLE_INTERFACE */

float_type PyEnsemblePotProxy::callCyEnergyMaybeDerivs(DerivList* derivList, bool calcDerivs, int i) {
	if (derivList == 0){
		calcDerivs = false;
	}
	PyGILState_STATE gstate;
	gstate = PyGILState_Ensure();


	/* evaluate result or handle exception */
	int error;
	float result = 0.0;
	cy_call_calc_energy_and_derivs_maybe(m_obj,i,derivList, esim, calcDerivs, &result, &error);

	/* Release the thread. No Python API allowed beyond this point. */
	PyGILState_Release(gstate);
	return result;
}

float_type PyEnsemblePotProxy::rms(){
	return -1.0f;
}

#if OLD_ENSEMBLE_INTERFACE > 0

int PyEnsemblePotProxy::violations(){
	return -1;
}

#else /* OLD_ENSEMBLE_INTERFACE */

float_type PyEnsemblePotProxy::violations(){
	return -1;
}

#endif /* OLD_ENSEMBLE_INTERFACE */

int PyEnsemblePotProxy::numRestraints(){
	return -1;
}

