#include "Python.h"
#include "PyEnsemblePotProxy.hh"
#include "stdio.h"

#include "ensembleSimulation.hh"
#include "pyEnsemblePot_api.h"
#include "derivList.hh"

PyEnsemblePotProxy::PyEnsemblePotProxy(const String& potName, const String& instanceName, Simulation* simulation, PyObject *obj):
  EnsemblePot(potName, instanceName, simulation),  m_obj(obj)
{

  printf("initted\n");

  if (import_pyEnsemblePot()) {
    printf("[c+]  error in import_pyEnsemblePot!\n");
  } else {
    Py_XINCREF(this->m_obj);
  }

}

PyEnsemblePotProxy::~PyEnsemblePotProxy()
{
  Py_XDECREF(this->m_obj);
}

float_type PyEnsemblePotProxy::energyMaybeDerivs0(DerivList& derivList, bool calcDerivs){return callCyEnergyMaybeDerivs(derivList, calcDerivs, 0);}
float_type PyEnsemblePotProxy::energyMaybeDerivs1(DerivList &derivList, bool calcDerivs){return callCyEnergyMaybeDerivs(derivList, calcDerivs, 1);}
float_type PyEnsemblePotProxy::energyMaybeDerivs2(DerivList &derivList, bool calcDerivs){return callCyEnergyMaybeDerivs(derivList, calcDerivs, 2);}
float_type PyEnsemblePotProxy::energyMaybeDerivs3(DerivList &derivList, bool calcDerivs){return callCyEnergyMaybeDerivs(derivList, calcDerivs, 3);}
float_type PyEnsemblePotProxy::energyMaybeDerivs4(DerivList &derivList, bool calcDerivs){return callCyEnergyMaybeDerivs(derivList, calcDerivs, 4);}


float_type PyEnsemblePotProxy::callCyEnergyMaybeDerivs(DerivList& derivList, bool calcDerivs, int i) {
	printf("proxy called\n");
	return 0.0;
}

float_type PyEnsemblePotProxy::rms(){
	return -1.0f;
}
int PyEnsemblePotProxy::violations(){
	return -1;
}
int PyEnsemblePotProxy::numRestraints(){
	return -1;
}

