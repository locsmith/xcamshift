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

float_type PyEnsemblePotProxy::rms(){
	return -1.0f;
}
int PyEnsemblePotProxy::violations(){
	return -1;
}
int PyEnsemblePotProxy::numRestraints(){
	return -1;
}

