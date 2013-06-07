#include "Python.h"
#include "PyEnsemblePotProxy.hh"
#include "stdio.h"

#include "pyEnsemblePot_api.h"


PyEnsemblePotProxy::PyEnsemblePotProxy(const String& potName, const String& instanceName, Simulation* simulation, PyObject *obj) :
  EnsemblePot(potName, instanceName, simulation),  m_obj(obj)
{
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

void PyEnsemblePotProxy::doRun(){
	printf("[c+]  PyEnsemblePotProxy::doRun() call run()\n");
	run();
	printf("[c+]  PyEnsemblePotProxy::doRun() run() called\n");
}

void
PyEnsemblePotProxy::run()
{
  if (this->m_obj) {
    printf("[c+]  PyEnsemblePotProxy::run()\n");
    int error;
    //cy_call_fct(this->m_obj, "run", &error);
    cy_call_run(this->m_obj, &error);
    if (error)
    	printf("** you must override 'run', it is a pure virtual function\n");
    return;
  }
  printf("** invalid cy-state: obj [%p]\n", this->m_obj);
  return;
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
