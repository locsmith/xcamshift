#include "Python.h"
#include "lib.h"
#include "stdio.h"

#include "cyfwk_api.h"


CyAlgBase::CyAlgBase(const String& potName, const String& instanceName, Simulation* simulation, PyObject *obj) :
  EnsemblePot(potName, instanceName, simulation),  m_obj(obj)
{
  if (import_cyfwk()) {
    printf("[c+]  error in import_cyfwk!\n");
  } else {
    Py_XINCREF(this->m_obj);
  }
}

CyAlgBase::~CyAlgBase()
{
  Py_XDECREF(this->m_obj);
}

void CyAlgBase::doRun(){
	printf("[c+]  CyAlgBase::doRun() call run()\n");
	run();
	printf("[c+]  CyAlgBase::doRun() run() called\n");
}

void
CyAlgBase::run()
{
  if (this->m_obj) {
    printf("[c+]  CyAlgBase::run()\n");
    int error;
    //cy_call_fct(this->m_obj, "run", &error);
    cy_call_run(this->m_obj, &error);
    if (error)
    	printf("** you must override 'run', it is a pure virtual function\n");
    return;
  }
  printf("** invalid cy-state: obj [%d]\n", this->m_obj);
  return;
}

float_type CyAlgBase::rms(){
	return 0.0f;
}
int CyAlgBase::violations(){
	return 0;
}
int CyAlgBase::numRestraints(){
	return 0;
}
