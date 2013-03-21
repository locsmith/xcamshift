#include "Python.h"
#include "lib.h"
#include "stdio.h"

#include "cyfwk_api.h"

#
IAlg::~IAlg()
{
}

IAlg::IAlg()
{
    printf("[c+]  IAlg::IAlg()\n");
}

void
IAlg::doRun()
{
	printf("[c+]  IAlg::doRun()\n");
  run();
}

void
IAlg::doJump(int i)
{
	printf("[c+]  IAlg::doJump( %i)\n",i);
  jump(i);
}

void
IAlg::jump(int i)
{
	printf("[c+]  IAlg::jump(%i)\n",i);
}

class Alg : public IAlg
{
public:
  Alg();
  ~Alg();

  void run();
  void jump(int i);
};

Alg::Alg()
{
	printf("[c+]  Alg::Alg()\n");
}

Alg::~Alg()
{
	printf("[c+]  Alg::~Alg()\n");
}

void
Alg::run()
{
	printf("[c+]  Alg::run()\n");
}

void
Alg::jump(int i)
{
  printf("[c+]  Alg::jump( %i\n",i);
}

IAlg *create_alg(int algid)
{
  switch (algid) {
  case 0:
    return new Alg;
  case 1:
    return 0;
  default:
    return 0;
  }

  return 0;
}

CyAlgBase::CyAlgBase(PyObject *obj) :
  m_obj(obj)
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

void
CyAlgBase::jump(int i)
{
  printf("[c+]  CyAlgBase::jump( %d)\n",i);
  if (this->m_obj) {
    int error;
    // call a virtual overload, if it exists
    cy_call_jump(this->m_obj, i, &error);
    if (error)
      // call parent method
      IAlg::jump(i);
    return;
  }
  printf("** invalid cy-state: obj [%d]\n", this->m_obj);
  return;
}
