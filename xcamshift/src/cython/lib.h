// -*- c++ -*-
#ifndef FWK_LIB_H
#define FWK_LIB_H 1

#include "ensemblePot.hh"

struct _object;
typedef _object PyObject;



class CyAlgBase : public EnsemblePot {
public:
  PyObject *m_obj;

  CyAlgBase(const String&, const String&, Simulation*, PyObject *obj);
  virtual ~CyAlgBase();
  void doRun();
  virtual void run();

  float_type rms();
  int violations();
  int numRestraints();
};

#endif // !FWK_LIB_H
