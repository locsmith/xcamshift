// -*- c++ -*-
#ifndef PY_ENSEMBLE_POT_PROXY_H
#define PY_ENSEMBLE_POT_PROXY_H 1

#include "ensemblePot.hh"

struct _object;
typedef _object PyObject;

class PyEnsemblePotProxy : public EnsemblePot  {

public:
  PyObject *m_obj;


public:
  PyEnsemblePotProxy(const String&, const String&, Simulation*, PyObject *obj);
  ~PyEnsemblePotProxy();
  
//  float_type energyMaybeDerivs0(DerivList&,bool calcDerivs);
//  float_type energyMaybeDerivs1(DerivList&,bool calcDerivs);
//  float_type energyMaybeDerivs2(DerivList&,bool calcDerivs);
//  float_type energyMaybeDerivs3(DerivList&,bool calcDerivs);
//  float_type energyMaybeDerivs4(DerivList&,bool calcDerivs);


  float_type energyMaybeDerivs0(DerivList*);
  float_type energyMaybeDerivs1(DerivList*);
  float_type energyMaybeDerivs2(DerivList*);
  float_type energyMaybeDerivs3(DerivList*);
  float_type energyMaybeDerivs4(DerivList*);

  
  float_type callCyEnergyMaybeDerivs(DerivList*, bool, int);
  
  EnsembleSimulation* ensembleSimulation() {return esim;};
  
  
  float_type rms();
  float_type violations();
  int numRestraints();
};

#endif // !PY_ENSEMBLE_POT_PROXY_H
