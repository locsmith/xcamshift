// -*- c++ -*-
#ifndef PY_ENSEMBLE_POT_PROXY_H
#define PY_ENSEMBLE_POT_PROXY_H 1

#include "ensemblePot.hh"


class PyEnsemblePotProxy : public EnsemblePot  {

public:
	PyEnsemblePotProxy(const String&, const String&, Simulation*);


public:
  float_type rms();
  int violations();
  int numRestraints();
};

#endif // !PY_ENSEMBLE_POT_PROXY_H
