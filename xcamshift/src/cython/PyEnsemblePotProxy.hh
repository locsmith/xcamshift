/*-------------------------------------------------------------------------------
 Copyright (c) 2013-2015 Gary Thompson & The University of Leeds.
 All rights reserved. This program and the accompanying materials
 are made available under the terms of the GNU Lesser Public License v3.0
 which accompanies this distribution, and is available at
 http://www.gnu.org/licenses/lgpl-3.0.html

  Contributors:
      gary thompson - initial API and implementation
-------------------------------------------------------------------------------*/
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
