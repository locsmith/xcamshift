#include "PyEnsemblePotProxy.hh"
#include "stdio.h"


PyEnsemblePotProxy::PyEnsemblePotProxy(const String& potName, const String& instanceName, Simulation* simulation) :
EnsemblePot(potName, instanceName, simulation){

	printf("initted\n");
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

