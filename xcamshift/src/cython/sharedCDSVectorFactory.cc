/*
 * SharedCDSVectorFactory.cc
 *
 *  Created on: 26 Jun 2013
 *      Author: garyt
 */

#include "sharedCDSVectorFactory.hh"
#include "ensembleSimulation.hh"

SharedVec* SharedCDSVectorFactory::createSharedVec(const int size,const float& i, const EnsembleSimulation* simulation) {

	return new SharedVec(size,0.0,EnsembleSimulation::SharedAlloc(simulation));
}

//TODO: this need to be moved into cython space but currentlly there are type problems
void SharedCDSVectorFactory::clearVector(FloatVec* vector) {
	vector->set(0.0);
}

