/*
 * SharedCDSVectorFactory.cc
 *
 *  Created on: 26 Jun 2013
 *      Author: garyt
 */

#include "sharedCDSVectorFactory.hh"
#include "ensembleSimulation.hh"

void* SharedCDSVectorFactory::createSharedVector(const int size,const float& i, const EnsembleSimulation* simulation) {

	return new SharedVec(size,0.0,EnsembleSimulation::SharedAlloc(simulation));
}

//TODO: this need to be moved into cython space but currentlly there are type problems
void SharedCDSVectorFactory::clearSharedVector(void* vector) {
	((SharedVec*)vector)->set(0.0);
}

//TODO: this need to be moved into cython space but currentlly there are type problems
void SharedCDSVectorFactory::resizeSharedVector(void* vector, int size) {
	((SharedVec*)vector)->resize(size);
}

