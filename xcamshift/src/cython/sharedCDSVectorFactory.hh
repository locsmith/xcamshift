/*
 * SharedCDSVectorFactory.hh
 *
 *  Created on: 26 Jun 2013
 *      Author: garyt
 */


#include <cdsVector.hh>
#include <ensembleSimulation.hh>


#ifndef SHARED_CDS_VECTOR_FACTORY_HH_
#define SHARED_CDS_VECTOR_FACTORY_HH_

//TODO: this file should disappear and be done in pure cython...

//TODO: doubles should  float_type
typedef CDSVector<double,0,EnsembleSimulation::SharedAlloc> SharedVec;
//typedef CDSVector<double> FloatVec;

class SharedCDSVectorFactory {
public:

	static void* createSharedVector(const int size,const float& i, const EnsembleSimulation* simulation);
	static void clearSharedVector(void*);
	static void resizeSharedVector(void*, int size);
	static inline void setSharedVectorValue(void* sharedVec, int offset ,double value) {(*((SharedVec*)(sharedVec)))[offset]=value;};
	static inline double getSharedVectorValue(void* sharedVec, int offset) {return (*((SharedVec*)(sharedVec)))[offset];};
	static inline void addToSharedVectorValue(void* sharedVec, int offset ,double value) {(*((SharedVec*)(sharedVec)))[offset]+=value;};

};

//template <class T>
//class SharedCDSVector : public CDSVector<T, 0,EnsembleSimulation::SharedAlloc> {
//
//	public:
//	SharedCDSVector(const int size, const float& i, const EnsembleSimulation* simulation): CDSVector<T,0,EnsembleSimulation::SharedAlloc>(size,i,EnsembleSimulation::SharedAlloc(simulation)) {}
//
//};

//template <class T> CDSVector<T,>* SharedCDSVectorFactory<T,0,EnsembleSimulation::SharedAlloc>::createSharedVec (const int size,const float& i, const EnsembleSimulation* simulation){
//	return new CDSVector<T,0,EnsembleSimulation::SharedAlloc>(size,i,EnsembleSimulation::SharedAlloc(simulation));
//}



#endif /* SHARED_CDS_VECTOR_FACTORY_HH_ */
