#include "cdsMath.hh"
#include "cdsVector.hh"
#include "cdsVector.cc"
#include "vec3.hh"
#include "ensembleSimulation.hh"
#include "cdsString.hh"
#include "cdsString.cc"

template class CDSVector<Vec3,0,CDS::DefaultAlloc>;
template class  CDSVectorBase<Vec3,CDS::DefaultAlloc>;
template class CDSVector<double,0,EnsembleSimulation::SharedAlloc>;
template class  CDSVectorBase<double,EnsembleSimulation::SharedAlloc>;
template class CDSVector<float,0,CDS::DefaultAlloc>;
template class  CDSVectorBase<float,CDS::DefaultAlloc>;
//template class CDSString<char>;

