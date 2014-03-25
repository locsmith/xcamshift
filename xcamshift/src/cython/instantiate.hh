#include "cdsMath.hh"
#include "cdsVector.hh"
#include "cdsVector.cc"
#include "cdsString.hh"
#include "cdsString.cc"
#include "cdsSStream.hh"
#include "cdsSStream.cc"
#include "vec3.hh"
#include "fixedVector.cc"
#include "ensembleSimulation.hh"
#include "cdsList.hh"
#include "cdsList.cc"
template class CDSVector<Vec3,0,CDS::DefaultAlloc>;
template class  CDSVectorBase<Vec3,CDS::DefaultAlloc>;
template class CDSVector<double,0,EnsembleSimulation::SharedAlloc>;
template class  CDSVectorBase<double,EnsembleSimulation::SharedAlloc>;
template class CDSVector<float,0,CDS::DefaultAlloc>;
template class  CDSVectorBase<float,CDS::DefaultAlloc>;
template class CDSString<char>;
template class CDSStringStreamBuf<char>;
template class CDSList<CDSString<char>, 10>;
