#include "cdsMath.hh"
#include "cdsVector.hh"
#include "vec3.hh"
#include "cdsVector.cc"
#include <set>
#include <vector>
template class CDSVector<Vec3,0,CDS::DefaultAlloc>;
template class  CDSVectorBase<Vec3,CDS::DefaultAlloc>;
template class CDSVector<int,0,CDS::DefaultAlloc>;
template class  CDSVectorBase<int,CDS::DefaultAlloc>;
template class std::set<int>;
template class std::vector<int>;
