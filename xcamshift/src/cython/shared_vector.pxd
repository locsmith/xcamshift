from xplor_access cimport CDSVector, EnsembleSimulation

cdef extern from "sharedCDSVectorFactory.hh" namespace "SharedCDSVectorFactory":
    void clearVector(CDSVector[double]*) nogil
    CDSVector[double]* createSharedVec (int size, float& i,  EnsembleSimulation* simulation) nogil

