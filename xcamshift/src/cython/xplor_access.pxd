#-------------------------------------------------------------------------------
# Copyright (c) 2013-2015 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
#
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 31 Jul 2012

@author: garyt
'''
from libc.string cimport const_char

cdef extern from *:
    ctypedef int FIVE "5"


cdef extern from "cdsString.hh":
     cdef cppclass String:
          String(char* , int) nogil

cdef extern from "fixedVector.hh":
     cdef cppclass FixedVector[T,S]:
         T& operator[](long) nogil
         T* pointer() nogil

cdef extern from "fixedVector.hh":
    cdef float norm(Vec3) nogil
    cdef float dot(Vec3&, Vec3&) nogil
    cdef Vec3 cross(Vec3&, Vec3&) nogil


cdef extern from "vec3.hh":
    cdef cppclass Vec3:
        Vec3() nogil
        Vec3(float,float,float) nogil
        Vec3(Vec3&) nogil
        float x() nogil
        float y() nogil
        float z() nogil
        Vec3& operator-(Vec3&) nogil
        Vec3& operator+(Vec3&) nogil
        Vec3& operator*(Vec3&) nogil
        float norm() nogil
        float& operator[](long) nogil
        bint operator==(Vec3&) nogil
        bint operator!=(Vec3&) nogil

cdef extern from "atom.hh":
    cdef cppclass Atom:
        const_char segmentName() nogil
        const_char residueName() nogil
        const_char atomName() nogil
        int   residueNum() nogil
        int index() nogil
        Vec3& pos() nogil

cdef extern from "dihedral.hh":
    cdef cppclass Dihedral:
        Dihedral() nogil
        Dihedral(Atom&, Atom&,
                 Atom&, Atom&) nogil
        float value() nogil

cdef extern from "bondAngle.hh":
     cdef cppclass BondAngle:
          BondAngle() nogil
          BondAngle(Atom&, Atom&,
                    Atom&,) nogil
          float value() nogil

cdef extern from "cdsVector.hh":
    cdef cppclass CDSVector[T]:
        T& data(int i) nogil
        int size() nogil
        CDSVector[T]& resize(int) nogil
        T& operator[] (int) nogil
        void set(T&)

cdef extern from "cdsList.hh":
    cdef cppclass CDSList[T]:
        int size() nogil
        CDSList[T]&  resize(int) nogil
        void append(T) nogil
        T& operator[] (int) nogil


cdef extern from "sharedCDSVectorFactory.hh" namespace "SharedCDSVectorFactory":
    void clearSharedVector(void *) nogil
    void deleteSharedVector(void *) nogil
    void resizeSharedVector(void *, int)  nogil
    void * createSharedVector (int size, float& i,  EnsembleSimulation* simulation) nogil
    void setSharedVectorValue(void * sharedVec, int offset ,double value) nogil
    void addToSharedVectorValue(void * sharedVec, int offset ,double value) nogil
    double getSharedVectorValue(void * sharedVec, int offset) nogil



cdef extern from 'simulation.hh':
    cdef cppclass Simulation:
        int id() nogil
        int numAtoms() nogil
        int size() nogil
        Atom atomByID(int index) nogil
#        CDSVector[Vec3] atomPosArr() nogil
        Vec3 atomPos(int index) nogil

cdef extern from "simulation.hh" namespace "Simulation":
    Simulation* currentSimulation() nogil
    int numSimulations() nogil

cdef extern from 'ensembleSimulation.hh':
    cdef cppclass EnsembleMemberSimulation:
        int id() nogil
        int memberIndex() nogil

    cdef cppclass EnsembleSimulation(Simulation):
         String name() nogil
         EnsembleMemberSimulation* member() nogil
         int rawID() nogil
         CDSVector[int] residueNumArr() nogil
         CDSVector[FixedVector[char,FIVE]] segmentNameArr() nogil
         CDSVector[Vec3] atomPosArr() nogil


cdef extern from "derivList.hh":
    cdef cppclass DerivList:
        void clear() nogil
        CDSVector[Vec3]& operator[](const Simulation* defaultSimulation) nogil

cdef extern from "sthead.hh":
    ctypedef double float_type
