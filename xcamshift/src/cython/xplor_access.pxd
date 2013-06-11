'''
Created on 31 Jul 2012

@author: garyt
'''
from libc.string cimport const_char

cdef extern from "cdsString.hh":
     cdef cppclass String:
          String(char* , int)
                  

cdef extern from "fixedVector.hh":
    cdef float norm(Vec3)
    cdef float dot(Vec3&, Vec3&)
    cdef Vec3 cross(Vec3&, Vec3&)

    
cdef extern from "vec3.hh":
    cdef cppclass Vec3:
        Vec3()
        Vec3(float,float,float)
        Vec3(Vec3&)
        float x()
        float y()
        float z() 
        Vec3& operator-(Vec3&)
        Vec3& operator+(Vec3&)
        Vec3& operator*(Vec3&)
        float norm()
        float& operator[](long)
        bint operator==(Vec3&)
        bint operator!=(Vec3&)
        
cdef extern from "atom.hh":
    cdef cppclass Atom:
        const_char segmentName()
        const_char residueName()
        const_char atomName()
        int   residueNum()
        int index()
        Vec3& pos()
        
cdef extern from "dihedral.hh":
    cdef cppclass Dihedral:
        Dihedral()
        Dihedral(Atom&, Atom&,
                 Atom&, Atom&)
        float value()
        
cdef extern from "cdsVector.hh":
    cdef cppclass CDSVector[T]:
        T& data(int i)
        int size()
        CDSVector[T]& resize(int)
        T& operator[] (int) 


cdef extern from 'simulation.hh':
    cdef cppclass Simulation:
        int id()
        int numAtoms()
        Atom atomByID(int index)
#        CDSVector[Vec3] atomPosArr()
        Vec3 atomPos(int index)

cdef extern from "simulation.hh" namespace "Simulation":
    Simulation* currentSimulation()
    

cdef extern from "pot.hh":
    cdef cppclass DerivList:
        pass

cdef extern from "sthead.hh":
    ctypedef double float_type
