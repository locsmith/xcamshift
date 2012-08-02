'''
Created on 31 Jul 2012

@author: garyt
'''

cdef extern from "fixedVector.hh":
    cdef float norm(Vec3)
    cdef float dot(Vec3&, Vec3&)
    
cdef extern from "vec3.hh":
    cdef cppclass Vec3:
        Vec3()
        Vec3(float,float,float)
        float x()
        float y()
        float z() 
        Vec3& operator-(Vec3&)
        Vec3& operator+(Vec3&)
        float norm()
        
cdef extern from "atom.hh":
    cdef cppclass Atom:
        pass
        
cdef extern from "dihedral.hh":
    cdef cppclass Dihedral:
        Dihedral()
        Dihedral(Atom&, Atom&,
                 Atom&, Atom&)
        float value()
        
cdef extern from "cdsVector.hh" namespace "CDS":
    cdef cppclass CDSVector[T]:
        T& data(int i)
        int size()


cdef extern from 'simulation.hh':
    cdef cppclass Simulation:
        int id()
        int numAtoms()
        Atom atomByID(int index)
        CDSVector[Vec3] atomPosArr()

cdef extern from "simulation.hh" namespace "Simulation":
    Simulation* currentSimulation()
    

     