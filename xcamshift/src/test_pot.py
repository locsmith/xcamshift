'''
Created on 27 Dec 2011

@author: garyt
'''



from pyPot import PyPot

from vec3 import norm


def mult(x,y):
    print "x",x
    print "y",y
    print 
    return x*y


class BondPot(PyPot):
    '''example class to evaluate energy, derivs of a single bond'''
    def __init__(self,name,atom1,atom2,length,forcec=1):
        '''constructor - force constant is optional.'''
        PyPot.__init__(self,name) #first call base class constructor 
        self.a1 = atom1 
        self.a2 = atom2
        self.length = length
        self.forcec = forcec
        
    def calcEnergy(self):
        self.q1 = self.a1.pos() 
        self.q2 = self.a2.pos()
        self.dist = norm(self.q1-self.q2)
        return 0.5 * self.forcec * (self.dist-self.length)**2
    
    
    def calcEnergyAndDerivs(self,derivs):
        energy = self.calcEnergy()
        delta = [self.forcec * (self.dist-self.length) / self.dist]*3
        dists = (self.q1[0]-self.q2[0],
                 self.q1[1]-self.q2[1],
                 self.q1[2]-self.q2[2] )
        print 'deltas ',delta
        print 'dists ', dists
        deriv1 = map(mult, delta, dists)
        print deriv1
        print self.a1.index()
        print self.a2.index()
        derivs[self.a1.index()] = deriv1
        derivs[self.a2.index()] = [-elem for elem in deriv1]
        print derivs
        return energy