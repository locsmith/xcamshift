'''
Created on 8 Jun 2013

@author: garyt
'''
import unittest
import protocol
from ensembleSimulation import EnsembleSimulation
import ivm

class Test(unittest.TestCase):


    def setUp(self):
        ensembleSize=3
     
        esim = EnsembleSimulation("ensemble",ensembleSize)
        print "ensemble size", esim.size(),"running on",esim.numThreads(),"processors."
        
        dyn  = ivm.IVM(esim) 
        protocol.torsionTopology(dyn)
        
        potList = PotList("test pot list")
 #       pot  = XCamshiftEnsemble('test')
#        potList.append(pot)
        
        protocol.initDynamics(dyn,
                      bathTemp=400,
                      initVelocities=1,
                      potList=potList,
                      finalTime=1,
                      numSteps=10000,
                      printInterval=0)
        print 'start'
        dyn.run()
 
    def tearDown(self):
        pass


    def testEnsemblePot(self):
        pass


if __name__ == "__main__":
    unittest.main()