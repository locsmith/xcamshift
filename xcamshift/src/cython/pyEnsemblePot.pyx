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


# from pyEnsemblePotProxy cimport PyEnsemblePotProxy  as  CPPPyEnsemblePotProxy
from pyPot import PyPot
cimport cpython.ref as cpy_ref
from  xplor_access cimport Simulation, String, DerivList, EnsembleSimulation, EnsembleMemberSimulation, numSimulations, currentSimulation
from libc.stdio cimport printf
from ensembleSimulation import fromSimulation
from simulation import   Simulation_simulationByID
#

cdef extern from "instantiate.hh":
    pass

cdef class PyEnsemblePotData:
    cdef PyEnsemblePotProxy* ensemblePotProxy
    cdef String* _instance_name
    cdef String* _potential_name
    cdef Simulation *_simulation
#
    def __init__(self, instance_name, potential_name='test'):
        self._instance_name = new String(instance_name, len(instance_name))
        self._potential_name = new String(potential_name, len(potential_name))
        self._simulation = <Simulation*><size_t>currentSimulation()
        self.ensemblePotProxy = new PyEnsemblePotProxy(self._instance_name[0], self._potential_name[0],self._simulation, <cpy_ref.PyObject*>self)

    def calcEnergyAndDerivList(self,derivList):
        pointer = int(derivList.this)
        result = self.ensemblePotProxy[0].calcEnergyAndDerivs((<DerivList*><size_t>pointer)[0])
        return result

    def calcEnergy(self):
        return  self.ensemblePotProxy[0].calcEnergy()

    #TODO: should this be simulation? will it conflict with xPyPot
    def ensembleSimulation(self):
        ensemble_simulation_id = self.ensemblePotProxy[0].ensembleSimulation()[0].rawID()
        simulation = Simulation_simulationByID(ensemble_simulation_id)
        if simulation.this != <size_t>self.ensemblePotProxy[0].ensembleSimulation():
            my_this = <size_t>self.ensemblePotProxy[0].ensembleSimulation()
            lookup_this = int(Simulation_simulationByID(ensemble_simulation_id).this)
            msg = "bad ensembleSimulation lookup values of this for the same id differ id-1=%i [%i] vs id-2=%i [%i]"
            msg = msg %  (ensemble_simulation_id,my_this,simulation.rawID(), lookup_this)
            raise Exception(msg )

        return fromSimulation(simulation)


class PyEnsemblePot(PyPot,PyEnsemblePotData):

    def __init__(self,name):
        PyPot.__init__(self,name)
        PyEnsemblePotData.__init__(self, self.instanceName(),self.potName())

    def calcEnergy(self):
        return PyEnsemblePotData.calcEnergy(self)

    def calcEnergyAndDerivsMaybe0(self, Py_ssize_t derivListPtr, Py_ssize_t ensembleSimulationPtr, bint calcDerivatives):
        return 0.0

    def calcEnergyAndDerivsMaybe1(self, Py_ssize_t derivListPtr, Py_ssize_t ensembleSimulationPtr, bint calcDerivatives):
        return 0.0

    def calcEnergyAndDerivsMaybe2(self, Py_ssize_t derivListPtr, Py_ssize_t ensembleSimulationPtr, bint calcDerivatives):
        return 0.0

    def calcEnergyAndDerivsMaybe3(self, Py_ssize_t derivListPtr, Py_ssize_t ensembleSimulationPtr, bint calcDerivatives):
        return 0.0

    def calcEnergyAndDerivsMaybe4(self, Py_ssize_t derivListPtr, Py_ssize_t ensembleSimulationPtr, bint calcDerivatives):
        return 0.0
#
cdef public api double cy_call_calc_energy_and_derivs_maybe(object self, int i, DerivList* derivList, EnsembleSimulation* esim, bint calcDerivs, float* result,  int *error):
    try:
        if i == 0:
            result[0] = self.calcEnergyAndDerivsMaybe0(<Py_ssize_t> derivList, <Py_ssize_t> esim, calcDerivs)
        elif i == 1:
            result[0] = self.calcEnergyAndDerivsMaybe1(<Py_ssize_t> derivList, <Py_ssize_t> esim, calcDerivs)
        elif i == 2:
            result[0] = self.calcEnergyAndDerivsMaybe2(<Py_ssize_t> derivList, <Py_ssize_t> esim, calcDerivs)
        elif i == 3:
            result[0] = self.calcEnergyAndDerivsMaybe3(<Py_ssize_t> derivList, <Py_ssize_t> esim, calcDerivs)
        elif i == 4:
            result[0] = self.calcEnergyAndDerivsMaybe4(<Py_ssize_t> derivList, <Py_ssize_t> esim, calcDerivs)
    except Exception as  e:
        #TODO this needs to be propogated
        print e
        import traceback
        traceback.print_exc()


