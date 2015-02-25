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

from cpython.ref cimport PyObject
from xplor_access cimport String, Simulation, currentSimulation, float_type, DerivList, EnsembleSimulation



cdef extern from "PyEnsemblePotProxy.hh":


    cdef cppclass PyEnsemblePotProxy:
        PyEnsemblePotProxy(String& potName,  String& instanceName, Simulation* simulation, PyObject *object)
        float_type calcEnergyAndDerivs(DerivList&) nogil
        float calcEnergy() nogil
        EnsembleSimulation *ensembleSimulation() nogil


