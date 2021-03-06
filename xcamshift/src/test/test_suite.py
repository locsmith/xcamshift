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
Created on 25 Jan 2012

@author: garyt
'''
import unittest2
from test.test_segment_manager import Test_segment_manager
from test.test_observed_chemical_shifts import TestObservedShiftTable
from test.test_xcamshift import TestXcamshift
from test.Test_table_manager import Test_table_manager
from test.test_component_list import Test_component_list,Test_native_component_list
from test.test_xcamshift_afa import TestXcamshiftAFA
from test.test_xcamshift_a4 import TestXcamshiftA4
from test.test_python_utils import Test_python_utils
#TODO move importesr test to its own suite
from test.test_table_importers import Test_table_importers
from test.test_xcamshift_aga import TestXcamshiftAGA
from test.test_xcamshift_vin import TestXcamshiftVIN
from test.test_xcamshift_agfa import TestXcamshiftAGFA
from test.test_xcamshift_agaga import TestXcamshiftAGAGA
from test import test_xcamshift_gb3
from test.test_cython_non_bonded_list import Test_cython_non_boned_list
from test.test_xcamshift_acaggaca import TestXcamshiftACAGGACA
from test.test_xcamshift_hbond_ingktlkg import TestXcamshiftHBondINGKTLKG
from test.test_xcamshift_print_shifts import TestXcamshiftPrintShifts

def load_tests(loader, tests, pattern):

    test_list = (
                TestXcamshiftAFA,
                TestObservedShiftTable,
                Test_segment_manager,
                TestXcamshift,
                Test_table_manager,
                Test_component_list,
                Test_native_component_list,
                TestXcamshiftA4,
                Test_python_utils,
                # todo move to its own suite
                #Test_table_importers,
                TestXcamshiftAGA,
                TestXcamshiftVIN,
                TestXcamshiftAGFA,
                TestXcamshiftAGAGA,
                Test_cython_non_boned_list,
                TestXcamshiftACAGGACA,
                TestXcamshiftHBondINGKTLKG,
                TestXcamshiftPrintShifts
)

    suite = unittest2.TestSuite()

    test_list = [loader.loadTestsFromTestCase(test_case) for test_case in test_list]
    test_list.append(loader.loadTestsFromName('TestXcamshiftGB3.test_total_forces_and_energy_10_step', test_xcamshift_gb3))

    suite.addTests(test_list)

    return suite

def run_tests():
    unittest2.main(verbosity=10)

if __name__=="__main__":
    run_tests()


