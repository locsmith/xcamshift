#-------------------------------------------------------------------------------
# Copyright (c) 2013 gary thompson.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 23 Dec 2011

@author: garyt
''' 
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
 

ext_modules = [Extension("test_distance",  ["test_distance.pyx"],
                         define_macros = [('CPLUSPLUS', '1') ,
                                          ('USE_CDS_NAMESPACE', '1')],
                        
                         language="c++",
                         
                         include_dirs=['/home/garyt/programs/xplor-nih/2.31.0/common',
                                       '/home/garyt/programs/xplor-nih/2.31.0/CDSlib',
                                       '/home/garyt/programs/xplor-nih/2.31.0/arch/Linux_i686/include']),
               
               Extension("shift_calculators",  ["Cython_shift_calculator.pyx"],
                         define_macros = [('CPLUSPLUS', '1') ,
                                          ('USE_CDS_NAMESPACE', '1')],
                        
                         language="c++",
                         
                         include_dirs=['/home/garyt/programs/xplor-nih/2.31.0/common',
                                       '/home/garyt/programs/xplor-nih/2.31.0/CDSlib',
                                       '/home/garyt/programs/xplor-nih/2.31.0/arch/Linux_i686/include']),
  
               Extension("fast_segment_manager",  ["Cython_segment_manager.pyx"],
                         define_macros = [('CPLUSPLUS', '1') ,
                                          ('USE_CDS_NAMESPACE', '1')],
                        
                         language="c++",
                         
                         include_dirs=['/home/garyt/programs/xplor-nih/2.31.0/common',
                                       '/home/garyt/programs/xplor-nih/2.31.0/CDSlib',
                                       '/home/garyt/programs/xplor-nih/2.31.0/arch/Linux_i686/include'])
               ]

if False:
    for e in ext_modules:
        e.pyrex_directives = {"boundscheck": False}
setup(
    ext_modules=ext_modules,
    cmdclass = {'build_ext': build_ext},
)
      
