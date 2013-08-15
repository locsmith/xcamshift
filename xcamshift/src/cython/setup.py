#-------------------------------------------------------------------------------
# Copyright (c) 2013 Gary Thompson & The University of Leeds.
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
 
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

# change path and architecture to match local installation
xplor_path ='/home/garyt/programs/xplor-nih/2.31.0'
architecture = 'Linux_i686'


params = {'path'         : xplor_path,
          'architecture' : architecture}

templates=[ '{path}/common',
            '{path}/CDSlib',
            '{path}/arch/{architecture}/include']
include_dirs = [template.format(**params) for template in templates]
                                       
ext_modules = [Extension("shift_calculators",  ["Cython_shift_calculator.pyx", 'sharedCDSVectorFactory.cc'],
                         define_macros = [('CPLUSPLUS', '1') ,
                                          ('USE_CDS_NAMESPACE', '1')],
                        
                         language="c++",
                         extra_compile_args=["-O3"],
                         include_dirs=include_dirs),
				       
               Extension("pyEnsemblePot",  ["PyEnsemblePotProxy.cc",'pyEnsemblePot.pyx'],   
                         define_macros = [('CPLUSPLUS', '1') ,
                                           ('USE_CDS_NAMESPACE', '1')],
                        
                         language="c++",
                         extra_compile_args=["-O3"],
                         include_dirs=include_dirs),
               
               Extension("fast_segment_manager",  ["Cython_segment_manager.pyx"],
                         define_macros = [('CPLUSPLUS', '1') ,
                                          ('USE_CDS_NAMESPACE', '1')],
                        
                         language="c++",
                         extra_compile_args=["-O3"],
                         include_dirs=include_dirs),
               
               
               Extension("xCamshiftEnsemble",  ["xCamshiftEnsemble.pyx"],
                         define_macros = [('CPLUSPLUS', '1') ,
                                          ('USE_CDS_NAMESPACE', '1')],
                        
                         language="c++",
                         extra_compile_args=["-O3"],
                         include_dirs=include_dirs)
              ]

if False:
    for e in ext_modules:
        e.pyrex_directives = {"boundscheck": False}
setup(
    ext_modules=ext_modules,
    cmdclass = {'build_ext': build_ext},
)
      
