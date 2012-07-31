'''
Created on 23 Dec 2011

@author: garyt
''' 
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
 
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("test_distance",  ["test_distance.pyx"],
                             define_macros = [('CPLUSPLUS', '1') ,
                                              ('USE_CDS_NAMESPACE', '1')],
                            
                             language="c++",
                             
                             include_dirs=['/home/garyt/programs/xplor-nih/2.31.0/common',
                                           '/home/garyt/programs/xplor-nih/2.31.0/CDSlib',
                                           '/home/garyt/programs/xplor-nih/2.31.0/arch/Linux_i686/include'])
                   ]
      )