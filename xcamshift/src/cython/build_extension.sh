#!/bin/sh -x


rm _PyEnsemblePotProxy.so pyEnsemblePotProxy_wrap.cxx pyEnsemblePotProxy_wrap_new_2.cpp \
   pyEnsemblePotProxy_wrap_new_2.o pyEnsemblePotProxy_wrap_new.cpp



swig -classptr -python -c++  -keyword -shadow                                              \
-I${xc_xplor_root}/python/                                  \
-I${xc_xplor_root}/python/bin.${xc_platform}_${xc_architecture}/include/python${xc_python_version}                            \
-I${xc_xplor_root}/CDSlib                                   \
-I${xc_xplor_root}/common                                   \
-I${xc_xplor_root}/intVar                                   \
-I${xc_xplor_root}/nmrPot                                   \
-I${xc_xplor_root}/vmd                                      \
-I${xc_xplor_root}/surfD                                    \
-I${xc_xplor_root}/cminpack                                 \
-I${xc_xplor_root}/sparta                                   \
-I${xc_xplor_root}/devel                                    \
-I${xc_xplor_root}/fortlib                                  \
-I. -I${xc_xplor_root}/common                               \
pyEnsemblePotProxy.i

#sed replace
sed -e 's/SWIG_/SWIGPY_/g' pyEnsemblePotProxy_wrap.cxx > pyEnsemblePotProxy_wrap_new.cpp

#include
${xc_xplor_root}/bin/includeCC pyEnsemblePotProxy_wrap_new.cpp --template-dir  \
${xc_xplor_root}/CDSlib --cc 'c++' -DX_MMAP_FLAGS=0        \
-DFORTRAN_INIT -O3 -DLINUX -D_REENTRANT -DNDEBUG                                          \
-I${xc_xplor_root}/python/                                 \
-I${xc_xplor_root}/arch/${xc_platform}_${xc_architecture}/include                 \
-DSWIG_VERSION=20004 -I/usr/share/swig2.0/python -DCPLUSPLUS -DUSE_CDS_NAMESPACE          \
-I${xc_xplor_root}/python/                                 \
-I${xc_xplor_root}/arch/${xc_platform}_${xc_architecture}/include -DSWIGPY_GLOBAL \
-I. -I${xc_xplor_root}/python/bin.${xc_platform}_${xc_architecture}/include/python${xc_python_version}                      \
-I${xc_xplor_root}/CDSlib                                  \
-I${xc_xplor_root}/common                                  \
-I${xc_xplor_root}/intVar                                  \
-I${xc_xplor_root}/nmrPot                                  \
-I${xc_xplor_root}/vmd                                     \
-I${xc_xplor_root}/surfD                                   \
-I${xc_xplor_root}/cminpack                                \
-I${xc_xplor_root}/sparta                                  \
-I${xc_xplor_root}/devel                                   \
-I${xc_xplor_root}/fortlib                                 \
-DNIHXPLOR_VERSION='"${xc_xplor_version}"' -DPYTHON_VERSION='"${xc_python_version} "'                               \
-DSWIGPY_PYTHON_SILENT_MEMLEAK > pyEnsemblePotProxy_wrap_new_2.cpp

#compile
gcc -pthread -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall                     \
-Wstrict-prototypes -fPIC -DCPLUSPLUS=1 -DUSE_CDS_NAMESPACE=1  -DSWIGPY_NOINCLUDE          \
-fno-use-cxa-get-exception-ptr  -DX_MMAP_FLAGS=0 -DFORTRAN_INIT -O3 -DLINUX  -D_REENTRANT  \
-DNIHXPLOR_VERSION='"${xc_xplor_version}"' -DPYTHON_VERSION='"${xc_python_version} "'                                \
-DSWIGPY_PYTHON_SILENT_MEMLEAK -DSWIGPY_NOINCLUDE                                          \
-I${xc_xplor_root}/common                                             \
-I${xc_xplor_root}/CDSlib                                             \
-I${xc_xplor_root}/arch/${xc_platform}_${xc_architecture}/include                            \
-I${xc_xplor_root}/python/bin.${xc_platform}_${xc_architecture}/include/python${xc_python_version}               \
-I.                                                                                        \
-I${xc_xplor_root}/python                                             \
-I${xc_xplor_root}/intVar                                             \
-I${xc_xplor_root}/vmd                                      \
-c pyEnsemblePotProxy_wrap_new_2.cpp

#link
gcc -shared build/temp.${xc_lower_platform}-${xc_architecture}-2.6/PyEnsemblePotProxy.o pyEnsemblePotProxy_wrap_new_2.o \
    -lstdc++ -o _PyEnsemblePotProxy.so
