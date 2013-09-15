#!/bin/sh

xplor_version=$1
xplor_root=$2
architecture=$3
python_version=$4


rm _PyEnsemblePotProxy.so pyEnsemblePotProxy_wrap.cxx pyEnsemblePotProxy_wrap_new_2.cpp \
   pyEnsemblePotProxy_wrap_new_2.o pyEnsemblePotProxy_wrap_new.cpp




swig -classptr -python -c++  -keyword -shadow                                              \
-I${xplor_root}/python/                                  \
-I${xplor_root}/python/bin.${architecture}/include/python${python_version}                           \
-I${xplor_root}/CDSlib                                   \
-I${xplor_root}/common                                   \
-I${xplor_root}/intVar                                   \
-I${xplor_root}/nmrPot                                   \
-I${xplor_root}/vmd                                      \
-I${xplor_root}/surfD                                    \
-I${xplor_root}/cminpack                                 \
-I${xplor_root}/sparta                                   \
-I${xplor_root}/devel                                    \
-I${xplor_root}/fortlib                                  \
-I. -I${xplor_root}/common                               \
pyEnsemblePotProxy.i

#sed replace
sed -e 's/SWIG_/SWIGPY_/g' pyEnsemblePotProxy_wrap.cxx > pyEnsemblePotProxy_wrap_new.cpp

#include
${xplor_root}/bin/includeCC pyEnsemblePotProxy_wrap_new.cpp --template-dir  \
${xplor_root}/CDSlib --cc 'c++' -DX_MMAP_FLAGS=0        \
-DFORTRAN_INIT -O3 -DLINUX -D_REENTRANT -DNDEBUG                                          \
-I${xplor_root}/python/                                 \
-I${xplor_root}/arch/Linux_i686/include                 \
-DSWIG_VERSION=20004 -I/usr/share/swig2.0/python -DCPLUSPLUS -DUSE_CDS_NAMESPACE          \
-I${xplor_root}/python/                                 \
-I${xplor_root}/arch/Linux_i686/include -DSWIGPY_GLOBAL \
-I. -I${xplor_root}/python/bin.${architecture}/include/python${python_version}                     \
-I${xplor_root}/CDSlib                                  \
-I${xplor_root}/common                                  \
-I${xplor_root}/intVar                                  \
-I${xplor_root}/nmrPot                                  \
-I${xplor_root}/vmd                                     \
-I${xplor_root}/surfD                                   \
-I${xplor_root}/cminpack                                \
-I${xplor_root}/sparta                                  \
-I${xplor_root}/devel                                   \
-I${xplor_root}/fortlib                                 \
-DNIHXPLOR_VERSION='"2.31-custom"' -DPYTHON_VERSION='"${python_version}"'                               \
-DSWIGPY_PYTHON_SILENT_MEMLEAK > pyEnsemblePotProxy_wrap_new_2.cpp

#compile
gcc -pthread -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall                     \
-Wstrict-prototypes -fPIC -DCPLUSPLUS=1 -DUSE_CDS_NAMESPACE=1  -DSWIGPY_NOINCLUDE          \
-fno-use-cxa-get-exception-ptr  -DX_MMAP_FLAGS=0 -DFORTRAN_INIT -O3 -DLINUX  -D_REENTRANT  \
-DNIHXPLOR_VERSION='"2.31-custom"' -DPYTHON_VERSION='"${python_version}"'                                \
-DSWIGPY_PYTHON_SILENT_MEMLEAK -DSWIGPY_NOINCLUDE                                          \
-I${xplor_root}/common                                             \
-I${xplor_root}/CDSlib                                             \
-I${xplor_root}/arch/Linux_i686/include                            \
-I${xplor_root}/python/bin.${architecture}/include/python${python_version}              \
-I.                                                                                        \
-I${xplor_root}/python                                             \
-I${xplor_root}/intVar                                             \
-I${xplor_root}/vmd                                      \
-c pyEnsemblePotProxy_wrap_new_2.cpp

#link
gcc -shared build/temp.linux-i686-2.6/PyEnsemblePotProxy.o pyEnsemblePotProxy_wrap_new_2.o \
    -lstdc++ -o _PyEnsemblePotProxy.so
