#!/bin/bash -x
root='pyEnsemblePotProxy'

rm _${root}.so ${root}_wrap.cxx ${root}_wrap_new_2.cpp \
   ${root}_wrap_new_2.o ${root}_wrap_new.cpp
mkdir -p  build/temp.${xc_platform}-${xc_architecture}-${xc_python_version}




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
${root}.i

#sed replace
sed -e 's/SWIG_/SWIGPY_/g' ${root}_wrap.cxx > ${root}_wrap_new.cpp

#include
${xc_xplor_root}/bin/includeCC ${root}_wrap_new.cpp    \
--template-dir  ${xc_xplor_root}/CDSlib                    \
    --cc 'g++'                                           \
    -O3                                                  \
    -ffast-math                                          \
    -funroll-loops                                       \
    -DX_MMAP_FLAGS=0                                     \
    -DFORTRAN_INIT                                       \
    -fno-common                                          \
    -DDARWIN                                             \
    -D_REENTRANT                                         \
    -Wall                                                \
    -Wno-unused                                          \
    -DNDEBUG                                             \
-I${xc_xplor_root}/python/                                 \
-I${xc_xplor_root}/arch/${xc_platform}_${xc_architecture}/include                 \
-DSWIG_VERSION=20004                                      \
-I/usr/local/Cellar/swig/2.0.4/share/swig/2.0.4/python    \
-DCPLUSPLUS                                               \
-DUSE_CDS_NAMESPACE                                       \
-I${xc_xplor_root}/python/                                 \
-I${xc_xplor_root}/arch/${xc_platform}_${xc_architecture}/include \
-DSWIGPY_GLOBAL                                           \
-I.                                                       \
-I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 \
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
-DNIHXPLOR_VERSION="${xc_xplor_version}"                   \
-DPYTHON_VERSION="${xc_python_version}"                    \
-DSWIGPY_PYTHON_SILENT_MEMLEAK > ${root}_wrap_new_2.cpp

#compile

g++ -c  ${root}_wrap_new_2.cpp                         \
  -O3                                                    \
  -ffast-math                                            \
  -funroll-loops                                         \
  -DX_MMAP_FLAGS=0                                       \
  -DFORTRAN_INIT                                         \
  -fno-common                                            \
  -DDARWIN                                               \
  -D_REENTRANT                                           \
  -Wall                                                  \
  -Wno-unused                                            \
  -DNDEBUG                                               \
-I${xc_xplor_root}/common                                             \
-I${xc_xplor_root}/CDSlib                                             \
-I${xc_xplor_root}/python/                               \
-I${xc_xplor_root}/arch/${xc_platform}_${xc_architecture}/include                            \
-I${xc_xplor_root}/python/bin.${xc_platform}_${xc_architecture}/include/python${xc_python_version}               \
-DSWIG_VERSION=20004                                       \
-I/usr/local/Cellar/swig/2.0.4/share/swig/2.0.4/python     \
-DCPLUSPLUS                                                \
-DUSE_CDS_NAMESPACE                                        \
-I${xc_xplor_root}/python/                                 \
-I${xc_xplor_root}/arch/Darwin_13_x86_64/include           \
-DSWIGPY_GLOBAL                                            \
-I.                                                                                        \
-I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 \
-I${xc_xplor_root}/intVar                                             \
-I${xc_xplor_root}/vmd                                      \
-I${xc_xplor_root}/python/swig_wrappers                    \
  -DNIHXPLOR_VERSION='"2.35-custom"'                     \
  -DPYTHON_VERSION='"2.7"'                               \
  -DSWIGPY_PYTHON_SILENT_MEMLEAK                         \
  -DSWIGPY_NOINCLUDE                                     \
-o _${root}.o

g++                                                        \
-dynamiclib                                                \
-flat_namespace                                            \
-undefined suppress                                        \
-single_module                                             \
_${root}.o                                                  \
-o  lib${root}.dylib                                       \
-L${xc_xplor_root}/bin.Darwin_13_x86_64/                   \
-u _PyMac_Error                                            \
/System/Library/Frameworks/Python.framework/Versions/2.7/Python  \
-ldl                                                       \
-framework  CoreFoundation                                 \
-L${xc_xplor_root}/bin.Darwin_13_x86_64/                   \
-lswigpy-xplor

g++                                        \
-bundle                                    \
-flat_namespace                            \
-undefined suppress                        \
_${root}.o                                 \
-o _${root}.so                             \
-L${xc_xplor_root}/bin.Darwin_13_x86_64/   \
-lcommon                                   \
-lswigpy-xplor                            \
-lpy                                       \
-u _PyMac_Error /System/Library/Frameworks/Python.framework/Versions/2.7/Python \
-framework CoreFoundation                  



