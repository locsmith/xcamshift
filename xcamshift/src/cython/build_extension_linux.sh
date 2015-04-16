#!/bin/bash -x

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

root='pyEnsemblePotProxy'
xc_platform_base="`echo $xc_platform_base|tr '[a-z]' '[A-Z]'`"

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
${xc_xplor_root}/bin/includeCC ${root}_wrap_new.cpp      \
--template-dir  ${xc_xplor_root}/CDSlib                  \
    --cc $CXX                                            \
    -O3                                                  \
    -DX_MMAP_FLAGS=0                                     \
    -DFORTRAN_INIT                                       \
    -DNDEBUG                                             \
    -D_REENTRANT                                         \
    -D${xc_platform_base}                                \
    -DSWIGPY_GLOBAL                                      \
    -DSWIG_VERSION=20004                                 \
    -DCPLUSPLUS                                          \
    -DUSE_CDS_NAMESPACE                                  \
    ${EXTRA_CC_FLAG}                                     \
-I${xc_xplor_root}/python/                                         \
-I${xc_xplor_root}/arch/${xc_platform}_${xc_architecture}/include  \
-I${SWIG_INCLUDE_PATH}/python                                      \
-I${xc_xplor_root}/python/                                 \
-I${xc_xplor_root}/arch/${xc_platform}_${xc_architecture}/include \
-I.                                                       \
-I${xc_xplor_root}/python/bin.${xc_platform}_${xc_architecture}/include/python${xc_python_version}                      \
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

$CXX                                                     \
  -O3                                                    \
  -DX_MMAP_FLAGS=0                                       \
  -DFORTRAN_INIT                                         \
  -DNDEBUG                                               \
  -Wall                                                  \
  -D_REENTRANT                                           \
  -DCPLUSPLUS=1                                          \
  -DUSE_CDS_NAMESPACE=1                                  \
  -DSWIGPY_NOINCLUDE                                     \
  -D${xc_platform_base}                                   \
  -pthread                                               \
  -fno-strict-aliasing                                   \
  -g                                                     \
  -fwrapv                                                \
  -Wstrict-prototypes                                    \
  -fPIC                                                  \
  -DSWIG_VERSION=20004                                   \
  -DSWIGPY_GLOBAL                                        \
-I${xc_xplor_root}/common                                \
-I${xc_xplor_root}/CDSlib                                \
-I${xc_xplor_root}/arch/${xc_platform}_${xc_architecture}/include                                  \
-I${xc_xplor_root}/python/bin.${xc_platform}_${xc_architecture}/include/python${xc_python_version} \
-I.                                                                                        \
-I${xc_xplor_root}/python                                \
-I${xc_xplor_root}/intVar                                \
-I${xc_xplor_root}/vmd                                   \
  -DSWIGPY_PYTHON_SILENT_MEMLEAK                         \
  -DNIHXPLOR_VERSION='"${xc_xplor_version}"'             \
  -DPYTHON_VERSION='"${xc_python_version}"'              \
-c  ${root}_wrap_new_2.cpp                               \
-o _${root}.o

#link
$CXX                                                     \
-shared build/temp.${xc_lower_platform}-${xc_architecture}-2.6/PyEnsemblePotProxy.o
  pyEnsemblePotProxy_wrap_new_2.o                        \
    -lstdc++                                             \
-o _${root}.so
