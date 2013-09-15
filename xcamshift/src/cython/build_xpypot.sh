#!/bin/bash
root='xPyPot'
new_root='XPyPot'
xplor_version=$1
xplor_root=$2
architecture=$3
python_version=$4
rm _${root}.so ${root}_wrap.cxx ${root}_wrap_new_2.cpp \
   ${root}_wrap_new_2.o ${root}_wrap_new.cpp




swig -classptr -python -c++  -keyword -shadow                                              \
-I${xplor_root}/python/                                  \
-I${xplor_root}/python/bin.${architecture}/include/python${python_version}         \
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
${root}.i

#sed replace
sed -e 's/SWIG_/SWIGPY_/g' ${root}_wrap.cxx > ${new_root}_wrap_new.cpp

#include
${xplor_root}/bin/includeCC ${new_root}_wrap_new.cpp --template-dir  \
${xplor_root}/CDSlib --cc 'c++' -DX_MMAP_FLAGS=0        \
-DFORTRAN_INIT -O3 -DLINUX -D_REENTRANT -DNDEBUG                                          \
-I${xplor_root}/python/                                 \
-I${xplor_root}/arch/Linux_i686/include                 \
-DSWIG_VERSION=20004 -I/usr/share/swig2.0/python -DCPLUSPLUS -DUSE_CDS_NAMESPACE          \
-I${xplor_root}/python/                                 \
-I${xplor_root}/arch/Linux_i686/include -DSWIGPY_GLOBAL \
-I. -I${xplor_root}/python/bin.${architecture}/include/python${python_version}                    \
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
-DSWIGPY_PYTHON_SILENT_MEMLEAK > ${new_root}_wrap_new_2.cpp

#compile

gcc -pthread -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall                     \
-Wstrict-prototypes -fPIC -DCPLUSPLUS=1 -DUSE_CDS_NAMESPACE=1  -DSWIGPY_NOINCLUDE          \
-fno-use-cxa-get-exception-ptr  -DX_MMAP_FLAGS=0 -DFORTRAN_INIT -O3 -DLINUX  -D_REENTRANT  \
-DNIHXPLOR_VERSION='"2.31-custom"' -DPYTHON_VERSION='"${python_version}"'                                \
-DSWIGPY_PYTHON_SILENT_MEMLEAK -DSWIGPY_NOINCLUDE                                          \
-I${xplor_root}/common                                             \
-I${xplor_root}/CDSlib                                             \
-I${xplor_root}/arch/Linux_i686/include                            \
-I${xplor_root}/python/bin.Linux_i686/include/python${python_version}            \
-I.                                                                                        \
-I${xplor_root}/python                                             \
-I${xplor_root}/intVar                                             \
-I${xplor_root}/vmd                                      \
-I${xplor_root}/python/swig_wrappers                               \
-c ${root}.cc -o build/temp.linux-i686-${python_version}/${root}.o

gcc -pthread -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall                     \
-Wstrict-prototypes -fPIC -DCPLUSPLUS=1 -DUSE_CDS_NAMESPACE=1  -DSWIGPY_NOINCLUDE          \
-fno-use-cxa-get-exception-ptr  -DX_MMAP_FLAGS=0 -DFORTRAN_INIT -O3 -DLINUX  -D_REENTRANT  \
-DNIHXPLOR_VERSION='"2.31-custom"' -DPYTHON_VERSION='"${python_version}"'                                \
-DSWIGPY_PYTHON_SILENT_MEMLEAK -DSWIGPY_NOINCLUDE                                          \
-I${xplor_root}/common                                             \
-I${xplor_root}/CDSlib                                             \
-I${xplor_root}/arch/Linux_i686/include                            \
-I${xplor_root}/python/bin.Linux_i686/include/python${python_version}            \
-I.                                                                                        \
-I${xplor_root}/python                                             \
-I${xplor_root}/intVar                                             \
-I${xplor_root}/vmd                                      \
-c ${new_root}_wrap_new_2.cpp

#link
gcc -shared build/temp.linux-i686-${python_version}/${root}.o ${new_root}_wrap_new_2.o \
    -lstdc++ -o _${root}.so
