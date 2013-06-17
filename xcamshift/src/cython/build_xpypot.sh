#!/bin/bash
root='xPyPot'
new_root='XPyPot'

rm _${root}.so ${root}_wrap.cxx ${root}_wrap_new_2.cpp \
   ${root}_wrap_new_2.o ${root}_wrap_new.cpp




swig -classptr -python -c++  -keyword -shadow                                              \
-I/home/garyt/programs/xplor-nih/2.31.0/python/                                  \
-I/home/garyt/programs/python/2.6.2/installation/include/python2.6                         \
-I/home/garyt/programs/xplor-nih/2.31.0/CDSlib                                   \
-I/home/garyt/programs/xplor-nih/2.31.0/common                                   \
-I/home/garyt/programs/xplor-nih/2.31.0/intVar                                   \
-I/home/garyt/programs/xplor-nih/2.31.0/nmrPot                                   \
-I/home/garyt/programs/xplor-nih/2.31.0/vmd                                      \
-I/home/garyt/programs/xplor-nih/2.31.0/surfD                                    \
-I/home/garyt/programs/xplor-nih/2.31.0/cminpack                                 \
-I/home/garyt/programs/xplor-nih/2.31.0/sparta                                   \
-I/home/garyt/programs/xplor-nih/2.31.0/devel                                    \
-I/home/garyt/programs/xplor-nih/2.31.0/fortlib                                  \
-I. -I/home/garyt/programs/xplor-nih/2.31.0/common                               \
${root}.i

#sed replace
sed -e 's/SWIG_/SWIGPY_/g' ${root}_wrap.cxx > ${new_root}_wrap_new.cpp

#include
~/programs/xplor-nih/2.31.0/bin/includeCC ${new_root}_wrap_new.cpp --template-dir  \
/home/garyt/programs/xplor-nih/2.31.0/CDSlib --cc 'c++' -DX_MMAP_FLAGS=0        \
-DFORTRAN_INIT -O3 -DLINUX -D_REENTRANT -DNDEBUG                                          \
-I/home/garyt/programs/xplor-nih/2.31.0/python/                                 \
-I/home/garyt/programs/xplor-nih/2.31.0/arch/Linux_i686/include                 \
-DSWIG_VERSION=20004 -I/usr/share/swig2.0/python -DCPLUSPLUS -DUSE_CDS_NAMESPACE          \
-I/home/garyt/programs/xplor-nih/2.31.0/python/                                 \
-I/home/garyt/programs/xplor-nih/2.31.0/arch/Linux_i686/include -DSWIGPY_GLOBAL \
-I. -I/home/garyt/programs/python/2.6.2/installation/include/python2.6                    \
-I/home/garyt/programs/xplor-nih/2.31.0/CDSlib                                  \
-I/home/garyt/programs/xplor-nih/2.31.0/common                                  \
-I/home/garyt/programs/xplor-nih/2.31.0/intVar                                  \
-I/home/garyt/programs/xplor-nih/2.31.0/nmrPot                                  \
-I/home/garyt/programs/xplor-nih/2.31.0/vmd                                     \
-I/home/garyt/programs/xplor-nih/2.31.0/surfD                                   \
-I/home/garyt/programs/xplor-nih/2.31.0/cminpack                                \
-I/home/garyt/programs/xplor-nih/2.31.0/sparta                                  \
-I/home/garyt/programs/xplor-nih/2.31.0/devel                                   \
-I/home/garyt/programs/xplor-nih/2.31.0/fortlib                                 \
-DNIHXPLOR_VERSION='"2.31-custom"' -DPYTHON_VERSION='"2.6"'                               \
-DSWIGPY_PYTHON_SILENT_MEMLEAK > ${new_root}_wrap_new_2.cpp

#compile

gcc -pthread -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall                     \
-Wstrict-prototypes -fPIC -DCPLUSPLUS=1 -DUSE_CDS_NAMESPACE=1  -DSWIGPY_NOINCLUDE          \
-fno-use-cxa-get-exception-ptr  -DX_MMAP_FLAGS=0 -DFORTRAN_INIT -O3 -DLINUX  -D_REENTRANT  \
-DNIHXPLOR_VERSION='"2.31-custom"' -DPYTHON_VERSION='"2.6"'                                \
-DSWIGPY_PYTHON_SILENT_MEMLEAK -DSWIGPY_NOINCLUDE                                          \
-I/home/garyt/programs/xplor-nih/2.31.0/common                                             \
-I/home/garyt/programs/xplor-nih/2.31.0/CDSlib                                             \
-I/home/garyt/programs/xplor-nih/2.31.0/arch/Linux_i686/include                            \
-I/home/garyt/programs/xplor-nih/2.31.0/python/bin.Linux_i686/include/python2.6            \
-I.                                                                                        \
-I/home/garyt/programs/xplor-nih/2.31.0/python                                             \
-I/home/garyt/programs/xplor-nih/2.31.0/intVar                                             \
-I/home/garyt/programs/xplor-nih/2.31.0/vmd                                      \
-I/home/garyt/programs/xplor-nih/2.31.0/python/swig_wrappers                               \
-c ${root}.cc -o build/temp.linux-i686-2.6/${root}.o

gcc -pthread -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall                     \
-Wstrict-prototypes -fPIC -DCPLUSPLUS=1 -DUSE_CDS_NAMESPACE=1  -DSWIGPY_NOINCLUDE          \
-fno-use-cxa-get-exception-ptr  -DX_MMAP_FLAGS=0 -DFORTRAN_INIT -O3 -DLINUX  -D_REENTRANT  \
-DNIHXPLOR_VERSION='"2.31-custom"' -DPYTHON_VERSION='"2.6"'                                \
-DSWIGPY_PYTHON_SILENT_MEMLEAK -DSWIGPY_NOINCLUDE                                          \
-I/home/garyt/programs/xplor-nih/2.31.0/common                                             \
-I/home/garyt/programs/xplor-nih/2.31.0/CDSlib                                             \
-I/home/garyt/programs/xplor-nih/2.31.0/arch/Linux_i686/include                            \
-I/home/garyt/programs/xplor-nih/2.31.0/python/bin.Linux_i686/include/python2.6            \
-I.                                                                                        \
-I/home/garyt/programs/xplor-nih/2.31.0/python                                             \
-I/home/garyt/programs/xplor-nih/2.31.0/intVar                                             \
-I/home/garyt/programs/xplor-nih/2.31.0/vmd                                      \
-c ${new_root}_wrap_new_2.cpp

#link
gcc -shared build/temp.linux-i686-2.6/${root}.o ${new_root}_wrap_new_2.o \
    -lstdc++ -o _${root}.so
