#!/bin/sh

unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
    export CFLAGS=-Qunused-arguments
    export CPPFLAGS=-Qunused-arguments
elif [[ "$unamestr" == 'Darwin' ]]; then
    export LDSHARED="/usr/local/Cellar/gcc48/4.8.2/bin/g++-4.8 -Wl,-F. -bundle -undefined dynamic_lookup" 
    export CC=/usr/local/Cellar/gcc48/4.8.2/bin/g++-4.8 
    export LINKCC=/usr/local/Cellar/gcc48/4.8.2/bin/g++-4.8
    export LDCXXSHARED="/usr/local/Cellar/gcc48/4.8.2/bin/g++-4.8 -bundle -undefined dynamic_lookup"
    export CXX=/usr/local/Cellar/gcc48/4.8.2/bin/g++-4.8 
fi

# this is the version of xplor  used
export xc_xplor_version='2.35.0'

# i name my xplor directories after the version of the programs used
# you may want a different path here
export xc_xplor_root="/Users/garythompson/programs/xplor-nih/${xc_xplor_version}"

# this is the operating system you are using:
# for linux this should be Linux 
# for OSX this should be Darwin_XX where XX is the version of darwin  underlying the 
# OSX build to identify the Darwin version have a look at bin.Darwin_XX_x86_64 directory
# in your xplor nih directory XX_Darwin be the platform and x86_64 is the architecture
export xc_platform='Darwin_10'

# this is the processor  architecture you are using
# for linux with a 64 bit processor use x86_64, with a 32 bit processor use i686 
# for OSX with a 64 bit processor use x86_64 for a 32 bit processsor use ?not currently know?
export xc_architecture=x86_64

# this is the version of python you are compiling against
# find it by running pyXplor
# and type 
# python> import sys
# python> print sys.version 
# version will be in the first two digits ie xx.yy in xx.yy.zz
export xc_python_version='2.7'

export xc_lower_platform=$(echo "${xc_platform}" | tr '[:upper:]' '[:lower:]')

# this is where to find the cython executable
export PYTHONPATH=/Users/garythompson/programs/cython/0.20.1

# the interface used by xplor to call the ensemble code changed at some point
# this controls whether the old or the new calling convention is used
export xc_old_ensemble_interface=0

${xc_xplor_root}/bin/pyXplor setup.py build_ext --inplace 
./build_extension.sh 
