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

export xc_xplor_version='2.35.0'
export xc_xplor_root="/Users/garythompson/programs/xplor-nih/${xc_xplor_version}"
export xc_platform='Darwin_10'
export xc_architecture='x86_64'
export xc_python_version='2.7'
export xc_lower_platform=$(echo "${xc_platform}" | tr '[:upper:]' '[:lower:]')

export PYTHONPATH=/Users/garythompson/programs/cython/0.20.1

${xc_xplor_root}/bin/pyXplor setup.py build_ext --inplace 
#./build_extension.sh 
