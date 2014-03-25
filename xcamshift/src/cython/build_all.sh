#!/bin/sh  
export CFLAGS=-Qunused-arguments
export CPPFLAGS=-Qunused-arguments

export xc_xplor_version='2.35.0'
export xc_xplor_root="/Users/garythompson/programs/xplor-nih/${xc_xplor_version}"
export xc_platform='Darwin_13'
export xc_architecture='x86_64'
export xc_python_version='2.7'
export xc_lower_platform=$(echo "${xc_platform}" | tr '[:upper:]' '[:lower:]')

export PYTHONPATH=/Users/garythompson/programs/cython/0.20.1

${xc_xplor_root}/bin/pyXplor setup.py build_ext --inplace 
#./build_extension.sh 
