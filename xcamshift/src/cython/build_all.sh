#!/bin/sh  

export xc_xplor_version='2.31.0'
export xc_xplor_root="/home/garyt/programs/xplor-nih/${xc_xplor_version}"
export xc_architecture='Linux_i686'
export xc_python_version='2.6'

export PYTHONPATH=/home/garyt/programs/cython/0.18.0

./build_xpypot.sh $xc_xplor_version $xc_xplor_root $xc_architecture $xc_python_version
${xc_xplor_root}/bin/pyXplor setup.py build_ext --inplace 
./build_extension.sh $xc_xplor_version $xc_xplor_root $xc_architecture $xc_python_version
