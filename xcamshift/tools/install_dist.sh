#!/bin/bash

dist_dir=$1
install_dir=$2
if [ $# != 2 ] ; then
    echo 'ERROR: $* expects two parameters <dist_dir> and <install_dir>' > /dev/stderr
    echo '       dist_dir: the distribution directory' > /dev/stderr
    echo '       install_dir: the installation directory' > /dev/stderr
    echo 'exiting...' > /dev/stderr
    exit -1
fi
if [ !  -r ${dist_dir}/python/xcamshift.py ] ; then
      echo ERROR: the distribution directory doesn\'t appear to be correct\; it should contain the file \<dist_dir\>/python/xcamshift.py > /dev/stderr
      echo "       your dist_dir was: $dist_dir" > /dev/stderr
      echo exiting... > /dev/stderr
      exit -1
fi
if [ !  -r ${install_dir}/bin/xplor ] ; then
      echo ERROR: the installation directory doesn\'t appear to be correct\; it should contain the file   \<install_dir\>/bin/xplor > /dev/stderr
      echo "       your install_dir was: $install_dir" > /dev/stderr
      echo exiting... > /dev/stderr
      exit -1
fi


echo  installing xcamshift:
echo  
echo    distribution directory: ${dist_dir}
echo    installation directory: ${install_dir}
echo
cp -r ${dist_dir}/modules/* ${install_dir}/python/bin.Darwin_*_x86_64/
cp -r ${dist_dir}/database/XCamShift ${install_dir}/databases/
cp -r ${dist_dir}/python/ ${install_dir}/python/
cp -r ${dist_dir}/test_data  ${install_dir}
echo testing distribution...
cd $install_dir
./bin/pyXplor ./python/test/test_suite.py

