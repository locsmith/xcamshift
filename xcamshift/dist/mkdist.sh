#!/bin/bash

rm -rf python
cp ../../xcamshift/tools/install_dist.sh .
cp -r ../../xcamshift/licenses .
cp -r ../../xcamshift/src python
cp -r ../../xcamshift/lib/PyYAML/3.10/installation/lib/python2.7/site-packages/yaml python
cp -r ../../xcamshift/lib/unittest2/0.5.1/installation/lib/python2.7/site-packages/unittest2-0.5.1-py2.7.egg/unittest2 python
unzip  -d python ../../xcamshift/lib/nanotime/0.5.2/installation/lib/python2.7/site-packages/nanotime-0.5.2-py2.7.egg nanotime/__init__.py
find python -name '*.pyc' -exec rm {} \;
find python -name '*.pyo' -exec rm {} \;

cd python
rm -rf cython test_pot.py cython_backup.tgz main.py cython_working.tgz case.py tools pyEnsemblePot.so camshift2shifts.py test_pypot.py
cd ..

find . -name '*.pickle' -exec rm {} \;
find . -name 'forces*.yaml' -exec rm {} \;
rmdir python/test/pickle
find . -name '*.marshal' -exec rm {} \;

pandoc  -f markdown -t plain  -o  README ../../xcamshift/web\ site/xcamshift.md
pandoc  -f markdown -t plain  -o  INSTALL ../../xcamshift/web\ site/install.md
