#!/bin/bash


XCAMSHIFT_HOME=$HOME/Dropbox/git/xcamshift/xcamshift
export PYTHONPATH=$XCAMSHIFT_HOME/src/:$XCAMSHIFT_HOME/lib/PyYAML/3.10/installation/lib/python2.7/site-packages/

python2.7 $XCAMSHIFT_HOME/src/table_builders/xcamshift/camshift_table_builder.py -t $XCAMSHIFT_HOME/data/xcamshift_templates/ $XCAMSHIFT_HOME/data/camshift_1_35_0_almost_1_0_4/  

