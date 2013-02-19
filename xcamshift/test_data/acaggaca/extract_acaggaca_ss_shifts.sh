#!/bin/bash

cat <<LICENSE
#-------------------------------------------------------------------------------
# Copyright (c) 2013 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
LICENSE

echo 
echo acaggaca_shifts = {
echo
~/programs/camshift/1.35.0-instrumented/bin/camshift --pdb acaggaca.pdb --data ~/programs/camshift/1.35.0-vanilla/data/ | grep CSS_ | awk '$7 != 0 {printf "(@@,%i,@%s@)  %+4.8f\n",$3,$6,$7}' | sort -g -k 1 | awk '{printf "    %-15s  : %+8.4f,\n",$1,$2}' | tr '@+' "' " 
echo }
