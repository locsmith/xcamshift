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
echo acaggaca_ss_shifts = {
echo
~/programs/camshift/1.35.0-instrumented/bin/camshift --pdb acaggaca.pdb --data ~/programs/camshift/1.35.0-vanilla/data/ | grep CSS_ | awk '$7 != 0 {printf "(@@,%i,@%s@)  %+4.8f\n",$3,$6,$7}' | sort -g -k 1 | sed -e 's/@H@/@HN@/' | awk '{printf "    %-15s  : %+8.4f,\n",$1,$2}' |  tr '@+' "' " 
echo }

echo
echo acaggaca_shifts = {
echo
~/programs/camshift/1.35.0-instrumented/bin/camshift --pdb acaggaca.pdb --data ~/programs/camshift/1.35.0-vanilla/data/ | egrep 'call|Calculated chemical shift' | tr ':' ' ' | awk ' /Calculated chemical shift at residue/ {printf "(@@,%i,@%s@) %9.5f\n", $6,$9,$10}'  | sed 's/@H@/@HN@/g' | awk '{printf "    %-15s  : %9.5f,\n",$1,$2}' | tr "@" "'"
echo }

echo
echo  acaggaca_hbond_shifts = {
echo
~/programs/camshift/1.35.0-instrumented/bin/camshift --pdb acaggaca.pdb --data ~/programs/camshift/1.35.0-vanilla/data/ | grep CHB_TOTAL| awk '$7 != 0 {printf "(@@,%i,@%s@)  %+4.8f\n",$3,$6,$7}' | sort -g -k 1 | sed -e 's/@H@/@HN@/' | awk '{printf "    %-15s  : %+8.4f,\n",$1,$2}' |  tr '@+' "' " 
echo }

echo
echo \#
echo \# corrections added because camshift adds HG protons back in for 
echo \# disulphide bonded cystines
echo \#  
echo acaggaca_ss_nbond_corrections =  {
echo
~/programs/camshift/1.35.0-instrumented/bin/camshift --data /home/garyt/programs/camshift/1.35.0-vanilla/data --pdb acaggaca.pdb | grep '^NB_COMPONENT' | egrep  ':2@HG1|7@HG1' | tr ':@[]' '    ' | awk '{printf "((@@,%i,@%s@),(@@,%i,@%s@),%i) %14.10f\n",$3,$4,$7,$8,$9,$10}' |  awk '{printf "    %-30s  : %+6.10f,\n",$1,$2}' |  tr '@+' "' " 
echo "}"
echo

echo
echo \#
echo \# corrections added because camshift adds HG protons back in for 
echo \# disulphide bonded cystines which causes errors in the sidechain shift
echo \#  
echo acaggaca_ss_sidechain_corrections =  {
echo
~/programs/camshift/1.35.0-instrumented/bin/camshift --data /home/garyt/programs/camshift/1.35.0-vanilla/data --pdb acaggaca.pdb | egrep 'CSC_TOTAL residue: [2|7] atom name: CB' | grep  'CSC_TOTAL residue: [2|7] atom name: CB' | tr ':@[]' '    ' | awk '{printf "(@@,%i,@%s@) %14.10f\n",$3,$6,$7}' |  awk '{printf "    %-20s  : %+6.10f,\n",$1,$2}' |  tr '@+' "' " 
echo "}"
echo
