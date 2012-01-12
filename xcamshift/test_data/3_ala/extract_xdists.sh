#!/bin/sh

touch xdists.py
rm xdists.py
echo "xdists_ala_3 = {" >> xdists.py
~/programs/camshift/1.35.0/bin/camshift --pdb 3ala.pdb --data ~/programs/camshift/1.35.0/data/ | grep ^XDIST | grep -v incomplete | grep -v complete | tr ',()' '    ' | awk '{print $4,$6,$10,$12,$14,$16}' | tr ':@' '  ' | tr -s ' '  | awk '{printf "          (%i, \"%s\", %i, \"%s\", %i, \"%s\") : (% 4.7f,% 4.7f),\n",$1,$2,$3,$4,$5,$6,$7,$8}' | sed -e 's/\(\".\",\)/\1 /g' |  sed -e 's/\(\".\"\))/\1 )/g'>> xdists.py 
echo "}" >> xdists.py
