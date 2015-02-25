

from pdbTool import PDBTool
from protocol import initStruct
from test_pot import BondPot

initStruct("../test_data/3_ala/3ala.psf")
PDBTool("../test_data/3_ala/3ala.pdb").read()

pot = BondPot('test','(resid 2 and name CA)','(resid 3 and name CA)',1)
print int(pot.simulation().this)
