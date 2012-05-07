#
# mutate a single amino acid - to derivatized version.
# read in coordinates of the original protein and build in the positions of
# the new atoms.
#
#
xplor.requireVersion("2.12")
from pdbTool import PDBTool

xplor.parseArguments() # check for typos on the command-line

seq="""
ALA GLY ALA
"""

simWorld.setRandomSeed(5521)

import protocol
protocol.initTopology(('protein'))
protocol.initParams(('protein'))

import psfGen
psfGen.seqToPSF(seq,seqType='prot')


notKnown=AtomSel("not known")

protocol.addUnknownAtoms(dyn_stepsize=0.01)
protocol.fixupCovalentGeom(sel=notKnown)


xplor.command("write psf output=AGA.psf end")
PDBTool("AGA.pdb").write()
