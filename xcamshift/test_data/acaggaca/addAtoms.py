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
ALA CYS ALA GLY GLY ALA CYS  ALA
"""

simWorld.setRandomSeed(5521)

import protocol
protocol.initTopology(('protein'))
protocol.initParams(('protein'))

import psfGen
psfGen.seqToPSF(seq,seqType='prot')
protocol.addDisulfideBond('resid 2','resid 7')


notKnown=AtomSel("not known")

protocol.addUnknownAtoms(dyn_stepsize=0.01)
protocol.fixupCovalentGeom(sel=notKnown)


xplor.command("write psf output=acaggaca.psf end")
PDBTool("acaggaca.pdb").write()
