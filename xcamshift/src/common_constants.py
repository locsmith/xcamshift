#-------------------------------------------------------------------------------
# Copyright (c) 2013-2015 Gary Thompson & The University of Leeds.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v3.0
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/lgpl-3.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 6 Apr 2012

@author: garyt
'''

BASE = ''
GLY = 'gly'
PRO = 'pro'
CAMSHIFT_RESIDUE_TYPES = (BASE, GLY, PRO)

SPHERE_1 = 'sphere_1'
SPHERE_2 = 'sphere_2'

SP3 = 'SP3'
SP2 = 'SP2'

N = 'N'
H = 'H'
CA = 'CA'
CB = 'CB'
CG = 'CG'
HA = 'HA'
C = 'C'
HN = 'HN'

h_keys = (HA,CA,H,N,C,CB)

O = 'O'

DATA = 'data'

BACK_BONE = "BB  "
XTRA = "EXTRA"
RANDOM_COIL = "RC  "
DIHEDRAL = "DIHEDRAL"
SIDE_CHAIN = "SIDECHAIN"       
RING = 'RING'
NON_BONDED = "nb"   
DISULPHIDE = 'DISULPHIDE'
HBOND = 'HBOND'

CAMSHIFT_SUB_POTENTIALS =  RANDOM_COIL, BACK_BONE, XTRA, SIDE_CHAIN, DIHEDRAL, RING, NON_BONDED, DISULPHIDE, HBOND    

AXIS_NAMES = ('X','Y','Z')

DATA_TABLES_CHANGED = "data tables"
STRUCTURE_CHANGED = "structure"
SHIFT_DATA_CHANGED = "shift data"
ROUND_CHANGED =  "round"
TARGET_ATOM_IDS_CHANGED = "target atoms"
