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
'''
Created on 31 Dec 2011

@author: garyt
'''
import sys
epsilon =  sys.float_info.epsilon
from protocol import initStruct
from pdbTool import PDBTool
import unittest2
from xcamshift import Hbond_backbone_donor_and_acceptor_indexer, Hbond_backbone_donor_indexer, Hbond_backbone_acceptor_indexer, Hydrogen_bond_context, Hbond_atom_type_indexer,\
     Hydrogen_bond_donor_context, Hydrogen_bond_acceptor_context, Hydrogen_bond_donor_component_factory, DONOR,ACCEPTOR,\
     Hbond_donor_atom_type_indexer, Hbond_acceptor_atom_type_indexer, Hydrogen_bond_acceptor_component_factory, Xcamshift, Hydrogen_bond_parameter_factory, \
     Hydrogen_bond_donor_lookup_factory, Hydrogen_bond_acceptor_lookup_factory,\
    Hydrogen_bond_potential, Hydrogen_bond_component_factory
from cython.shift_calculators import Fast_hydrogen_bond_calculator, allocate_array
from cython.fast_segment_manager import Segment_Manager
from utils import Atom_utils
from table_manager import Table_manager
from atomSel import AtomSel
from component_list import Native_component_list

BACKBONE=1
SIDE_CHAIN=0
EXPECTED_DONOR_ACCEPTORS_BASE =     ((ACCEPTOR, (7,  'O'),   (7,  'C'),   BACKBONE),
                                     (ACCEPTOR, (8,  'N'),   (7,  'C'),   BACKBONE),
                                     (DONOR,    (8,  'HN'),  (8,  'N'),   BACKBONE),
                                     (ACCEPTOR, (8,  'O'),   (8,  'C'),   BACKBONE),
                                     (ACCEPTOR, (9,  'N'),   (8,  'C'),   BACKBONE),
                                     (DONOR,    (9,  'HN'),  (9,  'N'),   BACKBONE),
                                     (ACCEPTOR, (9,  'O'),   (9,  'C'),   BACKBONE),
                                     (ACCEPTOR, (10, 'N'),   (9,  'C'),   BACKBONE),
                                     (DONOR,    (10, 'HN'),  (10, 'N'),   BACKBONE),
                                     (ACCEPTOR, (10, 'O'),   (10, 'C'),   BACKBONE),
                                     (DONOR,    (10, 'HZ1'), (10, 'NZ'),  SIDE_CHAIN),
                                     (DONOR,    (10, 'HZ2'), (10, 'NZ'),  SIDE_CHAIN),
                                     (DONOR,    (10, 'HZ3'), (10, 'NZ'),  SIDE_CHAIN),
                                     (ACCEPTOR, (10, 'NZ'),  (10, 'CE'),  SIDE_CHAIN),
                                     (DONOR,    (11, 'HG1'), (11, 'OG1'), SIDE_CHAIN),
                                     (ACCEPTOR, (11, 'N'),   (10, 'C'),   BACKBONE),
                                     (DONOR,    (11, 'HN'),  (11, 'N'),   BACKBONE),
                                     (ACCEPTOR, (11, 'O'),   (11, 'C'),   BACKBONE),
                                     (ACCEPTOR, (11, 'OG1'), (11, 'CB'),  SIDE_CHAIN),
                                     (ACCEPTOR, (12, 'N'),   (11, 'C'),   BACKBONE),
                                     (DONOR,    (12, 'HN'),  (12, 'N'),   BACKBONE),
                                     (ACCEPTOR, (12, 'O'),   (12, 'C'),   BACKBONE),
                                     (DONOR,    (13, 'HZ1'), (13, 'NZ'),  SIDE_CHAIN),
                                     (DONOR,    (13, 'HZ2'), (13, 'NZ'),  SIDE_CHAIN),
                                     (DONOR,    (13, 'HZ3'), (13, 'NZ'),  SIDE_CHAIN),
                                     (ACCEPTOR, (13, 'N'),   (12, 'C'),   BACKBONE),
                                     (DONOR,    (13, 'HN'),  (13, 'N'),   BACKBONE),
                                     (ACCEPTOR, (13, 'O'),   (13, 'C'),   BACKBONE),
                                     (ACCEPTOR, (13, 'NZ'),  (13, 'CE'),  SIDE_CHAIN),
                                     (ACCEPTOR, (14, 'N'),   (13, 'C'),   BACKBONE),
                                     (DONOR,    (14, 'HN'),  (14, 'N'),   BACKBONE))
   
EXPECTED_DONORS = [(elem[1][0],elem[1][1],elem[2][1],elem[3]) for elem in  EXPECTED_DONOR_ACCEPTORS_BASE if elem[0] == DONOR]

DIST=0
ANG_1=1
ANG_2=2

EXPECTED_DONOR_ACCEPTOR_ENERGIES = {
        ('',  7,  'O', 1, DIST):     3.27874313127,
        ('',  7,  'O', 1, ANG_1):  696.529729434,
        ('',  7,  'O', 1, ANG_2):  6683.4703906,
        ('',  9, 'HN', 0, DIST):     3.64372,
        ('',  9, 'HN', 0, ANG_1):  701.205,
        ('',  9, 'HN', 0, ANG_2): 6548.52,
        ('', 10,  'N', 1, DIST):     4.45286,
        ('', 10,  'N', 1, ANG_1):  672.789,
        ('', 10,  'N', 1, ANG_2): 6252.28,
        ('', 11, 'HN', 0, DIST):     4.45286,
        ('', 11, 'HN', 0, ANG_1):  672.789,
        ('', 11, 'HN', 0, ANG_2): 6252.28,
        ('', 11,  'N', 1, DIST):     4.66541,
        ('', 11,  'N', 1, ANG_1):  695.025,
        ('', 11,  'N', 1, ANG_2): 6238.2,
        ('', 12, 'HN', 0, DIST):     4.66541,
        ('', 12, 'HN', 0, ANG_1):  695.025,
        ('', 12, 'HN', 0, ANG_2): 6238.2,
        ('', 12,  'O', 1, DIST):     3.64372,
        ('', 12,  'O', 1, ANG_1):  701.205,
        ('', 12,  'O', 1, ANG_2): 6548.52,
        ('', 14, 'HN', 0, DIST):     3.27874313127,
        ('', 14, 'HN', 0, ANG_1):  696.529729434,
        ('', 14, 'HN', 0, ANG_2):  6683.4703906
}

EXPECTED_HBOND_SHIFT_DATA = (
         (ACCEPTOR, (7,  'O'),   (7,  'C'),   BACKBONE),
         (ACCEPTOR, (8,  'N'),   (7,  'C'),   BACKBONE),
         (DONOR,    (8,  'HN'),  (8,  'N'),   BACKBONE),
         (ACCEPTOR, (8,  'O'),   (8,  'C'),   BACKBONE),
         (ACCEPTOR, (9,  'N'),   (8,  'C'),   BACKBONE),
         (DONOR,    (9,  'HN'),  (9,  'N'),   BACKBONE),
         (ACCEPTOR, (9,  'O'),   (9,  'C'),   BACKBONE),
         (ACCEPTOR, (10, 'N'),   (9,  'C'),   BACKBONE),
         (DONOR,    (10, 'HN'),  (10, 'N'),   BACKBONE),
         (ACCEPTOR, (10, 'O'),   (10, 'C'),   BACKBONE),
         (ACCEPTOR, (11, 'N'),   (10, 'C'),   BACKBONE),
         (DONOR,    (11, 'HN'),  (11, 'N'),   BACKBONE),
         (ACCEPTOR, (11, 'O'),   (11, 'C'),   BACKBONE),
         (ACCEPTOR, (12, 'N'),   (11, 'C'),   BACKBONE),
         (DONOR,    (12, 'HN'),  (12, 'N'),   BACKBONE),
         (ACCEPTOR, (12, 'O'),   (12, 'C'),   BACKBONE),
         (ACCEPTOR, (13, 'N'),   (12, 'C'),   BACKBONE),
         (DONOR,    (13, 'HN'),  (13, 'N'),   BACKBONE),
         (ACCEPTOR, (13, 'O'),   (13, 'C'),   BACKBONE),
         (ACCEPTOR, (14, 'N'),   (13, 'C'),   BACKBONE),
         (DONOR,    (14, 'HN'),  (14, 'N'),   BACKBONE))




EXPECTED_SHIFT_OFFSET_DATA = {
            (8, 'N', 7, 'O', -1)      : (  -0.0000752700,   -0.0180410700,    0.0018754600),
            (8, 'N', 8, 'HN', 0)      : (  -0.0134281400,    0.0003982300,   -0.0001003400),
            (8, 'N', 8, 'O', 0)       : (  -0.0000844000,    0.0065627300,   -0.0006374100),
            (8, 'N', 9, 'HN', 1)      : (  -0.0155729300,    0.0117727900,   -0.0012193100),
            (8, 'HN', 7, 'O', -1)     : (  -0.0000001000,   -0.0012311600,    0.0001250700),
            (8, 'HN', 8, 'HN', 0)     : (  -0.0142871900,    0.0044821200,   -0.0005105900),
            (8, 'HN', 8, 'O', 0)      : (  -0.0000054700,    0.0002329400,   -0.0000156700),
            (8, 'HN', 9, 'HN', 1)     : (   0.0011321300,    0.0001720200,   -0.0000099700),
            (8, 'CA', 7, 'O', -1)     : (  -0.0000051500,    0.0014532100,   -0.0001394000),
            (8, 'CA', 8, 'HN', 0)     : (   0.0030282700,    0.0010660900,   -0.0001326300),
            (8, 'CA', 8, 'O', 0)      : (  -0.0000256900,   -0.0002918000,    0.0000327700),
            (8, 'CA', 9, 'HN', 1)     : (  -0.0052441200,    0.0034480700,   -0.0003803300),
            (8, 'HA', 7, 'O', -1)     : (  -0.0000000300,   -0.0005796100,    0.0000686600),
            (8, 'HA', 8, 'HN', 0)     : (  -0.0013889800,    0.0001070500,   -0.0000005600),
            (8, 'HA', 8, 'O', 0)      : (   0.0000011400,    0.0003084000,   -0.0000280800),
            (8, 'HA', 9, 'HN', 1)     : (  -0.0017768900,    0.0008554600,   -0.0000814700),
            (8, 'CB', 7, 'O', -1)     : (   0.0000330200,   -0.0016325700,    0.0001679500),
            (8, 'CB', 8, 'HN', 0)     : (  -0.0066043000,    0.0007734100,   -0.0000318000),
            (8, 'CB', 8, 'O', 0)      : (  -0.0000179100,    0.0011680200,   -0.0001242800),
            (8, 'CB', 9, 'HN', 1)     : (   0.0017154800,   -0.0008403200,    0.0001072600),
            (8, 'C', 7, 'O', -1)      : (   0.0000084800,    0.0039887400,   -0.0003917000),
            (8, 'C', 8, 'HN', 0)      : (  -0.0043335900,   -0.0004836800,    0.0000546400),
            (8, 'C', 8, 'O', 0)       : (  -0.0000169300,   -0.0113520400,    0.0011652700),
            (8, 'C', 9, 'HN', 1)      : (  -0.0077359500,    0.0057129800,   -0.0005442800),
            (9, 'N', 8, 'O', -1)      : (  -0.0000841600,   -0.0165728000,    0.0017066600),
            (9, 'N', 9, 'HN', 0)      : (  -0.0126117400,    0.0018564600,   -0.0002635900),
            (9, 'N', 9, 'O', 0)       : (  -0.0000794900,    0.0055338900,   -0.0005243700),
            (9, 'N', 10, 'HN', 1)     : (  -0.0132877500,    0.0128629400,   -0.0013345400),
            (9, 'HN', 8, 'O', -1)     : (  -0.0000007300,   -0.0012923600,    0.0001313100),
            (9, 'HN', 9, 'HN', 0)     : (  -0.0148477300,    0.0051241700,   -0.0005738600),
            (9, 'HN', 9, 'O', 0)      : (  -0.0000060100,    0.0000222400,    0.0000062900),
            (9, 'HN', 10, 'HN', 1)    : (   0.0009498300,    0.0000924900,   -0.0000008900),
            (9, 'CA', 8, 'O', -1)     : (  -0.0000070200,    0.0006044400,   -0.0000456600),
            (9, 'CA', 9, 'HN', 0)     : (   0.0030458700,    0.0013572500,   -0.0001635600),
            (9, 'CA', 9, 'O', 0)      : (  -0.0000247200,   -0.0002636900,    0.0000302800),
            (9, 'CA', 10, 'HN', 1)    : (  -0.0050053100,    0.0030799600,   -0.0003454400),
            (9, 'HA', 8, 'O', -1)     : (  -0.0000002900,   -0.0004180100,    0.0000524600),
            (9, 'HA', 9, 'HN', 0)     : (  -0.0012570400,    0.0001397800,   -0.0000042600),
            (9, 'HA', 9, 'O', 0)      : (   0.0000014700,    0.0001921100,   -0.0000161200),
            (9, 'HA', 10, 'HN', 1)    : (  -0.0018399300,    0.0009026500,   -0.0000829600),
            (9, 'C', 8, 'O', -1)      : (   0.0000098900,    0.0046344400,   -0.0004593400),
            (9, 'C', 9, 'HN', 0)      : (  -0.0049752900,   -0.0006142000,    0.0000670400),
            (9, 'C', 9, 'O', 0)       : (  -0.0000165000,   -0.0107069100,    0.0010977700),
            (9, 'C', 10, 'HN', 1)     : (  -0.0099267500,    0.0060277800,   -0.0005737400),
            (10, 'N', 9, 'O', -1)     : (  -0.0000752700,   -0.0180410700,    0.0018754600),
            (10, 'N', 10, 'HN', 0)    : (  -0.0134281400,    0.0003982300,   -0.0001003400),
            (10, 'N', 10, 'O', 0)     : (  -0.0000844000,    0.0065627300,   -0.0006374100),
            (10, 'N', 11, 'HN', 1)    : (  -0.0155729300,    0.0117727900,   -0.0012193100),
            (10, 'HN', 9, 'O', -1)    : (  -0.0000001000,   -0.0012311600,    0.0001250700),
            (10, 'HN', 10, 'HN', 0)   : (  -0.0142871900,    0.0044821200,   -0.0005105900),
            (10, 'HN', 10, 'O', 0)    : (  -0.0000054700,    0.0002329400,   -0.0000156700),
            (10, 'HN', 11, 'HN', 1)   : (   0.0011321300,    0.0001720200,   -0.0000099700),
            (10, 'CA', 9, 'O', -1)    : (  -0.0000051500,    0.0014532100,   -0.0001394000),
            (10, 'CA', 10, 'HN', 0)   : (   0.0030282700,    0.0010660900,   -0.0001326300),
            (10, 'CA', 10, 'O', 0)    : (  -0.0000256900,   -0.0002918000,    0.0000327700),
            (10, 'CA', 11, 'HN', 1)   : (  -0.0052441200,    0.0034480700,   -0.0003803300),
            (10, 'HA', 9, 'O', -1)    : (  -0.0000000300,   -0.0005796100,    0.0000686600),
            (10, 'HA', 10, 'HN', 0)   : (  -0.0013889800,    0.0001070500,   -0.0000005600),
            (10, 'HA', 10, 'O', 0)    : (   0.0000011400,    0.0003084000,   -0.0000280800),
            (10, 'HA', 11, 'HN', 1)   : (  -0.0017768900,    0.0008554600,   -0.0000814700),
            (10, 'CB', 9, 'O', -1)    : (   0.0000330200,   -0.0016325700,    0.0001679500),
            (10, 'CB', 10, 'HN', 0)   : (  -0.0066043000,    0.0007734100,   -0.0000318000),
            (10, 'CB', 10, 'O', 0)    : (  -0.0000179100,    0.0011680200,   -0.0001242800),
            (10, 'CB', 11, 'HN', 1)   : (   0.0017154800,   -0.0008403200,    0.0001072600),
            (10, 'C', 9, 'O', -1)     : (   0.0000084800,    0.0039887400,   -0.0003917000),
            (10, 'C', 10, 'HN', 0)    : (  -0.0043335900,   -0.0004836800,    0.0000546400),
            (10, 'C', 10, 'O', 0)     : (  -0.0000169300,   -0.0113520400,    0.0011652700),
            (10, 'C', 11, 'HN', 1)    : (  -0.0077359500,    0.0057129800,   -0.0005442800),
            (11, 'N', 10, 'O', -1)    : (  -0.0000752700,   -0.0180410700,    0.0018754600),
            (11, 'N', 11, 'HN', 0)    : (  -0.0134281400,    0.0003982300,   -0.0001003400),
            (11, 'N', 11, 'O', 0)     : (  -0.0000844000,    0.0065627300,   -0.0006374100),
            (11, 'N', 12, 'HN', 1)    : (  -0.0155729300,    0.0117727900,   -0.0012193100),
            (11, 'HN', 10, 'O', -1)   : (  -0.0000001000,   -0.0012311600,    0.0001250700),
            (11, 'HN', 11, 'HN', 0)   : (  -0.0142871900,    0.0044821200,   -0.0005105900),
            (11, 'HN', 11, 'O', 0)    : (  -0.0000054700,    0.0002329400,   -0.0000156700),
            (11, 'HN', 12, 'HN', 1)   : (   0.0011321300,    0.0001720200,   -0.0000099700),
            (11, 'CA', 10, 'O', -1)   : (  -0.0000051500,    0.0014532100,   -0.0001394000),
            (11, 'CA', 11, 'HN', 0)   : (   0.0030282700,    0.0010660900,   -0.0001326300),
            (11, 'CA', 11, 'O', 0)    : (  -0.0000256900,   -0.0002918000,    0.0000327700),
            (11, 'CA', 12, 'HN', 1)   : (  -0.0052441200,    0.0034480700,   -0.0003803300),
            (11, 'HA', 10, 'O', -1)   : (  -0.0000000300,   -0.0005796100,    0.0000686600),
            (11, 'HA', 11, 'HN', 0)   : (  -0.0013889800,    0.0001070500,   -0.0000005600),
            (11, 'HA', 11, 'O', 0)    : (   0.0000011400,    0.0003084000,   -0.0000280800),
            (11, 'HA', 12, 'HN', 1)   : (  -0.0017768900,    0.0008554600,   -0.0000814700),
            (11, 'CB', 10, 'O', -1)   : (   0.0000330200,   -0.0016325700,    0.0001679500),
            (11, 'CB', 11, 'HN', 0)   : (  -0.0066043000,    0.0007734100,   -0.0000318000),
            (11, 'CB', 11, 'O', 0)    : (  -0.0000179100,    0.0011680200,   -0.0001242800),
            (11, 'CB', 12, 'HN', 1)   : (   0.0017154800,   -0.0008403200,    0.0001072600),
            (11, 'C', 10, 'O', -1)    : (   0.0000084800,    0.0039887400,   -0.0003917000),
            (11, 'C', 11, 'HN', 0)    : (  -0.0043335900,   -0.0004836800,    0.0000546400),
            (11, 'C', 11, 'O', 0)     : (  -0.0000169300,   -0.0113520400,    0.0011652700),
            (11, 'C', 12, 'HN', 1)    : (  -0.0077359500,    0.0057129800,   -0.0005442800),
            (12, 'N', 11, 'O', -1)    : (  -0.0000752700,   -0.0180410700,    0.0018754600),
            (12, 'N', 12, 'HN', 0)    : (  -0.0134281400,    0.0003982300,   -0.0001003400),
            (12, 'N', 12, 'O', 0)     : (  -0.0000844000,    0.0065627300,   -0.0006374100),
            (12, 'N', 13, 'HN', 1)    : (  -0.0155729300,    0.0117727900,   -0.0012193100),
            (12, 'HN', 11, 'O', -1)   : (  -0.0000001000,   -0.0012311600,    0.0001250700),
            (12, 'HN', 12, 'HN', 0)   : (  -0.0142871900,    0.0044821200,   -0.0005105900),
            (12, 'HN', 12, 'O', 0)    : (  -0.0000054700,    0.0002329400,   -0.0000156700),
            (12, 'HN', 13, 'HN', 1)   : (   0.0011321300,    0.0001720200,   -0.0000099700),
            (12, 'CA', 11, 'O', -1)   : (  -0.0000051500,    0.0014532100,   -0.0001394000),
            (12, 'CA', 12, 'HN', 0)   : (   0.0030282700,    0.0010660900,   -0.0001326300),
            (12, 'CA', 12, 'O', 0)    : (  -0.0000256900,   -0.0002918000,    0.0000327700),
            (12, 'CA', 13, 'HN', 1)   : (  -0.0052441200,    0.0034480700,   -0.0003803300),
            (12, 'HA', 11, 'O', -1)   : (  -0.0000000300,   -0.0005796100,    0.0000686600),
            (12, 'HA', 12, 'HN', 0)   : (  -0.0013889800,    0.0001070500,   -0.0000005600),
            (12, 'HA', 12, 'O', 0)    : (   0.0000011400,    0.0003084000,   -0.0000280800),
            (12, 'HA', 13, 'HN', 1)   : (  -0.0017768900,    0.0008554600,   -0.0000814700),
            (12, 'CB', 11, 'O', -1)   : (   0.0000330200,   -0.0016325700,    0.0001679500),
            (12, 'CB', 12, 'HN', 0)   : (  -0.0066043000,    0.0007734100,   -0.0000318000),
            (12, 'CB', 12, 'O', 0)    : (  -0.0000179100,    0.0011680200,   -0.0001242800),
            (12, 'CB', 13, 'HN', 1)   : (   0.0017154800,   -0.0008403200,    0.0001072600),
            (12, 'C', 11, 'O', -1)    : (   0.0000084800,    0.0039887400,   -0.0003917000),
            (12, 'C', 12, 'HN', 0)    : (  -0.0043335900,   -0.0004836800,    0.0000546400),
            (12, 'C', 12, 'O', 0)     : (  -0.0000169300,   -0.0113520400,    0.0011652700),
            (12, 'C', 13, 'HN', 1)    : (  -0.0077359500,    0.0057129800,   -0.0005442800),
            (13, 'N', 12, 'O', -1)    : (  -0.0000752700,   -0.0180410700,    0.0018754600),
            (13, 'N', 13, 'HN', 0)    : (  -0.0134281400,    0.0003982300,   -0.0001003400),
            (13, 'N', 13, 'O', 0)     : (  -0.0000844000,    0.0065627300,   -0.0006374100),
            (13, 'N', 14, 'HN', 1)    : (  -0.0155729300,    0.0117727900,   -0.0012193100),
            (13, 'HN', 12, 'O', -1)   : (  -0.0000001000,   -0.0012311600,    0.0001250700),
            (13, 'HN', 13, 'HN', 0)   : (  -0.0142871900,    0.0044821200,   -0.0005105900),
            (13, 'HN', 13, 'O', 0)    : (  -0.0000054700,    0.0002329400,   -0.0000156700),
            (13, 'HN', 14, 'HN', 1)   : (   0.0011321300,    0.0001720200,   -0.0000099700),
            (13, 'CA', 12, 'O', -1)   : (  -0.0000051500,    0.0014532100,   -0.0001394000),
            (13, 'CA', 13, 'HN', 0)   : (   0.0030282700,    0.0010660900,   -0.0001326300),
            (13, 'CA', 13, 'O', 0)    : (  -0.0000256900,   -0.0002918000,    0.0000327700),
            (13, 'CA', 14, 'HN', 1)   : (  -0.0052441200,    0.0034480700,   -0.0003803300),
            (13, 'HA', 12, 'O', -1)   : (  -0.0000000300,   -0.0005796100,    0.0000686600),
            (13, 'HA', 13, 'HN', 0)   : (  -0.0013889800,    0.0001070500,   -0.0000005600),
            (13, 'HA', 13, 'O', 0)    : (   0.0000011400,    0.0003084000,   -0.0000280800),
            (13, 'HA', 14, 'HN', 1)   : (  -0.0017768900,    0.0008554600,   -0.0000814700),
            (13, 'CB', 12, 'O', -1)   : (   0.0000330200,   -0.0016325700,    0.0001679500),
            (13, 'CB', 13, 'HN', 0)   : (  -0.0066043000,    0.0007734100,   -0.0000318000),
            (13, 'CB', 13, 'O', 0)    : (  -0.0000179100,    0.0011680200,   -0.0001242800),
            (13, 'CB', 14, 'HN', 1)   : (   0.0017154800,   -0.0008403200,    0.0001072600),
            (13, 'C', 12, 'O', -1)    : (   0.0000084800,    0.0039887400,   -0.0003917000),
            (13, 'C', 13, 'HN', 0)    : (  -0.0043335900,   -0.0004836800,    0.0000546400),
            (13, 'C', 13, 'O', 0)     : (  -0.0000169300,   -0.0113520400,    0.0011652700),
            (13, 'C', 14, 'HN', 1)    : (  -0.0077359500,    0.0057129800,   -0.0005442800),
}

EXPECTED_INDIRECT_DONORS = {}
for elem in EXPECTED_DONORS:
    EXPECTED_INDIRECT_DONORS['',elem[0],elem[1]] = '',elem[0],elem[2]
    
EXPECTED_DIRECT_DONORS = [('',elem[0],elem[1]) for elem in EXPECTED_DONORS]

EXPECTED_DONOR_TYPES = sorted(['HON',])
EXPECTED_ACCEPTOR_TYPES = sorted(['ON',])

EXPECTED_INDIRECT_ACCEPTORS = {}
for elem in EXPECTED_DONOR_ACCEPTORS_BASE:
    if elem[0] ==  ACCEPTOR:
        EXPECTED_INDIRECT_ACCEPTORS['',elem[1][0],elem[1][1]] = '',elem[2][0],elem[2][1]

    
EXPECTED_DIRECT_ACCEPTORS = [('',elem[1][0],elem[1][1]) for elem in EXPECTED_DONOR_ACCEPTORS_BASE if elem[0] ==  ACCEPTOR]


EXPECTED_BACK_BONE_DONORS = [('',elem[0],elem[1]) for elem in EXPECTED_DONORS if elem[-1] == 1]
EXPECTED_BACK_BONE_ACCEPTORS = [('',elem[1][0],elem[1][1]) for elem in EXPECTED_DONOR_ACCEPTORS_BASE if elem[-1] == 1 and elem[0] == ACCEPTOR]

EXPECTED_BACKBONE_DONORS_AND_ACCEPTORS = [('',elem[1][0],elem[1][1]) for elem in EXPECTED_DONOR_ACCEPTORS_BASE if elem[-1] == 1]

class TestXcamshiftHBondINGKTLKG(unittest2.TestCase):

    def __init__(self,*args,**kwargs):
        super(TestXcamshiftHBondINGKTLKG, self).__init__(*args,**kwargs)
        self._esim = None
        
        
    def assertSequenceContains(self,expected,sequence):
        if expected not in sequence:
            sequence_strings = [`elem` for elem in sequence]
            sequence_string = '\n'.join(sequence_strings)
            raise AssertionError("element %s not found in the sequence:\n%s" % (`expected`,sequence_string))

        
    def assertSequenceDoesntContain(self,expected,sequence):
        if expected in sequence:
            sequence_strings = [`elem` for elem in sequence]
            sequence_string = '\n'.join(sequence_strings)
            raise AssertionError("element %s found in the sequence:\n%s" % (`expected`,sequence_string))        
        
    def get_single_member_ensemble_simulation(self):
        if self._esim.__class__ ==  None.__class__:
            #TODO note EnsembleSimulation can't have a single member that causes a crash!
            # therefore a hack
            self._esim =  Xcamshift().ensembleSimulation()
        return self._esim
    
    def assertEmpty(self, expected_keys, msg=""):
        return self.assertEqual(len(expected_keys), 0, msg)
    
    def assertLength(self,sequence,length,msg='bad length expected length of %i'):
        msg = msg % length
        self.assertEqual(len(sequence), length, msg)
        
    
    def check_almost_equal(self, list_1, list_2, delta = 1e-7):
        difference_offset = -1
        for i, (elem_1, elem_2) in enumerate(zip(list_1, list_2)):
            diff = abs(elem_1 - elem_2)
            if diff > delta:
                difference_offset = i
        
        return difference_offset
    
    def are_almost_equal_sequences(self, list_1, list_2, delta =  1e-7):
        result = True
        if self.check_almost_equal(list_1, list_2, delta) > 0:
            result = False
        return result
        
    def assertSequenceAlmostEqual(self,result,expected, delta = 1e-7, msg=""):
        len_result = len(result)
        len_expected = len(expected)
        if len_result != len_expected:
            raise AssertionError("the two lists are of different length %i and %i" % (len_result,len_expected))
        
        difference_offset = self.check_almost_equal(result, expected, delta)
        
            
        if difference_offset > 0:
            if msg != "":
                msg = msg + " "
                
            template = "%slists differ at item %i: %s - %s > %s"
            elem_1 = result[difference_offset]
            elem_2 = expected[difference_offset]
            message = template % (msg,difference_offset, `elem_1`,`elem_2`,delta)
            raise AssertionError(message)            
             
    def setUp(self):
        initStruct("test_data/ingktlkg_hbond/INGKTLKG.psf")
        PDBTool("test_data/ingktlkg_hbond/INGKTLKG.pdb").read()
        Atom_utils.clear_cache()

        table_manager =  Table_manager.get_default_table_manager()
        self.donor_and_acceptor_indexer = Hbond_backbone_donor_and_acceptor_indexer(table_manager)
        self.donor_indexer  = Hbond_backbone_donor_indexer(table_manager)
        self.acceptor_indexer  = Hbond_backbone_acceptor_indexer(table_manager)
        self.donor_atom_type_indexer =  Hbond_donor_atom_type_indexer(table_manager)
        self.acceptor_atom_type_indexer =  Hbond_acceptor_atom_type_indexer(table_manager)
        
        Segment_Manager.reset_segment_manager()
#         print "In method", self._testMethodName

    def test_backbone_donor_and_acceptor_indexers(self):
         
 
        donors = [donor for donor in self.donor_indexer.iter_keys()]
        self.assertSequenceEqual(sorted(donors), sorted(EXPECTED_BACK_BONE_DONORS)) 
         
        acceptors = [acceptor for acceptor in self.acceptor_indexer.iter_keys()]
        self.assertSequenceEqual(sorted(acceptors), sorted(EXPECTED_BACK_BONE_ACCEPTORS))
        
        donors_and_acceptors = [donor_or_acceptor for donor_or_acceptor in self.donor_and_acceptor_indexer.iter_keys()]
        self.assertEqual(sorted(EXPECTED_BACKBONE_DONORS_AND_ACCEPTORS), sorted(donors_and_acceptors))

    def test_backbone_donor_and_acceptor_indexers_get_max_index(self):
        self.assertEqual(self.donor_indexer.get_max_index(), len(EXPECTED_BACK_BONE_DONORS))
        self.assertEqual(self.acceptor_indexer.get_max_index(), len(EXPECTED_BACK_BONE_ACCEPTORS))
        self.assertEqual(self.donor_and_acceptor_indexer.get_max_index(), len(EXPECTED_BACKBONE_DONORS_AND_ACCEPTORS), )
 
    def test_backbone_donor_and_acceptor_indexers_get_name(self):    
        self.assertTrue('donor' in self.donor_indexer.get_name().lower())
        self.assertTrue('acceptor' in self.acceptor_indexer.get_name().lower())
        
        self.assertTrue('acceptor' in self.donor_and_acceptor_indexer.get_name().lower())
        self.assertTrue('donor' in self.donor_and_acceptor_indexer.get_name().lower())
        
    def test_backbone_donor_and_acceptor_indexers_get_index_for_key(self,):

        for i,acceptor in enumerate(EXPECTED_BACK_BONE_ACCEPTORS):
            self.assertEqual(i,self.acceptor_indexer.get_index_for_key(acceptor))
        self.assertEqual(i+1, self.acceptor_indexer.get_max_index())
 
        for j,donor in enumerate(EXPECTED_BACK_BONE_DONORS):
            self.assertEqual(j,self.donor_indexer.get_index_for_key(donor))
        self.assertEqual(j+1, self.donor_indexer.get_max_index())
        
        for k,donor_acceptor in enumerate(EXPECTED_BACKBONE_DONORS_AND_ACCEPTORS):
            self.assertEqual(k,self.donor_and_acceptor_indexer.get_index_for_key(donor_acceptor))
        self.assertEqual(k+1, self.donor_and_acceptor_indexer.get_max_index())
                        
      
    def test_backbone_donor_and_acceptor_indexers_get_index_get_key_for_index(self):
        for i,acceptor in enumerate(EXPECTED_BACK_BONE_ACCEPTORS):
            self.assertEqual(acceptor,self.acceptor_indexer.get_key_for_index(i))
        self.assertEqual(i+1, self.acceptor_indexer.get_max_index())
 
        for j,donor in enumerate(EXPECTED_BACK_BONE_DONORS):
            self.assertEqual(donor,self.donor_indexer.get_key_for_index(j))
        self.assertEqual(j+1, self.donor_indexer.get_max_index())
    
        for k,donor_acceptor in enumerate(EXPECTED_BACKBONE_DONORS_AND_ACCEPTORS):
            self.assertEqual(donor_acceptor,self.donor_and_acceptor_indexer.get_key_for_index(k))
        self.assertEqual(k+1, self.donor_and_acceptor_indexer.get_max_index())
         
        
    def test_atom_indexer_indexers(self):
         
 
        donor_atom_indices = [index for index in self.donor_atom_type_indexer.iter_keys()]
        self.assertEqual(len(EXPECTED_DONOR_TYPES), len(donor_atom_indices))
        
        self.assertEqual(sorted(EXPECTED_DONOR_TYPES), donor_atom_indices)

        acceptor_atom_indices = [index for index in self.acceptor_atom_type_indexer.iter_keys()]
        self.assertEqual(len(EXPECTED_ACCEPTOR_TYPES), len(acceptor_atom_indices))
        
        self.assertEqual(sorted(EXPECTED_ACCEPTOR_TYPES), acceptor_atom_indices)        
      
    def test_get_max_index(self):
        self.assertEqual(self.donor_atom_type_indexer.get_max_index(), len(EXPECTED_DONOR_TYPES))
        self.assertEqual(self.acceptor_atom_type_indexer.get_max_index(), len(EXPECTED_ACCEPTOR_TYPES))
        
    def test_get_name(self):    
        for elem in 'hydrogen', 'bond', 'atom','type':
            self.assertTrue(elem in self.donor_atom_type_indexer.get_name().lower(), elem)
            self.assertTrue(elem in self.acceptor_atom_type_indexer.get_name().lower(), elem)
        
 
 
    def test_get_index_for_key(self,):
        for i,atom_name in enumerate(sorted(EXPECTED_DONOR_TYPES)):
            self.assertEqual(i,self.donor_atom_type_indexer.get_index_for_key(atom_name))
        self.assertEqual(i+1, self.donor_atom_type_indexer.get_max_index())
        
        for i,atom_name in enumerate(sorted(EXPECTED_ACCEPTOR_TYPES)):
            self.assertEqual(i,self.acceptor_atom_type_indexer.get_index_for_key(atom_name))
        self.assertEqual(i+1, self.acceptor_atom_type_indexer.get_max_index()) 
        
                        
      
    def test_get_key_for_index(self):
        for i,atom_name in enumerate(EXPECTED_DONOR_TYPES):
            self.assertEqual(atom_name,self.donor_atom_type_indexer.get_key_for_index(i))
        self.assertEqual(i+1, self.donor_atom_type_indexer.get_max_index())

        for i,atom_name in enumerate(EXPECTED_ACCEPTOR_TYPES):
            self.assertEqual(atom_name,self.acceptor_atom_type_indexer.get_key_for_index(i))
        self.assertEqual(i+1, self.acceptor_atom_type_indexer.get_max_index())
     
          
    def test_hbond_context(self):
        atom = Atom_utils.find_atom('', 10, 'HN')[0]
        offset_data_0 = (0, 'O')
        table  =  Table_manager.get_default_table_manager().get_hydrogen_bond_table('LYS')
        indexer  = Hbond_backbone_donor_and_acceptor_indexer(Table_manager.get_default_table_manager())
        
        hbond_context_0 = Hydrogen_bond_context(atom,offset_data_0,table, indexer)
        
        self.assertTrue(hbond_context_0.complete)        
        self.assertEqual(Atom_utils.find_atom_ids('   ',10,'HN')[0], hbond_context_0.target_atom_index)
        self.assertEqual(hbond_context_0.hbond_index, EXPECTED_BACKBONE_DONORS_AND_ACCEPTORS.index(('',10,'O')))
        
        EXPECTED_COEFFS_0 = [-0.00000547, 0.00023294, -0.00001567]
        self.assertAlmostEqual(hbond_context_0.coeffs, EXPECTED_COEFFS_0)
        self.assertSequenceAlmostEqual(EXPECTED_COEFFS_0, hbond_context_0.coeffs)
        
        offset_data_1 = (-1, 'O')
        hbond_context_1 = Hydrogen_bond_context(atom,offset_data_1,table,indexer)

        EXPECTED_COEFFS_1 = [-0.00000010, -0.00123116, 0.00012507]
        self.assertTrue(hbond_context_1.complete)
        self.assertEqual(Atom_utils.find_atom_ids('   ',10,'HN')[0], hbond_context_1.target_atom_index)
        self.assertSequenceAlmostEqual(EXPECTED_COEFFS_1, hbond_context_1.coeffs)
        self.assertEqual(hbond_context_1.hbond_index, EXPECTED_BACKBONE_DONORS_AND_ACCEPTORS.index(('',9,'O')))

        offset_data_2 = (-3, 'HN')
        hbond_context_2 = Hydrogen_bond_context(atom,offset_data_2,table,indexer)

        self.assertFalse(hbond_context_2.complete)
        
    
    def test_hydrogen_bond_donor_context(self):

        atom_0 = Atom_utils.find_atom('', 10, 'HN')[0]
        table  =  Table_manager.get_default_table_manager().get_hydrogen_bond_table('LYS')
        
        hbond_donor_context_0 = Hydrogen_bond_donor_context(atom_0,table)
        self.assertTrue(hbond_donor_context_0.complete)
        self.assertEqual(hbond_donor_context_0.atom_type_id,0)
        self.assertEqual(hbond_donor_context_0.direct_atom_id, atom_0.index())
        self.assertEqual(hbond_donor_context_0.indirect_atom_id, Atom_utils.find_atom('', 10, 'N')[0].index())

    def test_hydrogen_bond_acceptor_context(self):
        
        atom_0 = Atom_utils.find_atom('', 10, 'O')[0]
        table  =  Table_manager.get_default_table_manager().get_hydrogen_bond_table('LYS')
        
        hbond_acceptor_context_1 = Hydrogen_bond_acceptor_context(atom_0,table)
        self.assertTrue(hbond_acceptor_context_1.complete)
        self.assertEqual(hbond_acceptor_context_1.atom_type_id,0)
        self.assertEqual(hbond_acceptor_context_1.direct_atom_id, atom_0.index())
        self.assertEqual(hbond_acceptor_context_1.indirect_atom_id, Atom_utils.find_atom('', 10, 'C')[0].index())
        
    


    def _build_component_list(self, factory, format):
        component_list = Native_component_list(format)
        table_provider = Table_manager.get_default_table_manager().get_hydrogen_bond_table
        segment = '    '
        #TODO note an oddity here as this takes a residu butr builds a list for all residues...
        target_residue_number = 10
        selected_atoms = AtomSel('(all)')
        factory.create_components(component_list, table_provider, segment, target_residue_number, selected_atoms)
        return component_list

    def _do_test_donor_acceptor_components(self, factory,expected_direct_donors_or_acceptors, expected_indirect_donors_or_acceptors, donor_or_acceptor, expected_backbone):

        component_list = self._build_component_list(factory, 'i' * 5)
        for i,component in enumerate(component_list):
            INDEX = 0
            DIRECT_ATOM_ID = 1
            INDIRECT_ATOM_ID = 2
            DONOR_OR_ACCEPTOR = 3
            ATOM_TYPE = 4
            BACKBONE = 5
            
            self.assertEqual(i,component[INDEX])
            
            atom_key = Atom_utils._get_atom_info_from_index(component[DIRECT_ATOM_ID])
            indirect_atom_key = Atom_utils._get_atom_info_from_index(component[INDIRECT_ATOM_ID])
            
            if component[BACKBONE] > -1:
                self.assertSequenceContains(atom_key, expected_backbone)
                self.assertEqual(component[BACKBONE], expected_backbone.index(atom_key))
            else:
                self.assertSequenceDoesntContain(atom_key, expected_backbone)
            
            self.assertIn(atom_key, expected_direct_donors_or_acceptors)
            if atom_key in expected_direct_donors_or_acceptors:
                expected_direct_donors_or_acceptors.remove(atom_key)
            
            self.assertEqual(indirect_atom_key,expected_indirect_donors_or_acceptors[atom_key])
            del expected_indirect_donors_or_acceptors[atom_key]
            
            self.assertEqual(component[DONOR_OR_ACCEPTOR], donor_or_acceptor)
            
            if donor_or_acceptor ==  DONOR:
                self.assertEqual(component[ATOM_TYPE],0)
            elif donor_or_acceptor ==  ACCEPTOR:
                self.assertEqual(component[ATOM_TYPE],0)
        
        self.assertEmpty(expected_direct_donors_or_acceptors)
        self.assertEmpty(expected_indirect_donors_or_acceptors)

    def test_hydrogen_bond_donor_components(self):
        self._do_test_donor_acceptor_components(Hydrogen_bond_donor_component_factory(), set(EXPECTED_DIRECT_DONORS), dict(EXPECTED_INDIRECT_DONORS), DONOR, EXPECTED_BACKBONE_DONORS_AND_ACCEPTORS)
        
    def test_hydrogen_bond_acceptor_components(self):
        self._do_test_donor_acceptor_components(Hydrogen_bond_acceptor_component_factory(), set(EXPECTED_DIRECT_ACCEPTORS), dict(EXPECTED_INDIRECT_ACCEPTORS), ACCEPTOR, EXPECTED_BACKBONE_DONORS_AND_ACCEPTORS)
    
    def test_fast_hydrogen_bond_calculator(self):
        test = Fast_hydrogen_bond_calculator(self.get_single_member_ensemble_simulation())
        format = 'i' * 6
        donor_components = self._build_component_list(Hydrogen_bond_donor_component_factory(),format)
        acceptor_components = self._build_component_list(Hydrogen_bond_acceptor_component_factory(),format)
        parameter_format =  ('i'*4) +('f'*8)
        parameter_components = self._build_component_list(Hydrogen_bond_parameter_factory(), parameter_format)
        donor_index_format = 'i'*3
        donor_index_components = self._build_component_list(Hydrogen_bond_donor_lookup_factory(), donor_index_format)
        acceptor_index_format = 'i'*6
        acceptor_index_components = self._build_component_list(Hydrogen_bond_acceptor_lookup_factory(), acceptor_index_format)
        components = {'DONR' : donor_components.get_native_components(), 
                      'ACCP' : acceptor_components.get_native_components(), 
                      'PARA' : parameter_components.get_native_components(),
                      'DIDX' : donor_index_components.get_native_components(),
                      'AIDX' : acceptor_index_components.get_native_components()
                      }
         
        num_donors_and_acceptors = self.donor_and_acceptor_indexer.get_max_index()
         
        energies = allocate_array(num_donors_and_acceptors*3, type='f')
        test(components, energies)
         
        for donor_acceptor_selector in EXPECTED_DONOR_ACCEPTOR_ENERGIES.keys():
            offset =  self.donor_and_acceptor_indexer.get_index_for_key(donor_acceptor_selector[:3])*3+donor_acceptor_selector[4]
            self.assertAlmostEqual(energies[offset]/EXPECTED_DONOR_ACCEPTOR_ENERGIES[donor_acceptor_selector],1.0,places=4)
            del EXPECTED_DONOR_ACCEPTOR_ENERGIES[donor_acceptor_selector]
            energies[offset] = 0.0
  
      
        self.assertEmpty(EXPECTED_DONOR_ACCEPTOR_ENERGIES)

        self.assertSequenceAlmostEqual(energies, [0.0]* len(energies))
        
    def test_parameter_components(self):
        factory = Hydrogen_bond_parameter_factory()
        format = ('i'*4) +('f'*8)
        component_list = self._build_component_list(factory, format)
        self.assertLength(component_list, 3)
        EXPECTED = [(0, 0, 0, 0, -0.40353099999999997, 13.9148, -0.42172500000000002, -3.0801599999999998, 0.98725099999999999, 0.0, 3.0, -16.5), 
                    (1, 1, 0, 0, 12.151400000000001, -1.1225799999999999, -1.1225799999999999, 1.3095300000000001, 1.3095300000000001, 0.0, 15.5, 3.0), 
                    (2, 2, 0, 0, 15.0486, -1.7946899999999999, -1.7946899999999999, 1.4939800000000001, 1.4939800000000001, 0.0, 15.5, 2.5)]
        
        for i in range(3):
            self.assertSequenceEqual(component_list[i][:4], EXPECTED[i][:4])
            self.assertSequenceAlmostEqual(component_list[i][4:], EXPECTED[i][4:])
        component_list.get_native_components()

    def test_donor_index_components(self):
        factory = Hydrogen_bond_donor_lookup_factory()
        format = ('ii')
        component_list = self._build_component_list(factory, format)
        EXPECTED = ((0,1,0),)
        self.assertLength(component_list, 1)
        self.assertSequenceEqual(component_list[0], EXPECTED[0])

    def test_acceptor_index_components(self):
        factory = Hydrogen_bond_acceptor_lookup_factory()
        format = ('ii')
        component_list = self._build_component_list(factory, format)
        EXPECTED = ((0,0,0,0,1,2),)
        self.assertLength(component_list, 1)
        self.assertSequenceEqual(component_list[0], EXPECTED[0])
        
    def test_hydrogen_bond_components(self):
        
        potential = Hydrogen_bond_potential(self.get_single_member_ensemble_simulation())
        hydrogen_bond_component_factory = potential._get_component_list('ATOM')
        expected_shift_offset_data_copy = set(EXPECTED_SHIFT_OFFSET_DATA)

        for component in hydrogen_bond_component_factory.get_all_components():
            target_segid, target_res_num, target_atom =  Atom_utils._get_atom_info_from_index(component[0])
            
            if target_atom == 'HA1' and target_res_num in (9,14):
                target_atom = 'HA' 
                
            hbond_residue_num,hbond_atom = EXPECTED_HBOND_SHIFT_DATA[component[1]][1]
            key = target_res_num, target_atom,hbond_residue_num,hbond_atom,hbond_residue_num-target_res_num

            expected_dist,expected_ang_1,expected_ang_2 = EXPECTED_SHIFT_OFFSET_DATA[key]
            
            self.assertAlmostEqual(((component[2]+epsilon)/(expected_dist+epsilon)),1.0,msg=key)
            self.assertAlmostEqual(((component[3]+epsilon)/(expected_ang_1+epsilon)),1.0,msg=key)
            self.assertAlmostEqual(((component[4]+epsilon)/(expected_ang_2+epsilon)),1.0,msg=key)
            
            expected_shift_offset_data_copy.remove(key)
            
        self.assertEmpty(expected_shift_offset_data_copy)
        
def run_tests():
    unittest2.main(module='test.test_xcamshift_hbond_ingktlkg')
    
if __name__ == "__main__":
    run_tests()
