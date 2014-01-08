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
from common_constants import ROUND_CHANGED, TARGET_ATOM_IDS_CHANGED
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
from cython.shift_calculators import Fast_hydrogen_bond_calculator, allocate_array, Fast_hydrogen_bond_shift_calculator, CDSSharedVectorFloat
from cython.fast_segment_manager import Segment_Manager
from utils import Atom_utils
from table_manager import Table_manager
from atomSel import AtomSel
from component_list import Native_component_list
from array import array

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
            (8, 'N', 7, 'O', -1)        :  ( -0.00007527458,  -0.01804107276,   0.00187546307),
            (8, 'N', 8, 'HN', 0)        :  ( -0.01342813965,   0.00039822518,  -0.00010034386),
            (8, 'N', 8, 'O', 0)         :  ( -0.00008440091,   0.00656272665,  -0.00063741094),
            (8, 'N', 9, 'HN', 1)        :  ( -0.01557293172,   0.01177279380,  -0.00121930908),
            (8, 'HN', 7, 'O', -1)       :  ( -0.00000010025,  -0.00123116240,   0.00012506817),
            (8, 'HN', 8, 'HN', 0)       :  ( -0.01428719399,   0.00448211681,  -0.00051058950),
            (8, 'HN', 8, 'O', 0)        :  ( -0.00000547163,   0.00023294214,  -0.00001567124),
            (8, 'HN', 9, 'HN', 1)       :  (  0.00113213021,   0.00017202282,  -0.00000996565),
            (8, 'CA', 7, 'O', -1)       :  ( -0.00000514598,   0.00145321195,  -0.00013940254),
            (8, 'CA', 8, 'HN', 0)       :  (  0.00302826943,   0.00106609006,  -0.00013263144),
            (8, 'CA', 8, 'O', 0)        :  ( -0.00002568669,  -0.00029179647,   0.00003277032),
            (8, 'CA', 9, 'HN', 1)       :  ( -0.00524412125,   0.00344807243,  -0.00038032800),
            (8, 'HA', 7, 'O', -1)       :  ( -0.00000002820,  -0.00057960923,   0.00006865534),
            (8, 'HA', 8, 'HN', 0)       :  ( -0.00138898021,   0.00010705011,  -0.00000055516),
            (8, 'HA', 8, 'O', 0)        :  (  0.00000114176,   0.00030839949,  -0.00002808213),
            (8, 'HA', 9, 'HN', 1)       :  ( -0.00177689196,   0.00085546479,  -0.00008146647),
            (8, 'CB', 7, 'O', -1)       :  (  0.00003302363,  -0.00163257205,   0.00016794747),
            (8, 'CB', 8, 'HN', 0)       :  ( -0.00660430277,   0.00077341075,  -0.00003180200),
            (8, 'CB', 8, 'O', 0)        :  ( -0.00001790638,   0.00116802227,  -0.00012428080),
            (8, 'CB', 9, 'HN', 1)       :  (  0.00171548083,  -0.00084031610,   0.00010726395),
            (8, 'C', 7, 'O', -1)        :  (  0.00000848419,   0.00398873528,  -0.00039169941),
            (8, 'C', 8, 'HN', 0)        :  ( -0.00433359160,  -0.00048367656,   0.00005464226),
            (8, 'C', 8, 'O', 0)         :  ( -0.00001692652,  -0.01135204363,   0.00116527444),
            (8, 'C', 9, 'HN', 1)        :  ( -0.00773595447,   0.00571298451,  -0.00054428381),
            (9, 'N', 8, 'O', -1)        :  ( -0.00008415758,  -0.01657280308,   0.00170666381),
            (9, 'N', 9, 'HN', 0)        :  ( -0.01261174271,   0.00185645693,  -0.00026359190),
            (9, 'N', 9, 'O', 0)         :  ( -0.00007948963,   0.00553389452,  -0.00052436707),
            (9, 'N', 10, 'HN', 1)       :  ( -0.01328775361,   0.01286293772,  -0.00133454024),
            (9, 'HN', 8, 'O', -1)       :  ( -0.00000072737,  -0.00129235758,   0.00013130866),
            (9, 'HN', 9, 'HN', 0)       :  ( -0.01484773310,   0.00512416682,  -0.00057385688),
            (9, 'HN', 9, 'O', 0)        :  ( -0.00000600656,   0.00002223877,   0.00000629368),
            (9, 'HN', 10, 'HN', 1)      :  (  0.00094982520,   0.00009249361,  -0.00000088809),
            (9, 'CA', 8, 'O', -1)       :  ( -0.00000701771,   0.00060444293,  -0.00004566089),
            (9, 'CA', 9, 'HN', 0)       :  (  0.00304587135,   0.00135725064,  -0.00016355607),
            (9, 'CA', 9, 'O', 0)        :  ( -0.00002471841,  -0.00026369124,   0.00003028417),
            (9, 'CA', 10, 'HN', 1)      :  ( -0.00500531320,   0.00307995577,  -0.00034543584),
            (9, 'HA', 8, 'O', -1)       :  ( -0.00000029080,  -0.00041801001,   0.00005246284),
            (9, 'HA', 9, 'HN', 0)       :  ( -0.00125703883,   0.00013978457,  -0.00000425579),
            (9, 'HA', 9, 'O', 0)        :  (  0.00000147256,   0.00019211417,  -0.00001611504),
            (9, 'HA', 10, 'HN', 1)      :  ( -0.00183993298,   0.00090265292,  -0.00008295994),
            (9, 'C', 8, 'O', -1)        :  (  0.00000988690,   0.00463443653,  -0.00045934464),
            (9, 'C', 9, 'HN', 0)        :  ( -0.00497529123,  -0.00061420316,   0.00006703670),
            (9, 'C', 9, 'O', 0)         :  ( -0.00001649747,  -0.01070691070,   0.00109776839),
            (9, 'C', 10, 'HN', 1)       :  ( -0.00992675184,   0.00602777674,  -0.00057374321),
            (10, 'N', 9, 'O', -1)       :  ( -0.00007527458,  -0.01804107276,   0.00187546307),
            (10, 'N', 10, 'HN', 0)      :  ( -0.01342813965,   0.00039822518,  -0.00010034386),
            (10, 'N', 10, 'O', 0)       :  ( -0.00008440091,   0.00656272665,  -0.00063741094),
            (10, 'N', 11, 'HN', 1)      :  ( -0.01557293172,   0.01177279380,  -0.00121930908),
            (10, 'HN', 9, 'O', -1)      :  ( -0.00000010025,  -0.00123116240,   0.00012506817),
            (10, 'HN', 10, 'HN', 0)     :  ( -0.01428719399,   0.00448211681,  -0.00051058950),
            (10, 'HN', 10, 'O', 0)      :  ( -0.00000547163,   0.00023294214,  -0.00001567124),
            (10, 'HN', 11, 'HN', 1)     :  (  0.00113213021,   0.00017202282,  -0.00000996565),
            (10, 'CA', 9, 'O', -1)      :  ( -0.00000514598,   0.00145321195,  -0.00013940254),
            (10, 'CA', 10, 'HN', 0)     :  (  0.00302826943,   0.00106609006,  -0.00013263144),
            (10, 'CA', 10, 'O', 0)      :  ( -0.00002568669,  -0.00029179647,   0.00003277032),
            (10, 'CA', 11, 'HN', 1)     :  ( -0.00524412125,   0.00344807243,  -0.00038032800),
            (10, 'HA', 9, 'O', -1)      :  ( -0.00000002820,  -0.00057960923,   0.00006865534),
            (10, 'HA', 10, 'HN', 0)     :  ( -0.00138898021,   0.00010705011,  -0.00000055516),
            (10, 'HA', 10, 'O', 0)      :  (  0.00000114176,   0.00030839949,  -0.00002808213),
            (10, 'HA', 11, 'HN', 1)     :  ( -0.00177689196,   0.00085546479,  -0.00008146647),
            (10, 'CB', 9, 'O', -1)      :  (  0.00003302363,  -0.00163257205,   0.00016794747),
            (10, 'CB', 10, 'HN', 0)     :  ( -0.00660430277,   0.00077341075,  -0.00003180200),
            (10, 'CB', 10, 'O', 0)      :  ( -0.00001790638,   0.00116802227,  -0.00012428080),
            (10, 'CB', 11, 'HN', 1)     :  (  0.00171548083,  -0.00084031610,   0.00010726395),
            (10, 'C', 9, 'O', -1)       :  (  0.00000848419,   0.00398873528,  -0.00039169941),
            (10, 'C', 10, 'HN', 0)      :  ( -0.00433359160,  -0.00048367656,   0.00005464226),
            (10, 'C', 10, 'O', 0)       :  ( -0.00001692652,  -0.01135204363,   0.00116527444),
            (10, 'C', 11, 'HN', 1)      :  ( -0.00773595447,   0.00571298451,  -0.00054428381),
            (11, 'N', 10, 'O', -1)      :  ( -0.00007527458,  -0.01804107276,   0.00187546307),
            (11, 'N', 11, 'HN', 0)      :  ( -0.01342813965,   0.00039822518,  -0.00010034386),
            (11, 'N', 11, 'O', 0)       :  ( -0.00008440091,   0.00656272665,  -0.00063741094),
            (11, 'N', 12, 'HN', 1)      :  ( -0.01557293172,   0.01177279380,  -0.00121930908),
            (11, 'HN', 10, 'O', -1)     :  ( -0.00000010025,  -0.00123116240,   0.00012506817),
            (11, 'HN', 11, 'HN', 0)     :  ( -0.01428719399,   0.00448211681,  -0.00051058950),
            (11, 'HN', 11, 'O', 0)      :  ( -0.00000547163,   0.00023294214,  -0.00001567124),
            (11, 'HN', 12, 'HN', 1)     :  (  0.00113213021,   0.00017202282,  -0.00000996565),
            (11, 'CA', 10, 'O', -1)     :  ( -0.00000514598,   0.00145321195,  -0.00013940254),
            (11, 'CA', 11, 'HN', 0)     :  (  0.00302826943,   0.00106609006,  -0.00013263144),
            (11, 'CA', 11, 'O', 0)      :  ( -0.00002568669,  -0.00029179647,   0.00003277032),
            (11, 'CA', 12, 'HN', 1)     :  ( -0.00524412125,   0.00344807243,  -0.00038032800),
            (11, 'HA', 10, 'O', -1)     :  ( -0.00000002820,  -0.00057960923,   0.00006865534),
            (11, 'HA', 11, 'HN', 0)     :  ( -0.00138898021,   0.00010705011,  -0.00000055516),
            (11, 'HA', 11, 'O', 0)      :  (  0.00000114176,   0.00030839949,  -0.00002808213),
            (11, 'HA', 12, 'HN', 1)     :  ( -0.00177689196,   0.00085546479,  -0.00008146647),
            (11, 'CB', 10, 'O', -1)     :  (  0.00003302363,  -0.00163257205,   0.00016794747),
            (11, 'CB', 11, 'HN', 0)     :  ( -0.00660430277,   0.00077341075,  -0.00003180200),
            (11, 'CB', 11, 'O', 0)      :  ( -0.00001790638,   0.00116802227,  -0.00012428080),
            (11, 'CB', 12, 'HN', 1)     :  (  0.00171548083,  -0.00084031610,   0.00010726395),
            (11, 'C', 10, 'O', -1)      :  (  0.00000848419,   0.00398873528,  -0.00039169941),
            (11, 'C', 11, 'HN', 0)      :  ( -0.00433359160,  -0.00048367656,   0.00005464226),
            (11, 'C', 11, 'O', 0)       :  ( -0.00001692652,  -0.01135204363,   0.00116527444),
            (11, 'C', 12, 'HN', 1)      :  ( -0.00773595447,   0.00571298451,  -0.00054428381),
            (12, 'N', 11, 'O', -1)      :  ( -0.00007527458,  -0.01804107276,   0.00187546307),
            (12, 'N', 12, 'HN', 0)      :  ( -0.01342813965,   0.00039822518,  -0.00010034386),
            (12, 'N', 12, 'O', 0)       :  ( -0.00008440091,   0.00656272665,  -0.00063741094),
            (12, 'N', 13, 'HN', 1)      :  ( -0.01557293172,   0.01177279380,  -0.00121930908),
            (12, 'HN', 11, 'O', -1)     :  ( -0.00000010025,  -0.00123116240,   0.00012506817),
            (12, 'HN', 12, 'HN', 0)     :  ( -0.01428719399,   0.00448211681,  -0.00051058950),
            (12, 'HN', 12, 'O', 0)      :  ( -0.00000547163,   0.00023294214,  -0.00001567124),
            (12, 'HN', 13, 'HN', 1)     :  (  0.00113213021,   0.00017202282,  -0.00000996565),
            (12, 'CA', 11, 'O', -1)     :  ( -0.00000514598,   0.00145321195,  -0.00013940254),
            (12, 'CA', 12, 'HN', 0)     :  (  0.00302826943,   0.00106609006,  -0.00013263144),
            (12, 'CA', 12, 'O', 0)      :  ( -0.00002568669,  -0.00029179647,   0.00003277032),
            (12, 'CA', 13, 'HN', 1)     :  ( -0.00524412125,   0.00344807243,  -0.00038032800),
            (12, 'HA', 11, 'O', -1)     :  ( -0.00000002820,  -0.00057960923,   0.00006865534),
            (12, 'HA', 12, 'HN', 0)     :  ( -0.00138898021,   0.00010705011,  -0.00000055516),
            (12, 'HA', 12, 'O', 0)      :  (  0.00000114176,   0.00030839949,  -0.00002808213),
            (12, 'HA', 13, 'HN', 1)     :  ( -0.00177689196,   0.00085546479,  -0.00008146647),
            (12, 'CB', 11, 'O', -1)     :  (  0.00003302363,  -0.00163257205,   0.00016794747),
            (12, 'CB', 12, 'HN', 0)     :  ( -0.00660430277,   0.00077341075,  -0.00003180200),
            (12, 'CB', 12, 'O', 0)      :  ( -0.00001790638,   0.00116802227,  -0.00012428080),
            (12, 'CB', 13, 'HN', 1)     :  (  0.00171548083,  -0.00084031610,   0.00010726395),
            (12, 'C', 11, 'O', -1)      :  (  0.00000848419,   0.00398873528,  -0.00039169941),
            (12, 'C', 12, 'HN', 0)      :  ( -0.00433359160,  -0.00048367656,   0.00005464226),
            (12, 'C', 12, 'O', 0)       :  ( -0.00001692652,  -0.01135204363,   0.00116527444),
            (12, 'C', 13, 'HN', 1)      :  ( -0.00773595447,   0.00571298451,  -0.00054428381),
            (13, 'N', 12, 'O', -1)      :  ( -0.00007527458,  -0.01804107276,   0.00187546307),
            (13, 'N', 13, 'HN', 0)      :  ( -0.01342813965,   0.00039822518,  -0.00010034386),
            (13, 'N', 13, 'O', 0)       :  ( -0.00008440091,   0.00656272665,  -0.00063741094),
            (13, 'N', 14, 'HN', 1)      :  ( -0.01557293172,   0.01177279380,  -0.00121930908),
            (13, 'HN', 12, 'O', -1)     :  ( -0.00000010025,  -0.00123116240,   0.00012506817),
            (13, 'HN', 13, 'HN', 0)     :  ( -0.01428719399,   0.00448211681,  -0.00051058950),
            (13, 'HN', 13, 'O', 0)      :  ( -0.00000547163,   0.00023294214,  -0.00001567124),
            (13, 'HN', 14, 'HN', 1)     :  (  0.00113213021,   0.00017202282,  -0.00000996565),
            (13, 'CA', 12, 'O', -1)     :  ( -0.00000514598,   0.00145321195,  -0.00013940254),
            (13, 'CA', 13, 'HN', 0)     :  (  0.00302826943,   0.00106609006,  -0.00013263144),
            (13, 'CA', 13, 'O', 0)      :  ( -0.00002568669,  -0.00029179647,   0.00003277032),
            (13, 'CA', 14, 'HN', 1)     :  ( -0.00524412125,   0.00344807243,  -0.00038032800),
            (13, 'HA', 12, 'O', -1)     :  ( -0.00000002820,  -0.00057960923,   0.00006865534),
            (13, 'HA', 13, 'HN', 0)     :  ( -0.00138898021,   0.00010705011,  -0.00000055516),
            (13, 'HA', 13, 'O', 0)      :  (  0.00000114176,   0.00030839949,  -0.00002808213),
            (13, 'HA', 14, 'HN', 1)     :  ( -0.00177689196,   0.00085546479,  -0.00008146647),
            (13, 'CB', 12, 'O', -1)     :  (  0.00003302363,  -0.00163257205,   0.00016794747),
            (13, 'CB', 13, 'HN', 0)     :  ( -0.00660430277,   0.00077341075,  -0.00003180200),
            (13, 'CB', 13, 'O', 0)      :  ( -0.00001790638,   0.00116802227,  -0.00012428080),
            (13, 'CB', 14, 'HN', 1)     :  (  0.00171548083,  -0.00084031610,   0.00010726395),
            (13, 'C', 12, 'O', -1)      :  (  0.00000848419,   0.00398873528,  -0.00039169941),
            (13, 'C', 13, 'HN', 0)      :  ( -0.00433359160,  -0.00048367656,   0.00005464226),
            (13, 'C', 13, 'O', 0)       :  ( -0.00001692652,  -0.01135204363,   0.00116527444),
            (13, 'C', 14, 'HN', 1)      :  ( -0.00773595447,   0.00571298451,  -0.00054428381)
}

CAMSHIFT_EXPECTED_SHIFTS = {
            # correction for an extra h bond at '',7,'O' - '',14,'HN'
#             ('', 8,'N',('O',-1))      :      (  -0.000246806004409, 3.27874302864, -12.56614343240),
#             ('', 8,'HN',('O',-1))     :      (  -3.28693988621e-07, 3.27874302864,  -0.85754120682),
#             ('', 8,'CA',('O',-1))     :      (  -1.68723460505e-05, 3.27874302864,   1.01220531862),
#             ('', 8,'HA',('O',-1))     :      (  -9.24605534077e-08, 3.27874302864,  -0.40371505707),
#             ('', 8,'CB',('O',-1))     :      (   0.000108275996643, 3.27874302864,  -1.13713495959),
#             ('', 8,'C',('O',-1))      :      (   2.78174788162e-05, 3.27874302864,   2.77827268417),
            
            ('',8,'HA',('H',1))       :      (  -0.006474500,    0.599857000,   -0.533485000),
            ('',9,'HA',('H',0))       :      (  -0.004580300,    0.098017700,   -0.027869100),
            ('',9,'HA',('O',0))       :      (   0.000006546,    0.135957000,   -0.104492000),
            ('',10,'HA',('O',-1))     :      (  -0.000000125,   -0.410183000,    0.445172000),
            ('',10,'HA',('H',1))      :      (  -0.007912250,    0.575547000,   -0.509351000),
            ('',11,'HA',('H',0))      :      (  -0.006184930,    0.072022100,   -0.003471000),
            ('',11,'HA',('H',1))      :      (  -0.008289930,    0.594569000,   -0.508204000),
            ('',12,'HA',('H',0))      :      (  -0.006480160,    0.074402500,   -0.003463190),
            ('',12,'HA',('O',0))      :      (   0.000004160,    0.216251000,   -0.183896000),
            ('',13,'HA',('O',-1))     :      (  -0.000000103,   -0.406425000,    0.449591000),
            ('',8,'CA',('H',1))       :      (  -0.019108100,    2.417810000,   -2.490590000),
            ('',9,'CA',('H',0))       :      (   0.011098300,    0.951711000,   -1.071050000),
            ('',9,'CA',('O',0))       :      (  -0.000109881,   -0.186612000,    0.196367000),
            ('',10,'CA',('O',-1))     :      (  -0.000022875,    1.028420000,   -0.903907000),
            ('',10,'CA',('H',1))      :      (  -0.023351300,    2.319830000,   -2.377920000),
            ('',11,'CA',('H',0))      :      (   0.013484500,    0.717254000,   -0.829249000),
            ('',11,'CA',('H',1))      :      (  -0.024466000,    2.396500000,   -2.372560000),
            ('',12,'CA',('H',0))      :      (   0.014128100,    0.740959000,   -0.827382000),
            ('',12,'CA',('O',0))      :      (  -0.000093595,   -0.204609000,    0.214597000),
            ('',13,'CA',('O',-1))     :      (  -0.000018750,    1.019000000,   -0.912880000),
            ('',8,'HN',('H',1))       :      (   0.004125170,    0.120623000,   -0.065260200),
            ('',9,'HN',('H',0))       :      (  -0.054101000,    3.593090000,   -3.757910000),
            ('',9,'HN',('O',0))       :      (  -0.000026701,    0.015738100,    0.040809200),
            ('',10,'HN',('O',-1))     :      (  -0.000000446,   -0.871281000,    0.810961000),
            ('',10,'HN',('H',1))      :      (   0.005041210,    0.115735000,   -0.062308000),
            ('',11,'HN',('H',0))      :      (  -0.063618800,    3.015520000,   -3.192350000),
            ('',11,'HN',('H',1))      :      (   0.005281850,    0.119560000,   -0.062167700),
            ('',12,'HN',('H',0))      :      (  -0.066655600,    3.115180000,   -3.185160000),
            ('',12,'HN',('O',0))      :      (  -0.000019937,    0.163340000,   -0.102623000),
            ('',13,'HN',('O',-1))     :      (  -0.000000365,   -0.863298000,    0.819011000),
            ('',8,'N',('H',1))        :      (  -0.056743400,    8.255150000,   -7.984670000),
            ('',9,'N',('H',0))        :      (  -0.045953700,    1.301760000,   -1.726140000),
            ('',9,'N',('O',0))        :      (  -0.000353355,    3.916280000,   -3.400070000),
            ('',10,'N',('O',-1))      :      (  -0.000334618,  -12.767500000,   12.160800000),
            ('',10,'N',('H',1))       :      (  -0.069344000,    7.920610000,   -7.623460000),
            ('',11,'N',('H',0))       :      (  -0.059793600,    0.267922000,   -0.627378000),
            ('',11,'N',('H',1))       :      (  -0.072654100,    8.182390000,   -7.606300000),
            ('',12,'N',('H',0))       :      (  -0.062647700,    0.276776000,   -0.625965000),
            ('',12,'N',('O',0))       :      (  -0.000307533,    4.601820000,   -4.174100000),
            ('',13,'N',('O',-1))      :      (  -0.000274280,  -12.650500000,   12.281500000),
            ('',8,'C',('H',1))        :      (  -0.028187700,    4.005980000,   -3.564250000),
            ('',9,'C',('H',0))        :      (  -0.018128600,   -0.430683000,    0.438991000),
            ('',9,'C',('O',0))        :      (  -0.000073336,   -7.577170000,    7.118090000),
            ('',10,'C',('O',-1))      :      (   0.000037715,    2.822790000,   -2.539840000),
            ('',10,'C',('H',1))       :      (  -0.034447100,    3.843630000,   -3.403010000),
            ('',11,'C',('H',0))       :      (  -0.019296900,   -0.325412000,    0.341639000),
            ('',11,'C',('H',1))       :      (  -0.036091400,    3.970670000,   -3.395350000),
            ('',12,'C',('H',0))       :      (  -0.020218000,   -0.336167000,    0.340870000),
            ('',12,'C',('O',0))       :      (  -0.000061676,   -7.960110000,    7.630820000),
            ('',13,'C',('O',-1))      :      (   0.000030914,    2.796920000,   -2.565050000),
            ('',8,'CB',('H',1))       :      (   0.006250740,   -0.589234000,    0.702420000),
            ('',10,'CB',('O',-1))     :      (   0.000146800,   -1.155350000,    1.089000000),
            ('',10,'CB',('H',1))      :      (   0.007638790,   -0.565355000,    0.670644000),
            ('',11,'CB',('H',0))      :      (  -0.029408000,    0.520342000,   -0.198835000),
            ('',11,'CB',('H',1))      :      (   0.008003420,   -0.584041000,    0.669134000),
            ('',12,'CB',('H',0))      :      (  -0.030811800,    0.537540000,   -0.198387000),
            ('',12,'CB',('O',0))      :      (  -0.000065246,    0.819023000,   -0.813855000),
            ('',13,'CB',('O',-1))     :      (   0.000120329,   -1.144770000,    1.099810000),
}

XCAMSHIFT_EXPECTED_SHIFTS = {
            ('', 8, 'N', -1, 'O')     :      (  -0.0002468060, -12.5661434324,  12.5346015674),
            ('', 8, 'N', 1, 'HN')     :      (  -0.0567437096,   8.2551514189,  -7.9846776601),
            ('', 8, 'HN', -1, 'O')    :      (  -0.0000003287,  -0.8575412068,   0.8358893890),
            ('', 8, 'HN', 1, 'HN')    :      (   0.0041251878,   0.1206234009,  -0.0652603218),
            ('', 8, 'CA', -1, 'O')    :      (  -0.0000168723,   1.0122053186,  -0.9316927240),
            ('', 8, 'CA', 1, 'HN')    :      (  -0.0191082128,   2.4178084231,  -2.4905879362),
            ('', 8, 'HA', -1, 'O')    :      (  -0.0000000925,  -0.4037150571,   0.4588559200),
            ('', 8, 'HA', 1, 'HN')    :      (  -0.0064745318,   0.5998568815,  -0.5334853268),
            ('', 8, 'CB', -1, 'O')    :      (   0.0001082760,  -1.1371349596,   1.1224719134),
            ('', 8, 'CB', 1, 'HN')    :      (   0.0062507656,  -0.5892345320,   0.7024208048),
            ('', 8, 'C', -1, 'O')     :      (   0.0000278175,   2.7782726842,  -2.6179113399),
            ('', 8, 'C', 1, 'HN')     :      (  -0.0281878044,   4.0059779340,  -3.5642568810),
            ('', 9, 'N', 0, 'HN')     :      (  -0.0459539076,   1.3017583863,  -1.7261385073),
            ('', 9, 'HN', 0, 'HN')    :      (  -0.0541012746,   3.5930955484,  -3.7579169097),
            ('', 9, 'CA', 0, 'HN')    :      (   0.0110983624,   0.9517120351,  -1.0710512369),
            ('', 9, 'HA1', 0, 'HN')   :      (  -0.0045803223,   0.0980177527,  -0.0278691530),
            ('', 9, 'C', 0, 'HN')     :      (  -0.0181286662,  -0.4306828246,   0.4389915975),
            ('', 10, 'N', 1, 'HN')    :      (  -0.0693441126,   7.9206033109,  -7.6234668710),
            ('', 10, 'HN', 1, 'HN')   :      (   0.0050412193,   0.1157350193,  -0.0623080758),
            ('', 10, 'CA', 1, 'HN')   :      (  -0.0233513471,   2.3198243653,  -2.3779187375),
            ('', 10, 'HA', 1, 'HN')   :      (  -0.0079122543,   0.5755470930,  -0.5093515216),
            ('', 10, 'CB', 1, 'HN')   :      (   0.0076387990,  -0.5653552247,   0.6706446976),
            ('', 10, 'C', 1, 'HN')    :      (  -0.0344471360,   3.8436317491,  -3.4030170545),
            ('', 11, 'N', 0, 'HN')    :      (  -0.0597936499,   0.2679214240,  -0.6273783284),
            ('', 11, 'N', 1, 'HN')    :      (  -0.0726542086,   8.1823920467,  -7.6063024761),
            ('', 11, 'HN', 0, 'HN')   :      (  -0.0636189001,   3.0155177988,  -3.1923506532),
            ('', 11, 'HN', 1, 'HN')   :      (   0.0052818587,   0.1195602487,  -0.0621677879),
            ('', 11, 'CA', 0, 'HN')   :      (   0.0134844652,   0.7172534067,  -0.8292494540),
            ('', 11, 'CA', 1, 'HN')   :      (  -0.0244660084,   2.3964983085,  -2.3725648038),
            ('', 11, 'HA', 0, 'HN')   :      (  -0.0061849369,   0.0720221105,  -0.0034710181),
            ('', 11, 'HA', 1, 'HN')   :      (  -0.0082899406,   0.5945698543,  -0.5082047060),
            ('', 11, 'CB', 0, 'HN')   :      (  -0.0294080474,   0.5203420574,  -0.1988351415),
            ('', 11, 'CB', 1, 'HN')   :      (   0.0080034321,  -0.5840411282,   0.6691347271),
            ('', 11, 'C', 0, 'HN')    :      (  -0.0192968844,  -0.3254121517,   0.3416389377),
            ('', 11, 'C', 1, 'HN')    :      (  -0.0360914476,   3.9706699881,  -3.3953550905),
            ('', 12, 'N', 0, 'HN')    :      (  -0.0626478608,   0.2767766599,  -0.6259657730),
            ('', 12, 'N', 0, 'O')     :      (  -0.0003075349,   4.6018220600,  -4.1741023473),
            ('', 12, 'HN', 0, 'HN')   :      (  -0.0666557069,   3.1151855338,  -3.1851630090),
            ('', 12, 'HN', 0, 'O')    :      (  -0.0000199372,   0.1633403821,  -0.1026235283),
            ('', 12, 'CA', 0, 'HN')   :      (   0.0141281374,   0.7409597905,  -0.8273823816),
            ('', 12, 'CA', 0, 'O')    :      (  -0.0000935956,  -0.2046093803,   0.2145973046),
            ('', 12, 'HA', 0, 'HN')   :      (  -0.0064801708,   0.0744025576,  -0.0034632030),
            ('', 12, 'HA', 0, 'O')    :      (   0.0000041603,   0.2162515144,  -0.1838965688),
            ('', 12, 'CB', 0, 'HN')   :      (  -0.0308118214,   0.5375402030,  -0.1983874600),
            ('', 12, 'CB', 0, 'O')    :      (  -0.0000652462,   0.8190240026,  -0.8138560957),
            ('', 12, 'C', 0, 'HN')    :      (  -0.0202180086,  -0.3361675491,   0.3408697305),
            ('', 12, 'C', 0, 'O')     :      (  -0.0000616758,  -7.9601189549,   7.6308303953),
            ('', 13, 'N', -1, 'O')    :      (  -0.0002742810, -12.6505050478,  12.2815193646),
            ('', 13, 'N', 1, 'HN')    :      (  -0.0510596413,   8.2001008176,  -8.1492159189),
            ('', 13, 'HN', -1, 'O')   :      (  -0.0000003653,  -0.8632982286,   0.8190122089),
            ('', 13, 'HN', 1, 'HN')   :      (   0.0037119640,   0.1198190074,  -0.0666051249),
            ('', 13, 'CA', -1, 'O')   :      (  -0.0000187506,   1.0190006633,  -0.9128812088),
            ('', 13, 'CA', 1, 'HN')   :      (  -0.0171941260,   2.4016849384,  -2.5419108599),
            ('', 13, 'HA', -1, 'O')   :      (  -0.0000001028,  -0.4064253599,   0.4495913042),
            ('', 13, 'HA', 1, 'HN')   :      (  -0.0058259721,   0.5958566542,  -0.5444787258),
            ('', 13, 'CB', -1, 'O')   :      (   0.0001203295,  -1.1447690076,   1.0998084356),
            ('', 13, 'CB', 1, 'HN')   :      (   0.0056246208,  -0.5853051413,   0.7168954150),
            ('', 13, 'C', -1, 'O')    :      (   0.0000309142,   2.7969243551,  -2.5650539144),
            ('', 13, 'C', 1, 'HN')    :      (  -0.0253642068,   3.9792635247,  -3.6377046326)
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


    def _calculate_energies(self):
        test = Fast_hydrogen_bond_calculator(self.get_single_member_ensemble_simulation())
        format = 'i' * 6
        donor_components = self._build_component_list(Hydrogen_bond_donor_component_factory(), format)
        acceptor_components = self._build_component_list(Hydrogen_bond_acceptor_component_factory(), format)
        parameter_format = ('i' * 4) + ('f' * 8)
        parameter_components = self._build_component_list(Hydrogen_bond_parameter_factory(), parameter_format)
        donor_index_format = 'i' * 3
        donor_index_components = self._build_component_list(Hydrogen_bond_donor_lookup_factory(), donor_index_format)
        acceptor_index_format = 'i' * 6
        acceptor_index_components = self._build_component_list(Hydrogen_bond_acceptor_lookup_factory(), acceptor_index_format)
        components = {'DONR':donor_components.get_native_components(), 'ACCP':acceptor_components.get_native_components(), 
            'PARA':parameter_components.get_native_components(), 
            'DIDX':donor_index_components.get_native_components(), 
            'AIDX':acceptor_index_components.get_native_components()}
        num_donors_and_acceptors = self.donor_and_acceptor_indexer.get_max_index()
        energies = allocate_array(num_donors_and_acceptors * 3, type='f')
        test(components, energies)
        return energies
    
    def test_fast_hydrogen_bond_calculator(self):
        energies = self._calculate_energies()
         
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
            self.assertAlmostEqual(((component[2]+epsilon)/(expected_dist+epsilon)),1.0,places=4,msg=key)
            self.assertAlmostEqual(((component[3]+epsilon)/(expected_ang_1+epsilon)),1.0,places=4,msg=key)
            self.assertAlmostEqual(((component[4]+epsilon)/(expected_ang_2+epsilon)),1.0,places=4,msg=key)
            
            expected_shift_offset_data_copy.remove(key)
            
        self.assertEmpty(expected_shift_offset_data_copy)

    def _get_total_shifts(self,data):
        result  =  {}
        for full_key,values in data.items():
            key = full_key[0],full_key[1],full_key[2]
            for i,value in enumerate(values):
                current_value = result.setdefault(key,0.0)
                result[key] =  current_value + value
        return result
    
    def test_shift_calculator(self):
        expected_shifts = self._get_total_shifts(CAMSHIFT_EXPECTED_SHIFTS)
        simulation = self.get_single_member_ensemble_simulation()
        
        potential = Hydrogen_bond_potential(simulation)
        components = potential._get_components()
        calculator  = Fast_hydrogen_bond_shift_calculator(components['SIMU'],"test")
        energies = self._calculate_energies()

        
        #modify the energies for testing to allow for some camshift idiosyncrasies
        energies[0] = 0.0       # camshift doesn't see the hbond 7 O - 14 HN
        energies[1]= 0.0
        energies[2]= 0.0
        energies[60] = 0.0      # camshift doesn't see the hbond 7 O - 14 HN
        energies[61]= 0.0
        energies[62]= 0.0
        energies[18] = 4.4453   # camshift includes 9  O  (GLY)     - [ HE2 LYS 13
        energies[19] = 707.69   # as a hydrogen bond... so we add an extra energy
        energies[20] = 6484.15
        

        components['HBLT'] =  energies
        atom_components = potential._get_component_list('ATOM')
        
        target_atom_ids  = []
        for i,component in enumerate(atom_components):
            if component[0] not in target_atom_ids:
                target_atom_ids.append(component[0])
             
        num_targets = len(target_atom_ids)

        shift_cache = CDSSharedVectorFloat(num_targets, simulation)
        component_to_target = [target_atom_ids.index(atom_component[0]) for atom_component in atom_components]
        component_to_target =  array('i',component_to_target)

        calculator(components, shift_cache, component_to_target, None)
        for index,atom_id in enumerate(target_atom_ids):
            atom_key = Atom_utils._get_atom_info_from_index(atom_id)
            if atom_key == ('',9,'HA1'):
                atom_key = ('',9,'HA')
            
            if abs(shift_cache[index]-expected_shifts[atom_key]) > 1.0e-4:
                print 'bad',atom_key,shift_cache[index],expected_shifts[atom_key],shift_cache[index]-expected_shifts[atom_key],index

    def test_potential_shifts(self):
        expected_shifts = self._get_total_shifts(XCAMSHIFT_EXPECTED_SHIFTS)
        
        simulation = self.get_single_member_ensemble_simulation()
        
        potential = Hydrogen_bond_potential(simulation)
        target_atom_ids = potential.get_target_atom_ids()
        
        potential._prepare(TARGET_ATOM_IDS_CHANGED,target_atom_ids)
        potential._prepare(ROUND_CHANGED,None)
        result = CDSSharedVectorFloat(len(target_atom_ids),simulation)
        potential.calc_shifts(target_atom_ids,result)

        for i,target_atom_id in enumerate(target_atom_ids):
            key = Atom_utils._get_atom_info_from_index(target_atom_id)
            self.assertAlmostEqual(result[i],expected_shifts[key])
            del expected_shifts[key]
        self.assertEmpty(expected_shifts)
                      
def run_tests():
    unittest2.main(module='test.test_xcamshift_hbond_ingktlkg')
    
if __name__ == "__main__":
    run_tests()
