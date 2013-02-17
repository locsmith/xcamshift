#-------------------------------------------------------------------------------
# Copyright (c) 2013 gary thompson.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the GNU Lesser Public License v2.1
# which accompanies this distribution, and is available at
# http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
# 
# Contributors:
#     gary thompson - initial API and implementation
#-------------------------------------------------------------------------------
'''
Created on 29 Mar 2012

@author: garyt
'''
 
from common_constants import BACK_BONE, DIHEDRAL, NON_BONDED, SIDE_CHAIN, XTRA, RANDOM_COIL, RING

agfa_shifts  =   {
    ('',2,'HA')   :    3.682520,
    ('',3,'HA')   :    4.620800,
    ('',2,'CA')   :   46.114000,
    ('',3,'CA')   :   56.881600,
    ('',2,'H')    :    8.748080,
    ('',3,'H')    :    9.180120,
    ('',2,'N')    :  110.690000,
    ('',3,'N')    :  122.149000,
    ('',2,'C')    :  174.633000,
    ('',3,'C')    :  174.222000,
    ('',3,'CB')   :   40.769300,
}
agfa_subpotential_shifts  =   {
    ('',2,'HA',RANDOM_COIL)    :   3.97000,
    ('',2,'HA',BACK_BONE)      :   0.85809,
    ('',2,'HA',SIDE_CHAIN)     :   0.10556,
    ('',2,'HA',NON_BONDED)     :  -0.00086,
    ('',2,'HA',XTRA)           :  -1.06630,
    ('',2,'HA',RING)           :   0.05199,
    ('',2,'HA',DIHEDRAL)       :  -0.23596,
    ('',3,'HA',RANDOM_COIL)    :   4.43024,
    ('',3,'HA',BACK_BONE)      :   0.68642,
    ('',3,'HA',SIDE_CHAIN)     :   0.00000,
    ('',3,'HA',NON_BONDED)     :   0.08199,
    ('',3,'HA',XTRA)           :  -0.99700,
    ('',3,'HA',RING)           :   0.26451,
    ('',3,'HA',DIHEDRAL)       :   0.15463,
    ('',2,'CA',RANDOM_COIL)    :  45.48650,
    ('',2,'CA',BACK_BONE)      :  33.26630,
    ('',2,'CA',SIDE_CHAIN)     :  -3.81784,
    ('',2,'CA',NON_BONDED)     :  -0.00043,
    ('',2,'CA',XTRA)           : -29.81120,
    ('',2,'CA',RING)           :   0.03352,
    ('',2,'CA',DIHEDRAL)       :   0.95716,
    ('',3,'CA',RANDOM_COIL)    :  57.64620,
    ('',3,'CA',BACK_BONE)      :  42.70070,
    ('',3,'CA',SIDE_CHAIN)     :  -1.49001,
    ('',3,'CA',NON_BONDED)     :  -0.05797,
    ('',3,'CA',XTRA)           : -40.98990,
    ('',3,'CA',RING)           :   0.17457,
    ('',3,'CA',DIHEDRAL)       :  -1.10204,
    ('',2,'HN',RANDOM_COIL)    :   8.33000,
    ('',2,'HN',BACK_BONE)      :  -7.47515,
    ('',2,'HN',SIDE_CHAIN)     :   0.07021,
    ('',2,'HN',NON_BONDED)     :   0.00000,
    ('',2,'HN',XTRA)           :   7.77572,
    ('',2,'HN',RING)           :   0.03383,
    ('',2,'HN',DIHEDRAL)       :   0.01348,
    ('',3,'HN',RANDOM_COIL)    :   8.30000,
    ('',3,'HN',BACK_BONE)      :  -4.49337,
    ('',3,'HN',SIDE_CHAIN)     :  -0.27149,
    ('',3,'HN',NON_BONDED)     :   0.34671,
    ('',3,'HN',XTRA)           :   5.17056,
    ('',3,'HN',RING)           :   0.14539,
    ('',3,'HN',DIHEDRAL)       :  -0.01768,
    ('',2,'N',RANDOM_COIL)     : 108.60000,
    ('',2,'N',BACK_BONE)       :  50.83020,
    ('',2,'N',SIDE_CHAIN)      :   2.87169,
    ('',2,'N',NON_BONDED)      :   0.00000,
    ('',2,'N',XTRA)            : -50.89940,
    ('',2,'N',RING)            :  -0.02586,
    ('',2,'N',DIHEDRAL)        :  -0.68631,
    ('',3,'N',RANDOM_COIL)     : 121.37000,
    ('',3,'N',BACK_BONE)       :  -9.82297,
    ('',3,'N',SIDE_CHAIN)      : -16.63280,
    ('',3,'N',NON_BONDED)      :   1.33040,
    ('',3,'N',XTRA)            :  26.52700,
    ('',3,'N',RING)            :  -0.12315,
    ('',3,'N',DIHEDRAL)        :  -0.49977,
    ('',2,'C',RANDOM_COIL)     : 173.60000,
    ('',2,'C',BACK_BONE)       : -11.61210,
    ('',2,'C',SIDE_CHAIN)      :   0.99617,
    ('',2,'C',NON_BONDED)      :   0.26502,
    ('',2,'C',XTRA)            :  11.15990,
    ('',2,'C',RING)            :   0.28962,
    ('',2,'C',DIHEDRAL)        :  -0.06535,
    ('',3,'C',RANDOM_COIL)     : 175.80000,
    ('',3,'C',BACK_BONE)       : -13.20110,
    ('',3,'C',SIDE_CHAIN)      :  -0.56996,
    ('',3,'C',NON_BONDED)      :  -0.11972,
    ('',3,'C',XTRA)            :  11.68820,
    ('',3,'C',RING)            :   0.17460,
    ('',3,'C',DIHEDRAL)        :   0.44968,
    ('',3,'CB',RANDOM_COIL)    :  39.30000,
    ('',3,'CB',BACK_BONE)      :  -5.37522,
    ('',3,'CB',SIDE_CHAIN)     :   0.00000,
    ('',3,'CB',NON_BONDED)     :   0.02288,
    ('',3,'CB',XTRA)           :   4.42699,
    ('',3,'CB',RING)           :   0.92892,
    ('',3,'CB',DIHEDRAL)       :   1.46570,
}
agfa_component_shifts_bb =  {
    (('',2,'HA'),('',1,'N'))        :   0.34649,
    (('',2,'HA'),('',1,'CA'))       :   0.00000,
    (('',2,'HA'),('',1,'HA'))       :   0.37253,
    (('',2,'HA'),('',1,'C'))        :   0.00000,
    (('',2,'HA'),('',1,'O'))        :  -0.87926,
    (('',2,'HA'),('',3,'N'))        :   0.00000,
    (('',2,'HA'),('',3,'HN'))       :   0.30064,
    (('',2,'HA'),('',3,'CA'))       :   0.25651,
    (('',2,'HA'),('',3,'HA'))       :   0.10066,
    (('',2,'HA'),('',3,'C'))        :   0.49062,
    (('',2,'HA'),('',2,'N'))        :   0.00000,
    (('',2,'HA'),('',2,'HN'))       :  -0.44828,
    (('',2,'HA'),('',2,'CA'))       :   0.00000,
    (('',2,'HA'),('',2,'HA'))       :   0.00000,
    (('',2,'HA'),('',2,'C'))        :   0.00000,
    (('',2,'HA'),('',2,'O'))        :   0.31820,
    (('',3,'HA'),('',2,'N'))        :   0.49225,
    (('',3,'HA'),('',2,'CA'))       :   0.00000,
    (('',3,'HA'),('',2,'HA'))       :   0.39453,
    (('',3,'HA'),('',2,'C'))        :   0.00000,
    (('',3,'HA'),('',2,'O'))        :  -1.14327,
    (('',3,'HA'),('',4,'N'))        :   0.00000,
    (('',3,'HA'),('',4,'HN'))       :   0.52482,
    (('',3,'HA'),('',4,'CA'))       :  -0.24179,
    (('',3,'HA'),('',4,'HA'))       :   0.18341,
    (('',3,'HA'),('',4,'C'))        :   0.64053,
    (('',3,'HA'),('',3,'N'))        :   0.00000,
    (('',3,'HA'),('',3,'HN'))       :  -0.60734,
    (('',3,'HA'),('',3,'CA'))       :   0.00000,
    (('',3,'HA'),('',3,'HA'))       :  -0.00000,
    (('',3,'HA'),('',3,'C'))        :   0.00000,
    (('',3,'HA'),('',3,'O'))        :   0.44329,
    (('',2,'CA'),('',1,'N'))        :  -2.37688,
    (('',2,'CA'),('',1,'CA'))       :  -2.02046,
    (('',2,'CA'),('',1,'HA'))       :  -6.89848,
    (('',2,'CA'),('',1,'C'))        :   9.42163,
    (('',2,'CA'),('',1,'O'))        :  -2.91953,
    (('',2,'CA'),('',3,'N'))        :  -3.64285,
    (('',2,'CA'),('',3,'HN'))       :   2.02788,
    (('',2,'CA'),('',3,'CA'))       :   3.98649,
    (('',2,'CA'),('',3,'HA'))       :  -0.01275,
    (('',2,'CA'),('',3,'C'))        :   2.62992,
    (('',2,'CA'),('',2,'N'))        :   0.00000,
    (('',2,'CA'),('',2,'HN'))       :   5.96767,
    (('',2,'CA'),('',2,'CA'))       :   0.00000,
    (('',2,'CA'),('',2,'HA'))       :   7.15802,
    (('',2,'CA'),('',2,'C'))        :  10.34850,
    (('',2,'CA'),('',2,'O'))        :   9.59711,
    (('',3,'CA'),('',2,'N'))        :  -2.49116,
    (('',3,'CA'),('',2,'CA'))       :  -1.20547,
    (('',3,'CA'),('',2,'HA'))       :  -9.10961,
    (('',3,'CA'),('',2,'C'))        :  12.23470,
    (('',3,'CA'),('',2,'O'))        :  -2.61998,
    (('',3,'CA'),('',4,'N'))        :   2.15857,
    (('',3,'CA'),('',4,'HN'))       :   2.68380,
    (('',3,'CA'),('',4,'CA'))       :   3.20494,
    (('',3,'CA'),('',4,'HA'))       :  -0.02794,
    (('',3,'CA'),('',4,'C'))        :   2.70097,
    (('',3,'CA'),('',3,'N'))        :   0.00000,
    (('',3,'CA'),('',3,'HN'))       :   8.31036,
    (('',3,'CA'),('',3,'CA'))       :  -0.00000,
    (('',3,'CA'),('',3,'HA'))       :   8.51123,
    (('',3,'CA'),('',3,'C'))        :   9.47428,
    (('',3,'CA'),('',3,'O'))        :   8.87605,
    (('',2,'HN'),('',1,'N'))        :  -2.68066,
    (('',2,'HN'),('',1,'CA'))       :   1.77948,
    (('',2,'HN'),('',1,'HA'))       :  -2.54176,
    (('',2,'HN'),('',1,'C'))        :   0.00000,
    (('',2,'HN'),('',1,'O'))        :   0.00000,
    (('',2,'HN'),('',3,'N'))        :  -3.77156,
    (('',2,'HN'),('',3,'HN'))       :   0.00000,
    (('',2,'HN'),('',3,'CA'))       :   2.42980,
    (('',2,'HN'),('',3,'HA'))       :  -0.31707,
    (('',2,'HN'),('',3,'C'))        :  -0.48296,
    (('',2,'HN'),('',2,'N'))        :   0.00000,
    (('',2,'HN'),('',2,'HN'))       :   0.00000,
    (('',2,'HN'),('',2,'CA'))       :   0.00000,
    (('',2,'HN'),('',2,'HA'))       :   0.03740,
    (('',2,'HN'),('',2,'C'))        :   1.66323,
    (('',2,'HN'),('',2,'O'))        :  -3.59105,
    (('',3,'HN'),('',2,'N'))        :  -2.32535,
    (('',3,'HN'),('',2,'CA'))       :   1.82209,
    (('',3,'HN'),('',2,'HA'))       :  -2.80768,
    (('',3,'HN'),('',2,'C'))        :   0.00000,
    (('',3,'HN'),('',2,'O'))        :   0.00000,
    (('',3,'HN'),('',4,'N'))        :  -1.46861,
    (('',3,'HN'),('',4,'HN'))       :   0.00000,
    (('',3,'HN'),('',4,'CA'))       :   1.40793,
    (('',3,'HN'),('',4,'HA'))       :  -0.20708,
    (('',3,'HN'),('',4,'C'))        :  -0.19421,
    (('',3,'HN'),('',3,'N'))        :   0.00000,
    (('',3,'HN'),('',3,'HN'))       :   0.00000,
    (('',3,'HN'),('',3,'CA'))       :   0.00000,
    (('',3,'HN'),('',3,'HA'))       :   0.16901,
    (('',3,'HN'),('',3,'C'))        :   0.93848,
    (('',3,'HN'),('',3,'O'))        :  -1.82794,
    (('',2,'N'),('',1,'N'))         :  -2.73641,
    (('',2,'N'),('',1,'CA'))        :   0.00000,
    (('',2,'N'),('',1,'HA'))        :  -9.85348,
    (('',2,'N'),('',1,'C'))         :  15.50190,
    (('',2,'N'),('',1,'O'))         :  -4.78287,
    (('',2,'N'),('',3,'N'))         :   3.60263,
    (('',2,'N'),('',3,'HN'))        :  12.81990,
    (('',2,'N'),('',3,'CA'))        : -12.99700,
    (('',2,'N'),('',3,'HA'))        :   3.11644,
    (('',2,'N'),('',3,'C'))         :  -0.39693,
    (('',2,'N'),('',2,'N'))         :  -0.00000,
    (('',2,'N'),('',2,'HN'))        :  14.93430,
    (('',2,'N'),('',2,'CA'))        : -26.30710,
    (('',2,'N'),('',2,'HA'))        :  17.42810,
    (('',2,'N'),('',2,'C'))         :  29.50830,
    (('',2,'N'),('',2,'O'))         :  10.99260,
    (('',3,'N'),('',2,'N'))         :  -3.17937,
    (('',3,'N'),('',2,'CA'))        :   0.00000,
    (('',3,'N'),('',2,'HA'))        :  -8.41672,
    (('',3,'N'),('',2,'C'))         :  14.27860,
    (('',3,'N'),('',2,'O'))         :  -3.67960,
    (('',3,'N'),('',4,'N'))         :   3.28718,
    (('',3,'N'),('',4,'HN'))        :  12.87450,
    (('',3,'N'),('',4,'CA'))        : -15.64570,
    (('',3,'N'),('',4,'HA'))        :   3.57662,
    (('',3,'N'),('',4,'C'))         :   0.25904,
    (('',3,'N'),('',3,'N'))         :  -0.00000,
    (('',3,'N'),('',3,'HN'))        :  -3.46839,
    (('',3,'N'),('',3,'CA'))        : -23.19440,
    (('',3,'N'),('',3,'HA'))        :  -4.54475,
    (('',3,'N'),('',3,'C'))         :   4.38028,
    (('',3,'N'),('',3,'O'))         :  13.64970,
    (('',2,'C'),('',1,'N'))         :  -0.18869,
    (('',2,'C'),('',1,'CA'))        :  -4.09823,
    (('',2,'C'),('',1,'HA'))        :  -1.37988,
    (('',2,'C'),('',1,'C'))         :   0.00000,
    (('',2,'C'),('',1,'O'))         :  -0.17330,
    (('',2,'C'),('',3,'N'))         :   0.00000,
    (('',2,'C'),('',3,'HN'))        :   0.00000,
    (('',2,'C'),('',3,'CA'))        :  -2.95131,
    (('',2,'C'),('',3,'HA'))        :  -0.66979,
    (('',2,'C'),('',3,'C'))         :  -3.39004,
    (('',2,'C'),('',2,'N'))         :   0.00000,
    (('',2,'C'),('',2,'HN'))        :   1.23911,
    (('',2,'C'),('',2,'CA'))        :   0.00000,
    (('',2,'C'),('',2,'HA'))        :   0.00000,
    (('',2,'C'),('',2,'C'))         :   0.00000,
    (('',2,'C'),('',2,'O'))         :   0.00000,
    (('',3,'C'),('',2,'N'))         :  -0.18922,
    (('',3,'C'),('',2,'CA'))        :  -1.54869,
    (('',3,'C'),('',2,'HA'))        :  -1.65937,
    (('',3,'C'),('',2,'C'))         :   0.00000,
    (('',3,'C'),('',2,'O'))         :  -1.64398,
    (('',3,'C'),('',4,'N'))         :   0.00000,
    (('',3,'C'),('',4,'HN'))        :   0.00000,
    (('',3,'C'),('',4,'CA'))        :  -3.62975,
    (('',3,'C'),('',4,'HA'))        :  -0.50979,
    (('',3,'C'),('',4,'C'))         :  -3.11858,
    (('',3,'C'),('',3,'N'))         :   0.00000,
    (('',3,'C'),('',3,'HN'))        :  -0.90172,
    (('',3,'C'),('',3,'CA'))        :   0.00000,
    (('',3,'C'),('',3,'HA'))        :   0.00000,
    (('',3,'C'),('',3,'C'))         :  -0.00000,
    (('',3,'C'),('',3,'O'))         :   0.00000,
    (('',3,'CB'),('',2,'N'))        :  -1.50426,
    (('',3,'CB'),('',2,'CA'))       :  -0.69935,
    (('',3,'CB'),('',2,'HA'))       :   1.91689,
    (('',3,'CB'),('',2,'C'))        :  -3.28512,
    (('',3,'CB'),('',2,'O'))        :   2.35177,
    (('',3,'CB'),('',4,'N'))        :   3.16920,
    (('',3,'CB'),('',4,'HN'))       :   2.28277,
    (('',3,'CB'),('',4,'CA'))       : -10.15570,
    (('',3,'CB'),('',4,'HA'))       :   0.61830,
    (('',3,'CB'),('',4,'C'))        :   0.54843,
    (('',3,'CB'),('',3,'N'))        :   8.82469,
    (('',3,'CB'),('',3,'HN'))       :  -2.56107,
    (('',3,'CB'),('',3,'CA'))       :   0.00000,
    (('',3,'CB'),('',3,'HA'))       :   3.25154,
    (('',3,'CB'),('',3,'C'))        : -14.45600,
    (('',3,'CB'),('',3,'O'))        :   4.32263,
}
agfa_component_shifts_sc =  {
    (('',2,'HA'),('',2,'HA2'))      :   0.10556,
    (('',3,'HA'),('',3,'CB'))       :   0.00000,
    (('',3,'HA'),('',3,'CG'))       :   0.00000,
    (('',3,'HA'),('',3,'CD1'))      :   0.00000,
    (('',3,'HA'),('',3,'CD2'))      :   0.00000,
    (('',3,'HA'),('',3,'CE1'))      :   0.00000,
    (('',3,'HA'),('',3,'CE2'))      :   0.00000,
    (('',3,'HA'),('',3,'CZ'))       :   0.00000,
    (('',3,'HA'),('',3,'HB1'))      :   0.00000,
    (('',3,'HA'),('',3,'HB2'))      :   0.00000,
    (('',3,'HA'),('',3,'HD1'))      :   0.00000,
    (('',3,'HA'),('',3,'HD2'))      :   0.00000,
    (('',3,'HA'),('',3,'HE1'))      :   0.00000,
    (('',3,'HA'),('',3,'HE2'))      :   0.00000,
    (('',3,'HA'),('',3,'HZ'))       :   0.00000,
    (('',2,'CA'),('',2,'HA2'))      :  -3.81784,
    (('',3,'CA'),('',3,'CB'))       :   0.00000,
    (('',3,'CA'),('',3,'CG'))       :   0.00000,
    (('',3,'CA'),('',3,'CD1'))      :   0.00000,
    (('',3,'CA'),('',3,'CD2'))      :   0.00000,
    (('',3,'CA'),('',3,'CE1'))      :   0.00000,
    (('',3,'CA'),('',3,'CE2'))      :   0.00000,
    (('',3,'CA'),('',3,'CZ'))       :   0.00000,
    (('',3,'CA'),('',3,'HB1'))      :  -0.74505,
    (('',3,'CA'),('',3,'HB2'))      :  -0.74496,
    (('',3,'CA'),('',3,'HD1'))      :   0.00000,
    (('',3,'CA'),('',3,'HD2'))      :   0.00000,
    (('',3,'CA'),('',3,'HE1'))      :   0.00000,
    (('',3,'CA'),('',3,'HE2'))      :   0.00000,
    (('',3,'CA'),('',3,'HZ'))       :   0.00000,
    (('',2,'HN'),('',2,'HA2'))      :   0.07021,
    (('',3,'HN'),('',3,'CB'))       :  -0.03620,
    (('',3,'HN'),('',3,'CG'))       :   0.00000,
    (('',3,'HN'),('',3,'CD1'))      :   0.00000,
    (('',3,'HN'),('',3,'CD2'))      :   0.00000,
    (('',3,'HN'),('',3,'CE1'))      :   0.00000,
    (('',3,'HN'),('',3,'CE2'))      :   0.00000,
    (('',3,'HN'),('',3,'CZ'))       :   0.00000,
    (('',3,'HN'),('',3,'HB1'))      :  -0.12771,
    (('',3,'HN'),('',3,'HB2'))      :  -0.10758,
    (('',3,'HN'),('',3,'HD1'))      :   0.00000,
    (('',3,'HN'),('',3,'HD2'))      :   0.00000,
    (('',3,'HN'),('',3,'HE1'))      :   0.00000,
    (('',3,'HN'),('',3,'HE2'))      :   0.00000,
    (('',3,'HN'),('',3,'HZ'))       :   0.00000,
    (('',2,'N'),('',2,'HA2'))       :   2.87169,
    (('',3,'N'),('',3,'CB'))        : -17.47630,
    (('',3,'N'),('',3,'CG'))        :   0.00000,
    (('',3,'N'),('',3,'CD1'))       :   0.00000,
    (('',3,'N'),('',3,'CD2'))       :   0.00000,
    (('',3,'N'),('',3,'CE1'))       :   0.00000,
    (('',3,'N'),('',3,'CE2'))       :   0.00000,
    (('',3,'N'),('',3,'CZ'))        :   5.35450,
    (('',3,'N'),('',3,'HB1'))       :  -2.57200,
    (('',3,'N'),('',3,'HB2'))       :  -1.93905,
    (('',3,'N'),('',3,'HD1'))       :   0.00000,
    (('',3,'N'),('',3,'HD2'))       :   0.00000,
    (('',3,'N'),('',3,'HE1'))       :   0.00000,
    (('',3,'N'),('',3,'HE2'))       :   0.00000,
    (('',3,'N'),('',3,'HZ'))        :   0.00000,
    (('',2,'C'),('',2,'HA2'))       :   0.99617,
    (('',3,'C'),('',3,'CB'))        :   4.06446,
    (('',3,'C'),('',3,'CG'))        :   0.00000,
    (('',3,'C'),('',3,'CD1'))       :   0.00000,
    (('',3,'C'),('',3,'CD2'))       :   0.00000,
    (('',3,'C'),('',3,'CE1'))       :   0.00000,
    (('',3,'C'),('',3,'CE2'))       :   0.00000,
    (('',3,'C'),('',3,'CZ'))        :   0.00000,
    (('',3,'C'),('',3,'HB1'))       :  -2.01179,
    (('',3,'C'),('',3,'HB2'))       :  -2.62263,
    (('',3,'C'),('',3,'HD1'))       :   0.00000,
    (('',3,'C'),('',3,'HD2'))       :   0.00000,
    (('',3,'C'),('',3,'HE1'))       :   0.00000,
    (('',3,'C'),('',3,'HE2'))       :   0.00000,
    (('',3,'C'),('',3,'HZ'))        :   0.00000,
    (('',3,'CB'),('',3,'CB'))       :   0.00000,
    (('',3,'CB'),('',3,'CG'))       :   0.00000,
    (('',3,'CB'),('',3,'CD1'))      :   0.00000,
    (('',3,'CB'),('',3,'CD2'))      :   0.00000,
    (('',3,'CB'),('',3,'CE1'))      :   0.00000,
    (('',3,'CB'),('',3,'CE2'))      :   0.00000,
    (('',3,'CB'),('',3,'CZ'))       :   0.00000,
    (('',3,'CB'),('',3,'HB1'))      :   0.00000,
    (('',3,'CB'),('',3,'HB2'))      :   0.00000,
    (('',3,'CB'),('',3,'HD1'))      :   0.00000,
    (('',3,'CB'),('',3,'HD2'))      :   0.00000,
    (('',3,'CB'),('',3,'HE1'))      :   0.00000,
    (('',3,'CB'),('',3,'HE2'))      :   0.00000,
    (('',3,'CB'),('',3,'HZ'))       :   0.00000,
}
agfa_component_shifts_ring =  {
    (('',2,'HA'),(3,'PHE',6))    :   0.05200,
    (('',3,'HA'),(3,'PHE',6))    :   0.26450,
    (('',2,'CA'),(3,'PHE',6))    :   0.03350,
    (('',3,'CA'),(3,'PHE',6))    :   0.17460,
    (('',2,'HN'),(3,'PHE',6))    :   0.03380,
    (('',3,'HN'),(3,'PHE',6))    :   0.14540,
    (('',2,'N'),(3,'PHE',6))     :  -0.02590,
    (('',3,'N'),(3,'PHE',6))     :  -0.12310,
    (('',2,'C'),(3,'PHE',6))     :   0.28960,
    (('',3,'C'),(3,'PHE',6))     :   0.17460,
    (('',3,'CB'),(3,'PHE',6))    :   0.92890,
}
