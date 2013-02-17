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
from common_constants import RANDOM_COIL, BACK_BONE,SIDE_CHAIN,NON_BONDED, XTRA, DIHEDRAL
vin_shifts  =   {
    ('',2,"HA")   :    3.530820,
    ('',2,"CA")   :   61.569600,
    ('',2,"HN")   :    8.664200,
    ('',2,"N")    :  126.150000,
    ('',2,"C")    :  174.367000,
    ('',2,"CB")   :   35.959500,
}
vin_subpotential_shifts  =   {
    ('',2,'HA',RANDOM_COIL)    :   4.19333,
    ('',2,'HA',BACK_BONE)      :   0.45119,
    ('',2,'HA',SIDE_CHAIN)     :  -0.07469,
    ('',2,'HA',NON_BONDED)     :   0.00000,
    ('',2,'HA',XTRA)           :  -0.99118,
    ('',2,'HA',DIHEDRAL)       :  -0.04784,
    ('',2,'CA',RANDOM_COIL)    :  61.20620,
    ('',2,'CA',BACK_BONE)      :  42.37060,
    ('',2,'CA',SIDE_CHAIN)     :  -1.37429,
    ('',2,'CA',NON_BONDED)     :   0.00000,
    ('',2,'CA',XTRA)           : -40.87570,
    ('',2,'CA',DIHEDRAL)       :   0.24273,
    ('',2,'HN',RANDOM_COIL)    :   8.05000,
    ('',2,'HN',BACK_BONE)      :  -4.80116,
    ('',2,'HN',SIDE_CHAIN)     :  -0.00656,
    ('',2,'HN',NON_BONDED)     :   0.00000,
    ('',2,'HN',XTRA)           :   5.34235,
    ('',2,'HN',DIHEDRAL)       :   0.07957,
    ('',2,'N',RANDOM_COIL)     : 123.73000,
    ('',2,'N',BACK_BONE)       :  -9.28369,
    ('',2,'N',SIDE_CHAIN)      : -16.75470,
    ('',2,'N',NON_BONDED)      :   0.00000,
    ('',2,'N',XTRA)            :  27.59400,
    ('',2,'N',DIHEDRAL)        :   0.86422,
    ('',2,'C',RANDOM_COIL)     : 176.80000,
    ('',2,'C',BACK_BONE)       : -12.22480,
    ('',2,'C',SIDE_CHAIN)      :  -1.23565,
    ('',2,'C',NON_BONDED)      :   0.00000,
    ('',2,'C',XTRA)            :  10.93780,
    ('',2,'C',DIHEDRAL)        :   0.08967,
    ('',2,'CB',RANDOM_COIL)    :  37.50000,
    ('',2,'CB',BACK_BONE)      :  -5.35914,
    ('',2,'CB',SIDE_CHAIN)     :  -0.04282,
    ('',2,'CB',NON_BONDED)     :   0.00000,
    ('',2,'CB',XTRA)           :   4.19188,
    ('',2,'CB',DIHEDRAL)       :  -0.33047,
}
vin_component_shifts_bb =  {
    (('',2,'HA1'),('',1,'N'))  :   0.55894,
    (('',2,'HA1'),('',1,'CA'))  :   0.00000,
    (('',2,'HA1'),('',1,'HA'))  :   0.38139,
    (('',2,'HA1'),('',1,'C'))  :   0.00000,
    (('',2,'HA1'),('',1,'O'))  :  -1.32886,
    (('',2,'HA1'),('',3,'N'))  :   0.00000,
    (('',2,'HA1'),('',3,'HN'))  :   0.32594,
    (('',2,'HA1'),('',3,'CA'))  :  -0.21019,
    (('',2,'HA1'),('',3,'HA'))  :   0.17353,
    (('',2,'HA1'),('',3,'C'))  :   0.52655,
    (('',2,'HA1'),('',2,'N'))  :   0.00000,
    (('',2,'HA1'),('',2,'HN'))  :  -0.54141,
    (('',2,'HA1'),('',2,'CA'))  :   0.00000,
    (('',2,'HA1'),('',2,'HA'))  :  -0.00000,
    (('',2,'HA1'),('',2,'C'))  :   0.00000,
    (('',2,'HA1'),('',2,'O'))  :   0.56530,
    (('',2,'CA'),('',1,'N'))   :  -2.51890,
    (('',2,'CA'),('',1,'CA'))  :  -1.20519,
    (('',2,'CA'),('',1,'HA'))  :  -9.07210,
    (('',2,'CA'),('',1,'C'))   :  12.23200,
    (('',2,'CA'),('',1,'O'))   :  -2.61911,
    (('',2,'CA'),('',3,'N'))   :   2.15853,
    (('',2,'CA'),('',3,'HN'))  :   2.68320,
    (('',2,'CA'),('',3,'CA'))  :   3.20547,
    (('',2,'CA'),('',3,'HA'))  :  -0.03019,
    (('',2,'CA'),('',3,'C'))   :   2.35717,
    (('',2,'CA'),('',2,'N'))   :   0.00000,
    (('',2,'CA'),('',2,'HN'))  :   8.30688,
    (('',2,'CA'),('',2,'CA'))  :  -0.00000,
    (('',2,'CA'),('',2,'HA1'))  :   8.52168,
    (('',2,'CA'),('',2,'C'))   :   9.47409,
    (('',2,'CA'),('',2,'O'))   :   8.87714,
    (('',2,'HN'),('',1,'N'))   :  -2.41978,
    (('',2,'HN'),('',1,'CA'))  :   1.82247,
    (('',2,'HN'),('',1,'HA'))  :  -2.74229,
    (('',2,'HN'),('',1,'C'))   :   0.00000,
    (('',2,'HN'),('',1,'O'))   :   0.00000,
    (('',2,'HN'),('',3,'N'))   :  -1.55503,
    (('',2,'HN'),('',3,'HN'))  :   0.00000,
    (('',2,'HN'),('',3,'CA'))  :   1.52055,
    (('',2,'HN'),('',3,'HA'))  :  -0.23847,
    (('',2,'HN'),('',3,'C'))   :  -0.19474,
    (('',2,'HN'),('',2,'N'))   :   0.00000,
    (('',2,'HN'),('',2,'HN'))  :   0.00000,
    (('',2,'HN'),('',2,'CA'))  :   0.00000,
    (('',2,'HN'),('',2,'HA1'))  :   0.15067,
    (('',2,'HN'),('',2,'C'))   :   1.07973,
    (('',2,'HN'),('',2,'O'))   :  -2.22427,
    (('',2,'N'),('',1,'N'))    :  -3.25342,
    (('',2,'N'),('',1,'CA'))   :   0.00000,
    (('',2,'N'),('',1,'HA'))   :  -8.33440,
    (('',2,'N'),('',1,'C'))    :  14.28050,
    (('',2,'N'),('',1,'O'))    :  -3.67888,
    (('',2,'N'),('',3,'N'))    :   3.45606,
    (('',2,'N'),('',3,'HN'))   :  13.99460,
    (('',2,'N'),('',3,'CA'))   : -16.09120,
    (('',2,'N'),('',3,'HA'))   :   3.88272,
    (('',2,'N'),('',3,'C'))    :   0.23974,
    (('',2,'N'),('',2,'N'))    :  -0.00000,
    (('',2,'N'),('',2,'HN'))   :  -3.46956,
    (('',2,'N'),('',2,'CA'))   : -23.17940,
    (('',2,'N'),('',2,'HA1'))  :  -4.54552,
    (('',2,'N'),('',2,'C'))    :   4.37924,
    (('',2,'N'),('',2,'O'))    :  13.03580,
    (('',2,'C'),('',1,'N'))    :  -0.16194,
    (('',2,'C'),('',1,'CA'))   :  -1.45642,
    (('',2,'C'),('',1,'HA'))   :  -1.71409,
    (('',2,'C'),('',1,'C'))    :   0.00000,
    (('',2,'C'),('',1,'O'))    :  -1.32308,
    (('',2,'C'),('',3,'N'))    :   0.00000,
    (('',2,'C'),('',3,'HN'))   :   0.00000,
    (('',2,'C'),('',3,'CA'))   :  -3.63116,
    (('',2,'C'),('',3,'HA'))   :  -0.61594,
    (('',2,'C'),('',3,'C'))    :  -2.28470,
    (('',2,'C'),('',2,'N'))    :   0.00000,
    (('',2,'C'),('',2,'HN'))   :  -1.03744,
    (('',2,'C'),('',2,'CA'))   :   0.00000,
    (('',2,'C'),('',2,'HA1'))  :   0.00000,
    (('',2,'C'),('',2,'C'))    :  -0.00000,
    (('',2,'C'),('',2,'O'))    :   0.00000,
    (('',2,'CB'),('',1,'N'))   :  -1.54240,
    (('',2,'CB'),('',1,'CA'))  :  -0.70697,
    (('',2,'CB'),('',1,'HA'))  :   1.91644,
    (('',2,'CB'),('',1,'C'))   :  -3.38257,
    (('',2,'CB'),('',1,'O'))   :   2.52904,
    (('',2,'CB'),('',3,'N'))   :   4.10823,
    (('',2,'CB'),('',3,'HN'))  :   3.61312,
    (('',2,'CB'),('',3,'CA'))  : -11.71740,
    (('',2,'CB'),('',3,'HA'))  :   0.72442,
    (('',2,'CB'),('',3,'C'))   :   0.55343,
    (('',2,'CB'),('',2,'N'))   :   8.82175,
    (('',2,'CB'),('',2,'HN'))  :  -2.52134,
    (('',2,'CB'),('',2,'CA'))  :   0.00000,
    (('',2,'CB'),('',2,'HA1'))  :   3.25134,
    (('',2,'CB'),('',2,'C'))   : -14.45290,
    (('',2,'CB'),('',2,'O'))   :   3.44668,
}
