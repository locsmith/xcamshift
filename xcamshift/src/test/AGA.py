from common_constants import BACK_BONE, SIDE_CHAIN, NON_BONDED, XTRA, DIHEDRAL, RANDOM_COIL


aga_shifts  =   {
    ('',2,"CA")   :   43.653700,
    ('',2,"HA")   :    4.187840,
    ('',2,"N")    :  112.700000,
    ('',2,"HN")   :    8.421560,
    ('',2,"C")    :  173.555000,

}

aga_subpotential_shifts = {
    ('',2,'CA',RANDOM_COIL)    :  45.48650,
    ('',2,'CA',BACK_BONE)      :  33.34750,
    ('',2,'CA',SIDE_CHAIN)     :  -3.81418,
    ('',2,'CA',NON_BONDED)     :   0.00000,
    ('',2,'CA',XTRA)           : -30.86970,
    ('',2,'CA',DIHEDRAL)       :  -0.49650,
    ('',2,'HA',RANDOM_COIL)    :   3.97000,
    ('',2,'HA',BACK_BONE)      :   0.94405,
    ('',2,'HA',SIDE_CHAIN)     :   0.10554,
    ('',2,'HA',NON_BONDED)     :   0.00000,
    ('',2,'HA',XTRA)           :  -1.03976,
    ('',2,'HA',DIHEDRAL)       :   0.20800,
    ('',2,'N',RANDOM_COIL)     : 108.60000,
    ('',2,'N',BACK_BONE)       :  52.33050,
    ('',2,'N',SIDE_CHAIN)      :   2.87098,
    ('',2,'N',NON_BONDED)      :   0.00000,
    ('',2,'N',XTRA)            : -51.63270,
    ('',2,'N',DIHEDRAL)        :   0.53134,
    ('',2,'HN',RANDOM_COIL)    :   8.33000,
    ('',2,'HN',BACK_BONE)      :  -6.22866,
    ('',2,'HN',SIDE_CHAIN)     :   0.07337,
    ('',2,'HN',NON_BONDED)     :   0.00000,
    ('',2,'HN',XTRA)           :   6.11590,
    ('',2,'HN',DIHEDRAL)       :   0.13095,
    ('',2,'C',RANDOM_COIL)     : 173.60000,
    ('',2,'C',BACK_BONE)       : -11.86160,
    ('',2,'C',SIDE_CHAIN)      :   0.99559,
    ('',2,'C',NON_BONDED)      :   0.00000,
    ('',2,'C',XTRA)            :  11.14800,
    ('',2,'C',DIHEDRAL)        :  -0.32696,
}

aga_component_shifts_bb =  {
    (('',2,'HA1'),('',1,'N'))  :   0.31777,
    (('',2,'HA1'),('',1,'CA'))  :   0.00000,
    (('',2,'HA1'),('',1,'HA'))  :   0.35426,
    (('',2,'HA1'),('',1,'C'))  :   0.00000,
    (('',2,'HA1'),('',1,'O'))  :  -0.70639,
    (('',2,'HA1'),('',3,'N'))  :   0.00000,
    (('',2,'HA1'),('',3,'HN'))  :   0.28724,
    (('',2,'HA1'),('',3,'CA'))  :   0.25412,
    (('',2,'HA1'),('',3,'HA'))  :   0.09700,
    (('',2,'HA1'),('',3,'C'))  :   0.49932,
    (('',2,'HA1'),('',2,'N'))  :   0.00000,
    (('',2,'HA1'),('',2,'HN'))  :  -0.48105,
    (('',2,'HA1'),('',2,'CA'))  :   0.00000,
    (('',2,'HA1'),('',2,'HA'))  :   0.00000,
    (('',2,'HA1'),('',2,'C'))  :   0.00000,
    (('',2,'HA1'),('',2,'O'))  :   0.32180,
    (('',2,'CA'),('',1,'N'))   :  -2.15920,
    (('',2,'CA'),('',1,'CA'))  :  -2.02084,
    (('',2,'CA'),('',1,'HA'))  :  -6.87327,
    (('',2,'CA'),('',1,'C'))   :   9.42452,
    (('',2,'CA'),('',1,'O'))   :  -2.92020,
    (('',2,'CA'),('',3,'N'))   :  -3.64148,
    (('',2,'CA'),('',3,'HN'))  :   2.02805,
    (('',2,'CA'),('',3,'CA'))  :   3.98611,
    (('',2,'CA'),('',3,'HA'))  :  -0.01331,
    (('',2,'CA'),('',3,'C'))   :   2.47396,
    (('',2,'CA'),('',2,'N'))   :   0.00000,
    (('',2,'CA'),('',2,'HN'))  :   5.96857,
    (('',2,'CA'),('',2,'CA'))  :   0.00000,
    (('',2,'CA'),('',2,'HA1'))  :   7.16081,
    (('',2,'CA'),('',2,'C'))   :  10.33980,
    (('',2,'CA'),('',2,'O'))   :   9.59391,
    (('',2,'HN'),('',1,'N'))   :  -2.06826,
    (('',2,'HN'),('',1,'CA'))  :   1.77909,
    (('',2,'HN'),('',1,'HA'))  :  -2.49955,
    (('',2,'HN'),('',1,'C'))   :   0.00000,
    (('',2,'HN'),('',1,'O'))   :   0.00000,
    (('',2,'HN'),('',3,'N'))   :  -3.72381,
    (('',2,'HN'),('',3,'HN'))  :   0.00000,
    (('',2,'HN'),('',3,'CA'))  :   2.33299,
    (('',2,'HN'),('',3,'HA'))  :  -0.33928,
    (('',2,'HN'),('',3,'C'))   :  -0.43901,
    (('',2,'HN'),('',2,'N'))   :   0.00000,
    (('',2,'HN'),('',2,'HN'))  :   0.00000,
    (('',2,'HN'),('',2,'CA'))  :   0.00000,
    (('',2,'HN'),('',2,'HA1'))  :   0.04013,
    (('',2,'HN'),('',2,'C'))   :   1.47328,
    (('',2,'HN'),('',2,'O'))   :  -2.78423,
    (('',2,'N'),('',1,'N'))    :  -2.32136,
    (('',2,'N'),('',1,'CA'))   :   0.00000,
    (('',2,'N'),('',1,'HA'))   :  -9.77234,
    (('',2,'N'),('',1,'C'))    :  15.50180,
    (('',2,'N'),('',1,'O'))    :  -4.78167,
    (('',2,'N'),('',3,'N'))    :   3.92362,
    (('',2,'N'),('',3,'HN'))   :  14.68060,
    (('',2,'N'),('',3,'CA'))   : -13.65860,
    (('',2,'N'),('',3,'HA'))   :   3.59346,
    (('',2,'N'),('',3,'C'))    :  -0.39677,
    (('',2,'N'),('',2,'N'))    :  -0.00000,
    (('',2,'N'),('',2,'HN'))   :  14.93850,
    (('',2,'N'),('',2,'CA'))   : -26.31550,
    (('',2,'N'),('',2,'HA1'))  :  17.42100,
    (('',2,'N'),('',2,'C'))    :  29.50300,
    (('',2,'N'),('',2,'O'))    :  10.01460,
    (('',2,'C'),('',1,'N'))    :  -0.17996,
    (('',2,'C'),('',1,'CA'))   :  -4.37931,
    (('',2,'C'),('',1,'HA'))   :  -1.49720,
    (('',2,'C'),('',1,'C'))    :   0.00000,
    (('',2,'C'),('',1,'O'))    :  -0.23275,
    (('',2,'C'),('',3,'N'))    :   0.00000,
    (('',2,'C'),('',3,'HN'))   :   0.00000,
    (('',2,'C'),('',3,'CA'))   :  -2.95250,
    (('',2,'C'),('',3,'HA'))   :  -0.73699,
    (('',2,'C'),('',3,'C'))    :  -2.98048,
    (('',2,'C'),('',2,'N'))    :   0.00000,
    (('',2,'C'),('',2,'HN'))   :   1.09760,
    (('',2,'C'),('',2,'CA'))   :   0.00000,
    (('',2,'C'),('',2,'HA1'))  :   0.00000,
    (('',2,'C'),('',2,'C'))    :   0.00000,
    (('',2,'C'),('',2,'O'))    :   0.00000,
}
