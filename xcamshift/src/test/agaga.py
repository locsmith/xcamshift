from common_constants import BACK_BONE, SIDE_CHAIN, NON_BONDED, XTRA, DIHEDRAL, RANDOM_COIL

agaga_shifts  =   {
    ('',2,"HA")   :    4.078030,
    ('',3,"HA")   :    4.278500,
    ('',4,"HA")   :    3.777630,
    ('',2,"CA")   :   43.433800,
    ('',3,"CA")   :   50.457200,
    ('',4,"CA")   :   45.060900,
    ('',2,"HN")   :    8.191160,
    ('',3,"HN")   :    8.446220,
    ('',4,"HN")   :    7.898600,
    ('',2,"N")    :  104.702000,
    ('',3,"N")    :  121.445000,
    ('',4,"N")    :  103.538000,
    ('',2,"C")    :  173.575000,
    ('',3,"C")    :  175.985000,
    ('',4,"C")    :  174.564000,
    ('',3,"CB")   :   19.149800,
}
agaga_subpotential_shifts  =   {
    ('',2,'HA',RANDOM_COIL)    :   3.97000,
    ('',2,'HA',BACK_BONE)      :   1.00570,
    ('',2,'HA',SIDE_CHAIN)     :   0.10550,
    ('',2,'HA',NON_BONDED)     :  -0.01306,
    ('',2,'HA',XTRA)           :  -1.17221,
    ('',2,'HA',DIHEDRAL)       :   0.18208,
    ('',3,'HA',RANDOM_COIL)    :   4.43285,
    ('',3,'HA',BACK_BONE)      :   0.62145,
    ('',3,'HA',SIDE_CHAIN)     :   0.00000,
    ('',3,'HA',NON_BONDED)     :  -0.00476,
    ('',3,'HA',XTRA)           :  -0.89550,
    ('',3,'HA',DIHEDRAL)       :   0.12446,
    ('',4,'HA',RANDOM_COIL)    :   3.97000,
    ('',4,'HA',BACK_BONE)      :   0.77690,
    ('',4,'HA',SIDE_CHAIN)     :   0.10559,
    ('',4,'HA',NON_BONDED)     :   0.00000,
    ('',4,'HA',XTRA)           :  -1.07290,
    ('',4,'HA',DIHEDRAL)       :  -0.00195,
    ('',2,'CA',RANDOM_COIL)    :  45.48650,
    ('',2,'CA',BACK_BONE)      :  32.65910,
    ('',2,'CA',SIDE_CHAIN)     :  -3.81404,
    ('',2,'CA',NON_BONDED)     :   0.05281,
    ('',2,'CA',XTRA)           : -29.81260,
    ('',2,'CA',DIHEDRAL)       :  -1.13805,
    ('',3,'CA',RANDOM_COIL)    :  52.26640,
    ('',3,'CA',BACK_BONE)      :  42.38790,
    ('',3,'CA',SIDE_CHAIN)     :  -5.31042,
    ('',3,'CA',NON_BONDED)     :   0.00546,
    ('',3,'CA',XTRA)           : -37.69700,
    ('',3,'CA',DIHEDRAL)       :  -1.19519,
    ('',4,'CA',RANDOM_COIL)    :  45.48650,
    ('',4,'CA',BACK_BONE)      :  32.36670,
    ('',4,'CA',SIDE_CHAIN)     :  -3.81529,
    ('',4,'CA',NON_BONDED)     :   0.00000,
    ('',4,'CA',XTRA)           : -29.36660,
    ('',4,'CA',DIHEDRAL)       :   0.38953,
    ('',2,'HN',RANDOM_COIL)    :   8.33000,
    ('',2,'HN',BACK_BONE)      :  -7.17339,
    ('',2,'HN',SIDE_CHAIN)     :   0.07708,
    ('',2,'HN',NON_BONDED)     :   0.00000,
    ('',2,'HN',XTRA)           :   6.87995,
    ('',2,'HN',DIHEDRAL)       :   0.07752,
    ('',3,'HN',RANDOM_COIL)    :   8.24000,
    ('',3,'HN',BACK_BONE)      :  -4.64577,
    ('',3,'HN',SIDE_CHAIN)     :   0.00000,
    ('',3,'HN',NON_BONDED)     :  -0.19851,
    ('',3,'HN',XTRA)           :   5.05573,
    ('',3,'HN',DIHEDRAL)       :  -0.00523,
    ('',4,'HN',RANDOM_COIL)    :   8.33000,
    ('',4,'HN',BACK_BONE)      :  -7.42168,
    ('',4,'HN',SIDE_CHAIN)     :   0.08835,
    ('',4,'HN',NON_BONDED)     :  -0.18214,
    ('',4,'HN',XTRA)           :   7.11641,
    ('',4,'HN',DIHEDRAL)       :  -0.03233,
    ('',2,'N',RANDOM_COIL)     : 108.60000,
    ('',2,'N',BACK_BONE)       :  48.06570,
    ('',2,'N',SIDE_CHAIN)      :   2.87069,
    ('',2,'N',NON_BONDED)      :   0.00000,
    ('',2,'N',XTRA)            : -54.60950,
    ('',2,'N',DIHEDRAL)        :  -0.22485,
    ('',3,'N',RANDOM_COIL)     : 124.87000,
    ('',3,'N',BACK_BONE)       : -10.44140,
    ('',3,'N',SIDE_CHAIN)      :   0.00000,
    ('',3,'N',NON_BONDED)      :  -0.93914,
    ('',3,'N',XTRA)            :   8.29230,
    ('',3,'N',DIHEDRAL)        :  -0.33665,
    ('',4,'N',RANDOM_COIL)     : 108.60000,
    ('',4,'N',BACK_BONE)       :  50.15900,
    ('',4,'N',SIDE_CHAIN)      :   2.87194,
    ('',4,'N',NON_BONDED)      :  -0.86192,
    ('',4,'N',XTRA)            : -57.47540,
    ('',4,'N',DIHEDRAL)        :   0.24424,
    ('',2,'C',RANDOM_COIL)     : 173.60000,
    ('',2,'C',BACK_BONE)       : -12.57090,
    ('',2,'C',SIDE_CHAIN)      :   0.99587,
    ('',2,'C',NON_BONDED)      :   0.43548,
    ('',2,'C',XTRA)            :  11.51390,
    ('',2,'C',DIHEDRAL)        :  -0.39972,
    ('',3,'C',RANDOM_COIL)     : 177.10000,
    ('',3,'C',BACK_BONE)       : -12.89250,
    ('',3,'C',SIDE_CHAIN)      :   0.00000,
    ('',3,'C',NON_BONDED)      :   0.16855,
    ('',3,'C',XTRA)            :  11.20160,
    ('',3,'C',DIHEDRAL)        :   0.40769,
    ('',4,'C',RANDOM_COIL)     : 173.60000,
    ('',4,'C',BACK_BONE)       : -11.38630,
    ('',4,'C',SIDE_CHAIN)      :   0.99575,
    ('',4,'C',NON_BONDED)      :   0.00000,
    ('',4,'C',XTRA)            :  11.33910,
    ('',4,'C',DIHEDRAL)        :   0.01501,
    ('',3,'CB',RANDOM_COIL)    :  19.00000,
    ('',3,'CB',BACK_BONE)      :  -5.30961,
    ('',3,'CB',SIDE_CHAIN)     :   0.00000,
    ('',3,'CB',NON_BONDED)     :   0.00000,
    ('',3,'CB',XTRA)           :   5.15571,
    ('',3,'CB',DIHEDRAL)       :   0.30374,
}

agaga_component_shifts_bb =  {
    (('',2,'HA1'),('',1,'N'))  :   0.31119,
    (('',2,'HA1'),('',1,'CA'))  :   0.00000,
    (('',2,'HA1'),('',1,'HA'))  :   0.36887,
    (('',2,'HA1'),('',1,'C'))  :   0.00000,
    (('',2,'HA1'),('',1,'O'))  :  -0.70004,
    (('',2,'HA1'),('',3,'N'))  :   0.00000,
    (('',2,'HA1'),('',3,'HN'))  :   0.33533,
    (('',2,'HA1'),('',3,'CA'))  :   0.26311,
    (('',2,'HA1'),('',3,'HA'))  :   0.10469,
    (('',2,'HA1'),('',3,'C'))  :   0.49683,
    (('',2,'HA1'),('',2,'N'))  :   0.00000,
    (('',2,'HA1'),('',2,'HN'))  :  -0.48220,
    (('',2,'HA1'),('',2,'CA'))  :   0.00000,
    (('',2,'HA1'),('',2,'HA'))  :   0.00000,
    (('',2,'HA1'),('',2,'C'))  :   0.00000,
    (('',2,'HA1'),('',2,'O'))  :   0.30794,
    (('',3,'HA'),('',2,'N'))   :   0.47693,
    (('',3,'HA'),('',2,'CA'))  :   0.00000,
    (('',3,'HA'),('',2,'HA1'))  :   0.41033,
    (('',3,'HA'),('',2,'C'))   :   0.00000,
    (('',3,'HA'),('',2,'O'))   :  -1.18222,
    (('',3,'HA'),('',4,'N'))   :   0.00000,
    (('',3,'HA'),('',4,'HN'))  :   0.53124,
    (('',3,'HA'),('',4,'CA'))  :  -0.24300,
    (('',3,'HA'),('',4,'HA1'))  :   0.20777,
    (('',3,'HA'),('',4,'C'))   :   0.57754,
    (('',3,'HA'),('',3,'N'))   :   0.00000,
    (('',3,'HA'),('',3,'HN'))  :  -0.59474,
    (('',3,'HA'),('',3,'CA'))  :   0.00000,
    (('',3,'HA'),('',3,'HA'))  :  -0.00000,
    (('',3,'HA'),('',3,'C'))   :   0.00000,
    (('',3,'HA'),('',3,'O'))   :   0.43758,
    (('',4,'HA1'),('',3,'N'))  :   0.32155,
    (('',4,'HA1'),('',3,'CA'))  :   0.00000,
    (('',4,'HA1'),('',3,'HA'))  :   0.43105,
    (('',4,'HA1'),('',3,'C'))  :   0.00000,
    (('',4,'HA1'),('',3,'O'))  :  -1.13062,
    (('',4,'HA1'),('',5,'N'))  :   0.00000,
    (('',4,'HA1'),('',5,'HN'))  :   0.34572,
    (('',4,'HA1'),('',5,'CA'))  :   0.26515,
    (('',4,'HA1'),('',5,'HA'))  :   0.10479,
    (('',4,'HA1'),('',5,'C'))  :   0.51490,
    (('',4,'HA1'),('',4,'N'))  :   0.00000,
    (('',4,'HA1'),('',4,'HN'))  :  -0.38024,
    (('',4,'HA1'),('',4,'CA'))  :   0.00000,
    (('',4,'HA1'),('',4,'HA'))  :   0.00000,
    (('',4,'HA1'),('',4,'C'))  :   0.00000,
    (('',4,'HA1'),('',4,'O'))  :   0.30459,
    (('',2,'CA'),('',1,'N'))   :  -2.04981,
    (('',2,'CA'),('',1,'CA'))  :  -2.02096,
    (('',2,'CA'),('',1,'HA'))  :  -7.81181,
    (('',2,'CA'),('',1,'C'))   :   9.42477,
    (('',2,'CA'),('',1,'O'))   :  -2.92058,
    (('',2,'CA'),('',3,'N'))   :  -3.64166,
    (('',2,'CA'),('',3,'HN'))  :   2.02799,
    (('',2,'CA'),('',3,'CA'))  :   3.98572,
    (('',2,'CA'),('',3,'HA'))  :  -0.01286,
    (('',2,'CA'),('',3,'C'))   :   2.61190,
    (('',2,'CA'),('',2,'N'))   :   0.00000,
    (('',2,'CA'),('',2,'HN'))  :   5.97060,
    (('',2,'CA'),('',2,'CA'))  :   0.00000,
    (('',2,'CA'),('',2,'HA1'))  :   7.15899,
    (('',2,'CA'),('',2,'C'))   :  10.34460,
    (('',2,'CA'),('',2,'O'))   :   9.59220,
    (('',3,'CA'),('',2,'N'))   :  -2.39056,
    (('',3,'CA'),('',2,'CA'))  :  -1.20523,
    (('',3,'CA'),('',2,'HA1'))  :  -9.34390,
    (('',3,'CA'),('',2,'C'))   :  12.23680,
    (('',3,'CA'),('',2,'O'))   :  -2.62001,
    (('',3,'CA'),('',4,'N'))   :   2.15881,
    (('',3,'CA'),('',4,'HN'))  :   2.68350,
    (('',3,'CA'),('',4,'CA'))  :   3.20499,
    (('',3,'CA'),('',4,'HA1'))  :  -0.03067,
    (('',3,'CA'),('',4,'C'))   :   2.52293,
    (('',3,'CA'),('',3,'N'))   :   0.00000,
    (('',3,'CA'),('',3,'HN'))  :   8.31120,
    (('',3,'CA'),('',3,'CA'))  :  -0.00000,
    (('',3,'CA'),('',3,'HA'))  :   8.50898,
    (('',3,'CA'),('',3,'C'))   :   9.47732,
    (('',3,'CA'),('',3,'O'))   :   8.87380,
    (('',4,'CA'),('',3,'N'))   :  -2.07937,
    (('',4,'CA'),('',3,'CA'))  :  -2.02125,
    (('',4,'CA'),('',3,'HA'))  :  -7.81370,
    (('',4,'CA'),('',3,'C'))   :   9.42142,
    (('',4,'CA'),('',3,'O'))   :  -2.91992,
    (('',4,'CA'),('',5,'N'))   :  -3.64170,
    (('',4,'CA'),('',5,'HN'))  :   2.02816,
    (('',4,'CA'),('',5,'CA'))  :   3.98595,
    (('',4,'CA'),('',5,'HA'))  :  -0.01330,
    (('',4,'CA'),('',5,'C'))   :   2.35261,
    (('',4,'CA'),('',4,'N'))   :   0.00000,
    (('',4,'CA'),('',4,'HN'))  :   5.97017,
    (('',4,'CA'),('',4,'CA'))  :   0.00000,
    (('',4,'CA'),('',4,'HA1'))  :   7.16121,
    (('',4,'CA'),('',4,'C'))   :  10.34080,
    (('',4,'CA'),('',4,'O'))   :   9.59564,
    (('',2,'HN'),('',1,'N'))   :  -1.70860,
    (('',2,'HN'),('',1,'CA'))  :   1.77954,
    (('',2,'HN'),('',1,'HA'))  :  -3.77587,
    (('',2,'HN'),('',1,'C'))   :   0.00000,
    (('',2,'HN'),('',1,'O'))   :   0.00000,
    (('',2,'HN'),('',3,'N'))   :  -2.96674,
    (('',2,'HN'),('',3,'HN'))  :   0.00000,
    (('',2,'HN'),('',3,'CA'))  :   1.98586,
    (('',2,'HN'),('',3,'HA'))  :  -0.25263,
    (('',2,'HN'),('',3,'C'))   :  -0.41420,
    (('',2,'HN'),('',2,'N'))   :   0.00000,
    (('',2,'HN'),('',2,'HN'))  :   0.00000,
    (('',2,'HN'),('',2,'CA'))  :   0.00000,
    (('',2,'HN'),('',2,'HA1'))  :   0.04022,
    (('',2,'HN'),('',2,'C'))   :   1.39278,
    (('',2,'HN'),('',2,'O'))   :  -3.25376,
    (('',3,'HN'),('',2,'N'))   :  -2.03219,
    (('',3,'HN'),('',2,'CA'))  :   1.82219,
    (('',3,'HN'),('',2,'HA1'))  :  -3.13166,
    (('',3,'HN'),('',2,'C'))   :   0.00000,
    (('',3,'HN'),('',2,'O'))   :   0.00000,
    (('',3,'HN'),('',4,'N'))   :  -1.34160,
    (('',3,'HN'),('',4,'HN'))  :   0.00000,
    (('',3,'HN'),('',4,'CA'))  :   1.34143,
    (('',3,'HN'),('',4,'HA1'))  :  -0.23444,
    (('',3,'HN'),('',4,'C'))   :  -0.15887,
    (('',3,'HN'),('',3,'N'))   :   0.00000,
    (('',3,'HN'),('',3,'HN'))  :   0.00000,
    (('',3,'HN'),('',3,'CA'))  :   0.00000,
    (('',3,'HN'),('',3,'HA'))  :   0.16551,
    (('',3,'HN'),('',3,'C'))   :   0.95603,
    (('',3,'HN'),('',3,'O'))   :  -2.03215,
    (('',4,'HN'),('',3,'N'))   :  -1.79991,
    (('',4,'HN'),('',3,'CA'))  :   1.77996,
    (('',4,'HN'),('',3,'HA'))  :  -3.78731,
    (('',4,'HN'),('',3,'C'))   :   0.00000,
    (('',4,'HN'),('',3,'O'))   :   0.00000,
    (('',4,'HN'),('',5,'N'))   :  -3.85426,
    (('',4,'HN'),('',5,'HN'))  :   0.00000,
    (('',4,'HN'),('',5,'CA'))  :   2.41911,
    (('',4,'HN'),('',5,'HA'))  :  -0.35571,
    (('',4,'HN'),('',5,'C'))   :  -0.41220,
    (('',4,'HN'),('',4,'N'))   :   0.00000,
    (('',4,'HN'),('',4,'HN'))  :   0.00000,
    (('',4,'HN'),('',4,'CA'))  :   0.00000,
    (('',4,'HN'),('',4,'HA1'))  :   0.03172,
    (('',4,'HN'),('',4,'C'))   :   1.55154,
    (('',4,'HN'),('',4,'O'))   :  -2.99463,
    (('',2,'N'),('',1,'N'))    :  -2.09956,
    (('',2,'N'),('',1,'CA'))   :   0.00000,
    (('',2,'N'),('',1,'HA'))   : -12.49740,
    (('',2,'N'),('',1,'C'))    :  15.49900,
    (('',2,'N'),('',1,'O'))    :  -4.78339,
    (('',2,'N'),('',3,'N'))    :   3.33781,
    (('',2,'N'),('',3,'HN'))   :  11.20370,
    (('',2,'N'),('',3,'CA'))   : -12.47220,
    (('',2,'N'),('',3,'HA'))   :   3.01943,
    (('',2,'N'),('',3,'C'))    :  -0.37800,
    (('',2,'N'),('',2,'N'))    :  -0.00000,
    (('',2,'N'),('',2,'HN'))   :  14.94420,
    (('',2,'N'),('',2,'CA'))   : -26.32170,
    (('',2,'N'),('',2,'HA1'))  :  17.42420,
    (('',2,'N'),('',2,'C'))    :  29.51070,
    (('',2,'N'),('',2,'O'))    :  11.67890,
    (('',3,'N'),('',2,'N'))    :  -2.94567,
    (('',3,'N'),('',2,'CA'))   :   0.00000,
    (('',3,'N'),('',2,'HA1'))  :  -8.89112,
    (('',3,'N'),('',2,'C'))    :  14.28670,
    (('',3,'N'),('',2,'O'))    :  -3.67899,
    (('',3,'N'),('',4,'N'))    :   3.02648,
    (('',3,'N'),('',4,'HN'))   :  11.05000,
    (('',3,'N'),('',4,'CA'))   : -14.97930,
    (('',3,'N'),('',4,'HA1'))  :   3.80192,
    (('',3,'N'),('',4,'C'))    :   0.22971,
    (('',3,'N'),('',3,'N'))    :  -0.00000,
    (('',3,'N'),('',3,'HN'))   :  -3.47078,
    (('',3,'N'),('',3,'CA'))   : -23.19550,
    (('',3,'N'),('',3,'HA'))   :  -4.54382,
    (('',3,'N'),('',3,'C'))    :   4.38087,
    (('',3,'N'),('',3,'O'))    :  14.48810,
    (('',4,'N'),('',3,'N'))    :  -2.15418,
    (('',4,'N'),('',3,'CA'))   :   0.00000,
    (('',4,'N'),('',3,'HA'))   : -12.51230,
    (('',4,'N'),('',3,'C'))    :  15.49200,
    (('',4,'N'),('',3,'O'))    :  -4.78208,
    (('',4,'N'),('',5,'N'))    :   4.11336,
    (('',4,'N'),('',5,'HN'))   :  15.74000,
    (('',4,'N'),('',5,'CA'))   : -14.05750,
    (('',4,'N'),('',5,'HA'))   :   3.78105,
    (('',4,'N'),('',5,'C'))    :  -0.35501,
    (('',4,'N'),('',4,'N'))    :  -0.00000,
    (('',4,'N'),('',4,'HN'))   :  14.93850,
    (('',4,'N'),('',4,'CA'))   : -26.32480,
    (('',4,'N'),('',4,'HA1'))  :  17.42490,
    (('',4,'N'),('',4,'C'))    :  29.50670,
    (('',4,'N'),('',4,'O'))    :   9.34839,
    (('',2,'C'),('',1,'N'))    :  -0.19314,
    (('',2,'C'),('',1,'CA'))   :  -4.48533,
    (('',2,'C'),('',1,'HA'))   :  -1.69687,
    (('',2,'C'),('',1,'C'))    :   0.00000,
    (('',2,'C'),('',1,'O'))    :  -0.25226,
    (('',2,'C'),('',3,'N'))    :   0.00000,
    (('',2,'C'),('',3,'HN'))   :   0.00000,
    (('',2,'C'),('',3,'CA'))   :  -2.95181,
    (('',2,'C'),('',3,'HA'))   :  -0.68359,
    (('',2,'C'),('',3,'C'))    :  -3.34554,
    (('',2,'C'),('',2,'N'))    :   0.00000,
    (('',2,'C'),('',2,'HN'))   :   1.03762,
    (('',2,'C'),('',2,'CA'))   :   0.00000,
    (('',2,'C'),('',2,'HA1'))  :   0.00000,
    (('',2,'C'),('',2,'C'))    :   0.00000,
    (('',2,'C'),('',2,'O'))    :   0.00000,
    (('',3,'C'),('',2,'N'))    :  -0.18020,
    (('',3,'C'),('',2,'CA'))   :  -1.53808,
    (('',3,'C'),('',2,'HA1'))  :  -1.68038,
    (('',3,'C'),('',2,'C'))    :   0.00000,
    (('',3,'C'),('',2,'O'))    :  -1.61015,
    (('',3,'C'),('',4,'N'))    :   0.00000,
    (('',3,'C'),('',4,'HN'))   :   0.00000,
    (('',3,'C'),('',4,'CA'))   :  -3.62910,
    (('',3,'C'),('',4,'HA1'))  :  -0.63630,
    (('',3,'C'),('',4,'C'))    :  -2.69972,
    (('',3,'C'),('',3,'N'))    :   0.00000,
    (('',3,'C'),('',3,'HN'))   :  -0.91858,
    (('',3,'C'),('',3,'CA'))   :   0.00000,
    (('',3,'C'),('',3,'HA'))   :   0.00000,
    (('',3,'C'),('',3,'C'))    :  -0.00000,
    (('',3,'C'),('',3,'O'))    :   0.00000,
    (('',4,'C'),('',3,'N'))    :  -0.17436,
    (('',4,'C'),('',3,'CA'))   :  -4.27081,
    (('',4,'C'),('',3,'HA'))   :  -1.55933,
    (('',4,'C'),('',3,'C'))    :   0.00000,
    (('',4,'C'),('',3,'O'))    :  -0.21115,
    (('',4,'C'),('',5,'N'))    :   0.00000,
    (('',4,'C'),('',5,'HN'))   :   0.00000,
    (('',4,'C'),('',5,'CA'))   :  -2.95212,
    (('',4,'C'),('',5,'HA'))   :  -0.73696,
    (('',4,'C'),('',5,'C'))    :  -2.63750,
    (('',4,'C'),('',4,'N'))    :   0.00000,
    (('',4,'C'),('',4,'HN'))   :   1.15591,
    (('',4,'C'),('',4,'CA'))   :   0.00000,
    (('',4,'C'),('',4,'HA1'))  :   0.00000,
    (('',4,'C'),('',4,'C'))    :   0.00000,
    (('',4,'C'),('',4,'O'))    :   0.00000,
    (('',3,'CB'),('',2,'N'))   :  -1.48109,
    (('',3,'CB'),('',2,'CA'))  :  -0.69635,
    (('',3,'CB'),('',2,'HA1'))  :   1.91500,
    (('',3,'CB'),('',2,'C'))   :  -3.24862,
    (('',3,'CB'),('',2,'O'))   :   2.28496,
    (('',3,'CB'),('',4,'N'))   :   3.43347,
    (('',3,'CB'),('',4,'HN'))  :   2.68574,
    (('',3,'CB'),('',4,'CA'))  : -10.57580,
    (('',3,'CB'),('',4,'HA1'))  :   0.61386,
    (('',3,'CB'),('',4,'C'))   :   0.60115,
    (('',3,'CB'),('',3,'N'))   :   8.82270,
    (('',3,'CB'),('',3,'HN'))  :  -2.57489,
    (('',3,'CB'),('',3,'CA'))  :   0.00000,
    (('',3,'CB'),('',3,'HA'))  :   3.25050,
    (('',3,'CB'),('',3,'C'))   : -14.45640,
    (('',3,'CB'),('',3,'O'))   :   4.11618,
}
