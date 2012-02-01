'''
Created on 22 Jan 2012

@author: garyt
'''


ala_3_total_shifts =  {
    1 :  {  "HA" :  0.0000,  "CA" :  0.0000, "HN" : 0.0000,  "N" :   0.0000,  "C" :   0.0000,   "CB" :  0.0000},
    2 :  {  "HA" :  4.2651,  "CA" : 52.6395, "HN" : 8.2373,  "N" : 120.2627,  "C" : 177.7649,   "CB" : 18.3505},
    1 :  {  "HA" :  0.0000,  "CA" :  0.0000, "HN" : 0.0000,  "N" :   0.0000,  "C" :   0.0000,   "CB" :  0.0000},
}
 
 
ala_3_test_shifts_harmonic = {
    (2, "CA") : 54.6395,
    (2, "CB") : 15.3505,
    (2, "C") : 179.7649,
    (2, "HA") :  5.2651,
    (2, "HN") :  9.2373,
    (2, "N") : 122.2627
}

ala_3_test_shifts_well = {
#      N        HN        CA        HA        CB         C
#    120.2627    8.2373   52.6395    4.2651   18.3505  177.7649
    (2, "CA") : 52.6395,
    (2, "CB") : 18.3505,
    (2, "C") : 177.7649,
    (2, "HA") :  4.2651,
    (2, "HN") :  8.2373,
    (2, "N") : 120.2627
}

ala_3_energies = {
    (2, "CA") :  0.425554 ,
    (2, "CB") :  1.68482,
    (2, "C")  :  0.583972,
    (2, "HA") :  3.59174,
    (2, "HN") :  1.65455,
    (2, "N")  :  0.0403878,
    
    'total' : 7.98102
}

ala_3_factors_harmonic =  {
     (2, "CA") :  0.552835,
     (2, "HA") :  8.064630,
     (2, "N")  :  0.096157,
     (2, "HN")  :  4.149330,
     (2, "C")  :  0.783775,
     (2, "CB") : -1.359170,
}

ala_3_distance_forces_harmonic = {
    ((2,"CA"),(1,"N"))   : -0.0746,
    ((2,"CA"),(1,"CA"))  : -0.0461,
    ((2,"CA"),(1,"HA"))  : -0.2823,
    ((2,"CA"),(1,"C"))   :  1.1473,
    ((2,"CA"),(1,"O"))   : -0.1899,
    ((2,"CA"),(3,"N"))   :  0.2029,
    ((2,"CA"),(3,"HN"))  :  0.2304,
    ((2,"CA"),(3,"CA"))  :  0.1226,
    ((2,"CA"),(3,"HA"))  : -0.0009,
    ((2,"CA"),(3,"C"))   :  0.0688,
    ((2,"CA"),(2,"N"))   :  0.0000,
    ((2,"CA"),(2,"HN"))  :  1.0222,
    ((2,"CA"),(2,"HA"))  :  4.0408,
    ((2,"CA"),(2,"C"))   :  2.2514,
    ((2,"CA"),(2,"O"))   :  0.8513,
    ((2,"HA"),(1,"N"))   :  0.1812,
    ((2,"HA"),(1,"CA"))  :  0.0000,
    ((2,"HA"),(1,"HA"))  :  0.1497,
    ((2,"HA"),(1,"C"))   :  0.0000,
    ((2,"HA"),(1,"O"))   : -1.1673,
    ((2,"HA"),(3,"N"))   :  0.0000,
    ((2,"HA"),(3,"HN"))  :  0.3368,
    ((2,"HA"),(3,"CA"))  : -0.0944,
    ((2,"HA"),(3,"HA"))  :  0.0666,
    ((2,"HA"),(3,"C"))   :  0.1780,
    ((2,"HA"),(2,"N"))   :  0.0000,
    ((2,"HA"),(2,"HN"))  : -0.6427,
    ((2,"HA"),(2,"CA"))  :  0.0000,
    ((2,"HA"),(2,"C"))   :  0.0000,
    ((2,"HA"),(2,"O"))   :  0.5590,
    ((2,"N"),(1,"N"))    : -0.0359,
    ((2,"N"),(1,"CA"))   :  0.0000,
    ((2,"N"),(1,"HA"))   : -0.0992,
    ((2,"N"),(1,"C"))    :  0.7777,
    ((2,"N"),(1,"O"))    : -0.0699,
    ((2,"N"),(3,"N"))    :  0.0358,
    ((2,"N"),(3,"HN"))   :  0.1578,
    ((2,"N"),(3,"CA"))   : -0.0801,
    ((2,"N"),(3,"HA"))   :  0.0157,
    ((2,"N"),(3,"C"))    :  0.0010,
    ((2,"N"),(2,"HN"))   : -0.3471,
    ((2,"N"),(2,"CA"))   : -1.0486,
    ((2,"N"),(2,"HA"))   : -0.1013,
    ((2,"N"),(2,"C"))    :  0.0700,
    ((2,"N"),(2,"O"))    :  0.1170,
    ((2,"HN"),(1,"N"))   : -1.3537,
    ((2,"HN"),(1,"CA"))  :  1.1754,
    ((2,"HN"),(1,"HA"))  : -1.5003,
    ((2,"HN"),(1,"C"))   :  0.0000,
    ((2,"HN"),(1,"O"))   :  0.0000,
    ((2,"HN"),(3,"N"))   : -0.6127,
    ((2,"HN"),(3,"HN"))  :  0.0000,
    ((2,"HN"),(3,"CA"))  :  0.3007,
    ((2,"HN"),(3,"HA"))  : -0.0392,
    ((2,"HN"),(3,"C"))   : -0.0310,
    ((2,"HN"),(2,"N"))   :  0.0000,
    ((2,"HN"),(2,"CA"))  :  0.0000,
    ((2,"HN"),(2,"HA"))  :  0.0920,
    ((2,"HN"),(2,"C"))   :  0.4950,
    ((2,"HN"),(2,"O"))   : -0.5832,
    ((2,"C"),(1,"N"))    : -0.0056,
    ((2,"C"),(1,"CA"))   : -0.0557,
    ((2,"C"),(1,"HA"))   : -0.0559,
    ((2,"C"),(1,"C"))    :  0.0000,
    ((2,"C"),(1,"O"))    : -0.0921,
    ((2,"C"),(3,"N"))    :  0.0000,
    ((2,"C"),(3,"HN"))   :  0.0000,
    ((2,"C"),(3,"CA"))   : -0.4826,
    ((2,"C"),(3,"HA"))   : -0.0583,
    ((2,"C"),(3,"C"))    : -0.2111,
    ((2,"C"),(2,"N"))    :  0.0000,
    ((2,"C"),(2,"HN"))   : -0.0898,
    ((2,"C"),(2,"CA"))   :  0.0000,
    ((2,"C"),(2,"HA"))   :  0.0000,
    ((2,"C"),(2,"O"))    :  0.0000,
    ((2,"CB"),(1,"N"))   :  0.0836,
    ((2,"CB"),(1,"CA"))  :  0.0467,
    ((2,"CB"),(1,"HA"))  : -0.0992,
    ((2,"CB"),(1,"C"))   :  0.4487,
    ((2,"CB"),(1,"O"))   : -0.3236,
    ((2,"CB"),(3,"N"))   : -0.5060,
    ((2,"CB"),(3,"HN"))  : -0.4337,
    ((2,"CB"),(3,"CA"))  :  0.7505,
    ((2,"CB"),(3,"HA"))  : -0.0385,
    ((2,"CB"),(3,"C"))   : -0.0277,
    ((2,"CB"),(2,"N"))   : -2.0020,
    ((2,"CB"),(2,"HN"))  :  0.4016,
    ((2,"CB"),(2,"CA"))  : -0.0000,
    ((2,"CB"),(2,"HA"))  : -0.9666,
    ((2,"CB"),(2,"C"))   :  3.1539,
    ((2,"CB"),(2,"O"))   : -0.4933
}

ala_3_distance_real_forces_harmonic = {
    ((2,"CA"),(1,"N"))   :   (-0.0699008,-0.2829600,-0.0893715),
    ((2,"CA"),(1,"CA"))  :   (-0.0130589,-0.1393100,-0.1056710),
    ((2,"CA"),(1,"HA"))  :   (-0.2551660,-0.8614670,-0.8959030),
    ((2,"CA"),(1,"C"))   :   ( 0.1239080, 1.8012600, 2.1213600),
    ((2,"CA"),(1,"O"))   :   ( 0.0127220,-0.1283600,-0.5083120),
    ((2,"CA"),(3,"N"))   :   ( 0.4735070, 0.0024345,-0.1336940),
    ((2,"CA"),(3,"HN"))  :   ( 0.4966180, 0.2129350,-0.2233050),
    ((2,"CA"),(3,"CA"))  :   ( 0.4503610,-0.0672108,-0.0997123),
    ((2,"CA"),(3,"HA"))  :   (-0.0033258, 0.0014571, 0.0010595),
    ((2,"CA"),(3,"C"))   :   ( 0.3042620,-0.0353585, 0.0353585),
    ((2,"CA"),(2,"N"))   :   ( 0.0000000, 0.0000000, 0.0000000),
    ((2,"CA"),(2,"HN"))  :   ( 0.3036060, 2.1446600,-0.0746237),
    ((2,"CA"),(2,"HA"))  :   (-2.5739700,-2.3072800, 2.6588300),
    ((2,"CA"),(2,"C"))   :   ( 3.0483600,-1.5624500,-0.2386450),
    ((2,"CA"),(2,"O"))   :   ( 1.2812100,-1.5706500, 0.2604980),
    ((2,"HA"),(1,"N"))   :   ( 0.2852490, 0.7908680, 0.0978618),
    ((2,"HA"),(1,"CA"))  :   ( 0.0000000, 0.0000000, 0.0000000),
    ((2,"HA"),(1,"HA"))  :   ( 0.2307500, 0.5425100, 0.3767470),
    ((2,"HA"),(1,"C"))   :   ( 0.0000000, 0.0000000, 0.0000000),
    ((2,"HA"),(1,"O"))   :   (-0.6653750,-1.4556500,-2.3568300),
    ((2,"HA"),(3,"N"))   :   ( 0.0000000, 0.0000000,-0.0000000),
    ((2,"HA"),(3,"HN"))  :   ( 0.9404370, 0.5035650,-0.5480270),
    ((2,"HA"),(3,"CA"))  :   (-0.4067640,-0.0021712, 0.1388610),
    ((2,"HA"),(3,"HA"))  :   ( 0.2813920,-0.0666683,-0.1199500),
    ((2,"HA"),(3,"C"))   :   ( 0.9007720, 0.0101470,-0.0256346),
    ((2,"HA"),(2,"N"))   :   ( 0.0000000, 0.0000000,-0.0000000),
    ((2,"HA"),(2,"HN"))  :   (-0.6002700,-1.7153300, 0.4698040),
    ((2,"HA"),(2,"CA"))  :   ( 0.0000000, 0.0000000,-0.0000000),
    ((2,"HA"),(2,"C"))   :   ( 0.0000000,-0.0000000,-0.0000000),
    ((2,"HA"),(2,"O"))   :   ( 1.1974600,-0.7122150,-0.1967810),
    ((2,"N"),(1,"N"))    :   (-0.0280608,-0.0878067,-0.0236113),
    ((2,"N"),(1,"CA"))   :   ( 0.0000000, 0.0000000, 0.0000000),
    ((2,"N"),(1,"HA"))   :   (-0.0743065,-0.1692480,-0.2613130),
    ((2,"N"),(1,"C"))    :   (-0.0365531, 0.1742110, 1.0180400),
    ((2,"N"),(1,"O"))    :   ( 0.0155099, 0.0468090,-0.1493000),
    ((2,"N"),(3,"N"))    :   ( 0.0779743,-0.0477365,-0.0429056),
    ((2,"N"),(3,"HN"))   :   ( 0.3156090,-0.0665935,-0.2381270),
    ((2,"N"),(3,"CA"))   :   (-0.2818630, 0.1517910, 0.1084330),
    ((2,"N"),(3,"HA"))   :   ( 0.0538178,-0.0457443,-0.0263837),
    ((2,"N"),(3,"C"))    :   ( 0.0043430,-0.0018927,-0.0000265),
    ((2,"N"),(2,"HN"))   :   (-0.0492836,-0.2609950, 0.2127530),
    ((2,"N"),(2,"CA"))   :   ( 0.1625320, 1.4114000, 0.5662390),
    ((2,"N"),(2,"HA"))   :   ( 0.0802140, 0.1941540,-0.0119511),
    ((2,"N"),(2,"C"))    :   ( 0.0839461,-0.1428270,-0.0452287),
    ((2,"N"),(2,"O"))    :   ( 0.1579780,-0.3734140,-0.0273829),
    ((2,"HN"),(1,"N"))   :   (-0.8663370,-2.2944400,-1.7204900),
    ((2,"HN"),(1,"CA"))  :   (-0.0164554, 1.0825300, 2.7774400),
    ((2,"HN"),(1,"HA"))  :   (-0.9106520,-1.4312400,-4.8713100),
    ((2,"HN"),(1,"C"))   :   (-0.0000000,-0.0000000, 0.0000000),
    ((2,"HN"),(1,"O"))   :   (-0.0000000,-0.0000000, 0.0000000),
    ((2,"HN"),(3,"N"))   :   (-1.2480900, 1.2781100, 0.3590480),
    ((2,"HN"),(3,"HN"))  :   ( 0.0000000,-0.0000000,-0.0000000),
    ((2,"HN"),(3,"CA"))  :   ( 1.0147000,-0.7955240,-0.2224820),
    ((2,"HN"),(3,"HA"))  :   (-0.1289740, 0.1438270, 0.0419333),
    ((2,"HN"),(3,"C"))   :   (-0.1279050, 0.0809711,-0.0181968),
    ((2,"HN"),(2,"N"))   :   (-0.0000000,-0.0000000, 0.0000000),
    ((2,"HN"),(2,"CA"))  :   (-0.0000000,-0.0000000, 0.0000000),
    ((2,"HN"),(2,"HA"))  :   (-0.0859479,-0.2456050, 0.0672676),
    ((2,"HN"),(2,"C"))   :   ( 0.5232180,-1.3820500,-0.0163351),
    ((2,"HN"),(2,"O"))   :   (-0.7044810, 2.2994800,-0.2210250),
    ((2,"C"),(1,"N"))    :   ( 0.0023357,-0.0251328,-0.0073040),
    ((2,"C"),(1,"CA"))   :   ( 0.0596069,-0.2066480,-0.1333500),
    ((2,"C"),(1,"HA"))   :   ( 0.0251522,-0.2093780,-0.1833310),
    ((2,"C"),(1,"C"))    :   (-0.0000000, 0.0000000, 0.0000000),
    ((2,"C"),(1,"O"))    :   ( 0.1308320,-0.1261360,-0.2562320),
    ((2,"C"),(3,"N"))    :   ( 0.0000000, 0.0000000,-0.0000000),
    ((2,"C"),(3,"HN"))   :   ( 0.0000000, 0.0000000,-0.0000000),
    ((2,"C"),(3,"CA"))   :   (-1.1185600,-0.0704527, 0.3411650),
    ((2,"C"),(3,"HA"))   :   (-0.1302610, 0.0511947, 0.0604657),
    ((2,"C"),(3,"C"))    :   (-0.6479500,-0.0380029,-0.1308990),
    ((2,"C"),(2,"N"))    :   (-0.0000000, 0.0000000, 0.0000000),
    ((2,"C"),(2,"HN"))   :   ( 0.0949609,-0.2508330,-0.0029647),
    ((2,"C"),(2,"CA"))   :   (-0.0000000, 0.0000000, 0.0000000),
    ((2,"C"),(2,"HA"))   :   (-0.0000000, 0.0000000, 0.0000000),
    ((2,"C"),(2,"O"))    :   ( 0.0000000,-0.0000000, 0.0000000),
    ((2,"CB"),(1,"N"))   :   ( 0.1328010, 0.3115680, 0.2157070),
    ((2,"CB"),(1,"CA"))  :   ( 0.0436672, 0.1379600, 0.1715400),
    ((2,"CB"),(1,"HA"))  :   (-0.1544180,-0.2964320,-0.4522390),
    ((2,"CB"),(1,"C"))   :   ( 0.3409950, 0.6752600, 1.4501300),
    ((2,"CB"),(1,"O"))   :   (-0.1892780,-0.1976900,-1.3136200),
    ((2,"CB"),(3,"N"))   :   (-1.5110500, 0.0268204,-0.3663770),
    ((2,"CB"),(3,"HN"))  :   (-1.2175000,-0.3725810,-0.1795680),
    ((2,"CB"),(3,"CA"))  :   ( 3.2450100,-0.4600350, 0.4277650),
    ((2,"CB"),(3,"HA"))  :   (-0.1634070, 0.0630888,-0.0092494),
    ((2,"CB"),(3,"C"))   :   (-0.1404630, 0.0160253,-0.0525042),
    ((2,"CB"),(2,"N"))   :   (-1.6156500,-2.5646200,-3.8499300),
    ((2,"CB"),(2,"HN"))  :   ( 0.3810730, 0.8163550, 0.5260330),
    ((2,"CB"),(2,"CA"))  :   (-0.0000000, 0.0000000,-0.0000000),
    ((2,"CB"),(2,"HA"))  :   (-0.0144993, 0.6147720,-1.9728800),
    ((2,"CB"),(2,"C"))   :   ( 6.3267000,-2.3938000, 4.0275100),
    ((2,"CB"),(2,"O"))   :   (-1.0640100, 0.9421670,-0.8331520),
}

