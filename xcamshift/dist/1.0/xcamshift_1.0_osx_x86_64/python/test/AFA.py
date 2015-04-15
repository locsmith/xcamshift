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
Created on 17 Feb 2012

@author: garyt
'''
afa_total_shifts = {
    1 :  {  "HA" :   0.0000,  "CA" :  0.000000,  "HN" : 0.000000,  "N" :   0.000000,  "C" :   0.00000,   "CB" :  0.000000 },
    2 :  {  "HA":   4.05589,  "CA" : 58.457000,  "HN" : 8.938928,  "N" : 120.772093,  "C" : 173.516971,  "CB" : 37.960958 },
    3 :  {  "HA" :  0.00000,  "CA" :  0.000000,  "HN" : 0.000000,  "N" :   0.000000,  "C" :   0.00000,   "CB" :  0.000000 }
}


expected_ring_centre = (0.3425, 3.63183, -2.0035)

expected_ring_normals = (-1.37097,0.844347,0.54263)

expected_ring_shifts = {
    ('',2,"CA") :  0.137933,
    ('',2,"HA") :  0.231032,
    ('',2,"N")  : -0.081760,
    ('',2,"HN") :  0.111872,
    ('',2,"C")  : -0.043368,
    ('',2,"CB") :  0.928512
}

expected_shifts = {
    ('',2,"CA")   :   58.457000,
    ('',2,"HA")   :    4.055890,
    ('',2,"N")    :  120.772000,
    ('',2,"HN")   :    8.938930,
    ('',2,"C")    :  173.517000,
    ('',2,"CB")   :   37.961000
}

force_factors_harmonic = {
    ('',2,"CA")  :     5.927050,
    ('',2,"HA")  :    34.031100,
    ('',2,"N")   :    -2.567150,
    ('',2,"HN")  :    51.320900,
    ('',2,"C")   :    29.743300,
    ('',2,"CB")  :    18.561900
}

target_forces_harmonic = {
    ('',2,"CA")  : ( 0.049920,-0.069289,-0.006512),
    ('',2,"HA")  : ( 0.405700,-0.387447,-0.157811),
    ('',2,"N")   : (-0.002675, 0.010851,-0.001887),
    ('',2,"HN")  : ( 0.040646,-0.136146, 0.002054),
    ('',2,"C")   : ( 0.218287, 0.031321,-0.208593),
    ('',2,"CB")  : (-0.162282,-0.578938, 0.490843)
}

ring_forces_harmonic = {
    (('',2,"CA"),('',2,"CG"))   :   ( 0.006661, 0.002331,-0.004828),
    (('',2,"CA"),('',2,"CD1"))  :   (-0.093080, 0.063740, 0.034617),
    (('',2,"CA"),('',2,"CE1"))  :   ( 0.061458,-0.031427,-0.026533),
    (('',2,"CA"),('',2,"CZ"))   :   (-0.023271, 0.020765, 0.007019),
    (('',2,"CA"),('',2,"CE2"))  :   ( 0.076409,-0.040644,-0.032467),
    (('',2,"CA"),('',2,"CD2"))  :   (-0.078098, 0.054523, 0.028704),
    
    (('',2,"HA"),('',2,"CG"))   :   ( 0.207681,-0.104944,-0.082529),
    (('',2,"HA"),('',2,"CD1"))  :   (-0.902432, 0.578688, 0.356590),
    (('',2,"HA"),('',2,"CE1"))  :   ( 0.491901,-0.280020,-0.195156),
    (('',2,"HA"),('',2,"CZ"))   :   (-0.342814, 0.234093, 0.135358),
    (('',2,"HA"),('',2,"CE2"))  :   ( 0.767098,-0.449539,-0.304212),
    (('',2,"HA"),('',2,"CD2"))  :   (-0.627134, 0.409169, 0.247759),
    
    (('',2,"N"),('',2,"CG"))    :   (-0.001391,-0.000678, 0.001040),
    (('',2,"N"),('',2,"CD1"))   :   ( 0.010145,-0.007781,-0.003522),
    (('',2,"N"),('',2,"CE1"))   :   (-0.007416, 0.003033, 0.003426),
    (('',2,"N"),('',2,"CZ"))    :   ( 0.002279,-0.002939,-0.000413),
    (('',2,"N"),('',2,"CE2"))   :   (-0.009249, 0.004164, 0.004154),
    (('',2,"N"),('',2,"CD2"))   :   ( 0.008308,-0.006650,-0.002797),
    
    (('',2,"HN"),('',2,"CG"))   :   ( 0.041877,-0.007260,-0.019567),
    (('',2,"HN"),('',2,"CD1"))  :   (-0.188895, 0.134843, 0.071710),
    (('',2,"HN"),('',2,"CE1"))  :   ( 0.126695,-0.059510,-0.053170),
    (('',2,"HN"),('',2,"CZ"))   :   (-0.055386, 0.052642, 0.018930),
    (('',2,"HN"),('',2,"CE2"))  :   ( 0.175306,-0.089461,-0.072442),
    (('',2,"HN"),('',2,"CD2"))  :   (-0.140244, 0.104892, 0.052485),
    
    (('',2,"C"),('',2,"CG"))    :   (-0.062485, 0.010933, 0.045178),
    (('',2,"C"),('',2,"CD1"))   :   (-0.354042, 0.190344, 0.160416),
    (('',2,"C"),('',2,"CE1"))   :   ( 0.307384,-0.216937,-0.101297),
    (('',2,"C"),('',2,"CZ"))    :   (-0.010030,-0.021373, 0.024416),
    (('',2,"C"),('',2,"CE2"))   :   ( 0.281033,-0.200784,-0.090948),
    (('',2,"C"),('',2,"CD2"))   :   (-0.380146, 0.206497, 0.170828),
    
    (('',2,"CB"),('',2,"CG"))   :   ( 0.027047, 0.096490,-0.081807),
    (('',2,"CB"),('',2,"CD1"))  :   ( 0.027049, 0.096489,-0.081808),
    (('',2,"CB"),('',2,"CE1"))  :   ( 0.027045, 0.096491,-0.081807),
    (('',2,"CB"),('',2,"CZ"))   :   ( 0.027047, 0.096490,-0.081807),
    (('',2,"CB"),('',2,"CE2"))  :   ( 0.027045, 0.096491,-0.081807),
    (('',2,"CB"),('',2,"CD2"))  :   ( 0.027049, 0.096489,-0.081808),
}
