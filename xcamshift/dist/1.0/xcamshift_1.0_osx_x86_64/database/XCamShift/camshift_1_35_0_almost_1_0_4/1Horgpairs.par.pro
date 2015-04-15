# Parameter set created by ./extractParams.pl
# Rmsd of fit: 00.000000000 
# -----------------------------------------------------------------
# define parameters for the fit: constants, coefficients, exponents
# -----------------------------------------------------------------
# weighting contributions of different atoms against each other
WEIGHT    1.0
# setting the limit of the flat bottom potential
FLATBTM   00.000000000 
# scaling of the harmonic potential (past flat bottom potential)
SCALEHARM 0.62
# maximum value that tanh approaches on top of harmonic potential
TANHAMPLI 20.0
# shift difference past the flat bottom at which harmonic potential turns into tanh
ENDHARMON 4.0
# maximum deviation allowed between random coil value and predicted shift, set
# to rmsd during fit, a multiple of this will be used in CamShift
# 0.0 = deactivated
# MAXRCDEVI 0.646
MAXRCDEVI 0.0
# random coil values for this atom type and all 20 residues
RANDCOIL  8.24 8.23 8.3 8.34 8.32 8.32 8.42 8.33 8.32 8.05 8.16 8.29 8.28 8.3 0.0 8.31 8.15 8.15 8.12 8.03
# parameters for the equations that were fitted to the dihedral angles phi, psi, and chi1
DIHEDPHI  -0.4 -0.3 0.2 3.6 -3.4
DIHEDPSI  -0.104686 0.174676 0.0639679 3.51401 -2.82811
DIHEDCHI1 0.0343048 0.0448763 0.0406634 3.42026 4.15534
# constants
# a common constant to all amino acids
CONST     0
# amino acid specific constants:
#	 	ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL
CONSTAA   8.24 8.23 8.3 8.34 8.32 8.32 8.42 8.33 8.32 8.05 8.16 8.29 8.28 8.3 0.0 8.31 8.15 8.15 8.12 8.03
# adjust for residue i-1
CONSTAA-1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
# adjust for residue i+1
CONSTAA+1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
# coefficients
# for atom distances along the backbone
COBB1     00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COBB2     00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
# for atom distances along the various amino acid side chains
#		CB,  HB1, HB2, HB3
COSCALA1  00.000000000 00.000000000 00.000000000 00.000000000 
COSCALA2  00.000000000 00.000000000 00.000000000 00.000000000 
#	   	CB,  CG,  CD,  NE,  CZ,  NH1, NH2, HB1, HB2, HG1, HG2, HD1, HD2, HE,  HH11, HH12, HH21, HH22
COSCARG1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCARG2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, CG, OD1, ND2, HB1, HB2, HD21, HD22
COSCASN1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCASN2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#	    	CB, CG, OD1, OD2, HB1, HB2
COSCASP1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCASP2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, SG, HB1, HB2, HG1
COSCCYS1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCCYS2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, CG, CD, OE1, NE2, HB1, HB2, HG1, HG2, HE21, HE22
COSCGLN1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCGLN2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, CG, CD, OE1, OE2, HB1, HB2, HG1, HG2
COSCGLU1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCGLU2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		HA2
COSCGLY1  00.000000000 
COSCGLY2  00.000000000 
#		CB, CG, ND1, CD2, CE1, NE2, HB1, HB2, HD1, HD2, HE1, HE2
COSCHIS1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCHIS2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, CG1, CG2, CD, HB, HG11, HG12, HG21, HG22, HG23, HD1, HD2, HD3
COSCILE1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCILE2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, CG, CD1, CD2, HB1, HB2, HG1, HD11, HD12, HD13, HD21, HD22, HD23
COSCLEU1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCLEU2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, CG, CD, CE, NZ, HB1, HB2, HG1, HG2, HD1, HD2, HE1, HE2, HZ1, HZ2, HZ3
COSCLYS1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCLYS2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, CG, SD, CE, HB1, HB2, HG1, HG2, HE1, HE2, HE3
COSCMET1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCMET2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, CG, CD1, CD2, CE1, CE2, CZ, HB1, HB2, HD1, HD2, HE1, HE2, HZ
COSCPHE1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCPHE2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, CG, CD, HB1, HB2, HG1, HG2, HD1, HD2
COSCPRO1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCPRO2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, OG, HB1, HB2, HG1
COSCSER1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCSER2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, OG1, CG2, HB, HG1, HG21, HG22, HG23
COSCTHR1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCTHR2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2, HB1, HB2, HD1, HE1, HE3, HZ2, HZ3, HH2
COSCTRP1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCTRP2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, CG, CD1, CD2, CE1, CE2, CZ, OH, HB1, HB2, HD1, HD2, HE1, HE2, HH
COSCTYR1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCTYR2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#		CB, CG1, CG2, HB, HG11, HG12, HG13, HG21, HG22, HG23
COSCVAL1  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
COSCVAL2  00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#                for atoms within a certain cut off distance (C, N, O in two different hybridization states)
#                Type C, H, N, O, S, C_2, N_2, O_2
SPHERE1   00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
SPHERE2   00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#                for ring current effects
#                Phe, Tyr, Trp_1, Trp_2, His
RINGS     00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#                dihedral angles
#                phi, psi, chi1
DIHEDRALS 00.000000000 00.000000000 00.000000000 
#                hydrogen bonds, O i-1, H i, O i, H i+1
#                length, angle1, angle2, length, angle1, angle2, length, angle1, angle2, length, angle1, angle2
HBONDS    00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 
#                Additional distances
XTRADISTS 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 00.000000000 