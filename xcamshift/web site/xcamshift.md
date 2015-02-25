XCamShift
=========


XCamShift (XCS) is a new implementation of the camshift chemical shift force field[^1] for XPLOR-NIH [^2]. The camshaft implementation has a number of interesting and useful features

* It is very thoroughly tested against the original camshaft implementation in almost [^3] and has a large and full test suite which checks that XCS produces results that are identical with camshaft in almost.
* XCS uses cython[^4] to give a forcefield with native performance written in python
* XCS provides extensions against the standard camshift forcefield including the ability to carry out  ensemble calculations, weighting of shift terms and variable well widths to allow for variation in statistical errors in measured chemical shift.
* XCS is modular and individual chemical shift terms can be re weighted during calculations and there is scope for whole components of the force field [e.g hydrogen bonding, aromatic ring currents etc] to be modified or replaced independently. 
* XCS provides simple human readable files defining the components of the force field.
* XCS is open source and is free software under the Library Gnu Public License (LGPL)

For information on how to install, compile and use XCS see the sections below.

The Details
-----

###Installation

XCamShift can be installed from pre compiled binaries or can be compiled on the users own computer targeting a particular installation of xplor-nih. The precompiled binary python modules include all the software required to run XCS and do not require any other software to be installed other than XPLOR-NIH    

download XCamShift binaries (xplor nih 2.35):

1. OSX x86_64 (64 bit)
2. OSX i686 (32 bit) - coming
3. Linux x86_64 (64 bit)
4. Linux x86_64 (32 bit) - coming
6. sources: see Downloading the source from Gihub [below](#download_source) 

To install the XCamShift (see ~~setting paths~~ below if you want to use XCamShifts without modifying your xplor-nih installation):

> note: in the following <code>xxx</code> refers to a specific version e.g. v1, <code>\<platform\></code> and may include a computer platform  e.g. OSX or Linux  and a cpu architecture <code>\<architecture\></code> e.g. x86_64 or i686


0. if md5sum is installed check the md5sum of the download is correct <code>md5sum XCamShift_xxx.tgz</code> 
1. extract the XCamShift_xxx.tgz into a convenient temporary directory using the command <code>tar -zxvf XCamShift_xxx.tgz</code> 
2. move the .so files in xcamshaft_xxx/modules into the folder <code>python/bin.\<platform>\_\<architecture></code> in you xplor installation; on OSX \<platform> will be Darwin_\<XX> (where \<XX> is the darwin version) and just Linux for linux. \<architecture> will be either X86_64 or i686 for 64 and 32 bit computers respectively.
3. Move the directory database/XCamShift which contains the forcefield definition into the database directory of your XPLOR-NIH installation.
4. test the installation by running the script <code>src/suite.py</code> with the command line <code>pyXplor src/suite.py</code>. You shouldn't see any errors of the form <code>fail</code> or <code>error</code> or <code>abnormal program termination</code>; any of these show a bug is present and should be reported as described in bug reporting below.


###Using XCamShift

XCS can only be accessed from the python interface of xplor NIH currently. Here are some snippets showing how to carry out the main operations required to run XCS.

######importing shifts

 
    csFile=open("test_data/alvin/cs.dat")
    cs=csFile.read()
    csFile.close()

the format of a chemical shift restraint file is detail below under restraints

######setting up a calculation

    from XCamShift import XCamShift
    XCamShift.addRestraints(cs)
    potList.append(XCamShift)
    
######setting up weights

    XCamShift.setScale(6) # this will change the weights of all terms


######Format of there camshaft restraint list

camshaft restraint files are typical xplor restraints files and accept acclamation marks (!) as comments that run to the end of lines. The format is

    <weight-statement>|<class-statement>| <assign-statement>

    weight-statement = WEIGht <float>
    class-statement = CLASs <four-letter-code> 
    assign-statement = ASSIgn <selection> <float>|
                       ASSIgn <selection> <float> <float>

the <code>weight</code> and <code>\<class></code> statement set weights for all following <code>\<assign></code> statements. <assign> statements can take either  one and two float values which are in order: the chemical shift and the error of the chemical shift. It is not possible to set individual weights without also setting an error without recourse to the <code>weight</code> statement. The <code>class</code> statement can be used to tag shifts for later manipulation (currently you can just use this to reset the weight of a set of restraints.

An example of an XCamShift restraint file would be 

    class TEST                                     ! comment
    assign ( resid 2 and name C ) 177.477          ! restraint with default error
    weight 2.0
    assign ( resid 3 and name C ) 175.002 0.1      ! error is 0.1 weight 2.0
    weight 0.5
    assign ( resid 4 and name C ) 176.647 1.0      ! error 1.0 weight 0.5

It should be noted that chemical shift restraint files can't currently take sophisticated expressions as provided by the classic xplor syntax.

notes:

1. restraints from the first and last residues in a segment are ignored and a warning issued as these shift cannot be used by camshift as it calculates shifts for the atoms in the central residue of a triple of residues
2. segments with less than 3 residues are ignored and a warning is issued again because camshift as it calculates shifts for the atoms in the central residue of a triple of residues
 
######Restraint objects

Restraint objects can be retrieved from the XCamShift forcefield using the method <code>getRestraints()</code> on the XCamShift potential object. Restraint objects have the following attributes

* obs 	    -  the observed shift (can be set)
* calcd    -  the last calculated shift, if shifts have been calculated
* diff      -  obs - diff
* err       -  the err on the observed shift this will be the 1/2 width of the square well where no forces will be applied (can be set)
* weight	- the weight of the shift	(can be set)
* comment   - any end of line comment present when the shift was read in
* name      - the name of the atom being restrained (segid resid restype atom-name)
* calcd2    - calcd**2
* obs2      - obs**2

######Restraint statistics

the XCamShift force field object provides the following methods to retrieve information about the number of restraint and violations

* rms - the rms of diff between obs and calcd (see below for how a threshold is set) 
* violations - the number of violations
* numRestraints - the number of restraints

As the range of chemical shifts is not uniformly distributed the value of the threshold [set using the attribute <code>threshold</code> of the forcefield] for restraint violations has to be treated differently and is used as a multiplier of the error. Therefore a restraint is violated when the following equality tests true

    |diff| > threshold * err + err

######Listing predicted shifts

The chemical shifts predicted by XCamShift can be printed using the methods <code>print_shifts</code> of the potential. The shifts are printed in the following format

    SEGID     RESID    RES        HA       CA       HN        N        C       CB
    |AGB3|    2        GLN    4.5501  54.2330   8.7379 121.0727 175.1000  29.2632
    |AGB3|    3        TYR    5.0569  56.5513   9.5938 123.5669 174.7284  40.6070
    |AGB3|    4        LYS    4.8969  54.5661   9.4434 121.8343 174.0538  35.2424
    |AGB3|    5        LEU    4.6559  52.7604   9.2490 124.7154 174.3615  43.3438
    |AGB3|    6        VAL    4.2273  60.9981   8.9896 124.2064 174.2570  32.7248
    |AGB3|    7        ILE    4.2590  59.6785   9.0951 126.6597 174.7850  37.6970
    |AGB3|    8        ASN    4.6711  51.7622   9.3382 125.9824 173.3525  38.0188
    |AGB3|    9        GLY    4.0863  45.1157   8.6680 110.2235 174.1825        .
    |AGB3|    10       LYS    3.9673  58.6702   9.3420 122.0362 176.9272  31.7829
    |AGB3|    11       THR    4.2896  63.8501   8.1219 113.1026 174.1789  69.9168
    ...
    |AGB3|    54       VAL    4.4473  59.3198   8.6110 121.8589 173.1164  32.8183
    |AGB3|    55       THR    4.8757  60.3260   8.6678 120.6908 173.4372  71.6735

note that missing chemical shifts which are mythically possible e.g. GLY CB are replaced by a <code>.</code> and that segids are wrapped with a  <code>||</code> pair to make spaces visible.


###Rebuilding XCamShift

You should only need to recompile XCamShift if you want to make changes to the code in XCamShift,  re-compile it on a system which isn't supported or use a different compiler. If you only need to recompile with a different compiler or different target version of xplor-nih you won't need cython, you can just recompile the C code produced by cython.

####Pre requisites

To rebuild the camshift tables (it is unlikely that you will need to do this!)  if you are running on an older version of python you may need to build argparse[^5] and add it to the <code>PYTHONPATH</code>

####[Downloading the source from Gihub](id:download_source)

source code can be downloaded from the github master branch using the following command <code></code>


####Rebuilding the C code

1. Firstly make sure you have a compiler that is compatible with GCC. On linux GCC ~~and ICC~~ have been tested. On OSX I have only used GCC from the homebrew project not Clang.
2. Make sure you have installed the source for your distribution of XPLOR-NIH
3. Download the source code for the version of python you want to target and build it. To test which version of python explore is using run <code>\<xplor-nih directory>/bin/pyXplor test.py</code> with test.py containing  <code>import sys ; print sys.version</code>. This will give you the version of python, typically the first line of output should be of the form <code>2.7.5 (default, Mar  9 2014, 22:15:05)</code> and in this case the version is 2.7.5. 
4. edit <code>build_all.py</code> and change the variables
5. run <code>build_all.py</code>

####Setting Paths
XCamShift can be tested without modifying xplor-nih by adding the following to you PYTHONPATH

<code>\<XCamShift-directory\>/src/cython:\<XCamShift-directory\>/XCamShift/src</code>

where <XCamShift-directory> is the directory that contains the XCamShift distribution

if they are not installed already you will also have to add python paths for

1. nanotime
2. unittest2
3. PyYAML

####Testing XCamShift

If XCamShift is installed or the paths are correctly set XCamShift can be tested using the following command

<code>pyXplor \<XCamShift-directory\>/src/test/test_XCamShift_print_shifts.py</code>


####Reporting bugs

Bugs happen! XCamShift is well test but new code, please report bugs (or make requests for enhancements or even submit enhancements!) on the XCamShift github site: 

The author can also be contacted at 'g dot s dot thompson at leeds dot ac dot uk'


###Tested versions

* __Targeted version of Camshift:__ XCamShift is currently tested against version 1.35.0 of camshift which is implemented in almost XXX.yyyy.
* __Supported versions of XPLOR-NIH:__ XCamShift currently supports version of XPLOR-NIH from 2.31 to 2.35 and the precompiled binaries are compiled against XPLOR-NIH 2.35.Please note that for some versions of xplor-nih before 2.35 you may need to recompile from source as the interface for the ensemble code has changed (this definitely included xplor-nih 2.33).  
* * __Compiler environment:__ XCamShift binaries are compiled with GCC version 4.8.2 on OSX and Cython 0.20.2 (Note Cython is only required if the python code is modified XCamShift can be re-compiled from the C++ source without Cython being present).
* __Targteted platforms:__ XCamShift can be used on OSX and Linux under either a 32Bit or 64Bit intel environment (however 32bit OSX is not tested). 

    
###Todo


1. list all restraints add and remove restraints **done restraints work**
3. list violations  **done**
5. setup of ensemble calculations and weighting
7. listing calculated shifts
12. fill in  refs
13. test release 
14. hold all c data structures in cython **done**
14. release
11. modified non bonded potential
12. re weighting individual forcefield components
13. hbonds
15. format of forcefield component files
15. list status of all forcefield components
15.  set weights for classes 
16.  properly use modified
17.  check restraint base class properly implemented
15. test atom subsets
16. clean code and make sure 
17. integration with ccpn
18. listing the shift table
19. add holds for all native structures
20. vectorised calculations...
21. make restraints true restraints
22. add signal catchers
23. propogate exception from cython
24. all calls should go through the kitchen handler
25. remove dependancy on unittest2
 






[^1]: camshift reference
[^2]: xplor nih
[^3]: almost
[^4]: cython
[^5]: argparse



