
%module xPyPot

%{
#include <xPyPot.hh>
%}

%include "cdsTypeMaps.i"
%include "help.i"
%include "accessor.i"
%include "enum.i"
%include "pot.i"

AddPot(PyPot)


//rewrite the PyPot constructor to remove need to supply the
// second self argument.

%pythoncode %{
PyPot = realPyPot #PyPot will never be created stand-alone,
                  #so it doesn't need to be PotProxy-wrapped.
		  #However, user-created Pots w/swig-wrapped components
		  #will have to be explicitly wrapped.
oldConstructor = PyPot.__init__
def construct(self, *args):
    if len(args)==1:
        args = list(args)
	args.append(self)
    oldConstructor(self,*args)
    _swig_setattr(self, PyPot, 'thisown', 1)
PyPot.__init__ = construct

import derivList
%}

//trying to *not* generate Python methods for these...
// unsuccessfully.
%feature("notabstract") PyPot;
%ignore Pot::calcEnergy;
%ignore PyPot::calcEnergy;
%ignore Pot::calcEnergyAndDerivs;
%ignore PyPot::calcEnergyAndDerivs;

%include "xPyPot.hh"

ClassHelp(PyPot,"nih-py-xPyPot")
ModuleHelp("nih-py-xPyPot")
