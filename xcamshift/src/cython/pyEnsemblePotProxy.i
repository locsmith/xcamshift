
%module  PyEnsemblePotProxy

%{
#include <PyEnsemblePotProxy.hh>
%}

%include "cdsTypeMaps.i"
%include "help.i"
%include "cdsList_template.i"
%include "accessor.i"
%include "enum.i"
%include "pot.i"
%include <rc_ptr.hh>
%import <sthead.hh>


AddPot(PyEnsemblePotProxy)




%import  <ensemblePot.hh>
%include <PyEnsemblePotProxy.hh>




ClassHelp(PyEnsemblePotProxy,"nih-py-py-ensemble-pot-proxy")
#ModuleHelp("nih-py-ensemble-pot-proxy")
