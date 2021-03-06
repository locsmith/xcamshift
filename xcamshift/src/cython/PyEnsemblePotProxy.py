# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.4
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_PyEnsemblePotProxy', [dirname(__file__)])
        except ImportError:
            import _PyEnsemblePotProxy
            return _PyEnsemblePotProxy
        if fp is not None:
            try:
                _mod = imp.load_module('_PyEnsemblePotProxy', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _PyEnsemblePotProxy = swig_import_helper()
    del swig_import_helper
else:
    import _PyEnsemblePotProxy
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class Modified(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Modified, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Modified, name)
    __repr__ = _swig_repr
    MOD_SELF = _PyEnsemblePotProxy.Modified_MOD_SELF
    MOD_SIMULATION = _PyEnsemblePotProxy.Modified_MOD_SIMULATION
    def __init__(self, *args, **kwargs): 
        this = _PyEnsemblePotProxy.new_Modified(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
    def set(self, *args, **kwargs): return _PyEnsemblePotProxy.Modified_set(self, *args, **kwargs)
    def clear(self, *args, **kwargs): return _PyEnsemblePotProxy.Modified_clear(self, *args, **kwargs)
    def update(self, *args, **kwargs): return _PyEnsemblePotProxy.Modified_update(self, *args, **kwargs)
    def value(self, *args, **kwargs): return _PyEnsemblePotProxy.Modified_value(self, *args, **kwargs)
    def __call__(self, *args, **kwargs): return _PyEnsemblePotProxy.Modified___call__(self, *args, **kwargs)
    __swig_destroy__ = _PyEnsemblePotProxy.delete_Modified
    __del__ = lambda self : None;

class ModifiedPtr(Modified):
    def __init__(self, this):
        try: self.this.append(this)
        except: self.this = this
        self.this.own(0)
        self.__class__ = Modified

Modified_swigregister = _PyEnsemblePotProxy.Modified_swigregister
Modified_swigregister(Modified)

class ModifiedBase(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ModifiedBase, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ModifiedBase, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_setmethods__["modified"] = _PyEnsemblePotProxy.ModifiedBase_modified_set
    __swig_getmethods__["modified"] = _PyEnsemblePotProxy.ModifiedBase_modified_get
    if _newclass:modified = _swig_property(_PyEnsemblePotProxy.ModifiedBase_modified_get, _PyEnsemblePotProxy.ModifiedBase_modified_set)
    __swig_setmethods__["registeredSimulations"] = _PyEnsemblePotProxy.ModifiedBase_registeredSimulations_set
    __swig_getmethods__["registeredSimulations"] = _PyEnsemblePotProxy.ModifiedBase_registeredSimulations_get
    if _newclass:registeredSimulations = _swig_property(_PyEnsemblePotProxy.ModifiedBase_registeredSimulations_get, _PyEnsemblePotProxy.ModifiedBase_registeredSimulations_set)
    __swig_destroy__ = _PyEnsemblePotProxy.delete_ModifiedBase
    __del__ = lambda self : None;
    def registerTo(self, *args, **kwargs): return _PyEnsemblePotProxy.ModifiedBase_registerTo(self, *args, **kwargs)
    def unRegister(self, *args, **kwargs): return _PyEnsemblePotProxy.ModifiedBase_unRegister(self, *args, **kwargs)
    def updateValues(self, *args, **kwargs): return _PyEnsemblePotProxy.ModifiedBase_updateValues(self, *args, **kwargs)

class ModifiedBasePtr(ModifiedBase):
    def __init__(self, this):
        try: self.this.append(this)
        except: self.this = this
        self.this.own(0)
        self.__class__ = ModifiedBase

ModifiedBase_swigregister = _PyEnsemblePotProxy.ModifiedBase_swigregister
ModifiedBase_swigregister(ModifiedBase)

class PyEnsemblePotProxy(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, PyEnsemblePotProxy, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, PyEnsemblePotProxy, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _PyEnsemblePotProxy.new_PyEnsemblePotProxy(*args)
        try: self.this.append(this)
        except: self.this = this
    def __deref__(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy___deref__(self, *args, **kwargs)
    def __ref__(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy___ref__(self, *args, **kwargs)
    def registerInstanceData(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_registerInstanceData(self, *args, **kwargs)
    def decrRefCnt(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_decrRefCnt(self, *args, **kwargs)
    def incrRefCnt(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_incrRefCnt(self, *args, **kwargs)
    def refCnt(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_refCnt(self, *args, **kwargs)
    def instanceData(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_instanceData(self, *args, **kwargs)
    def help(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_help(self, *args, **kwargs)
    __oldinit__=__init__
    def __init__(self, *args):
        self.__oldinit__(*args)
        self.registerInstanceData(self)

    __swig_destroy__ = _PyEnsemblePotProxy.delete_PyEnsemblePotProxy
    __del__ = lambda self : None;
    __swig_setmethods__["m_obj"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_m_obj_set
    __swig_getmethods__["m_obj"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_m_obj_get
    if _newclass:m_obj = _swig_property(_PyEnsemblePotProxy.PyEnsemblePotProxy_m_obj_get, _PyEnsemblePotProxy.PyEnsemblePotProxy_m_obj_set)
    def energyMaybeDerivs0(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_energyMaybeDerivs0(self, *args, **kwargs)
    def energyMaybeDerivs1(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_energyMaybeDerivs1(self, *args, **kwargs)
    def energyMaybeDerivs2(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_energyMaybeDerivs2(self, *args, **kwargs)
    def energyMaybeDerivs3(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_energyMaybeDerivs3(self, *args, **kwargs)
    def energyMaybeDerivs4(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_energyMaybeDerivs4(self, *args, **kwargs)
    def callCyEnergyMaybeDerivs(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_callCyEnergyMaybeDerivs(self, *args, **kwargs)
    def ensembleSimulation(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_ensembleSimulation(self, *args, **kwargs)
    def rms(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_rms(self, *args, **kwargs)
    def violations(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_violations(self, *args, **kwargs)
    def numRestraints(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_numRestraints(self, *args, **kwargs)
    def help(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_help(self, *args, **kwargs)
    def calcEnergy(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_calcEnergy(self, *args, **kwargs)
    def calcEnergyAndDerivs(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_calcEnergyAndDerivs(self, *args, **kwargs)
    def energyMaybeDerivsPre(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_energyMaybeDerivsPre(self, *args, **kwargs)
    def energyMaybeDerivsPost(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_energyMaybeDerivsPost(self, *args, **kwargs)
    def simulation(self, *args): return _PyEnsemblePotProxy.PyEnsemblePotProxy_simulation(self, *args)
    def ensWeight(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_ensWeight(self, *args, **kwargs)
    def ensWeights(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_ensWeights(self, *args, **kwargs)
    def setEnsWeights(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_setEnsWeights(self, *args, **kwargs)
    def addEnsWeights(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_addEnsWeights(self, *args, **kwargs)
    def getEnsWeights(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_getEnsWeights(self, *args, **kwargs)
    def updateEnsWeights(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_updateEnsWeights(self, *args, **kwargs)
    def useSimEnsWeights(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_useSimEnsWeights(self, *args, **kwargs)
    def setUseSimEnsWeights(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_setUseSimEnsWeights(self, *args, **kwargs)
    def calcWDerivs(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_calcWDerivs(self, *args, **kwargs)
    def setCalcWDerivs(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_setCalcWDerivs(self, *args, **kwargs)
    def ensWeightsInfo(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_ensWeightsInfo(self, *args, **kwargs)
    def potName(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_potName(self, *args, **kwargs)
    def instanceName(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_instanceName(self, *args, **kwargs)
    def resetPotName(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_resetPotName(self, *args, **kwargs)
    def scale(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_scale(self, *args, **kwargs)
    def setScale(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_setScale(self, *args, **kwargs)
    def threshold(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_threshold(self, *args, **kwargs)
    def setThreshold(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_setThreshold(self, *args, **kwargs)
    def updateValues(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_updateValues(self, *args, **kwargs)
    def updateDelta(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_updateDelta(self, *args, **kwargs)
    __swig_setmethods__["instanceData_"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_instanceData__set
    __swig_getmethods__["instanceData_"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_instanceData__get
    if _newclass:instanceData_ = _swig_property(_PyEnsemblePotProxy.PyEnsemblePotProxy_instanceData__get, _PyEnsemblePotProxy.PyEnsemblePotProxy_instanceData__set)
    __swig_setmethods__["instanceDataCreate"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_instanceDataCreate_set
    __swig_getmethods__["instanceDataCreate"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_instanceDataCreate_get
    if _newclass:instanceDataCreate = _swig_property(_PyEnsemblePotProxy.PyEnsemblePotProxy_instanceDataCreate_get, _PyEnsemblePotProxy.PyEnsemblePotProxy_instanceDataCreate_set)
    __swig_setmethods__["instanceDataCleanup"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_instanceDataCleanup_set
    __swig_getmethods__["instanceDataCleanup"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_instanceDataCleanup_get
    if _newclass:instanceDataCleanup = _swig_property(_PyEnsemblePotProxy.PyEnsemblePotProxy_instanceDataCleanup_get, _PyEnsemblePotProxy.PyEnsemblePotProxy_instanceDataCleanup_set)
    __swig_setmethods__["modified"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_modified_set
    __swig_getmethods__["modified"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_modified_get
    if _newclass:modified = _swig_property(_PyEnsemblePotProxy.PyEnsemblePotProxy_modified_get, _PyEnsemblePotProxy.PyEnsemblePotProxy_modified_set)
    __swig_setmethods__["registeredSimulations"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_registeredSimulations_set
    __swig_getmethods__["registeredSimulations"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_registeredSimulations_get
    if _newclass:registeredSimulations = _swig_property(_PyEnsemblePotProxy.PyEnsemblePotProxy_registeredSimulations_get, _PyEnsemblePotProxy.PyEnsemblePotProxy_registeredSimulations_set)
    def registerTo(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_registerTo(self, *args, **kwargs)
    def unRegister(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_unRegister(self, *args, **kwargs)

class PyEnsemblePotProxyPtr(PyEnsemblePotProxy):
    def __init__(self, this):
        try: self.this.append(this)
        except: self.this = this
        self.this.own(0)
        self.__class__ = PyEnsemblePotProxy

PyEnsemblePotProxy_swigregister = _PyEnsemblePotProxy.PyEnsemblePotProxy_swigregister
PyEnsemblePotProxy_swigregister(PyEnsemblePotProxy)

realPyEnsemblePotProxy = PyEnsemblePotProxy
def PyEnsemblePotProxy(*args):
    from potProxy import PotProxy
    return PotProxy( realPyEnsemblePotProxy(*args) )
#this next bit may be specific to a swig version (works w/ 1.3.40)
realPyEnsemblePotProxy.__setattr__ = lambda self, name, value: _swig_setattr(self, realPyEnsemblePotProxy, name, value)
realPyEnsemblePotProxy.__getattr__ = lambda self, name: _swig_getattr(self, realPyEnsemblePotProxy, name)

PY_ENSEMBLE_POT_PROXY_H = _PyEnsemblePotProxy.PY_ENSEMBLE_POT_PROXY_H
class PyEnsemblePotProxy_LetterClass(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, PyEnsemblePotProxy_LetterClass, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, PyEnsemblePotProxy_LetterClass, name)
    __repr__ = _swig_repr
    __swig_setmethods__["m_obj"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_m_obj_set
    __swig_getmethods__["m_obj"] = _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_m_obj_get
    if _newclass:m_obj = _swig_property(_PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_m_obj_get, _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_m_obj_set)
    def __init__(self, *args, **kwargs): 
        this = _PyEnsemblePotProxy.new_PyEnsemblePotProxy_LetterClass(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _PyEnsemblePotProxy.delete_PyEnsemblePotProxy_LetterClass
    __del__ = lambda self : None;
    def energyMaybeDerivs0(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_energyMaybeDerivs0(self, *args, **kwargs)
    def energyMaybeDerivs1(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_energyMaybeDerivs1(self, *args, **kwargs)
    def energyMaybeDerivs2(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_energyMaybeDerivs2(self, *args, **kwargs)
    def energyMaybeDerivs3(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_energyMaybeDerivs3(self, *args, **kwargs)
    def energyMaybeDerivs4(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_energyMaybeDerivs4(self, *args, **kwargs)
    def callCyEnergyMaybeDerivs(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_callCyEnergyMaybeDerivs(self, *args, **kwargs)
    def ensembleSimulation(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_ensembleSimulation(self, *args, **kwargs)
    def rms(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_rms(self, *args, **kwargs)
    def violations(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_violations(self, *args, **kwargs)
    def numRestraints(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_numRestraints(self, *args, **kwargs)
    def help(self, *args, **kwargs): return _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_help(self, *args, **kwargs)

class PyEnsemblePotProxy_LetterClassPtr(PyEnsemblePotProxy_LetterClass):
    def __init__(self, this):
        try: self.this.append(this)
        except: self.this = this
        self.this.own(0)
        self.__class__ = PyEnsemblePotProxy_LetterClass

PyEnsemblePotProxy_LetterClass_swigregister = _PyEnsemblePotProxy.PyEnsemblePotProxy_LetterClass_swigregister
PyEnsemblePotProxy_LetterClass_swigregister(PyEnsemblePotProxy_LetterClass)

# This file is compatible with both classic and new-style classes.


