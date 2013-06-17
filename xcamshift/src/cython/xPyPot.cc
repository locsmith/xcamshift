
#include "xPyPot.hh"

#include <pyConvert.hh>
#include <derivList.hh>

static swig_type_info* SWIGTYPE_PyPotPtr=0;
static swig_type_info* SWIGTYPE_CDSVectorTVec3Ptr=0;
static swig_type_info* SWIGTYPE_DerivListPtr=0;


static void
init()
{
 if ( SWIGTYPE_PyPotPtr ) return;

 // cache type info
 SWIGTYPE_PyPotPtr          = SWIGPY_TypeQuery("PyPot *");
 SWIGTYPE_CDSVectorTVec3Ptr = SWIGPY_TypeQuery("CDSVector<Vec3 > *");
 SWIGTYPE_DerivListPtr      = SWIGPY_TypeQuery("DerivList *");

 
} /* init */


PyPot::PyPot(const String&     instanceName,
		   PyObject*   pobj,
		   Simulation* sim)
  : Pot("PyPot",instanceName), pobj(pobj), sim(sim)
{
 init();

 //
 // avoidance of reference cycles: note that special care as been taken to 
 // avoid reference cycles in the constructor and destructor. In Python,
 // the destructor will only be called when the reference count hits zero,
 // and if this object has references to the Python object, a cycle (self-
 // reference) results. At destruction, it is assumed that this destructor
 // is called before that of the derived class has exited. Otherwise, a
 // memory leak will result.
 //

 String name = "PyPot: ";
 // disabled for python 2.2
 //if ( !PyInstance_Check(pobj) )
 //  throw PyConvert_Exception( name + "expected an instance" );

 energyMeth = PyObject_GetAttrString(pobj,"calcEnergy");
 if ( !energyMeth )
   throw PyConvert_Exception( name + 
			      "instance has no method named calcEnergy" );
 // Py_DECREF(energyMeth);

 if ( !PyCallable_Check(energyMeth) )
   throw PyConvert_Exception( name + "member calcEnergy is not callable" );
 
 bool ok=0;
 
 energyAndDerivsMeth    = PyObject_GetAttrString(pobj,"calcEnergyAndDerivs");
 energyAndDerivListMeth = PyObject_GetAttrString(pobj,
						 "calcEnergyAndDerivList");
 if ( energyAndDerivsMeth ) {
   if ( !PyCallable_Check(energyAndDerivsMeth) )
     throw PyConvert_Exception( name + 
			      "member calcEnergyAndDerivs is not callable" );
   ok=1;
   //   Py_DECREF(energyAndDerivsMeth);
 } else if ( energyAndDerivListMeth ) {

   if ( !PyCallable_Check(energyAndDerivListMeth) )
     throw PyConvert_Exception( name + 
				"member calcEnergyAndDerivList is not callable" );
 
   ok=1;
   //   Py_DECREF(energyAndDerivListMeth);
 }

 if ( !ok ) 
   throw PyConvert_Exception(name + 
			     "instance has no method named" +
			     "calcEnergyAndDerivs or calcEnergyAndDerivList" );


 PyObject* module = PyImport_ImportModule("cdsVector");
 if (!module)
   throw PyConvert_Exception(name + 
			     "unable to load module cdsVector");

  
 new_VecVec3 = PyObject_GetAttrString(module,"CDSVector_Vec3Ptr");
 if (!new_VecVec3)
   throw PyConvert_Exception(name + 
			     "module cdsVector does not contain the function "+
			     "CDSVector_Vec3Ptr");

 Py_DECREF(module);

} /* constructor */

PyPot::~PyPot()
{
 Py_DECREF(new_VecVec3);

 //  these are old defs
 Py_DECREF(energyAndDerivsMeth);
 Py_DECREF(energyMeth);
} /* destructor */



float_type
PyPot::calcEnergy()
{
 
 PyObject* pargs  = Py_BuildValue("()");
 PyObject* result = PyEval_CallObject(energyMeth,pargs);

 Py_DECREF(pargs);

 if ( !result ) {
   PyErr_Print();
   throw PyConvert_Exception( "PyPot::calcEnergy: python exception" );
 }
 
 PyObject* fresult = PyNumber_Float(result);
 Py_DECREF(result);
 String msg = "PyPot::calcEnergy: bad return type";
 if ( !fresult ) {
   cerr << msg << endl;
   throw PyConvert_TypeError(msg);
 }
 double pot=PyFloat_AsDouble(fresult);
 
 Py_DECREF(fresult);

 return pot;
} /* calcEnergy */


float_type
PyPot::calcEnergyAndDerivs(DerivList& dListArray)
{
 float_type energy=0.;
 
 if ( energyAndDerivListMeth ) {
   PyObject* dListPtr = 
     SWIGPY_NewPointerObj((void*)&dListArray,
			  SWIGTYPE_DerivListPtr,0);
   PyObject* pargs  = Py_BuildValue("(O)",dListPtr);
   PyObject* result = PyEval_CallObject(energyAndDerivListMeth,pargs);
  
  
   Py_DECREF(pargs);
   Py_DECREF(dListPtr);
  
   if ( !result ) {
     PyErr_Print();
     throw PyConvert_Exception(
                 "PyPot::calcEnergyAndDerivList: python exception" );
   }
  
  
   PyObject* fresult = PyNumber_Float(result);
   Py_DECREF(result);
   String msg = "PyPot::calcEnergyAndDerivList: bad return type";
   if ( !fresult ) {
     cerr << msg << endl;
     throw PyConvert_TypeError(msg);
   }
   energy=PyFloat_AsDouble(fresult);
   
   Py_DECREF(fresult);
 } else {
   DerivList::VectorVec3 dList(sim->numAtoms());
   for (int i=0 ; i<dList.size() ; i++)
     dList(i) = Vec3(0.,0.,0.);
  
   PyObject* dListPtr = 
     SWIGPY_NewPointerObj((void*)&dList,
  			SWIGTYPE_CDSVectorTVec3Ptr,0);
   PyObject* pargs  = Py_BuildValue("(O)",dListPtr);
  // PyObject* pDerivList = PyEval_CallObject(new_VecVec3,pargs1);
  // PyObject* pargs  = Py_BuildValue("(O)",pDerivList);
   PyObject* result = PyEval_CallObject(energyAndDerivsMeth,pargs);
  
  
   Py_DECREF(pargs);
   Py_DECREF(dListPtr);
  
   if ( !result ) {
     PyErr_Print();
     throw PyConvert_Exception( "PyPot::calcEnergyAndDerivs: python exception" );
   }
  
  
   PyObject* fresult = PyNumber_Float(result);
   Py_DECREF(result);
   String msg = "PyPot::calcEnergyAndDerivs: bad return type";
   if ( !fresult ) {
     cerr << msg << endl;
     throw PyConvert_TypeError(msg);
   }
   energy=PyFloat_AsDouble(fresult);
   
   for (int i=0 ; i<dList.size() ; i++)
     dListArray[sim](i) += dList(i);
  
   Py_DECREF(fresult);
 }

 return energy;
} /* calcEnergyAndDerivs */


void 
fromPy(PyPot*& ret, 
       PyObject*         obj)
{
 if ((SWIGPY_ConvertPtr(obj,
		      (void **) &ret,
		      SWIGTYPE_PyPotPtr,0)) != -1) 
   return;
 
 // check that object is an instance of a class derived from 
 // PyPot
 init();
 PyObject* thisObj = PyObject_GetAttrString(obj,"this");
 if ( !thisObj )
   throw PyConvert_Exception( String("instance has no member named this") );

 if ((SWIGPY_ConvertPtr(thisObj,
		      (void **) &ret,
		      SWIGTYPE_PyPotPtr,0)) != -1) 
   return;
 throw PyConvert_Exception("PyPot");
} /* from PyPot */

