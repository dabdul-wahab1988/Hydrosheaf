#include <Python.h>
#include "pytensor_mod_helper.h"
#include <math.h>
#include <numpy/arrayobject.h>
#include <numpy/arrayscalars.h>
#include <numpy/npy_math.h>
//////////////////////
////  Support Code
//////////////////////

        /** ParamsType _Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437 **/
        #ifndef _PARAMS_23422FF40AD38AEC9C495B47C60E1ACCBED4344350F555C9609B89DC84EC7EDD_740788EE32F5ABCBDC5270FA40DB0968AA376AF2929E8E215CC39A0377B73437
        #define _PARAMS_23422FF40AD38AEC9C495B47C60E1ACCBED4344350F555C9609B89DC84EC7EDD_740788EE32F5ABCBDC5270FA40DB0968AA376AF2929E8E215CC39A0377B73437
        struct _Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437 {
            /* Attributes, */
            int _Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437_error;
            
        PyArrayObject* _new_order;
        
            typedef npy_int64 dtype__new_order;
            

                typedef npy_int64 dtype_input_ndim;
            
        npy_int64 input_ndim;
        

            /* Constructor. */
            _Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437() {
                _Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437_error = 0;
                
        _new_order = NULL;
        

        input_ndim = 0;
        
            }

            /* Destructor. */
            ~_Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437() {
                // cleanup() is defined below.
                cleanup();
            }

            /* Cleanup method. */
            void cleanup() {
                
        if (_new_order) {
            Py_XDECREF(_new_order);
        }
        

            }

            /* Extraction methods. */
            
            void extract__new_order(PyObject* py__new_order) {
                
            _new_order = NULL;
            if (py__new_order == Py_None) {
                // We can either fail here or set _new_order to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {this->setErrorOccurred(); return;}
            }
            if (!PyArray_Check(py__new_order)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {this->setErrorOccurred(); return;}
            }
            // We expect NPY_INT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py__new_order)) {
                PyArrayObject * tmp = (PyArrayObject*) py__new_order;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_INT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_INT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py__new_order),
                             (long int) PyArray_NDIM(tmp),
                             (long int) (PyArray_NDIM(tmp) >= 3 ?
            PyArray_DIMS(tmp)[PyArray_NDIM(tmp)-3] : -1),
                             (long int) (PyArray_NDIM(tmp) >= 2 ?
            PyArray_DIMS(tmp)[PyArray_NDIM(tmp)-2] : -1),
                             (long int) (PyArray_NDIM(tmp) >= 1 ?
            PyArray_DIMS(tmp)[PyArray_NDIM(tmp)-1] : -1),
                             (long int) (PyArray_NDIM(tmp) >= 3 ?
            PyArray_STRIDES(tmp)[PyArray_NDIM(tmp)-3] : -1),
                             (long int) (PyArray_NDIM(tmp) >= 2 ?
            PyArray_STRIDES(tmp)[PyArray_NDIM(tmp)-2] : -1),
                             (long int) (PyArray_NDIM(tmp) >= 1 ?
            PyArray_STRIDES(tmp)[PyArray_NDIM(tmp)-1] : -1)
            );
                {this->setErrorOccurred(); return;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py__new_order) != NPY_INT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_INT64) got %d",
                             NPY_INT64, PyArray_TYPE((PyArrayObject*) py__new_order));
                {this->setErrorOccurred(); return;}
            }
            
        _new_order = (PyArrayObject*)(py__new_order);
        Py_XINCREF(_new_order);
        
            }
            


            void extract_input_ndim(PyObject* py_input_ndim) {
                
            if (!PyObject_TypeCheck(py_input_ndim, &PyInt64ArrType_Type))
            {
                PyErr_Format(PyExc_ValueError,
                    "Scalar check failed (npy_int64)");
                {this->setErrorOccurred(); return;}
            }
            
        PyArray_ScalarAsCtype(py_input_ndim, &input_ndim);
        
            }
            

            /* Extract method. */
            
        void extract(PyObject* object, int field_pos) {
            switch(field_pos) {
                // Extraction cases.
                case 0: extract__new_order(object); break;
case 1: extract_input_ndim(object); break;
                // Default case.
                default:
                    PyErr_Format(PyExc_TypeError, "ParamsType: no extraction defined for a field %d.", field_pos);
                    this->setErrorOccurred();
                    break;
            }
        }
        

            /* Other methods. */
            void setErrorOccurred() {
                ++_Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437_error;
            }
            int errorOccurred() {
                return _Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437_error;
            }
        };
        #endif
        /** End ParamsType _Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437 **/
        

#define APPLY_SPECIFIC(str) str##_node_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b_0
#define PARAMS_TYPE _Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437
#define DTYPE_PARAM__new_order npy_int64
#define DTYPE_PARAM_input_ndim npy_int64


int APPLY_SPECIFIC(cpu_dimshuffle)(PyArrayObject *input, PyArrayObject **res, PARAMS_TYPE *params) {
    npy_int64* new_order;
    npy_intp nd_in;
    npy_intp nd_out;
    npy_intp* dimensions;
    npy_intp* strides;

    if (!PyArray_IS_C_CONTIGUOUS(params->_new_order)) {
        PyErr_SetString(PyExc_RuntimeError, "DimShuffle: param _new_order must be C-contiguous.");
        return 1;
    }
    new_order = (npy_int64*) PyArray_DATA(params->_new_order);
    nd_in = (npy_intp)(params->input_ndim);
    nd_out = PyArray_SIZE(params->_new_order);

    if (PyArray_NDIM(input) != nd_in) {
        PyErr_SetString(PyExc_ValueError, "DimShuffle: Input has less dimensions than expected.");
        return 1;
    }

    // Compute new dimensions and strides
    dimensions = (npy_intp*) malloc(nd_out * sizeof(npy_intp));
    strides = (npy_intp*) malloc(nd_out * sizeof(npy_intp));
    if (dimensions == NULL || strides == NULL) {
        PyErr_NoMemory();
        free(dimensions);
        free(strides);
        return 1;
    };

    npy_intp original_size = PyArray_SIZE(input);
    npy_intp new_size = 1;
    for (npy_intp i = 0; i < nd_out; ++i) {
        // We set the strides of length 1 dimensions to PyArray_ITEMSIZE(input).
        // The value is arbitrary, because there is never a next element.
        // np.expand_dims(x, 0) and x[None] do different things here.
        // I would prefer zero, but there are some poorly implemented BLAS operations
        // That don't handle zero strides correctly. At least they won't fail because of DimShuffle.
        if (new_order[i] != -1) {
            dimensions[i] = PyArray_DIMS(input)[new_order[i]];
            strides[i] = PyArray_DIMS(input)[new_order[i]] == 1 ? PyArray_ITEMSIZE(input) : PyArray_STRIDES(input)[new_order[i]];
        } else {
            dimensions[i] = 1;
            strides[i] = PyArray_ITEMSIZE(input);
        }
        new_size *= dimensions[i];
    }

    if (original_size != new_size) {
        PyErr_SetString(PyExc_ValueError, "DimShuffle: Attempting to squeeze axes with size not equal to one.");
        free(dimensions);
        free(strides);
        return 1;
    }

    if (*res)
        Py_XDECREF(*res);

    // Create the new array.
    *res = (PyArrayObject*)PyArray_New(&PyArray_Type, nd_out, dimensions,
                                       PyArray_TYPE(input), strides,
                                       PyArray_DATA(input), PyArray_ITEMSIZE(input),
                                       // borrow only the writable flag from the base
                                       // the NPY_OWNDATA flag will default to 0.
                                       (NPY_ARRAY_WRITEABLE * PyArray_ISWRITEABLE(input)),
                                       NULL);

    if (*res == NULL) {
        free(dimensions);
        free(strides);
        return 1;
    }

    // Declare it a view of the original input
    Py_INCREF((PyObject*)input);
    PyArray_SetBaseObject(*res, (PyObject*)input);

    // recalculate flags: CONTIGUOUS, FORTRAN, ALIGNED
    PyArray_UpdateFlags(*res, NPY_ARRAY_UPDATE_ALL);

    free(strides);
    free(dimensions);
    return 0;
}
#undef APPLY_SPECIFIC
#undef PARAMS_TYPE
#undef DTYPE_PARAM__new_order
#undef DTYPE_PARAM_input_ndim

    namespace {
    struct __struct_compiled_op_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b {
        PyObject* __ERROR;

        PyObject* storage_V3;
PyObject* storage_V1;
PyObject* storage_V5;
        
    PyObject* py_V5;
    
        _Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437* V5;
        

        __struct_compiled_op_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b() {
            // This is only somewhat safe because we:
            //  1) Are not a virtual class
            //  2) Do not use any virtual classes in the members
            //  3) Deal with mostly POD and pointers

            // If this changes, we would have to revise this, but for
            // now I am tired of chasing segfaults because
            // initialization code had an error and some pointer has
            // a junk value.
            #ifndef PYTENSOR_DONT_MEMSET_STRUCT
            memset(this, 0, sizeof(*this));
            #endif
        }
        ~__struct_compiled_op_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b(void) {
            cleanup();
        }

        int init(PyObject* __ERROR, PyObject* storage_V3, PyObject* storage_V1, PyObject* storage_V5) {
            Py_XINCREF(storage_V3);
Py_XINCREF(storage_V1);
Py_XINCREF(storage_V5);
            this->storage_V3 = storage_V3;
this->storage_V1 = storage_V1;
this->storage_V5 = storage_V5;
            



    py_V5 = PyList_GET_ITEM(storage_V5, 0);
    {Py_XINCREF(py_V5);}
    
        /* Seems c_init() is not called for a op param. So I call `new` here. */
        V5 = new _Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437;

        { // This need a separate namespace for Clinker
        const char* fields[] = {"_new_order", "input_ndim"};
        if (py_V5 == Py_None) {
            PyErr_SetString(PyExc_ValueError, "ParamsType: expected an object, not None.");
            {
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
            }
        return 5;
}
        }
        for (int i = 0; i < 2; ++i) {
            PyObject* o = PyDict_GetItemString(py_V5, fields[i]);
            if (o == NULL) {
                PyErr_Format(PyExc_TypeError, "ParamsType: missing expected attribute \"%s\" in object.", fields[i]);
                {
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
            }
        return 5;
}
            }
            V5->extract(o, i);
            if (V5->errorOccurred()) {
                /* The extract code from attribute type should have already raised a Python exception,
                 * so we just print the attribute name in stderr. */
                fprintf(stderr, "\nParamsType: error when extracting value for attribute \"%s\".\n", fields[i]);
                {
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
            }
        return 5;
}
            }
        }
        }
        

            this->__ERROR = __ERROR;
            return 0;
        }
        void cleanup(void) {
            __label_1:

double __DUMMY_1;
__label_3:

double __DUMMY_3;
__label_5:

        delete V5;
        V5 = NULL;
        
    {Py_XDECREF(py_V5);}
    
double __DUMMY_5;
__label_8:

double __DUMMY_8;

            Py_XDECREF(this->storage_V3);
Py_XDECREF(this->storage_V1);
Py_XDECREF(this->storage_V5);
        }
        int run(void) {
            int __failure = 0;
            
    PyObject* py_V1;
    
        PyArrayObject* V1;
        
            typedef npy_bool dtype_V1;
            
    PyObject* py_V3;
    
        PyArrayObject* V3;
        
{

    py_V1 = PyList_GET_ITEM(storage_V1, 0);
    {Py_XINCREF(py_V1);}
    
        if (py_V1 == Py_None)
        {
            
        V1 = NULL;
        
        }
        else
        {
            
        V1 = (PyArrayObject*)(py_V1);
        Py_XINCREF(V1);
        
        }
        
{

    py_V3 = PyList_GET_ITEM(storage_V3, 0);
    {Py_XINCREF(py_V3);}
    
        V3 = (PyArrayObject*)(py_V3);
        Py_XINCREF(V3);
        
{

{
// Op class DimShuffle

                #define APPLY_SPECIFIC(str) str##_node_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b_0
#define PARAMS_TYPE _Params_23422ff40ad38aec9c495b47c60e1accbed4344350f555c9609b89dc84ec7edd_740788ee32f5abcbdc5270fa40db0968aa376af2929e8e215cc39a0377b73437
#define DTYPE_PARAM__new_order npy_int64
#define DTYPE_PARAM_input_ndim npy_int64
                {
                  if (APPLY_SPECIFIC(cpu_dimshuffle)(V3, &V1, V5) != 0) {
                    {
        __failure = 7;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_7;}
                  }
                }
                #undef APPLY_SPECIFIC
#undef PARAMS_TYPE
#undef DTYPE_PARAM__new_order
#undef DTYPE_PARAM_input_ndim
                __label_7:

double __DUMMY_7;

}
__label_6:

double __DUMMY_6;

}
__label_4:

        if (V3) {
            Py_XDECREF(V3);
        }
        
    {Py_XDECREF(py_V3);}
    
double __DUMMY_4;

}
__label_2:

    if (!__failure) {
      
        {Py_XDECREF(py_V1);}
        if (!V1) {
            Py_INCREF(Py_None);
            py_V1 = Py_None;
        }
        else if ((void*)py_V1 != (void*)V1) {
            py_V1 = (PyObject*)V1;
        }

        {Py_XINCREF(py_V1);}

        if (V1 && !PyArray_ISALIGNED((PyArrayObject*) py_V1)) {
            PyErr_Format(PyExc_NotImplementedError,
                         "c_sync: expected an aligned array, got non-aligned array of type %ld"
                         " with %ld dimensions, with 3 last dims "
                         "%ld, %ld, %ld"
                         " and 3 last strides %ld %ld, %ld.",
                         (long int) PyArray_TYPE((PyArrayObject*) py_V1),
                         (long int) PyArray_NDIM(V1),
                         (long int) (PyArray_NDIM(V1) >= 3 ?
        PyArray_DIMS(V1)[PyArray_NDIM(V1)-3] : -1),
                         (long int) (PyArray_NDIM(V1) >= 2 ?
        PyArray_DIMS(V1)[PyArray_NDIM(V1)-2] : -1),
                         (long int) (PyArray_NDIM(V1) >= 1 ?
        PyArray_DIMS(V1)[PyArray_NDIM(V1)-1] : -1),
                         (long int) (PyArray_NDIM(V1) >= 3 ?
        PyArray_STRIDES(V1)[PyArray_NDIM(V1)-3] : -1),
                         (long int) (PyArray_NDIM(V1) >= 2 ?
        PyArray_STRIDES(V1)[PyArray_NDIM(V1)-2] : -1),
                         (long int) (PyArray_NDIM(V1) >= 1 ?
        PyArray_STRIDES(V1)[PyArray_NDIM(V1)-1] : -1)
        );
            {
        __failure = 2;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_2;}
        }
        
      PyObject* old = PyList_GET_ITEM(storage_V1, 0);
      {Py_XINCREF(py_V1);}
      PyList_SET_ITEM(storage_V1, 0, py_V1);
      {Py_XDECREF(old);}
    }
    
        if (V1) {
            Py_XDECREF(V1);
        }
        
    {Py_XDECREF(py_V1);}
    
double __DUMMY_2;

}

            
        if (__failure) {
            // When there is a failure, this code puts the exception
            // in __ERROR.
            PyObject* err_type = NULL;
            PyObject* err_msg = NULL;
            PyObject* err_traceback = NULL;
            PyErr_Fetch(&err_type, &err_msg, &err_traceback);
            if (!err_type) {err_type = Py_None;Py_INCREF(Py_None);}
            if (!err_msg) {err_msg = Py_None; Py_INCREF(Py_None);}
            if (!err_traceback) {err_traceback = Py_None; Py_INCREF(Py_None);}
            PyObject* old_err_type = PyList_GET_ITEM(__ERROR, 0);
            PyObject* old_err_msg = PyList_GET_ITEM(__ERROR, 1);
            PyObject* old_err_traceback = PyList_GET_ITEM(__ERROR, 2);
            PyList_SET_ITEM(__ERROR, 0, err_type);
            PyList_SET_ITEM(__ERROR, 1, err_msg);
            PyList_SET_ITEM(__ERROR, 2, err_traceback);
            {Py_XDECREF(old_err_type);}
            {Py_XDECREF(old_err_msg);}
            {Py_XDECREF(old_err_traceback);}
        }
        // The failure code is returned to index what code block failed.
        return __failure;
        
        }
    };
    }
    

        static int __struct_compiled_op_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b_executor(__struct_compiled_op_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b *self) {
            return self->run();
        }

        static void __struct_compiled_op_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b_destructor(PyObject *capsule) {
            __struct_compiled_op_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b *self = (__struct_compiled_op_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b *)PyCapsule_GetContext(capsule);
            delete self;
        }
    
//////////////////////
////  Functions
//////////////////////
static PyObject * instantiate(PyObject * self, PyObject *argtuple) {
  assert(PyTuple_Check(argtuple));
  if (4 != PyTuple_Size(argtuple)){ 
     PyErr_Format(PyExc_TypeError, "Wrong number of arguments, expected 4, got %%i", (int)PyTuple_Size(argtuple));
     return NULL;
  }
  __struct_compiled_op_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b* struct_ptr = new __struct_compiled_op_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b();
  if (struct_ptr->init( PyTuple_GET_ITEM(argtuple, 0),PyTuple_GET_ITEM(argtuple, 1),PyTuple_GET_ITEM(argtuple, 2),PyTuple_GET_ITEM(argtuple, 3) ) != 0) {
    delete struct_ptr;
    return NULL;
  }
    PyObject* thunk = PyCapsule_New((void*)(&__struct_compiled_op_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b_executor), NULL, __struct_compiled_op_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b_destructor);
    if (thunk != NULL && PyCapsule_SetContext(thunk, struct_ptr) != 0) {
        PyErr_Clear();
        Py_DECREF(thunk);
        thunk = NULL;
    }

  return thunk; }

//////////////////////
////  Module init
//////////////////////
static PyMethodDef MyMethods[] = {
	{"instantiate", instantiate, METH_VARARGS, "undocumented"} ,
	{NULL, NULL, 0, NULL}
};
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b",
  NULL,
  -1,
  MyMethods,
};

PyMODINIT_FUNC PyInit_mbdf3bf8f9d988df04d522f1f0ba8154f278cbc7f96cb85fe9a37a82af513189b(void) {
   import_array();
   
    PyObject *m = PyModule_Create(&moduledef);
    return m;
}
