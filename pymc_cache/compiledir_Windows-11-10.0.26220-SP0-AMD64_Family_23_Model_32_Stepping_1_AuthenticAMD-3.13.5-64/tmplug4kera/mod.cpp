#include <Python.h>
#include "pytensor_mod_helper.h"
#include <math.h>
#include <numpy/arrayobject.h>
#include <numpy/arrayscalars.h>
#include <numpy/npy_math.h>
//////////////////////
////  Support Code
//////////////////////

        /** ParamsType _Params_5f3a79074738131db04c07d9c11b5a873d43f6e049ea4bda50e4a6d382dc15ba_50f3d776e2f7b727b9f80294a199c920d8d16d39c954f6da00d13e7cb320a19f **/
        #ifndef _PARAMS_5F3A79074738131DB04C07D9C11B5A873D43F6E049EA4BDA50E4A6D382DC15BA_50F3D776E2F7B727B9F80294A199C920D8D16D39C954F6DA00D13E7CB320A19F
        #define _PARAMS_5F3A79074738131DB04C07D9C11B5A873D43F6E049EA4BDA50E4A6D382DC15BA_50F3D776E2F7B727B9F80294A199C920D8D16D39C954F6DA00D13E7CB320A19F
        struct _Params_5f3a79074738131db04c07d9c11b5a873d43f6e049ea4bda50e4a6d382dc15ba_50f3d776e2f7b727b9f80294a199c920d8d16d39c954f6da00d13e7cb320a19f {
            /* Attributes, */
            int _Params_5f3a79074738131db04c07d9c11b5a873d43f6e049ea4bda50e4a6d382dc15ba_50f3d776e2f7b727b9f80294a199c920d8d16d39c954f6da00d13e7cb320a19f_error;
            
                typedef npy_bool dtype_inplace;
            
        npy_bool inplace;
        

                typedef npy_bool dtype_set_instead_of_inc;
            
        npy_bool set_instead_of_inc;
        

            /* Constructor. */
            _Params_5f3a79074738131db04c07d9c11b5a873d43f6e049ea4bda50e4a6d382dc15ba_50f3d776e2f7b727b9f80294a199c920d8d16d39c954f6da00d13e7cb320a19f() {
                _Params_5f3a79074738131db04c07d9c11b5a873d43f6e049ea4bda50e4a6d382dc15ba_50f3d776e2f7b727b9f80294a199c920d8d16d39c954f6da00d13e7cb320a19f_error = 0;
                
        inplace = 0;
        

        set_instead_of_inc = 0;
        
            }

            /* Destructor. */
            ~_Params_5f3a79074738131db04c07d9c11b5a873d43f6e049ea4bda50e4a6d382dc15ba_50f3d776e2f7b727b9f80294a199c920d8d16d39c954f6da00d13e7cb320a19f() {
                // cleanup() is defined below.
                cleanup();
            }

            /* Cleanup method. */
            void cleanup() {
                

            }

            /* Extraction methods. */
            
            void extract_inplace(PyObject* py_inplace) {
                
            if (!PyObject_TypeCheck(py_inplace, &PyBoolArrType_Type))
            {
                PyErr_Format(PyExc_ValueError,
                    "Scalar check failed (npy_bool)");
                {this->setErrorOccurred(); return;}
            }
            
        PyArray_ScalarAsCtype(py_inplace, &inplace);
        
            }
            


            void extract_set_instead_of_inc(PyObject* py_set_instead_of_inc) {
                
            if (!PyObject_TypeCheck(py_set_instead_of_inc, &PyBoolArrType_Type))
            {
                PyErr_Format(PyExc_ValueError,
                    "Scalar check failed (npy_bool)");
                {this->setErrorOccurred(); return;}
            }
            
        PyArray_ScalarAsCtype(py_set_instead_of_inc, &set_instead_of_inc);
        
            }
            

            /* Extract method. */
            
        void extract(PyObject* object, int field_pos) {
            switch(field_pos) {
                // Extraction cases.
                case 0: extract_inplace(object); break;
case 1: extract_set_instead_of_inc(object); break;
                // Default case.
                default:
                    PyErr_Format(PyExc_TypeError, "ParamsType: no extraction defined for a field %d.", field_pos);
                    this->setErrorOccurred();
                    break;
            }
        }
        

            /* Other methods. */
            void setErrorOccurred() {
                ++_Params_5f3a79074738131db04c07d9c11b5a873d43f6e049ea4bda50e4a6d382dc15ba_50f3d776e2f7b727b9f80294a199c920d8d16d39c954f6da00d13e7cb320a19f_error;
            }
            int errorOccurred() {
                return _Params_5f3a79074738131db04c07d9c11b5a873d43f6e049ea4bda50e4a6d382dc15ba_50f3d776e2f7b727b9f80294a199c920d8d16d39c954f6da00d13e7cb320a19f_error;
            }
        };
        #endif
        /** End ParamsType _Params_5f3a79074738131db04c07d9c11b5a873d43f6e049ea4bda50e4a6d382dc15ba_50f3d776e2f7b727b9f80294a199c920d8d16d39c954f6da00d13e7cb320a19f **/
        

    namespace {
    struct __struct_compiled_op_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402 {
        PyObject* __ERROR;

        PyObject* storage_V7;
PyObject* storage_V5;
PyObject* storage_V3;
PyObject* storage_V1;
PyObject* storage_V9;
        
    PyObject* py_V9;
    
        _Params_5f3a79074738131db04c07d9c11b5a873d43f6e049ea4bda50e4a6d382dc15ba_50f3d776e2f7b727b9f80294a199c920d8d16d39c954f6da00d13e7cb320a19f* V9;
        

        __struct_compiled_op_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402() {
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
        ~__struct_compiled_op_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402(void) {
            cleanup();
        }

        int init(PyObject* __ERROR, PyObject* storage_V7, PyObject* storage_V5, PyObject* storage_V3, PyObject* storage_V1, PyObject* storage_V9) {
            Py_XINCREF(storage_V7);
Py_XINCREF(storage_V5);
Py_XINCREF(storage_V3);
Py_XINCREF(storage_V1);
Py_XINCREF(storage_V9);
            this->storage_V7 = storage_V7;
this->storage_V5 = storage_V5;
this->storage_V3 = storage_V3;
this->storage_V1 = storage_V1;
this->storage_V9 = storage_V9;
            





    py_V9 = PyList_GET_ITEM(storage_V9, 0);
    {Py_XINCREF(py_V9);}
    
        /* Seems c_init() is not called for a op param. So I call `new` here. */
        V9 = new _Params_5f3a79074738131db04c07d9c11b5a873d43f6e049ea4bda50e4a6d382dc15ba_50f3d776e2f7b727b9f80294a199c920d8d16d39c954f6da00d13e7cb320a19f;

        { // This need a separate namespace for Clinker
        const char* fields[] = {"inplace", "set_instead_of_inc"};
        if (py_V9 == Py_None) {
            PyErr_SetString(PyExc_ValueError, "ParamsType: expected an object, not None.");
            {
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
            }
        return 9;
}
        }
        for (int i = 0; i < 2; ++i) {
            PyObject* o = PyDict_GetItemString(py_V9, fields[i]);
            if (o == NULL) {
                PyErr_Format(PyExc_TypeError, "ParamsType: missing expected attribute \"%s\" in object.", fields[i]);
                {
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
            }
        return 9;
}
            }
            V9->extract(o, i);
            if (V9->errorOccurred()) {
                /* The extract code from attribute type should have already raised a Python exception,
                 * so we just print the attribute name in stderr. */
                fprintf(stderr, "\nParamsType: error when extracting value for attribute \"%s\".\n", fields[i]);
                {
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
            }
        return 9;
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

double __DUMMY_5;
__label_7:

double __DUMMY_7;
__label_9:

        delete V9;
        V9 = NULL;
        
    {Py_XDECREF(py_V9);}
    
double __DUMMY_9;
__label_12:

double __DUMMY_12;

            Py_XDECREF(this->storage_V7);
Py_XDECREF(this->storage_V5);
Py_XDECREF(this->storage_V3);
Py_XDECREF(this->storage_V1);
Py_XDECREF(this->storage_V9);
        }
        int run(void) {
            int __failure = 0;
            
    PyObject* py_V1;
    
        PyArrayObject* V1;
        
            typedef npy_float64 dtype_V1;
            
    PyObject* py_V3;
    
        PyArrayObject* V3;
        
    PyObject* py_V5;
    
        PyArrayObject* V5;
        
    PyObject* py_V7;
    
        PyArrayObject* V7;
        
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

    py_V5 = PyList_GET_ITEM(storage_V5, 0);
    {Py_XINCREF(py_V5);}
    
        V5 = (PyArrayObject*)(py_V5);
        Py_XINCREF(V5);
        
{

    py_V7 = PyList_GET_ITEM(storage_V7, 0);
    {Py_XINCREF(py_V7);}
    
        V7 = (PyArrayObject*)(py_V7);
        Py_XINCREF(V7);
        
{

{
// Op class AdvancedIncSubtensor1

            if (V9->inplace)
            {
                if (V7 != V1)
                {
                    Py_XDECREF(V1);
                    Py_INCREF(V7);
                    V1 = V7;
                }
            }
            else
            {
                Py_XDECREF(V1);
                V1 = (PyArrayObject*)PyArray_FromAny(py_V7, NULL, 0, 0,
                NPY_ARRAY_ENSURECOPY, NULL);
                if (!V1) {
                    // Exception already set
                    {
        __failure = 11;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_11;}
                }
            }

            if (PyArray_NDIM(V1) != 1) {
                PyErr_Format(PyExc_ValueError, "AdvancedIncSubtensor1: first input (x) ndim should be 1, got %d", PyArray_NDIM(V1));
                {
        __failure = 11;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_11;}
            }
            if (PyArray_SHAPE(V1)[0] != 2) {
                PyErr_Format(PyExc_ValueError, "AdvancedIncSubtensor1: first input (x) shape should be 2, got %d", PyArray_SHAPE(V1)[0]);
                {
        __failure = 11;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_11;}
            }
            if (PyArray_NDIM(V3) != 1) {
                PyErr_Format(PyExc_ValueError, "AdvancedIncSubtensor1: indices ndim should be 1, got %d", PyArray_NDIM(V3));
                {
        __failure = 11;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_11;}
            }
            if (PyArray_NDIM(V5) != 1) {
                PyErr_Format(PyExc_ValueError, "AdvancedIncSubtensor1: second input (y) ndim should be 1, got %d", PyArray_NDIM(V5));
                {
        __failure = 11;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_11;}
            }
            if (PyArray_SHAPE(V5)[0] != PyArray_SHAPE(V3)[0]) {
                if ((PyArray_NDIM(V5) == 1) && (PyArray_SHAPE(V5)[0] == 1)){
                    PyErr_Format(PyExc_ValueError, "Runtime broadcasting not allowed. AdvancedIncSubtensor1 was asked to broadcast the second input (y) along a dimension that was not marked as broadcastable. If broadcasting was intended, use `specify_broadcastable` on the relevant dimension(s).");
                } else {
                    PyErr_Format(PyExc_ValueError,
                    "AdvancedIncSubtensor1: Shapes of second input (y) and indices do not match: %d, %d",
                    PyArray_SHAPE(V5)[0], PyArray_SHAPE(V3)[0]);
                }
                {
        __failure = 11;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_11;}
            }

            {
                npy_intp out_shape0 = PyArray_SHAPE(V1)[0];
                npy_float64* out_data = (npy_float64*)PyArray_DATA(V1);
                npy_float64* y_data = (npy_float64*)PyArray_DATA(V5);
                npy_int64* idx_data = (npy_int64*)PyArray_DATA(V3);
                npy_intp n = PyArray_SHAPE(V3)[0];
                npy_intp out_jump = PyArray_STRIDES(V1)[0] / PyArray_ITEMSIZE(V1);
                npy_intp y_jump = PyArray_STRIDES(V5)[0] / PyArray_ITEMSIZE(V5);
                npy_intp idx_jump = PyArray_STRIDES(V3)[0] / PyArray_ITEMSIZE(V3);

                for(int i = 0; i < n; i++){
                    npy_int64 idx = idx_data[i * idx_jump];
                    if (0){
                        if (idx < 0) {
                            idx += out_shape0;
                        }
                    }
                    if (0){
                        if ((idx < 0) || (idx >= out_shape0)) {
                            PyErr_Format(PyExc_IndexError,"index %d out of bounds for array with shape %d", idx_data[i * idx_jump], out_shape0);
                            {
        __failure = 11;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_11;}
                        }
                    }
                    out_data[idx * out_jump] += y_data[i * y_jump];
                }

            }
            __label_11:

double __DUMMY_11;

}
__label_10:

double __DUMMY_10;

}
__label_8:

        if (V7) {
            Py_XDECREF(V7);
        }
        
    {Py_XDECREF(py_V7);}
    
double __DUMMY_8;

}
__label_6:

        if (V5) {
            Py_XDECREF(V5);
        }
        
    {Py_XDECREF(py_V5);}
    
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
    

        static int __struct_compiled_op_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402_executor(__struct_compiled_op_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402 *self) {
            return self->run();
        }

        static void __struct_compiled_op_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402_destructor(PyObject *capsule) {
            __struct_compiled_op_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402 *self = (__struct_compiled_op_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402 *)PyCapsule_GetContext(capsule);
            delete self;
        }
    
//////////////////////
////  Functions
//////////////////////
static PyObject * instantiate(PyObject * self, PyObject *argtuple) {
  assert(PyTuple_Check(argtuple));
  if (6 != PyTuple_Size(argtuple)){ 
     PyErr_Format(PyExc_TypeError, "Wrong number of arguments, expected 6, got %%i", (int)PyTuple_Size(argtuple));
     return NULL;
  }
  __struct_compiled_op_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402* struct_ptr = new __struct_compiled_op_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402();
  if (struct_ptr->init( PyTuple_GET_ITEM(argtuple, 0),PyTuple_GET_ITEM(argtuple, 1),PyTuple_GET_ITEM(argtuple, 2),PyTuple_GET_ITEM(argtuple, 3),PyTuple_GET_ITEM(argtuple, 4),PyTuple_GET_ITEM(argtuple, 5) ) != 0) {
    delete struct_ptr;
    return NULL;
  }
    PyObject* thunk = PyCapsule_New((void*)(&__struct_compiled_op_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402_executor), NULL, __struct_compiled_op_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402_destructor);
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
  "m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402",
  NULL,
  -1,
  MyMethods,
};

PyMODINIT_FUNC PyInit_m9418323c09b5b828b2a70e1b581e0a9707e57fd05d31c0bcd06864d008613402(void) {
   import_array();
   
    PyObject *m = PyModule_Create(&moduledef);
    return m;
}
