#include <Python.h>
#include "pytensor_mod_helper.h"
#include <math.h>
#include <numpy/arrayobject.h>
#include <numpy/arrayscalars.h>
#include <numpy/npy_math.h>
#include <vector>
#include <algorithm>
//////////////////////
////  Support Code
//////////////////////

    namespace {
    struct __struct_compiled_op_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a {
        PyObject* __ERROR;

        PyObject* storage_V35;
PyObject* storage_V33;
PyObject* storage_V31;
PyObject* storage_V29;
PyObject* storage_V27;
PyObject* storage_V25;
PyObject* storage_V23;
PyObject* storage_V21;
PyObject* storage_V19;
PyObject* storage_V17;
PyObject* storage_V15;
PyObject* storage_V13;
PyObject* storage_V11;
PyObject* storage_V9;
PyObject* storage_V7;
PyObject* storage_V5;
PyObject* storage_V3;
PyObject* storage_V1;
        

        __struct_compiled_op_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a() {
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
        ~__struct_compiled_op_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a(void) {
            cleanup();
        }

        int init(PyObject* __ERROR, PyObject* storage_V35, PyObject* storage_V33, PyObject* storage_V31, PyObject* storage_V29, PyObject* storage_V27, PyObject* storage_V25, PyObject* storage_V23, PyObject* storage_V21, PyObject* storage_V19, PyObject* storage_V17, PyObject* storage_V15, PyObject* storage_V13, PyObject* storage_V11, PyObject* storage_V9, PyObject* storage_V7, PyObject* storage_V5, PyObject* storage_V3, PyObject* storage_V1) {
            Py_XINCREF(storage_V35);
Py_XINCREF(storage_V33);
Py_XINCREF(storage_V31);
Py_XINCREF(storage_V29);
Py_XINCREF(storage_V27);
Py_XINCREF(storage_V25);
Py_XINCREF(storage_V23);
Py_XINCREF(storage_V21);
Py_XINCREF(storage_V19);
Py_XINCREF(storage_V17);
Py_XINCREF(storage_V15);
Py_XINCREF(storage_V13);
Py_XINCREF(storage_V11);
Py_XINCREF(storage_V9);
Py_XINCREF(storage_V7);
Py_XINCREF(storage_V5);
Py_XINCREF(storage_V3);
Py_XINCREF(storage_V1);
            this->storage_V35 = storage_V35;
this->storage_V33 = storage_V33;
this->storage_V31 = storage_V31;
this->storage_V29 = storage_V29;
this->storage_V27 = storage_V27;
this->storage_V25 = storage_V25;
this->storage_V23 = storage_V23;
this->storage_V21 = storage_V21;
this->storage_V19 = storage_V19;
this->storage_V17 = storage_V17;
this->storage_V15 = storage_V15;
this->storage_V13 = storage_V13;
this->storage_V11 = storage_V11;
this->storage_V9 = storage_V9;
this->storage_V7 = storage_V7;
this->storage_V5 = storage_V5;
this->storage_V3 = storage_V3;
this->storage_V1 = storage_V1;
            



















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

double __DUMMY_9;
__label_11:

double __DUMMY_11;
__label_13:

double __DUMMY_13;
__label_15:

double __DUMMY_15;
__label_17:

double __DUMMY_17;
__label_19:

double __DUMMY_19;
__label_21:

double __DUMMY_21;
__label_23:

double __DUMMY_23;
__label_25:

double __DUMMY_25;
__label_27:

double __DUMMY_27;
__label_29:

double __DUMMY_29;
__label_31:

double __DUMMY_31;
__label_33:

double __DUMMY_33;
__label_35:

double __DUMMY_35;
__label_38:

double __DUMMY_38;

            Py_XDECREF(this->storage_V35);
Py_XDECREF(this->storage_V33);
Py_XDECREF(this->storage_V31);
Py_XDECREF(this->storage_V29);
Py_XDECREF(this->storage_V27);
Py_XDECREF(this->storage_V25);
Py_XDECREF(this->storage_V23);
Py_XDECREF(this->storage_V21);
Py_XDECREF(this->storage_V19);
Py_XDECREF(this->storage_V17);
Py_XDECREF(this->storage_V15);
Py_XDECREF(this->storage_V13);
Py_XDECREF(this->storage_V11);
Py_XDECREF(this->storage_V9);
Py_XDECREF(this->storage_V7);
Py_XDECREF(this->storage_V5);
Py_XDECREF(this->storage_V3);
Py_XDECREF(this->storage_V1);
        }
        int run(void) {
            int __failure = 0;
            
    PyObject* py_V1;
    
        PyArrayObject* V1;
        
            typedef npy_float64 dtype_V1;
            
    PyObject* py_V3;
    
        PyArrayObject* V3;
        
            typedef npy_float64 dtype_V3;
            
    PyObject* py_V5;
    
        PyArrayObject* V5;
        
            typedef npy_float64 dtype_V5;
            
    PyObject* py_V7;
    
        PyArrayObject* V7;
        
            typedef npy_bool dtype_V7;
            
    PyObject* py_V9;
    
        PyArrayObject* V9;
        
            typedef npy_float64 dtype_V9;
            
    PyObject* py_V11;
    
        PyArrayObject* V11;
        
            typedef npy_float64 dtype_V11;
            
    PyObject* py_V13;
    
        PyArrayObject* V13;
        
            typedef npy_float64 dtype_V13;
            
    PyObject* py_V15;
    
        PyArrayObject* V15;
        
            typedef npy_bool dtype_V15;
            
    PyObject* py_V17;
    
        PyArrayObject* V17;
        
            typedef npy_float64 dtype_V17;
            
    PyObject* py_V19;
    
        PyArrayObject* V19;
        
            typedef npy_float64 dtype_V19;
            
    PyObject* py_V21;
    
        PyArrayObject* V21;
        
            typedef npy_float64 dtype_V21;
            
    PyObject* py_V23;
    
        PyArrayObject* V23;
        
            typedef npy_float64 dtype_V23;
            
    PyObject* py_V25;
    
        PyArrayObject* V25;
        
            typedef npy_float64 dtype_V25;
            
    PyObject* py_V27;
    
        PyArrayObject* V27;
        
            typedef npy_float64 dtype_V27;
            
    PyObject* py_V29;
    
        PyArrayObject* V29;
        
            typedef npy_bool dtype_V29;
            
    PyObject* py_V31;
    
        PyArrayObject* V31;
        
            typedef npy_bool dtype_V31;
            
    PyObject* py_V33;
    
        PyArrayObject* V33;
        
            typedef npy_float64 dtype_V33;
            
    PyObject* py_V35;
    
        PyArrayObject* V35;
        
            typedef npy_float64 dtype_V35;
            
{

    py_V1 = PyList_GET_ITEM(storage_V1, 0);
    {Py_XINCREF(py_V1);}
    
        if (py_V1 == Py_None)
        {
            
        V1 = NULL;
        
        }
        else
        {
            
            V1 = NULL;
            if (py_V1 == Py_None) {
                // We can either fail here or set V1 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 2;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_2;}
            }
            if (!PyArray_Check(py_V1)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 2;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_2;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V1)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V1;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V1),
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
                {
        __failure = 2;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_2;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V1) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V1));
                {
        __failure = 2;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_2;}
            }
            
        V1 = (PyArrayObject*)(py_V1);
        Py_XINCREF(V1);
        
        }
        
{

    py_V3 = PyList_GET_ITEM(storage_V3, 0);
    {Py_XINCREF(py_V3);}
    
        if (py_V3 == Py_None)
        {
            
        V3 = NULL;
        
        }
        else
        {
            
            V3 = NULL;
            if (py_V3 == Py_None) {
                // We can either fail here or set V3 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 4;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_4;}
            }
            if (!PyArray_Check(py_V3)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 4;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_4;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V3)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V3;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V3),
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
                {
        __failure = 4;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_4;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V3) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V3));
                {
        __failure = 4;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_4;}
            }
            
        V3 = (PyArrayObject*)(py_V3);
        Py_XINCREF(V3);
        
        }
        
{

    py_V5 = PyList_GET_ITEM(storage_V5, 0);
    {Py_XINCREF(py_V5);}
    
        if (py_V5 == Py_None)
        {
            
        V5 = NULL;
        
        }
        else
        {
            
            V5 = NULL;
            if (py_V5 == Py_None) {
                // We can either fail here or set V5 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 6;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_6;}
            }
            if (!PyArray_Check(py_V5)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 6;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_6;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V5)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V5;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V5),
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
                {
        __failure = 6;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_6;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V5) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V5));
                {
        __failure = 6;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_6;}
            }
            
        V5 = (PyArrayObject*)(py_V5);
        Py_XINCREF(V5);
        
        }
        
{

    py_V7 = PyList_GET_ITEM(storage_V7, 0);
    {Py_XINCREF(py_V7);}
    
        if (py_V7 == Py_None)
        {
            
        V7 = NULL;
        
        }
        else
        {
            
            V7 = NULL;
            if (py_V7 == Py_None) {
                // We can either fail here or set V7 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 8;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_8;}
            }
            if (!PyArray_Check(py_V7)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 8;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_8;}
            }
            // We expect NPY_BOOL
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V7)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V7;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_BOOL), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_BOOL,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V7),
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
                {
        __failure = 8;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_8;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V7) != NPY_BOOL) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_BOOL) got %d",
                             NPY_BOOL, PyArray_TYPE((PyArrayObject*) py_V7));
                {
        __failure = 8;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_8;}
            }
            
        V7 = (PyArrayObject*)(py_V7);
        Py_XINCREF(V7);
        
        }
        
{

    py_V9 = PyList_GET_ITEM(storage_V9, 0);
    {Py_XINCREF(py_V9);}
    
        if (py_V9 == Py_None)
        {
            
        V9 = NULL;
        
        }
        else
        {
            
            V9 = NULL;
            if (py_V9 == Py_None) {
                // We can either fail here or set V9 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 10;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_10;}
            }
            if (!PyArray_Check(py_V9)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 10;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_10;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V9)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V9;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V9),
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
                {
        __failure = 10;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_10;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V9) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V9));
                {
        __failure = 10;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_10;}
            }
            
        V9 = (PyArrayObject*)(py_V9);
        Py_XINCREF(V9);
        
        }
        
{

    py_V11 = PyList_GET_ITEM(storage_V11, 0);
    {Py_XINCREF(py_V11);}
    
        if (py_V11 == Py_None)
        {
            
        V11 = NULL;
        
        }
        else
        {
            
            V11 = NULL;
            if (py_V11 == Py_None) {
                // We can either fail here or set V11 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 12;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_12;}
            }
            if (!PyArray_Check(py_V11)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 12;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_12;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V11)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V11;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V11),
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
                {
        __failure = 12;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_12;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V11) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V11));
                {
        __failure = 12;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_12;}
            }
            
        V11 = (PyArrayObject*)(py_V11);
        Py_XINCREF(V11);
        
        }
        
{

    py_V13 = PyList_GET_ITEM(storage_V13, 0);
    {Py_XINCREF(py_V13);}
    
        if (py_V13 == Py_None)
        {
            
        V13 = NULL;
        
        }
        else
        {
            
            V13 = NULL;
            if (py_V13 == Py_None) {
                // We can either fail here or set V13 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 14;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_14;}
            }
            if (!PyArray_Check(py_V13)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 14;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_14;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V13)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V13;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V13),
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
                {
        __failure = 14;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_14;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V13) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V13));
                {
        __failure = 14;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_14;}
            }
            
        V13 = (PyArrayObject*)(py_V13);
        Py_XINCREF(V13);
        
        }
        
{

    py_V15 = PyList_GET_ITEM(storage_V15, 0);
    {Py_XINCREF(py_V15);}
    
        if (py_V15 == Py_None)
        {
            
        V15 = NULL;
        
        }
        else
        {
            
            V15 = NULL;
            if (py_V15 == Py_None) {
                // We can either fail here or set V15 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 16;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_16;}
            }
            if (!PyArray_Check(py_V15)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 16;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_16;}
            }
            // We expect NPY_BOOL
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V15)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V15;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_BOOL), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_BOOL,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V15),
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
                {
        __failure = 16;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_16;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V15) != NPY_BOOL) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_BOOL) got %d",
                             NPY_BOOL, PyArray_TYPE((PyArrayObject*) py_V15));
                {
        __failure = 16;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_16;}
            }
            
        V15 = (PyArrayObject*)(py_V15);
        Py_XINCREF(V15);
        
        }
        
{

    py_V17 = PyList_GET_ITEM(storage_V17, 0);
    {Py_XINCREF(py_V17);}
    
        if (py_V17 == Py_None)
        {
            
        V17 = NULL;
        
        }
        else
        {
            
            V17 = NULL;
            if (py_V17 == Py_None) {
                // We can either fail here or set V17 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 18;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_18;}
            }
            if (!PyArray_Check(py_V17)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 18;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_18;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V17)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V17;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V17),
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
                {
        __failure = 18;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_18;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V17) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V17));
                {
        __failure = 18;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_18;}
            }
            
        V17 = (PyArrayObject*)(py_V17);
        Py_XINCREF(V17);
        
        }
        
{

    py_V19 = PyList_GET_ITEM(storage_V19, 0);
    {Py_XINCREF(py_V19);}
    
            V19 = NULL;
            if (py_V19 == Py_None) {
                // We can either fail here or set V19 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 20;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_20;}
            }
            if (!PyArray_Check(py_V19)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 20;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_20;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V19)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V19;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V19),
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
                {
        __failure = 20;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_20;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V19) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V19));
                {
        __failure = 20;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_20;}
            }
            
        V19 = (PyArrayObject*)(py_V19);
        Py_XINCREF(V19);
        
{

    py_V21 = PyList_GET_ITEM(storage_V21, 0);
    {Py_XINCREF(py_V21);}
    
            V21 = NULL;
            if (py_V21 == Py_None) {
                // We can either fail here or set V21 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 22;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_22;}
            }
            if (!PyArray_Check(py_V21)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 22;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_22;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V21)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V21;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V21),
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
                {
        __failure = 22;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_22;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V21) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V21));
                {
        __failure = 22;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_22;}
            }
            
        V21 = (PyArrayObject*)(py_V21);
        Py_XINCREF(V21);
        
{

    py_V23 = PyList_GET_ITEM(storage_V23, 0);
    {Py_XINCREF(py_V23);}
    
            V23 = NULL;
            if (py_V23 == Py_None) {
                // We can either fail here or set V23 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 24;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_24;}
            }
            if (!PyArray_Check(py_V23)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 24;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_24;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V23)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V23;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V23),
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
                {
        __failure = 24;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_24;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V23) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V23));
                {
        __failure = 24;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_24;}
            }
            
        V23 = (PyArrayObject*)(py_V23);
        Py_XINCREF(V23);
        
{

    py_V25 = PyList_GET_ITEM(storage_V25, 0);
    {Py_XINCREF(py_V25);}
    
            V25 = NULL;
            if (py_V25 == Py_None) {
                // We can either fail here or set V25 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 26;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_26;}
            }
            if (!PyArray_Check(py_V25)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 26;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_26;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V25)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V25;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V25),
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
                {
        __failure = 26;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_26;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V25) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V25));
                {
        __failure = 26;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_26;}
            }
            
        V25 = (PyArrayObject*)(py_V25);
        Py_XINCREF(V25);
        
{

    py_V27 = PyList_GET_ITEM(storage_V27, 0);
    {Py_XINCREF(py_V27);}
    
            V27 = NULL;
            if (py_V27 == Py_None) {
                // We can either fail here or set V27 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 28;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_28;}
            }
            if (!PyArray_Check(py_V27)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 28;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_28;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V27)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V27;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V27),
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
                {
        __failure = 28;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_28;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V27) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V27));
                {
        __failure = 28;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_28;}
            }
            
        V27 = (PyArrayObject*)(py_V27);
        Py_XINCREF(V27);
        
{

    py_V29 = PyList_GET_ITEM(storage_V29, 0);
    {Py_XINCREF(py_V29);}
    
            V29 = NULL;
            if (py_V29 == Py_None) {
                // We can either fail here or set V29 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 30;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_30;}
            }
            if (!PyArray_Check(py_V29)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 30;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_30;}
            }
            // We expect NPY_BOOL
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V29)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V29;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_BOOL), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_BOOL,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V29),
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
                {
        __failure = 30;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_30;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V29) != NPY_BOOL) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_BOOL) got %d",
                             NPY_BOOL, PyArray_TYPE((PyArrayObject*) py_V29));
                {
        __failure = 30;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_30;}
            }
            
        V29 = (PyArrayObject*)(py_V29);
        Py_XINCREF(V29);
        
{

    py_V31 = PyList_GET_ITEM(storage_V31, 0);
    {Py_XINCREF(py_V31);}
    
            V31 = NULL;
            if (py_V31 == Py_None) {
                // We can either fail here or set V31 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 32;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_32;}
            }
            if (!PyArray_Check(py_V31)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 32;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_32;}
            }
            // We expect NPY_BOOL
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V31)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V31;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_BOOL), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_BOOL,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V31),
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
                {
        __failure = 32;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_32;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V31) != NPY_BOOL) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_BOOL) got %d",
                             NPY_BOOL, PyArray_TYPE((PyArrayObject*) py_V31));
                {
        __failure = 32;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_32;}
            }
            
        V31 = (PyArrayObject*)(py_V31);
        Py_XINCREF(V31);
        
{

    py_V33 = PyList_GET_ITEM(storage_V33, 0);
    {Py_XINCREF(py_V33);}
    
            V33 = NULL;
            if (py_V33 == Py_None) {
                // We can either fail here or set V33 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 34;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_34;}
            }
            if (!PyArray_Check(py_V33)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 34;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_34;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V33)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V33;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V33),
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
                {
        __failure = 34;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_34;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V33) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V33));
                {
        __failure = 34;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_34;}
            }
            
        V33 = (PyArrayObject*)(py_V33);
        Py_XINCREF(V33);
        
{

    py_V35 = PyList_GET_ITEM(storage_V35, 0);
    {Py_XINCREF(py_V35);}
    
            V35 = NULL;
            if (py_V35 == Py_None) {
                // We can either fail here or set V35 to NULL and rely on Ops
                // using tensors to handle the NULL case, but if they fail to do so
                // they'll end up with nasty segfaults, so this is public service.
                PyErr_SetString(PyExc_ValueError, "expected an ndarray, not None");
                {
        __failure = 36;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_36;}
            }
            if (!PyArray_Check(py_V35)) {
                PyErr_SetString(PyExc_ValueError, "expected an ndarray");
                {
        __failure = 36;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_36;}
            }
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V35)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V35;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
                             (long int) PyArray_TYPE((PyArrayObject*) py_V35),
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
                {
        __failure = 36;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_36;}
            }
            // This is a TypeError to be consistent with DEBUG_MODE
            // Note: DEBUG_MODE also tells the name of the container
            if (PyArray_TYPE((PyArrayObject*) py_V35) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V35));
                {
        __failure = 36;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_36;}
            }
            
        V35 = (PyArrayObject*)(py_V35);
        Py_XINCREF(V35);
        
{
// Op class Elemwise
npy_float64* V35_iter;
npy_float64* V33_iter;
npy_bool* V31_iter;
npy_bool* V29_iter;
npy_float64* V27_iter;
npy_float64* V25_iter;
npy_float64* V23_iter;
npy_float64* V21_iter;
npy_float64* V19_iter;


npy_float64* V3_iter;

{
    npy_intp dims[0];

    if (!V3) {
        V3 = (PyArrayObject*)PyArray_EMPTY(0,
                                              dims,
                                              NPY_FLOAT64,
                                              0);
    }
    else {
        PyArray_Dims new_dims;
        new_dims.len = 0;
        new_dims.ptr = dims;
        PyObject* success = PyArray_Resize(V3, &new_dims, 0, NPY_CORDER);
        if (!success) {
            // If we can't resize the ndarray we have we can allocate a new one.
            PyErr_Clear();
            Py_XDECREF(V3);
            V3 = (PyArrayObject*)PyArray_EMPTY(0, dims, NPY_FLOAT64, 0);
        } else {
            Py_DECREF(success);
        }
    }
    if (!V3) {
        {
__failure = 37;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_37;}
    }
}
npy_float64* V1_iter;

{
    npy_intp dims[0];

    if (!V1) {
        V1 = (PyArrayObject*)PyArray_EMPTY(0,
                                              dims,
                                              NPY_FLOAT64,
                                              0);
    }
    else {
        PyArray_Dims new_dims;
        new_dims.len = 0;
        new_dims.ptr = dims;
        PyObject* success = PyArray_Resize(V1, &new_dims, 0, NPY_CORDER);
        if (!success) {
            // If we can't resize the ndarray we have we can allocate a new one.
            PyErr_Clear();
            Py_XDECREF(V1);
            V1 = (PyArrayObject*)PyArray_EMPTY(0, dims, NPY_FLOAT64, 0);
        } else {
            Py_DECREF(success);
        }
    }
    if (!V1) {
        {
__failure = 37;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_37;}
    }
}

            if (V17) {
                Py_XDECREF(V17);
            }
            V17 = V27;
            Py_XINCREF(V17);
            
            if (V15) {
                Py_XDECREF(V15);
            }
            V15 = V31;
            Py_XINCREF(V15);
            
            if (V13) {
                Py_XDECREF(V13);
            }
            V13 = V25;
            Py_XINCREF(V13);
            
            if (V11) {
                Py_XDECREF(V11);
            }
            V11 = V23;
            Py_XINCREF(V11);
            
            if (V9) {
                Py_XDECREF(V9);
            }
            V9 = V21;
            Py_XINCREF(V9);
            
            if (V7) {
                Py_XDECREF(V7);
            }
            V7 = V29;
            Py_XINCREF(V7);
            
            if (V5) {
                Py_XDECREF(V5);
            }
            V5 = V19;
            Py_XINCREF(V5);
            

                {
                  #define V17_i V27_i
#define V15_i V31_i
#define V13_i V25_i
#define V11_i V23_i
#define V9_i V21_i
#define V7_i V29_i
#define V5_i V19_i

                  V35_iter = (npy_float64*)(PyArray_DATA(V35));
V33_iter = (npy_float64*)(PyArray_DATA(V33));
V31_iter = (npy_bool*)(PyArray_DATA(V31));
V29_iter = (npy_bool*)(PyArray_DATA(V29));
V27_iter = (npy_float64*)(PyArray_DATA(V27));
V25_iter = (npy_float64*)(PyArray_DATA(V25));
V23_iter = (npy_float64*)(PyArray_DATA(V23));
V21_iter = (npy_float64*)(PyArray_DATA(V21));
V19_iter = (npy_float64*)(PyArray_DATA(V19));
V3_iter = (npy_float64*)(PyArray_DATA(V3));
V1_iter = (npy_float64*)(PyArray_DATA(V1));

                  npy_float64& V35_i = *V35_iter;
npy_float64& V33_i = *V33_iter;
npy_bool& V31_i = *V31_iter;
npy_bool& V29_i = *V29_iter;
npy_float64& V27_i = *V27_iter;
npy_float64& V25_i = *V25_iter;
npy_float64& V23_i = *V23_iter;
npy_float64& V21_i = *V21_iter;
npy_float64& V19_i = *V19_iter;
npy_float64& V3_i = *V3_iter;
npy_float64& V1_i = *V1_iter;

                  {
npy_float64 V37_tmp1;
V37_tmp1 = exp((npy_float64)V33_i);
npy_float64 V37_tmp2;
V37_tmp2 = (6.0) - V19_i;
npy_float64 V37_tmp3;
V37_tmp3 = V37_tmp2 / V37_tmp1;
npy_bool V37_tmp4;
V37_tmp4 = (V37_tmp1 >= (0.0));
npy_float64 V37_tmp5;
V37_tmp5 = exp((npy_float64)V35_i);
npy_float64 V37_tmp6;
V37_tmp6 = (14.0) - V21_i;
npy_float64 V37_tmp7;
V37_tmp7 = V37_tmp6 / V37_tmp5;
npy_bool V37_tmp8;
V37_tmp8 = (V37_tmp5 >= (0.0));
npy_float64 V37_tmp9;
V37_tmp9 = log((npy_float64)V23_i);
npy_float64 V37_tmp10;
V37_tmp10 = V25_i + V37_tmp9;
npy_float64 V37_tmp11;
V37_tmp11 = (3.0) * V37_tmp10;
npy_float64 V37_tmp12;
V37_tmp12 = V37_tmp3 * V37_tmp3;
npy_float64 V37_tmp13;
V37_tmp13 = (-0.5) * V37_tmp12;
npy_float64 V37_tmp14;
V37_tmp14 = (-0.9189385332046727) + V37_tmp13;
npy_float64 V37_tmp15;
V37_tmp15 = V37_tmp14 - V33_i;
npy_bool V37_tmp16;
V37_tmp16 = (V37_tmp1 > (0));
npy_float64 V37_tmp17;
V37_tmp17 = V37_tmp16 ? V37_tmp15 : (-INFINITY);
npy_float64 V37_tmp18;
V37_tmp18 = V37_tmp7 * V37_tmp7;
npy_float64 V37_tmp19;
V37_tmp19 = (-0.5) * V37_tmp18;
npy_float64 V37_tmp20;
V37_tmp20 = (-0.9189385332046727) + V37_tmp19;
npy_float64 V37_tmp21;
V37_tmp21 = V37_tmp20 - V35_i;
npy_bool V37_tmp22;
V37_tmp22 = (V37_tmp5 > (0));
npy_float64 V37_tmp23;
V37_tmp23 = V37_tmp22 ? V37_tmp21 : (-INFINITY);
npy_float64 V37_tmp24;
V37_tmp24 = (0.42857142857142855) * V37_tmp1;
npy_float64 V37_tmp25;
V37_tmp25 = V37_tmp24 * V37_tmp24;
npy_float64 V37_tmp26;
V37_tmp26 = (-0.5) * V37_tmp25;
npy_float64 V37_tmp27;
V37_tmp27 = (-1.0730892130319312) + V37_tmp26;
npy_float64 V37_tmp28;
V37_tmp28 = V37_tmp4 ? V37_tmp27 : (-INFINITY);
npy_float64 V37_tmp29;
V37_tmp29 = (0.42857142857142855) * V37_tmp5;
npy_float64 V37_tmp30;
V37_tmp30 = V37_tmp29 * V37_tmp29;
npy_float64 V37_tmp31;
V37_tmp31 = (-0.5) * V37_tmp30;
npy_float64 V37_tmp32;
V37_tmp32 = (-1.0730892130319312) + V37_tmp31;
npy_float64 V37_tmp33;
V37_tmp33 = V37_tmp8 ? V37_tmp32 : (-INFINITY);
npy_float64 V37_tmp34;
V37_tmp34 = (3.0) * V27_i;
npy_bool V37_tmp35;
V37_tmp35 = (V31_i | V29_i);
npy_float64 V37_tmp36;
V37_tmp36 = V37_tmp35 ? (-INFINITY) : (0.6931471805599453);
npy_float64 V37_tmp37;
V37_tmp37 = (1.0986122886681098) + V37_tmp36 + V37_tmp34 + V37_tmp33 + V35_i + V37_tmp28 + V33_i + V37_tmp23 + V37_tmp17;
npy_float64 V37_tmp38;
V37_tmp38 = V37_tmp37 - V37_tmp11;
V17_i = V37_tmp38;
V15_i = V37_tmp8;
V13_i = V37_tmp7;
V11_i = V37_tmp6;
V9_i = V37_tmp5;
V7_i = V37_tmp4;
V5_i = V37_tmp3;
V3_i = V37_tmp2;
V1_i = V37_tmp1;
}

                  #undef V17_i
#undef V15_i
#undef V13_i
#undef V11_i
#undef V9_i
#undef V7_i
#undef V5_i

                }
                
__label_37:

double __DUMMY_37;

}
__label_36:

        if (V35) {
            Py_XDECREF(V35);
        }
        
    {Py_XDECREF(py_V35);}
    
double __DUMMY_36;

}
__label_34:

        if (V33) {
            Py_XDECREF(V33);
        }
        
    {Py_XDECREF(py_V33);}
    
double __DUMMY_34;

}
__label_32:

        if (V31) {
            Py_XDECREF(V31);
        }
        
    {Py_XDECREF(py_V31);}
    
double __DUMMY_32;

}
__label_30:

        if (V29) {
            Py_XDECREF(V29);
        }
        
    {Py_XDECREF(py_V29);}
    
double __DUMMY_30;

}
__label_28:

        if (V27) {
            Py_XDECREF(V27);
        }
        
    {Py_XDECREF(py_V27);}
    
double __DUMMY_28;

}
__label_26:

        if (V25) {
            Py_XDECREF(V25);
        }
        
    {Py_XDECREF(py_V25);}
    
double __DUMMY_26;

}
__label_24:

        if (V23) {
            Py_XDECREF(V23);
        }
        
    {Py_XDECREF(py_V23);}
    
double __DUMMY_24;

}
__label_22:

        if (V21) {
            Py_XDECREF(V21);
        }
        
    {Py_XDECREF(py_V21);}
    
double __DUMMY_22;

}
__label_20:

        if (V19) {
            Py_XDECREF(V19);
        }
        
    {Py_XDECREF(py_V19);}
    
double __DUMMY_20;

}
__label_18:

    if (!__failure) {
      
        {Py_XDECREF(py_V17);}
        if (!V17) {
            Py_INCREF(Py_None);
            py_V17 = Py_None;
        }
        else if ((void*)py_V17 != (void*)V17) {
            py_V17 = (PyObject*)V17;
        }

        {Py_XINCREF(py_V17);}

        if (V17 && !PyArray_ISALIGNED((PyArrayObject*) py_V17)) {
            PyErr_Format(PyExc_NotImplementedError,
                         "c_sync: expected an aligned array, got non-aligned array of type %ld"
                         " with %ld dimensions, with 3 last dims "
                         "%ld, %ld, %ld"
                         " and 3 last strides %ld %ld, %ld.",
                         (long int) PyArray_TYPE((PyArrayObject*) py_V17),
                         (long int) PyArray_NDIM(V17),
                         (long int) (PyArray_NDIM(V17) >= 3 ?
        PyArray_DIMS(V17)[PyArray_NDIM(V17)-3] : -1),
                         (long int) (PyArray_NDIM(V17) >= 2 ?
        PyArray_DIMS(V17)[PyArray_NDIM(V17)-2] : -1),
                         (long int) (PyArray_NDIM(V17) >= 1 ?
        PyArray_DIMS(V17)[PyArray_NDIM(V17)-1] : -1),
                         (long int) (PyArray_NDIM(V17) >= 3 ?
        PyArray_STRIDES(V17)[PyArray_NDIM(V17)-3] : -1),
                         (long int) (PyArray_NDIM(V17) >= 2 ?
        PyArray_STRIDES(V17)[PyArray_NDIM(V17)-2] : -1),
                         (long int) (PyArray_NDIM(V17) >= 1 ?
        PyArray_STRIDES(V17)[PyArray_NDIM(V17)-1] : -1)
        );
            {
        __failure = 18;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_18;}
        }
        
      PyObject* old = PyList_GET_ITEM(storage_V17, 0);
      {Py_XINCREF(py_V17);}
      PyList_SET_ITEM(storage_V17, 0, py_V17);
      {Py_XDECREF(old);}
    }
    
        if (V17) {
            Py_XDECREF(V17);
        }
        
    {Py_XDECREF(py_V17);}
    
double __DUMMY_18;

}
__label_16:

    if (!__failure) {
      
        {Py_XDECREF(py_V15);}
        if (!V15) {
            Py_INCREF(Py_None);
            py_V15 = Py_None;
        }
        else if ((void*)py_V15 != (void*)V15) {
            py_V15 = (PyObject*)V15;
        }

        {Py_XINCREF(py_V15);}

        if (V15 && !PyArray_ISALIGNED((PyArrayObject*) py_V15)) {
            PyErr_Format(PyExc_NotImplementedError,
                         "c_sync: expected an aligned array, got non-aligned array of type %ld"
                         " with %ld dimensions, with 3 last dims "
                         "%ld, %ld, %ld"
                         " and 3 last strides %ld %ld, %ld.",
                         (long int) PyArray_TYPE((PyArrayObject*) py_V15),
                         (long int) PyArray_NDIM(V15),
                         (long int) (PyArray_NDIM(V15) >= 3 ?
        PyArray_DIMS(V15)[PyArray_NDIM(V15)-3] : -1),
                         (long int) (PyArray_NDIM(V15) >= 2 ?
        PyArray_DIMS(V15)[PyArray_NDIM(V15)-2] : -1),
                         (long int) (PyArray_NDIM(V15) >= 1 ?
        PyArray_DIMS(V15)[PyArray_NDIM(V15)-1] : -1),
                         (long int) (PyArray_NDIM(V15) >= 3 ?
        PyArray_STRIDES(V15)[PyArray_NDIM(V15)-3] : -1),
                         (long int) (PyArray_NDIM(V15) >= 2 ?
        PyArray_STRIDES(V15)[PyArray_NDIM(V15)-2] : -1),
                         (long int) (PyArray_NDIM(V15) >= 1 ?
        PyArray_STRIDES(V15)[PyArray_NDIM(V15)-1] : -1)
        );
            {
        __failure = 16;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_16;}
        }
        
      PyObject* old = PyList_GET_ITEM(storage_V15, 0);
      {Py_XINCREF(py_V15);}
      PyList_SET_ITEM(storage_V15, 0, py_V15);
      {Py_XDECREF(old);}
    }
    
        if (V15) {
            Py_XDECREF(V15);
        }
        
    {Py_XDECREF(py_V15);}
    
double __DUMMY_16;

}
__label_14:

    if (!__failure) {
      
        {Py_XDECREF(py_V13);}
        if (!V13) {
            Py_INCREF(Py_None);
            py_V13 = Py_None;
        }
        else if ((void*)py_V13 != (void*)V13) {
            py_V13 = (PyObject*)V13;
        }

        {Py_XINCREF(py_V13);}

        if (V13 && !PyArray_ISALIGNED((PyArrayObject*) py_V13)) {
            PyErr_Format(PyExc_NotImplementedError,
                         "c_sync: expected an aligned array, got non-aligned array of type %ld"
                         " with %ld dimensions, with 3 last dims "
                         "%ld, %ld, %ld"
                         " and 3 last strides %ld %ld, %ld.",
                         (long int) PyArray_TYPE((PyArrayObject*) py_V13),
                         (long int) PyArray_NDIM(V13),
                         (long int) (PyArray_NDIM(V13) >= 3 ?
        PyArray_DIMS(V13)[PyArray_NDIM(V13)-3] : -1),
                         (long int) (PyArray_NDIM(V13) >= 2 ?
        PyArray_DIMS(V13)[PyArray_NDIM(V13)-2] : -1),
                         (long int) (PyArray_NDIM(V13) >= 1 ?
        PyArray_DIMS(V13)[PyArray_NDIM(V13)-1] : -1),
                         (long int) (PyArray_NDIM(V13) >= 3 ?
        PyArray_STRIDES(V13)[PyArray_NDIM(V13)-3] : -1),
                         (long int) (PyArray_NDIM(V13) >= 2 ?
        PyArray_STRIDES(V13)[PyArray_NDIM(V13)-2] : -1),
                         (long int) (PyArray_NDIM(V13) >= 1 ?
        PyArray_STRIDES(V13)[PyArray_NDIM(V13)-1] : -1)
        );
            {
        __failure = 14;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_14;}
        }
        
      PyObject* old = PyList_GET_ITEM(storage_V13, 0);
      {Py_XINCREF(py_V13);}
      PyList_SET_ITEM(storage_V13, 0, py_V13);
      {Py_XDECREF(old);}
    }
    
        if (V13) {
            Py_XDECREF(V13);
        }
        
    {Py_XDECREF(py_V13);}
    
double __DUMMY_14;

}
__label_12:

    if (!__failure) {
      
        {Py_XDECREF(py_V11);}
        if (!V11) {
            Py_INCREF(Py_None);
            py_V11 = Py_None;
        }
        else if ((void*)py_V11 != (void*)V11) {
            py_V11 = (PyObject*)V11;
        }

        {Py_XINCREF(py_V11);}

        if (V11 && !PyArray_ISALIGNED((PyArrayObject*) py_V11)) {
            PyErr_Format(PyExc_NotImplementedError,
                         "c_sync: expected an aligned array, got non-aligned array of type %ld"
                         " with %ld dimensions, with 3 last dims "
                         "%ld, %ld, %ld"
                         " and 3 last strides %ld %ld, %ld.",
                         (long int) PyArray_TYPE((PyArrayObject*) py_V11),
                         (long int) PyArray_NDIM(V11),
                         (long int) (PyArray_NDIM(V11) >= 3 ?
        PyArray_DIMS(V11)[PyArray_NDIM(V11)-3] : -1),
                         (long int) (PyArray_NDIM(V11) >= 2 ?
        PyArray_DIMS(V11)[PyArray_NDIM(V11)-2] : -1),
                         (long int) (PyArray_NDIM(V11) >= 1 ?
        PyArray_DIMS(V11)[PyArray_NDIM(V11)-1] : -1),
                         (long int) (PyArray_NDIM(V11) >= 3 ?
        PyArray_STRIDES(V11)[PyArray_NDIM(V11)-3] : -1),
                         (long int) (PyArray_NDIM(V11) >= 2 ?
        PyArray_STRIDES(V11)[PyArray_NDIM(V11)-2] : -1),
                         (long int) (PyArray_NDIM(V11) >= 1 ?
        PyArray_STRIDES(V11)[PyArray_NDIM(V11)-1] : -1)
        );
            {
        __failure = 12;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_12;}
        }
        
      PyObject* old = PyList_GET_ITEM(storage_V11, 0);
      {Py_XINCREF(py_V11);}
      PyList_SET_ITEM(storage_V11, 0, py_V11);
      {Py_XDECREF(old);}
    }
    
        if (V11) {
            Py_XDECREF(V11);
        }
        
    {Py_XDECREF(py_V11);}
    
double __DUMMY_12;

}
__label_10:

    if (!__failure) {
      
        {Py_XDECREF(py_V9);}
        if (!V9) {
            Py_INCREF(Py_None);
            py_V9 = Py_None;
        }
        else if ((void*)py_V9 != (void*)V9) {
            py_V9 = (PyObject*)V9;
        }

        {Py_XINCREF(py_V9);}

        if (V9 && !PyArray_ISALIGNED((PyArrayObject*) py_V9)) {
            PyErr_Format(PyExc_NotImplementedError,
                         "c_sync: expected an aligned array, got non-aligned array of type %ld"
                         " with %ld dimensions, with 3 last dims "
                         "%ld, %ld, %ld"
                         " and 3 last strides %ld %ld, %ld.",
                         (long int) PyArray_TYPE((PyArrayObject*) py_V9),
                         (long int) PyArray_NDIM(V9),
                         (long int) (PyArray_NDIM(V9) >= 3 ?
        PyArray_DIMS(V9)[PyArray_NDIM(V9)-3] : -1),
                         (long int) (PyArray_NDIM(V9) >= 2 ?
        PyArray_DIMS(V9)[PyArray_NDIM(V9)-2] : -1),
                         (long int) (PyArray_NDIM(V9) >= 1 ?
        PyArray_DIMS(V9)[PyArray_NDIM(V9)-1] : -1),
                         (long int) (PyArray_NDIM(V9) >= 3 ?
        PyArray_STRIDES(V9)[PyArray_NDIM(V9)-3] : -1),
                         (long int) (PyArray_NDIM(V9) >= 2 ?
        PyArray_STRIDES(V9)[PyArray_NDIM(V9)-2] : -1),
                         (long int) (PyArray_NDIM(V9) >= 1 ?
        PyArray_STRIDES(V9)[PyArray_NDIM(V9)-1] : -1)
        );
            {
        __failure = 10;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_10;}
        }
        
      PyObject* old = PyList_GET_ITEM(storage_V9, 0);
      {Py_XINCREF(py_V9);}
      PyList_SET_ITEM(storage_V9, 0, py_V9);
      {Py_XDECREF(old);}
    }
    
        if (V9) {
            Py_XDECREF(V9);
        }
        
    {Py_XDECREF(py_V9);}
    
double __DUMMY_10;

}
__label_8:

    if (!__failure) {
      
        {Py_XDECREF(py_V7);}
        if (!V7) {
            Py_INCREF(Py_None);
            py_V7 = Py_None;
        }
        else if ((void*)py_V7 != (void*)V7) {
            py_V7 = (PyObject*)V7;
        }

        {Py_XINCREF(py_V7);}

        if (V7 && !PyArray_ISALIGNED((PyArrayObject*) py_V7)) {
            PyErr_Format(PyExc_NotImplementedError,
                         "c_sync: expected an aligned array, got non-aligned array of type %ld"
                         " with %ld dimensions, with 3 last dims "
                         "%ld, %ld, %ld"
                         " and 3 last strides %ld %ld, %ld.",
                         (long int) PyArray_TYPE((PyArrayObject*) py_V7),
                         (long int) PyArray_NDIM(V7),
                         (long int) (PyArray_NDIM(V7) >= 3 ?
        PyArray_DIMS(V7)[PyArray_NDIM(V7)-3] : -1),
                         (long int) (PyArray_NDIM(V7) >= 2 ?
        PyArray_DIMS(V7)[PyArray_NDIM(V7)-2] : -1),
                         (long int) (PyArray_NDIM(V7) >= 1 ?
        PyArray_DIMS(V7)[PyArray_NDIM(V7)-1] : -1),
                         (long int) (PyArray_NDIM(V7) >= 3 ?
        PyArray_STRIDES(V7)[PyArray_NDIM(V7)-3] : -1),
                         (long int) (PyArray_NDIM(V7) >= 2 ?
        PyArray_STRIDES(V7)[PyArray_NDIM(V7)-2] : -1),
                         (long int) (PyArray_NDIM(V7) >= 1 ?
        PyArray_STRIDES(V7)[PyArray_NDIM(V7)-1] : -1)
        );
            {
        __failure = 8;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_8;}
        }
        
      PyObject* old = PyList_GET_ITEM(storage_V7, 0);
      {Py_XINCREF(py_V7);}
      PyList_SET_ITEM(storage_V7, 0, py_V7);
      {Py_XDECREF(old);}
    }
    
        if (V7) {
            Py_XDECREF(V7);
        }
        
    {Py_XDECREF(py_V7);}
    
double __DUMMY_8;

}
__label_6:

    if (!__failure) {
      
        {Py_XDECREF(py_V5);}
        if (!V5) {
            Py_INCREF(Py_None);
            py_V5 = Py_None;
        }
        else if ((void*)py_V5 != (void*)V5) {
            py_V5 = (PyObject*)V5;
        }

        {Py_XINCREF(py_V5);}

        if (V5 && !PyArray_ISALIGNED((PyArrayObject*) py_V5)) {
            PyErr_Format(PyExc_NotImplementedError,
                         "c_sync: expected an aligned array, got non-aligned array of type %ld"
                         " with %ld dimensions, with 3 last dims "
                         "%ld, %ld, %ld"
                         " and 3 last strides %ld %ld, %ld.",
                         (long int) PyArray_TYPE((PyArrayObject*) py_V5),
                         (long int) PyArray_NDIM(V5),
                         (long int) (PyArray_NDIM(V5) >= 3 ?
        PyArray_DIMS(V5)[PyArray_NDIM(V5)-3] : -1),
                         (long int) (PyArray_NDIM(V5) >= 2 ?
        PyArray_DIMS(V5)[PyArray_NDIM(V5)-2] : -1),
                         (long int) (PyArray_NDIM(V5) >= 1 ?
        PyArray_DIMS(V5)[PyArray_NDIM(V5)-1] : -1),
                         (long int) (PyArray_NDIM(V5) >= 3 ?
        PyArray_STRIDES(V5)[PyArray_NDIM(V5)-3] : -1),
                         (long int) (PyArray_NDIM(V5) >= 2 ?
        PyArray_STRIDES(V5)[PyArray_NDIM(V5)-2] : -1),
                         (long int) (PyArray_NDIM(V5) >= 1 ?
        PyArray_STRIDES(V5)[PyArray_NDIM(V5)-1] : -1)
        );
            {
        __failure = 6;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_6;}
        }
        
      PyObject* old = PyList_GET_ITEM(storage_V5, 0);
      {Py_XINCREF(py_V5);}
      PyList_SET_ITEM(storage_V5, 0, py_V5);
      {Py_XDECREF(old);}
    }
    
        if (V5) {
            Py_XDECREF(V5);
        }
        
    {Py_XDECREF(py_V5);}
    
double __DUMMY_6;

}
__label_4:

    if (!__failure) {
      
        {Py_XDECREF(py_V3);}
        if (!V3) {
            Py_INCREF(Py_None);
            py_V3 = Py_None;
        }
        else if ((void*)py_V3 != (void*)V3) {
            py_V3 = (PyObject*)V3;
        }

        {Py_XINCREF(py_V3);}

        if (V3 && !PyArray_ISALIGNED((PyArrayObject*) py_V3)) {
            PyErr_Format(PyExc_NotImplementedError,
                         "c_sync: expected an aligned array, got non-aligned array of type %ld"
                         " with %ld dimensions, with 3 last dims "
                         "%ld, %ld, %ld"
                         " and 3 last strides %ld %ld, %ld.",
                         (long int) PyArray_TYPE((PyArrayObject*) py_V3),
                         (long int) PyArray_NDIM(V3),
                         (long int) (PyArray_NDIM(V3) >= 3 ?
        PyArray_DIMS(V3)[PyArray_NDIM(V3)-3] : -1),
                         (long int) (PyArray_NDIM(V3) >= 2 ?
        PyArray_DIMS(V3)[PyArray_NDIM(V3)-2] : -1),
                         (long int) (PyArray_NDIM(V3) >= 1 ?
        PyArray_DIMS(V3)[PyArray_NDIM(V3)-1] : -1),
                         (long int) (PyArray_NDIM(V3) >= 3 ?
        PyArray_STRIDES(V3)[PyArray_NDIM(V3)-3] : -1),
                         (long int) (PyArray_NDIM(V3) >= 2 ?
        PyArray_STRIDES(V3)[PyArray_NDIM(V3)-2] : -1),
                         (long int) (PyArray_NDIM(V3) >= 1 ?
        PyArray_STRIDES(V3)[PyArray_NDIM(V3)-1] : -1)
        );
            {
        __failure = 4;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_4;}
        }
        
      PyObject* old = PyList_GET_ITEM(storage_V3, 0);
      {Py_XINCREF(py_V3);}
      PyList_SET_ITEM(storage_V3, 0, py_V3);
      {Py_XDECREF(old);}
    }
    
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
    

        static int __struct_compiled_op_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a_executor(__struct_compiled_op_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a *self) {
            return self->run();
        }

        static void __struct_compiled_op_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a_destructor(PyObject *capsule) {
            __struct_compiled_op_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a *self = (__struct_compiled_op_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a *)PyCapsule_GetContext(capsule);
            delete self;
        }
    
//////////////////////
////  Functions
//////////////////////
static PyObject * instantiate(PyObject * self, PyObject *argtuple) {
  assert(PyTuple_Check(argtuple));
  if (19 != PyTuple_Size(argtuple)){ 
     PyErr_Format(PyExc_TypeError, "Wrong number of arguments, expected 19, got %%i", (int)PyTuple_Size(argtuple));
     return NULL;
  }
  __struct_compiled_op_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a* struct_ptr = new __struct_compiled_op_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a();
  if (struct_ptr->init( PyTuple_GET_ITEM(argtuple, 0),PyTuple_GET_ITEM(argtuple, 1),PyTuple_GET_ITEM(argtuple, 2),PyTuple_GET_ITEM(argtuple, 3),PyTuple_GET_ITEM(argtuple, 4),PyTuple_GET_ITEM(argtuple, 5),PyTuple_GET_ITEM(argtuple, 6),PyTuple_GET_ITEM(argtuple, 7),PyTuple_GET_ITEM(argtuple, 8),PyTuple_GET_ITEM(argtuple, 9),PyTuple_GET_ITEM(argtuple, 10),PyTuple_GET_ITEM(argtuple, 11),PyTuple_GET_ITEM(argtuple, 12),PyTuple_GET_ITEM(argtuple, 13),PyTuple_GET_ITEM(argtuple, 14),PyTuple_GET_ITEM(argtuple, 15),PyTuple_GET_ITEM(argtuple, 16),PyTuple_GET_ITEM(argtuple, 17),PyTuple_GET_ITEM(argtuple, 18) ) != 0) {
    delete struct_ptr;
    return NULL;
  }
    PyObject* thunk = PyCapsule_New((void*)(&__struct_compiled_op_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a_executor), NULL, __struct_compiled_op_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a_destructor);
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
  "ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a",
  NULL,
  -1,
  MyMethods,
};

PyMODINIT_FUNC PyInit_ma7da1fdbbf1c2b95f6e5e8051d6a3dbe6122b0429fb2d56e251f828cf1a1ef0a(void) {
   import_array();
   
    PyObject *m = PyModule_Create(&moduledef);
    return m;
}
