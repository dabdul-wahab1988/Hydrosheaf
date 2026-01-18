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
    struct __struct_compiled_op_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea {
        PyObject* __ERROR;

        PyObject* storage_V11;
PyObject* storage_V9;
PyObject* storage_V7;
PyObject* storage_V5;
PyObject* storage_V3;
PyObject* storage_V1;
        

        __struct_compiled_op_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea() {
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
        ~__struct_compiled_op_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea(void) {
            cleanup();
        }

        int init(PyObject* __ERROR, PyObject* storage_V11, PyObject* storage_V9, PyObject* storage_V7, PyObject* storage_V5, PyObject* storage_V3, PyObject* storage_V1) {
            Py_XINCREF(storage_V11);
Py_XINCREF(storage_V9);
Py_XINCREF(storage_V7);
Py_XINCREF(storage_V5);
Py_XINCREF(storage_V3);
Py_XINCREF(storage_V1);
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
__label_14:

double __DUMMY_14;

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
        
            typedef npy_float64 dtype_V7;
            
    PyObject* py_V9;
    
        PyArrayObject* V9;
        
            typedef npy_float64 dtype_V9;
            
    PyObject* py_V11;
    
        PyArrayObject* V11;
        
            typedef npy_float64 dtype_V11;
            
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
            // We expect NPY_FLOAT64
            if (!PyArray_ISALIGNED((PyArrayObject*) py_V7)) {
                PyArrayObject * tmp = (PyArrayObject*) py_V7;
                PyErr_Format(PyExc_NotImplementedError,
                             "expected an aligned array of type %ld "
                             "(NPY_FLOAT64), got non-aligned array of type %ld"
                             " with %ld dimensions, with 3 last dims "
                             "%ld, %ld, %ld"
                             " and 3 last strides %ld %ld, %ld.",
                             (long int) NPY_FLOAT64,
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
            if (PyArray_TYPE((PyArrayObject*) py_V7) != NPY_FLOAT64) {
                PyErr_Format(PyExc_TypeError,
                             "expected type_num %d (NPY_FLOAT64) got %d",
                             NPY_FLOAT64, PyArray_TYPE((PyArrayObject*) py_V7));
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
        
{

    py_V9 = PyList_GET_ITEM(storage_V9, 0);
    {Py_XINCREF(py_V9);}
    
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
        
{

    py_V11 = PyList_GET_ITEM(storage_V11, 0);
    {Py_XINCREF(py_V11);}
    
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
        
{
// Op class Elemwise
npy_float64* V11_iter;
npy_intp V11_n0;
npy_intp V11_stride0;
npy_intp V11_jump0_0;
npy_float64* V9_iter;
npy_intp V9_n0;
npy_intp V9_stride0;
npy_intp V9_jump0_0;
npy_float64* V7_iter;
npy_intp V7_jumpx_0;


if (PyArray_NDIM(V11) < 1) {
    PyErr_SetString(PyExc_ValueError, "Not enough dimensions on input.");
                {
    __failure = 13;
    if (!PyErr_Occurred()) {
        PyErr_SetString(PyExc_RuntimeError,
            "Unexpected error in an Op's C code. "
            "No Python exception was set.");
    }
    goto __label_13;}
}
V11_n0 = PyArray_DIMS(V11)[0];
V11_stride0 = PyArray_STRIDES(V11)[0] / sizeof(npy_float64);
V11_jump0_0 = (V11_stride0) - (0);

if (PyArray_NDIM(V9) < 1) {
    PyErr_SetString(PyExc_ValueError, "Not enough dimensions on input.");
                {
    __failure = 13;
    if (!PyErr_Occurred()) {
        PyErr_SetString(PyExc_RuntimeError,
            "Unexpected error in an Op's C code. "
            "No Python exception was set.");
    }
    goto __label_13;}
}
V9_n0 = PyArray_DIMS(V9)[0];
V9_stride0 = PyArray_STRIDES(V9)[0] / sizeof(npy_float64);
V9_jump0_0 = (V9_stride0) - (0);
V7_jumpx_0 = -(0);

        if (V11_n0 != V9_n0)
        {
            if (V11_n0 == 1 || V9_n0 == 1)
            {
                PyErr_Format(PyExc_ValueError, "Runtime broadcasting not allowed. One input had a distinct dimension length of 1, but was not marked as broadcastable: (input[%i].shape[%i] = %lld, input[%i].shape[%i] = %lld). If broadcasting was intended, use `specify_broadcastable` on the relevant input.",
               0,
               0,
               (long long int) V11_n0,
               1,
               0,
               (long long int) V9_n0
                );
            } else {
                PyErr_Format(PyExc_ValueError, "Input dimension mismatch: (input[%i].shape[%i] = %lld, input[%i].shape[%i] = %lld)",
                   0,
                   0,
                   (long long int) V11_n0,
                   1,
                   0,
                   (long long int) V9_n0
                );
            }
            {
__failure = 13;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_13;}
        }

npy_float64* V5_iter;
npy_intp V5_n0;
npy_intp V5_stride0;
npy_intp V5_jump0_0;

{
    npy_intp dims[1];
    dims[0] = V11_n0;

    if (!V5) {
        V5 = (PyArrayObject*)PyArray_EMPTY(1,
                                              dims,
                                              NPY_FLOAT64,
                                              PyArray_ISFORTRAN(V11) && PyArray_ISFORTRAN(V9));
    }
    else {
        PyArray_Dims new_dims;
        new_dims.len = 1;
        new_dims.ptr = dims;
        PyObject* success = PyArray_Resize(V5, &new_dims, 0, NPY_CORDER);
        if (!success) {
            // If we can't resize the ndarray we have we can allocate a new one.
            PyErr_Clear();
            Py_XDECREF(V5);
            V5 = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_FLOAT64, 0);
        } else {
            Py_DECREF(success);
        }
    }
    if (!V5) {
        {
__failure = 13;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_13;}
    }
}

if (PyArray_NDIM(V5) < 1) {
    PyErr_SetString(PyExc_ValueError, "Not enough dimensions on input.");
                {
    __failure = 13;
    if (!PyErr_Occurred()) {
        PyErr_SetString(PyExc_RuntimeError,
            "Unexpected error in an Op's C code. "
            "No Python exception was set.");
    }
    goto __label_13;}
}
V5_n0 = PyArray_DIMS(V5)[0];
V5_stride0 = PyArray_STRIDES(V5)[0] / sizeof(npy_float64);
V5_jump0_0 = (V5_stride0) - (0);
npy_float64* V3_iter;
npy_intp V3_n0;
npy_intp V3_stride0;
npy_intp V3_jump0_0;

{
    npy_intp dims[1];
    dims[0] = V11_n0;

    if (!V3) {
        V3 = (PyArrayObject*)PyArray_EMPTY(1,
                                              dims,
                                              NPY_FLOAT64,
                                              PyArray_ISFORTRAN(V11) && PyArray_ISFORTRAN(V9));
    }
    else {
        PyArray_Dims new_dims;
        new_dims.len = 1;
        new_dims.ptr = dims;
        PyObject* success = PyArray_Resize(V3, &new_dims, 0, NPY_CORDER);
        if (!success) {
            // If we can't resize the ndarray we have we can allocate a new one.
            PyErr_Clear();
            Py_XDECREF(V3);
            V3 = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_FLOAT64, 0);
        } else {
            Py_DECREF(success);
        }
    }
    if (!V3) {
        {
__failure = 13;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_13;}
    }
}

if (PyArray_NDIM(V3) < 1) {
    PyErr_SetString(PyExc_ValueError, "Not enough dimensions on input.");
                {
    __failure = 13;
    if (!PyErr_Occurred()) {
        PyErr_SetString(PyExc_RuntimeError,
            "Unexpected error in an Op's C code. "
            "No Python exception was set.");
    }
    goto __label_13;}
}
V3_n0 = PyArray_DIMS(V3)[0];
V3_stride0 = PyArray_STRIDES(V3)[0] / sizeof(npy_float64);
V3_jump0_0 = (V3_stride0) - (0);
npy_float64* V1_iter;
npy_intp V1_n0;
npy_intp V1_stride0;
npy_intp V1_jump0_0;

{
    npy_intp dims[1];
    dims[0] = V11_n0;

    if (!V1) {
        V1 = (PyArrayObject*)PyArray_EMPTY(1,
                                              dims,
                                              NPY_FLOAT64,
                                              PyArray_ISFORTRAN(V11) && PyArray_ISFORTRAN(V9));
    }
    else {
        PyArray_Dims new_dims;
        new_dims.len = 1;
        new_dims.ptr = dims;
        PyObject* success = PyArray_Resize(V1, &new_dims, 0, NPY_CORDER);
        if (!success) {
            // If we can't resize the ndarray we have we can allocate a new one.
            PyErr_Clear();
            Py_XDECREF(V1);
            V1 = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_FLOAT64, 0);
        } else {
            Py_DECREF(success);
        }
    }
    if (!V1) {
        {
__failure = 13;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_13;}
    }
}

if (PyArray_NDIM(V1) < 1) {
    PyErr_SetString(PyExc_ValueError, "Not enough dimensions on input.");
                {
    __failure = 13;
    if (!PyErr_Occurred()) {
        PyErr_SetString(PyExc_RuntimeError,
            "Unexpected error in an Op's C code. "
            "No Python exception was set.");
    }
    goto __label_13;}
}
V1_n0 = PyArray_DIMS(V1)[0];
V1_stride0 = PyArray_STRIDES(V1)[0] / sizeof(npy_float64);
V1_jump0_0 = (V1_stride0) - (0);

{
        V11_iter = (npy_float64*)(PyArray_DATA(V11));
V9_iter = (npy_float64*)(PyArray_DATA(V9));
V7_iter = (npy_float64*)(PyArray_DATA(V7));
V5_iter = (npy_float64*)(PyArray_DATA(V5));
V3_iter = (npy_float64*)(PyArray_DATA(V3));
V1_iter = (npy_float64*)(PyArray_DATA(V1));

        for (npy_intp ITER_0 = 0; ITER_0<V1_n0; ITER_0++) {
            npy_float64 &V11_i = * ( V11_iter + ITER_0 * V11_jump0_0 );
npy_float64 &V9_i = * ( V9_iter + ITER_0 * V9_jump0_0 );
npy_float64 &V7_i = * ( V7_iter + ITER_0 * V7_jumpx_0 );
npy_float64 &V5_i = * ( V5_iter + ITER_0 * V5_jump0_0 );
npy_float64 &V3_i = * ( V3_iter + ITER_0 * V3_jump0_0 );
npy_float64 &V1_i = * ( V1_iter + ITER_0 * V1_jump0_0 );

            
        {
            
            {
npy_float64 V13_tmp1;
V13_tmp1 = V9_i + V7_i;
npy_float64 V13_tmp2;
V13_tmp2 = V11_i - V13_tmp1;
npy_float64 V13_tmp3;
V13_tmp3 = (0.04000000000000001) * V13_tmp2;
npy_float64 V13_tmp4;
V13_tmp4 = (0.2) * V13_tmp2;
npy_float64 V13_tmp5;
V13_tmp5 = V13_tmp4 * V13_tmp4;
npy_float64 V13_tmp6;
V13_tmp6 = (-0.5) * V13_tmp5;
npy_float64 V13_tmp7;
V13_tmp7 = (-2.5283764757095555) + V13_tmp6;
V5_i = V13_tmp7;
V3_i = V13_tmp2;
V1_i = V13_tmp3;
}

            
        }
        
        }
        }
__label_13:

double __DUMMY_13;

}
__label_12:

        if (V11) {
            Py_XDECREF(V11);
        }
        
    {Py_XDECREF(py_V11);}
    
double __DUMMY_12;

}
__label_10:

        if (V9) {
            Py_XDECREF(V9);
        }
        
    {Py_XDECREF(py_V9);}
    
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
    

        static int __struct_compiled_op_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea_executor(__struct_compiled_op_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea *self) {
            return self->run();
        }

        static void __struct_compiled_op_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea_destructor(PyObject *capsule) {
            __struct_compiled_op_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea *self = (__struct_compiled_op_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea *)PyCapsule_GetContext(capsule);
            delete self;
        }
    
//////////////////////
////  Functions
//////////////////////
static PyObject * instantiate(PyObject * self, PyObject *argtuple) {
  assert(PyTuple_Check(argtuple));
  if (7 != PyTuple_Size(argtuple)){ 
     PyErr_Format(PyExc_TypeError, "Wrong number of arguments, expected 7, got %%i", (int)PyTuple_Size(argtuple));
     return NULL;
  }
  __struct_compiled_op_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea* struct_ptr = new __struct_compiled_op_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea();
  if (struct_ptr->init( PyTuple_GET_ITEM(argtuple, 0),PyTuple_GET_ITEM(argtuple, 1),PyTuple_GET_ITEM(argtuple, 2),PyTuple_GET_ITEM(argtuple, 3),PyTuple_GET_ITEM(argtuple, 4),PyTuple_GET_ITEM(argtuple, 5),PyTuple_GET_ITEM(argtuple, 6) ) != 0) {
    delete struct_ptr;
    return NULL;
  }
    PyObject* thunk = PyCapsule_New((void*)(&__struct_compiled_op_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea_executor), NULL, __struct_compiled_op_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea_destructor);
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
  "mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea",
  NULL,
  -1,
  MyMethods,
};

PyMODINIT_FUNC PyInit_mc9032a477d4a29bc8482480523f4d121361c68a592381e47740524ba4b4e7eea(void) {
   import_array();
   
    PyObject *m = PyModule_Create(&moduledef);
    return m;
}
