#include <Python.h>
#include "pytensor_mod_helper.h"
#include <math.h>
#include <numpy/arrayobject.h>
#include <numpy/arrayscalars.h>
#include <numpy/npy_math.h>
//////////////////////
////  Support Code
//////////////////////

    namespace {
    struct __struct_compiled_op_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b {
        PyObject* __ERROR;

        PyObject* storage_V9;
PyObject* storage_V7;
PyObject* storage_V5;
PyObject* storage_V3;
PyObject* storage_V1;
        

        __struct_compiled_op_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b() {
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
        ~__struct_compiled_op_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b(void) {
            cleanup();
        }

        int init(PyObject* __ERROR, PyObject* storage_V9, PyObject* storage_V7, PyObject* storage_V5, PyObject* storage_V3, PyObject* storage_V1) {
            Py_XINCREF(storage_V9);
Py_XINCREF(storage_V7);
Py_XINCREF(storage_V5);
Py_XINCREF(storage_V3);
Py_XINCREF(storage_V1);
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
__label_12:

double __DUMMY_12;

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
        
    PyObject* py_V5;
    
        PyArrayObject* V5;
        
    PyObject* py_V7;
    
        PyArrayObject* V7;
        
    PyObject* py_V9;
    
        PyArrayObject* V9;
        
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

    py_V9 = PyList_GET_ITEM(storage_V9, 0);
    {Py_XINCREF(py_V9);}
    
        V9 = (PyArrayObject*)(py_V9);
        Py_XINCREF(V9);
        
{
// Op class Join

        int axis = 0;
        PyArrayObject* arrays[3] = {V7,V5,V3};
        int out_is_valid = V1 != NULL;

        

        if (out_is_valid) {
            // Check if we can reuse output
            npy_intp join_size = 0;
            npy_intp out_shape[1];
            npy_intp *shape = PyArray_SHAPE(arrays[0]);

            for (int i = 0; i < 3; i++) {
                if (PyArray_NDIM(arrays[i]) != 1) {
                    PyErr_SetString(PyExc_ValueError, "Input to join has wrong ndim");
                    {
        __failure = 11;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_11;}
                }

                join_size += PyArray_SHAPE(arrays[i])[axis];

                if (i > 0){
                    for (int j = 0; j < 1; j++) {
                        if ((j != axis) && (PyArray_SHAPE(arrays[i])[j] != shape[j])) {
                            PyErr_SetString(PyExc_ValueError, "Arrays shape must match along non join axis");
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
                }
            }

            memcpy(out_shape, shape, 1 * sizeof(npy_intp));
            out_shape[axis] = join_size;

            for (int i = 0; i < 1; i++) {
                out_is_valid &= (PyArray_SHAPE(V1)[i] == out_shape[i]);
            }
        }

        if (!out_is_valid) {
            // Use PyArray_Concatenate
            Py_XDECREF(V1);
            PyObject* arrays_tuple = PyTuple_New(3);
            Py_INCREF(V7); PyTuple_SetItem(arrays_tuple, 0, (PyObject*)V7);
Py_INCREF(V5); PyTuple_SetItem(arrays_tuple, 1, (PyObject*)V5);
Py_INCREF(V3); PyTuple_SetItem(arrays_tuple, 2, (PyObject*)V3);
            V1 = (PyArrayObject *)PyArray_Concatenate(arrays_tuple, axis);
            Py_DECREF(arrays_tuple);
            if(!V1){
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
        else {
            // Copy the data to the pre-allocated output buffer

            // Create view into output buffer
            PyArrayObject_fields *view;

            // PyArray_NewFromDescr steals a reference to descr, so we need to increase it
            Py_INCREF(PyArray_DESCR(V1));
            view = (PyArrayObject_fields *)PyArray_NewFromDescr(&PyArray_Type,
                                                                  PyArray_DESCR(V1),
                                                                  1,
                                                                  PyArray_SHAPE(arrays[0]),
                                                                  PyArray_STRIDES(V1),
                                                                  PyArray_DATA(V1),
                                                                  NPY_ARRAY_WRITEABLE,
                                                                  NULL);
            if (view == NULL) {
                {
        __failure = 11;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_11;}
            }

            // Copy data into output buffer
            for (int i = 0; i < 3; i++) {
                view->dimensions[axis] = PyArray_SHAPE(arrays[i])[axis];

                if (PyArray_CopyInto((PyArrayObject*)view, arrays[i]) != 0) {
                    Py_DECREF(view);
                    {
        __failure = 11;
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError,
                "Unexpected error in an Op's C code. "
                "No Python exception was set.");
        }
        goto __label_11;}
                }

                view->data += (view->dimensions[axis] * view->strides[axis]);
            }

            Py_DECREF(view);
        }
        __label_11:

double __DUMMY_11;

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
    

        static int __struct_compiled_op_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b_executor(__struct_compiled_op_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b *self) {
            return self->run();
        }

        static void __struct_compiled_op_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b_destructor(PyObject *capsule) {
            __struct_compiled_op_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b *self = (__struct_compiled_op_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b *)PyCapsule_GetContext(capsule);
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
  __struct_compiled_op_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b* struct_ptr = new __struct_compiled_op_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b();
  if (struct_ptr->init( PyTuple_GET_ITEM(argtuple, 0),PyTuple_GET_ITEM(argtuple, 1),PyTuple_GET_ITEM(argtuple, 2),PyTuple_GET_ITEM(argtuple, 3),PyTuple_GET_ITEM(argtuple, 4),PyTuple_GET_ITEM(argtuple, 5) ) != 0) {
    delete struct_ptr;
    return NULL;
  }
    PyObject* thunk = PyCapsule_New((void*)(&__struct_compiled_op_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b_executor), NULL, __struct_compiled_op_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b_destructor);
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
  "mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b",
  NULL,
  -1,
  MyMethods,
};

PyMODINIT_FUNC PyInit_mb633900005814c1d6e63d89f920ece560feebaaf01e53f323fef5d5084fdbb5b(void) {
   import_array();
   
    PyObject *m = PyModule_Create(&moduledef);
    return m;
}
