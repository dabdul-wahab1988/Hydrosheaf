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
    struct __struct_compiled_op_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4 {
        PyObject* __ERROR;

        PyObject* storage_V5;
PyObject* storage_V3;
PyObject* storage_V1;
        

        __struct_compiled_op_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4() {
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
        ~__struct_compiled_op_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4(void) {
            cleanup();
        }

        int init(PyObject* __ERROR, PyObject* storage_V5, PyObject* storage_V3, PyObject* storage_V1) {
            Py_XINCREF(storage_V5);
Py_XINCREF(storage_V3);
Py_XINCREF(storage_V1);
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
__label_8:

double __DUMMY_8;

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
        
{

    py_V5 = PyList_GET_ITEM(storage_V5, 0);
    {Py_XINCREF(py_V5);}
    
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
        
{
// Op class SoftmaxGrad

    PyArrayObject* op[3];
    npy_uint32 op_flags[3];
    npy_uint32 iter_flags;
    NpyIter* iter;
    NpyIter_IterNextFunc* get_next;
    char** data_ptr;

    int sm_ndim = PyArray_NDIM(V3);
    int axis = 0;
    int iterate_axis = !(axis == NPY_RAVEL_AXIS || sm_ndim == 1);

    // Validate inputs
    if ((PyArray_TYPE(V5) != NPY_DOUBLE) &&
        (PyArray_TYPE(V5) != NPY_FLOAT))
    {
        PyErr_SetString(PyExc_TypeError, "types should be float or float64");
        {
__failure = 7;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_7;};
    }
    if ((PyArray_TYPE(V3) != NPY_DOUBLE) &&
        (PyArray_TYPE(V3) != NPY_FLOAT))
    {
        PyErr_SetString(PyExc_TypeError, "types should be float or float64");
        {
__failure = 7;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_7;};
    }

    if (iterate_axis)
    {
        if (axis < 0) axis = sm_ndim + axis;
        if ((axis < 0) || (iterate_axis && (axis > sm_ndim)))
        {
            PyErr_SetString(PyExc_ValueError, "invalid axis in SoftmaxGrad");
            {
__failure = 7;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_7;};
        }
    }
    if ((V1 == NULL)
        || !(PyArray_CompareLists(PyArray_DIMS(V1), PyArray_DIMS(V3), sm_ndim)))
    {
        Py_XDECREF(V1);
        V1 = (PyArrayObject*)PyArray_SimpleNew(sm_ndim,
                                                 PyArray_DIMS(V3),
                                                 PyArray_TYPE(V3));
        if (!V1)
        {
            PyErr_SetString(PyExc_MemoryError, "failed to alloc SoftMaxGrad dx output");
            {
__failure = 7;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_7;};
        }
    }

    // Create numpy iterator
    op[0] = V5;
    op[1] = V3;
    op[2] = V1;
    op_flags[0] = NPY_ITER_READONLY;
    op_flags[1] = NPY_ITER_READONLY;
    op_flags[2] = NPY_ITER_READWRITE;
    iter_flags = (iterate_axis)? NPY_ITER_MULTI_INDEX : 0;
    iter = NpyIter_MultiNew(
        3,
        op,
        iter_flags,
        NPY_KEEPORDER,
        NPY_NO_CASTING,
        op_flags,
        NULL
    );

    if (iter == NULL)
    {
        PyErr_SetString(PyExc_MemoryError, "failed to create softmax iterator");
        {
__failure = 7;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_7;};
    }

    // SoftmaxGrad is applied across the entire array
    if (!iterate_axis)
    {
        get_next = NpyIter_GetIterNext(iter, NULL);
        if (get_next == NULL)
        {
            NpyIter_Deallocate(iter);
            PyErr_SetString(PyExc_RuntimeError, "Failed to obtain SoftMaxGrad GetIterNext");
            {
__failure = 7;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_7;};
        }
        data_ptr = NpyIter_GetDataPtrArray(iter);

        // Compute and accumulate dy * sm
        dtype_V1 sum_dy_times_sm = 0.0;
        do
        {
            dtype_V5* dy_ptr = (dtype_V5*)data_ptr[0];
            dtype_V3* sm_ptr = (dtype_V3*)data_ptr[1];
            dtype_V1* dx_ptr = (dtype_V1*)data_ptr[2];

            *dx_ptr = (dtype_V1)((*dy_ptr) * (*sm_ptr));
            sum_dy_times_sm += *dx_ptr;
        } while(get_next(iter));

        // Reset Iterator
        if (NpyIter_GotoIterIndex(iter, 0) == NPY_FAIL)
        {
            PyErr_SetString(PyExc_RuntimeError, "Failed to reset softmax iterator");
            {
__failure = 7;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_7;};
        }

        // Subtract sum(dy*sm) * sm
        do
        {
            dtype_V3* sm_ptr = (dtype_V3*)data_ptr[1];
            dtype_V1* dx_ptr = (dtype_V1*)data_ptr[2];
            *dx_ptr -= sum_dy_times_sm * ((dtype_V1)(*sm_ptr));
        } while(get_next(iter));
    }

    // SoftmaxGrad is applied across a specific axis
    else {
        // Collect axis strides and remove it from iteration
        npy_intp axis_size = PyArray_DIM(V3, axis);
        npy_intp* axis_stride = NpyIter_GetAxisStrideArray(iter, axis);
        if  (axis_stride == NULL)
        {
            PyErr_SetString(PyExc_RuntimeError, "Failed to obtain softmax axis strides");
            {
__failure = 7;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_7;};
        }
        npy_intp dy_axis_stride = axis_stride[0] / sizeof(dtype_V5);
        npy_intp sm_axis_stride = axis_stride[1] / sizeof(dtype_V3);
        npy_intp dx_axis_stride = axis_stride[2] / sizeof(dtype_V1);

        if (NpyIter_RemoveAxis(iter, axis) == NPY_FAIL)
        {
            PyErr_SetString(PyExc_RuntimeError, "Failed to remove SoftmaxGrad axis from iterator");
            {
__failure = 7;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_7;};
        }

        // Iterate over remaining axes
        get_next = NpyIter_GetIterNext(iter, NULL);
        if (get_next == NULL)
        {
            NpyIter_Deallocate(iter);
            PyErr_SetString(PyExc_RuntimeError, "Failed to obtain SoftamGrad GetIterNext");
            {
__failure = 7;
if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_RuntimeError,
        "Unexpected error in an Op's C code. "
        "No Python exception was set.");
}
goto __label_7;};
        }

        data_ptr = NpyIter_GetDataPtrArray(iter);
        do
        {
            dtype_V5* dy_axis = (dtype_V5*)data_ptr[0];
            dtype_V3* sm_axis = (dtype_V3*)data_ptr[1];
            dtype_V1* dx_axis = (dtype_V1*)data_ptr[2];

            // Compute and accumulate dy * sm
            dtype_V1 sum_dy_times_sm = 0.0;
            for (npy_intp i = 0; i < axis_size; i++)
            {
                dx_axis[i * dx_axis_stride] = (dtype_V1)(dy_axis[i * dy_axis_stride] * sm_axis[i * sm_axis_stride]);
                sum_dy_times_sm += dx_axis[i * dx_axis_stride];
            }

            // Subtract sum(dy*sm) * sm
            for (npy_intp i = 0; i < axis_size; i++)
            {
                dx_axis[i * dx_axis_stride] -= sum_dy_times_sm * (dtype_V1)(sm_axis[i * sm_axis_stride]);
            }

        } while(get_next(iter));
    }
    NpyIter_Deallocate(iter);
__label_7:

double __DUMMY_7;

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
    

        static int __struct_compiled_op_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4_executor(__struct_compiled_op_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4 *self) {
            return self->run();
        }

        static void __struct_compiled_op_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4_destructor(PyObject *capsule) {
            __struct_compiled_op_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4 *self = (__struct_compiled_op_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4 *)PyCapsule_GetContext(capsule);
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
  __struct_compiled_op_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4* struct_ptr = new __struct_compiled_op_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4();
  if (struct_ptr->init( PyTuple_GET_ITEM(argtuple, 0),PyTuple_GET_ITEM(argtuple, 1),PyTuple_GET_ITEM(argtuple, 2),PyTuple_GET_ITEM(argtuple, 3) ) != 0) {
    delete struct_ptr;
    return NULL;
  }
    PyObject* thunk = PyCapsule_New((void*)(&__struct_compiled_op_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4_executor), NULL, __struct_compiled_op_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4_destructor);
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
  "ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4",
  NULL,
  -1,
  MyMethods,
};

PyMODINIT_FUNC PyInit_ma94ac5ab5b989a0f08f6a151b1eced200bd65e8169014a7955f58d8dcfcb0cc4(void) {
   import_array();
   
    PyObject *m = PyModule_Create(&moduledef);
    return m;
}
