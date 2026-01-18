
        #include "numpy/npy_3kcompat.h"
        #include "pytensor_mod_helper.h"

        extern "C"{
        static PyObject *
        run_cthunk(PyObject *self, PyObject *args)
        {
          PyObject *py_cthunk = NULL;
          if(!PyArg_ParseTuple(args,"O",&py_cthunk))
            return NULL;

          if (!NpyCapsule_Check(py_cthunk)) {
            PyErr_SetString(PyExc_ValueError,
                           "Argument to run_cthunk must be a NpyCapsule.");
            return NULL;
          }
          void * ptr_addr = NpyCapsule_AsVoidPtr(py_cthunk);
          int (*fn)(void*) = (int (*)(void*))(ptr_addr);
          void* it = NpyCapsule_GetDesc(py_cthunk);
          int failure = fn(it);

          return Py_BuildValue("i", failure);
         }
         static PyMethodDef CutilsExtMethods[] = {
            {"run_cthunk",  run_cthunk, METH_VARARGS|METH_KEYWORDS,
             "Run an pytensor cthunk."},
            {NULL, NULL, 0, NULL}        /* Sentinel */
        };
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "cutils_ext",
        NULL,
        -1,
        CutilsExtMethods,
    };

    PyMODINIT_FUNC
    PyInit_cutils_ext(void) {
        return PyModule_Create(&moduledef);
    }
    }
    
