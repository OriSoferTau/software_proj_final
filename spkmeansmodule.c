#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "spkmeans.h"

static double** makeVectors(double** vectors, PyObject* pyVectors){
    int i;
    int j;
    int size = PyList_Size(pyVectors);
    PyObject* vectorObject =  PyList_GetItem(pyVectors, 0);
    int vectorSize = PyList_Size(vectorObject);

    vectors =(double**) safe_malloc(size*(sizeof(double*)));

    for(i = 0; i < size; i++ ) {
        vectors[i] = (double *) safe_malloc((vectorSize) * sizeof(double));
    }
    for(i=0;i<size;i++) {
        vectorObject = PyList_GetItem(pyVectors, i);
        for (j = 0; j < vectorSize; j++) {
            vectors[i][j] = PyFloat_AsDouble(PyList_GetItem(vectorObject, j));
        }
    }
    return vectors;
}


static PyObject* mainPy_capi(PyObject *self, PyObject *args){
    PyObject* matrix_obj;
    int num_of_rows;
    int num_of_culs;
    int func_num;
    int k;
    double** matrix = NULL;
    if(!PyArg_ParseTuple(args, "Oiiii", &matrix_obj, &num_of_rows, &num_of_culs, &func_num, &k)){
        return NULL;
    }
    matrix = makeVectors(matrix, matrix_obj);

    return Py_BuildValue("O",  mainPy(matrix, num_of_rows, num_of_culs, func_num, k));


}


static PyObject* fit_capi(PyObject *self, PyObject *args){
    PyObject* vectors_obj;
    PyObject* centroids_obj;
    int k;
    int dim;
    int num_of_vectors;
    int max_iter;
    int eps;
    double **vectors = NULL;
    double **centroids = NULL;
    if(!PyArg_ParseTuple(args, "OOiiiid", &vectors_obj, &centroids_obj, &k, &dim, &num_of_vectors, &max_iter, &eps )){
        return NULL;
    }

    centroids=makeVectorsPlusone(centroids, centroids_obj);
    vectors=makeVectorsPlusone(vectors,vectors_obj);


    return Py_BuildValue("O", fit(vectors, centroids, k, dim, num_of_vectors, max_iter, eps));

}


static PyMethodDef capiMethods[] = {
        {"mainPy",
                (PyCFunction) mainPy_capi,
                METH_VARARGS,
                        PyDoc_STR("doing spkmeans")},
        {"fit",
                    (PyCFunction) fit_capi,
                    METH_VARARGS,
                        PyDoc_STR("doing kmeans")},
        {NULL, NULL, 0, NULL}
};


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "myspk",
        NULL,
        -1,
        capiMethods,
        NULL,NULL,NULL,NULL
};


PyMODINIT_FUNC
PyInit_myspk(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if(!m){
        return NULL;
    }
    return m;
}







