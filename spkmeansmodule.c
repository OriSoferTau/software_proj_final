#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "spkmeans.h"

void free_eigTup(eigenTuple** egarr,int num_of_rows){
    int i;
    for ( i = 0; i <num_of_rows ; ++i) {
        safe_free(egarr[i]);
    }
    safe_free(egarr);
}


static double** makeVectors(double** vectors, PyObject* pyVectors){
    int i;
    int j;
    int size = PyList_Size(pyVectors);

    PyObject* vectorObject =  PyList_GetItem(pyVectors, 0);
    int vectorSize = PyList_Size(vectorObject);

    vectors = (double**) safe_malloc(size*(sizeof(double*)));
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
PyObject* makePyMatrix(double** matrix, int num_of_rows, int num_of_culs){

    PyObject* centObject;
    PyObject* vectorObject;
    int i;
    int j;

    centObject=PyList_New(0);
    if(centObject==NULL){
        exit_func(0);
    }
    for(i=0; i < num_of_rows; i++){
        vectorObject= PyList_New(0);
        if(vectorObject==NULL){
            exit_func(0);
        }
        for(j=0; j < num_of_culs; j++){
            PyList_Append(vectorObject,PyFloat_FromDouble(matrix[i][j]));
        }
        PyList_Append(centObject,vectorObject);
    }
    return centObject;
}


PyObject* makePyMatrixFromJacobi(jacobiMatrix* jacobi, int num_of_rows, int num_of_culs){
    PyObject* centObject;
    PyObject* vectorObject;
    int i;
    int j;
    centObject=PyList_New(0);
    if(centObject==NULL){
        exit_func(0);
    }
    vectorObject= PyList_New(0);
    if(vectorObject==NULL){
        exit_func(0);
    }
    for (i = 0; i < num_of_culs; ++i) {
        PyList_Append(vectorObject,PyFloat_FromDouble(jacobi->eigenValues[i]));
    }
    PyList_Append(centObject,vectorObject);

    for(i=0; i < num_of_rows; i++){
        vectorObject= PyList_New(0);
        if(vectorObject==NULL){
            exit_func(0);
        }
        for(j=0; j < num_of_culs; j++){
            PyList_Append(vectorObject,PyFloat_FromDouble(jacobi->eigenVectors[i][j]));
        }
        PyList_Append(centObject,vectorObject);
    }
    return centObject;



}



static PyObject* mainPy(double** matrix, int num_of_rows, int num_of_culs, int func_num, int k){
    jacobiMatrix * matJacobi;
    PyObject* centObject;
    double ** adjMatrix;
    double ** diagMatrix;
    double ** lnormMatrix;
    double** T;
    eigenTuple ** egarr;

    if(func_num==0){
        adjMatrix=wam(matrix, num_of_rows, num_of_culs);
        centObject = makePyMatrix(adjMatrix, num_of_rows,num_of_rows);
        free_matrix((void**) adjMatrix, num_of_rows);
        return centObject;
    }
    else if (func_num==1){
        adjMatrix=wam(matrix, num_of_rows, num_of_culs);
        diagMatrix=ddg(adjMatrix, num_of_rows, 0);

        centObject = makePyMatrix(diagMatrix, num_of_rows,num_of_rows);

        free_matrix((void**) adjMatrix, num_of_rows);
        free_matrix((void**) diagMatrix, num_of_rows);
        return centObject;

    }
    else if (func_num==2){
        adjMatrix=wam(matrix, num_of_rows, num_of_culs);
        diagMatrix=ddg(adjMatrix, num_of_rows, 1);
        lnormMatrix = lnorm(diagMatrix, adjMatrix , num_of_rows);

        centObject = makePyMatrix(lnormMatrix, num_of_rows,num_of_rows);

        free_matrix((void**) adjMatrix, num_of_rows);
        free_matrix((void**) diagMatrix, num_of_rows);
        free_matrix((void**) lnormMatrix, num_of_rows);
        return  centObject;
    }
    else if (func_num==3){

        matJacobi = jacobi(matrix, num_of_rows, num_of_culs);

        centObject = makePyMatrixFromJacobi(matJacobi,num_of_rows,num_of_culs);

        free_jacobi(matJacobi, num_of_rows);
        return centObject;
    }
    else /*if (func_num==4)*/{

        adjMatrix=wam(matrix, num_of_rows, num_of_culs);
        /*printf("after wam\n");*/
        diagMatrix=ddg(adjMatrix, num_of_rows, 1);

        lnormMatrix = lnorm(diagMatrix, adjMatrix , num_of_rows);
        matJacobi = jacobi(lnormMatrix, num_of_rows, num_of_rows);  /* ori: changed to num_of_rows! instead num_of_culs*/
        egarr= buildEigenTuple(matJacobi,num_of_rows);
        k=getK(num_of_rows,egarr,k);
        T=getT(matJacobi,num_of_rows,egarr,k);
        centObject=makePyMatrix(T,num_of_rows,k);

        free_matrix((void**) adjMatrix, num_of_rows);
        free_matrix((void**) diagMatrix, num_of_rows);
        free_matrix((void**) lnormMatrix, num_of_rows);
        free_jacobi(matJacobi, num_of_rows);
        free_matrix((void**) egarr, num_of_rows);
        free_matrix((void**) T, num_of_rows);

        return centObject;

    }
}
static double** makeVectorsPlusone(double** vectors, PyObject* pyVectors){
    int i;
    int j;
    int size = PyList_Size(pyVectors);
    PyObject* vectorObject =  PyList_GetItem(pyVectors, 0);
    int vectorSize = PyList_Size(vectorObject);

    vectors =(double**) safe_malloc(size*(sizeof(double*)));

    for(i = 0; i < size; i++ ) {
        vectors[i] = (double *) safe_malloc((vectorSize + 1) * sizeof(double));

    }

    for(i=0;i<size;i++) {
        vectorObject = PyList_GetItem(pyVectors, i);
        for (j = 0; j < vectorSize; j++) {
            vectors[i][j] = PyFloat_AsDouble(PyList_GetItem(vectorObject, j));
        }

        vectors[i][vectorSize] = 0;


    }

    return vectors;
}

PyObject* fit(double** vector_array, double** centroids, int k, int dim, int num_of_vectors, int max_iter, int eps){
    PyObject* centObject;
    PyObject* vectorObject;
    int i;
    int ind;
    int is_converged;

    for ( ind = 0; ind < max_iter; ind++) {
        for (i = 0; i < num_of_vectors;i++ ) {
            calc_norm(vector_array[i],dim,centroids,k);

        }
        /*printf("centroids in c: \n");
        print_vector_array(centroids,k,dim);
        printf("end of centroids in c:\n");*/
        update_centroids(vector_array,k,dim,num_of_vectors,centroids);
        is_converged=check_convergence(vector_array,dim,num_of_vectors,centroids,eps);
        if(is_converged==1){
            break;
        }
    }

    centObject=PyList_New(0);
    if(centObject==NULL){
        exit_func(0);
    }
    for(ind=0;ind<k;ind++){
        vectorObject= PyList_New(0);
        for(i=0;i<dim;i++){
            PyList_Append(vectorObject,PyFloat_FromDouble(centroids[ind][i]));

        }
        PyList_Append(centObject,vectorObject);
    }

    free_matrix((void**) vector_array,num_of_vectors);
    free_matrix((void**) centroids,k);

    return centObject;
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







