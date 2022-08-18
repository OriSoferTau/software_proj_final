#ifndef SOFTWARE_PROJ_FINAL_SPKMEANS_H
#define SOFTWARE_PROJ_FINAL_SPKMEANS_H


typedef struct{
    double** eigenVectors;
    double* eigenValues;
}jacobiMatrix;


typedef struct { /* for pivot matrix */
    double c;
    double s;
    int pivotRow;
    int pivotCul;
}pivoter ;


typedef struct {
    double val;
    int index;
} eigenTuple;


void exit_func(int idx);
void free_matrix(double** vector_array,int num_of_lines);/*free matrix */
void free_mat(char** file_lines,int num_of_lines );
void free_jacobi(jacobiMatrix* jacobi,int num_of_lines);
void* safe_malloc(size_t size);
void* safe_calloc(size_t size);
void is_valid_filename (char* filename);
void check_input_legal (int argc, char** argv);
void print_matrix(double** matrix,int num_of_lines);/* print square matrix */
void print_jacobi(jacobiMatrix* jacobi,int num_of_lines);/* print eigen values and eigen vectors */
int get_num_lines (char *filename);
int get_dim (char *filename);
int* line_lengths (char *filename , int N);
char** file_to_lines (char *filename, int *lines_lens, int num_of_lines);
double** to_vector_array(char** file_lines, int dimention, int num_of_lines);
double** matBuild(int num_of_rows,int num_of_culs);
double calc_exp_norm(double* vector1,double* vector2,int dim); /* calc norm with exponent*/
double** wam(double** matrix,int num_of_rows, int num_of_culs);
double** ddg(double** adjMatrix,int num_of_rows, int flag);
double ** lnorm(double** diagMatrix,double** adjMatrix,int num_of_rows);
void PMatrix(double**P,double** A,double** V,int num_of_rows,int num_of_culs,pivoter* rotator);/* create rotation matrix */
int checkJacobiConverge(double** AtagMatrix,double **A ,int num_of_rows,int num_of_culs);
int buildATag(double** AtagMatrix,double **A,pivoter* rotator, int num_of_rows,int num_of_culs);
jacobiMatrix * jacobi(double** A,int num_of_rows,int num_of_culs);
int cmpfunc (const void * a, const void * b);/* help function for qsort */
int getK(int num_of_rows,eigenTuple** arr,int k); /* get number of clusters */
double** getT(jacobiMatrix* eigens,int num_of_rows,eigenTuple** arr,int k);/* build T matrix */
int check_convergence(double** vector_array,int dim,int num_of_vectors, double** centroids,double epsilon);
void update_centroids(double** vector_array,int k,int dim,int num_of_vectors, double** centroids);
void calc_norm(double* vector,int dim, double**centroids, int k);
eigenTuple** buildEigenTuple(jacobiMatrix* eigens,int num_of_rows);
int main(int argc, char** argv);
void print_vector_array(double** matrix,int num_of_lines,int num_of_culs);/*delete later*/







#endif /*SOFTWARE_PROJ_FINAL_SPKMEANS_H*/
