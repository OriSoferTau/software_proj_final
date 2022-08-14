/*#define PY_SSIZE_T_CLEAN*/
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>


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

void exit_func(int idx){
    if(idx == 1){
        printf("Invalid Input!\n");
    }

    if (idx == 0){
        printf("An Error Has Occurred\n");
    }
    exit(1);
}

void free_matrix(double** vector_array,int num_of_lines){/*free matrix */
    int i;
    for (i = 0; i < num_of_lines;i++) {
        free(vector_array[i]);

    }
    free(vector_array);
}
void free_mat(char** file_lines,int num_of_lines ){
    int i;
    for (i = 0; i < num_of_lines;i++) {
        free(file_lines[i]);
    }
    free(file_lines);
}

void* safe_malloc(size_t size) {
    void * ptr = malloc(size);
    if(ptr == NULL){
        exit_func(0);
    }
    return ptr;
}

void* safe_calloc(size_t size) {
    void * ptr = calloc(size,sizeof(double));
    if(ptr == NULL){
        exit_func(0);
    }
    return ptr;
}

void is_valid_filename (char* filename){
    FILE *fp = fopen(filename, "r");
    if (fp == NULL){
        exit_func(1);
    }
    fclose(fp);
}

void check_input_legal (int argc, char** argv){
    int flag;
    if(argc!=3){
    }
        exit_func(1);

    if((strcmp(argv[1],"wam"))!=0 && strcmp(argv[1],"ddg")!=0 && strcmp(argv[1],"lnorm")!=0 && strcmp(argv[1],"jacobi")!=0){
        exit_func(1);
    }
    is_valid_filename(argv[2]);
}
int get_num_lines (char *filename) {
    int count_lines;
    int chr;
    FILE *fp = fopen(filename, "r");
    if (fp == NULL){
        exit_func(1);
    }
    count_lines = 0;
    chr = fgetc(fp);
    while (!feof(fp)){
        if(chr == '\n'){
            count_lines++;
        }
        chr = fgetc(fp);
    }
    fclose(fp);
    return count_lines;
}
int get_dim (char *filename) {
    int count_commas;
    int chr;
    FILE *fp = fopen(filename, "r");
    if (fp == NULL){
        exit_func(1);
    }
    count_commas = 1;
    chr = fgetc(fp);
    while (chr != '\n'){
        if(chr == ','){
            count_commas++;
        }
        chr = fgetc(fp);
    }
    fclose(fp);
    return count_commas;
}

int* line_lengths (char *filename , int N){
    int* line_lens;
    int count;
    int index;
    char chr;
    FILE *fp = fopen(filename, "r");
    if (fp == NULL){
        exit_func(1);
    }
    line_lens = safe_malloc(N*sizeof(int));

    if(line_lens==NULL){
        exit_func(0);
    }
    index=0;
    chr = fgetc(fp);
    while (chr != EOF){
        count++;
        if(chr == '\n'){
            line_lens[index]=count;
            index++;
            count=0;
        }
        chr = fgetc(fp);
    }
    fclose(fp);
    return line_lens;
}

char** file_to_lines (char *filename, int *lines_lens, int num_of_lines){
    FILE *fp;
    int i;
    char** lines_array = safe_malloc(sizeof(char*) * num_of_lines);
    if(lines_array == NULL){
        exit_func(0);
    }
    fp = fopen(filename, "r");
    if (fp == NULL){
        exit_func(1);
    }
    for (i = 0; i < num_of_lines; ++i) {
        lines_array[i] = safe_malloc((lines_lens[i]+1)*sizeof(char));
        fscanf(fp, "%s", lines_array[i]);
    }
    return lines_array;
}
double** to_vector_array(char** file_lines, int dimention, int num_of_lines)
{
    int j;
    int i;
    char* ptr;
    double** vector_array = (double **) safe_malloc(num_of_lines * sizeof (double *));
    for (i = 0; i < num_of_lines; ++i) {
        vector_array[i] = (double *) safe_malloc((dimention) * sizeof (double ));
    }

    for (i = 0; i < num_of_lines; i++) {
        vector_array[i][0]=strtod(file_lines[i], &ptr);
        ptr++;

        for (j = 1; j < dimention; j++) {
            vector_array[i][j]=strtod(ptr, &ptr);
            ptr++;
        }
        vector_array[i][dimention]=0;
    }
    return vector_array;
}

void matBuild(double** matrix,int num_of_rows,int num_of_culs) {
    int i;
    matrix = (double **) safe_malloc(sizeof(double *) * num_of_rows);
    for (i = 0; i < num_of_rows; ++i) {
        matrix[i] = (double *) safe_calloc(num_of_culs);
    }
}

double calc_exp_norm(double* vector1,double* vector2,int dim){ /* calc norm with exponent*/
    int i;
    double norm=0;
    for (i = 0; i <dim ; ++i) {
        norm+=pow((vector1[i]-vector2[i]),2);
    }
    norm=sqrt(norm);
    norm=-norm/2;
    norm=exp(norm);
    return norm;
}

double** wam(double** matrix,int num_of_rows, int num_of_culs){
    int i;
    int j;
    double norm;
    double **adjMatrix = NULL;
    matBuild(adjMatrix,num_of_rows,num_of_rows);
    for(i=0;i<num_of_rows; ++i){
        matrix[i][i]=0;
        for ( j =i+1; j <num_of_rows ; ++j) {
            norm=calc_exp_norm(matrix[i],matrix[j],num_of_culs);
            adjMatrix[i][j]=norm;
            adjMatrix[j][i]=norm;

        }
    }
    return adjMatrix;

}


double** ddg(double** adjMatrix,int num_of_rows, int flag){
    int i;
    int j;
    double sum;
    double** diagMatrix=NULL;
    matBuild(diagMatrix,num_of_rows,num_of_rows);
    for ( i = 0; i <num_of_rows ; ++i) {
        sum=0;
        for ( j = 0; j <num_of_rows ; ++j) {
            sum+=adjMatrix[i][j];
            diagMatrix[i][j]=0;
        }
        if(flag==1){
            sum= 1/sqrt(sum);
        }
        diagMatrix[i][i]=sum;
    }
    return diagMatrix;
}

double ** lnorm(double** diagMatrix,double** adjMatrix,int num_of_rows){
   int i;
   int j;
   double** lnormMatrix=NULL;
   matBuild(lnormMatrix,num_of_rows,num_of_rows);
    for ( i = 0; i < num_of_rows; ++i) {

        for ( j = 0; j < num_of_rows; ++j) {
            lnormMatrix[i][j]=1-diagMatrix[i][i]*adjMatrix[i][j]*diagMatrix[j][j];/* each element multiplied by i diag and j diag*/
        }
    }
    return lnormMatrix;
}

void PMatrix(double**P,double** A,double** V,int num_of_rows,int num_of_culs,pivoter* rotator) {/* create rotation matrix */
        int i;
        int j;
        int k;
        double max;
        int pivotRow;
        int pivotCul;
        double theta;
        double c;
        double t;
        double s;
        int sign;
        double temp;
        double **tempMatrix = NULL;
        max = 0;
        for (i = 0; i < num_of_rows; ++i) { /* finding max element */
            for (j = 0; j < num_of_culs; ++j) {
                if (i == j) {
                    continue;
                }
                if (fabs(A[i][j]) >= max) {
                    max = A[i][j];
                    pivotRow = i;
                    pivotCul = j;
                }
            }
        }

        theta = (A[pivotCul][pivotCul] - A[pivotRow][pivotRow]) / (2 * A[pivotRow][pivotCul]);
        sign = theta >= 0 ? 1 : 0;
        t = sign / (fabs(theta) + sqrt(pow(theta, 2) + 1));
        c = 1 / (sqrt(pow(t, 2) + 1));
        s = t * c;
        P[pivotRow][pivotRow] = c;
        P[pivotCul][pivotCul] = c;
        P[pivotRow][pivotCul] = s;
        P[pivotCul][pivotRow] = -s;
        rotator->c = c;
        rotator->s = s;
        rotator->pivotRow = pivotRow;
        rotator->pivotCul = pivotCul;
        matBuild(tempMatrix, num_of_rows, num_of_culs);
        for (i = 0; i < num_of_rows; ++i) { /*V=V*P_s */
            for (j = 0; j < num_of_rows; ++j) {
                temp = 0;
                for (k = 0; k < num_of_rows; ++k) {
                    temp += V[i][k] * P[k][j];
                }
                tempMatrix[i][j] = temp;
            }

        }
        for (i = 0; i < num_of_rows; ++i) {
            for (j = 0; j < num_of_rows; ++j) {
                V[i][j] = tempMatrix[i][j];
            }

        }
        P[pivotRow][pivotRow] = 1;
        P[pivotCul][pivotCul] = 1;
        P[pivotRow][pivotCul] = 0;
        P[pivotCul][pivotRow] = 0;
        free_matrix(tempMatrix, num_of_rows);


    }

int checkJacobiConverge(double** AtagMatrix,double **A ,int num_of_rows,int num_of_culs){
    double fA;
    double fATag;
    int i;
    int j;
    double epsilon=0.00001;
    fA=0;
    fATag=0;
    for ( i = 0; i <num_of_rows ; ++i) {
        for (j = 0; j <num_of_culs ; ++j) {
            if(i==j){
                continue;
            }
            fA+=pow(A[i][j],2);
            fATag+=pow(AtagMatrix[i][j],2);
        }

    }
    if(fA-fATag<=epsilon){
        return 1;
    }
    return 0;
}
int buildATag(double** AtagMatrix,double **A,pivoter* rotator, int num_of_rows,int num_of_culs){
    int i;
    int j;
    int flag;
    for ( i = 0; i <num_of_rows; ++i) {
        for (j = 0; j <num_of_culs ; ++j) {
            if(i!=rotator->pivotRow && i!=rotator->pivotCul && j==rotator->pivotRow){
                AtagMatrix[i][j]=rotator->c*A[i][j]-rotator->s*A[i][rotator->pivotCul];
            }
            else if(i!=rotator->pivotRow && i!=rotator->pivotCul && j==rotator->pivotCul){
                AtagMatrix[i][j]=rotator->c*A[i][j]+rotator->s*A[i][rotator->pivotRow];

            }
            else if(i==rotator->pivotRow && j==rotator->pivotRow){
                AtagMatrix[i][j]=pow(rotator->c,2)*A[i][j]+ pow(rotator->s,2)*A[rotator->pivotCul][rotator->pivotCul]
                                 -(2*rotator->s*rotator->c*A[i][rotator->pivotCul]);
            }
            else if(i==rotator->pivotCul && j==rotator->pivotCul){
                AtagMatrix[i][j]=pow(rotator->s,2)*A[rotator->pivotRow][rotator->pivotRow]+pow(rotator->c,2)*A[i][j]
                                 +(2*rotator->s*rotator->c*A[rotator->pivotRow][j]);
            }
            else if(i==rotator->pivotRow && j==rotator->pivotCul){
                AtagMatrix[i][j]=0;
            }
            else{
                AtagMatrix[i][j]=A[i][j];
            }
        }
    }
    flag=checkJacobiConverge(AtagMatrix,A,num_of_rows,num_of_culs);
    for ( i = 0; i <num_of_rows; ++i) { /*change A to Atag*/
        for (j = 0; j <num_of_culs ; ++j) {
            A[i][j]=AtagMatrix[i][j];
        }

    }


return flag;
}


jacobiMatrix * jacobi(double** A,int num_of_rows,int num_of_culs){
    int pivotRow;/* sent to P as pointer in order to remember pivotRow,pivotCul */
    int pivotCul;
    int i;
    double** P;
    int flag;
    pivoter* rotator;
    double **ATag = NULL;
    jacobiMatrix* ret= (jacobiMatrix *) safe_malloc(sizeof(jacobiMatrix));
    rotator=(pivoter*) safe_malloc(sizeof(pivoter));
    ret->eigenVectors=(double**) safe_malloc(sizeof(double*)*num_of_rows);/*create V matrix*/
    for ( i = 0; i <num_of_rows; ++i) {
        ret->eigenVectors[i]=(double *) safe_calloc(num_of_rows);
        ret->eigenVectors[i][i]=1;
    }
    P=(double **) safe_malloc(sizeof(double*)*num_of_rows);
    for ( i = 0; i <num_of_rows ; ++i) {
        P[i]= (double*) safe_calloc(num_of_culs);
        P[i][i]=1;
    }
    ret->eigenValues=(double *) safe_calloc(num_of_rows);
    matBuild(ATag,num_of_rows,num_of_culs);
    for ( i = 0; i <100; ++i) {
        PMatrix(P,A,ret->eigenVectors,num_of_rows,num_of_culs,rotator);
        flag=buildATag(ATag,A,rotator,num_of_rows,num_of_culs);
        if(flag==1){
            break;
        }
    }
    for (i = 0; i <num_of_rows; ++i) {
        ret->eigenValues[i]=A[i][i];

    }
    return ret;
}
int cmpfunc (const void * a, const void * b) {/* help function for qsort */
    eigenTuple* one=(eigenTuple*) a;
    eigenTuple* two=(eigenTuple*) b;
    if (one->val>two->val){
        return -1;
    }
    if(one->val<two->val){
        return 1;
    }
        return one->index-two->index;


}

int getK(jacobiMatrix* eigens,int num_of_rows,int num_of_culs,eigenTuple** arr,k) /* get number of clusters */{
    int i;
    int j;
    double temp;
    double max;
    int indexMax;
    arr=(eigenTuple**) safe_malloc(sizeof(eigenTuple*)*num_of_rows);
    for ( i = 0; i <num_of_rows; ++i) {
        arr[i]=(eigenTuple*) safe_malloc(sizeof(eigenTuple));
        arr[i]->val=eigens->eigenValues[i];
        arr[i]->index=i;
    }
    qsort(arr, num_of_rows, sizeof(eigenTuple*), cmpfunc);
    if (k>0){
        return k;
    }
    max=-1;
    for ( i = 0; i <num_of_rows/2; ++i) {
        temp=arr[i]->val-arr[i+1]->val;
        if(temp>max){
            max=temp;
            indexMax=i;
        }

    }
    return indexMax;

}
double** getT(jacobiMatrix* eigens,int num_of_rows,int num_of_culs,eigenTuple** arr,int k){/* build T matrix */
    double** T=NULL;
    int i;
    int j;
    int index;
    int norm;
    matBuild(T,num_of_rows,k);
    for (j = 0; j < k; ++j) {
        index=arr[j]->index;
        for (i = 0; i <num_of_rows; ++i) {
            T[i][j]=eigens->eigenVectors[i][index];
        }

    }
    for (j = 0; j <k ; ++j) {
        norm=0;
        for (i = 0; i <num_of_rows; ++i) {
            norm+=pow(T[i][j],2);

        }
        norm=sqrt(norm);
        for (i = 0; i <num_of_rows; ++i) {
            T[i][j]=T[i][j]/norm;

        }
    }
    return T;
}
int main(int argc, char** argv){
    jacobiMatrix * matJacobi;
    double** vector_array;
    double ** matrix;
    double ** adjMatrix;
    double ** diagMatrix;
    double ** lnormMatrix;
    int num_of_lines;
    int dimention;
    int* lines_lens;
    char** file_lines;
    check_input_legal(argc,argv);
    num_of_lines = get_num_lines(argv[2]);
    dimention= get_dim(argv[2]);
    lines_lens= line_lengths(argv[2],num_of_lines);
    file_lines= file_to_lines(argv[2],lines_lens, num_of_lines);
    vector_array= to_vector_array(file_lines,dimention,num_of_lines);
    free(lines_lens);
    free_mat(file_lines,num_of_lines);/*free string matrix*/
    if(strcmp(argv[1],"wam")==0){
        adjMatrix=wam(vector_array,num_of_lines,dimention)
        print_matrix(matrix);/*to add */
        free_matrix(matrix)
    }
    else if (strcmp(argv[1],"ddg")==0){
        adjMatrix=wam(vector_array,num_of_lines,dimention);
        diagMatrix=ddg(matrix,num_of_lines,0);
        print_matrix(diagMatrix);
        free_matrix(matrix);
        free_matrix(diagMatrix);

    }
    else if (strcmp(argv[1],"lnorm")==0){
        adjMatrix=wam(vector_array,num_of_lines,dimention);
        diagMatrix=ddg(matrix,num_of_lines,1);
        matrix = lnorm(diagMatrix,matrix,num_of_lines);
        print_matrix(matrix);
        free(diagMatrix);
        free_matrix(matrix);
    }
    else if (strcmp(argv[1],"jacobi")==0){
        adjMatrix=wam(vector_array,num_of_lines,dimention);
        diagMatrix=ddg(matrix,num_of_lines,1);
        lnormMatrix = lnorm(diagMatrix,matrix,num_of_lines);
        matJacobi = jacobi(lnormMatrix,num_of_lines,dimention);
        print_jacobi(matJacobi);
        /*free everything */
    }



}









