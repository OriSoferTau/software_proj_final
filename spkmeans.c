
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "spkmeans.h"




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
void free_jacobi(jacobiMatrix* jacobi,int num_of_lines){
    free(jacobi->eigenValues);
    free_matrix(jacobi->eigenVectors,num_of_lines);
    free(jacobi);
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
    if(argc!=3){

        exit_func(1);
    }


    if((strcmp(argv[1],"wam"))!=0 && strcmp(argv[1],"ddg")!=0 && strcmp(argv[1],"lnorm")!=0 && strcmp(argv[1],"jacobi")!=0){
        exit_func(1);
    }
    is_valid_filename(argv[2]);
}
void print_matrix(double** matrix,int num_of_lines){/* print square matrix */
    int i;
    int j;
    for (i = 0; i < num_of_lines; ++i) {
        for (j = 0; j < num_of_lines ; ++j) {
            if(num_of_lines-1==j){
                printf("%.4f",matrix[i][j]);
            }
            else{
                printf("%.4f,",matrix[i][j]);
            }
        }
        printf("\n");
    }

}
void print_vector_array(double** matrix,int num_of_lines,int num_of_culs){/* delete later */
    int i;
    int j;
    for (i = 0; i < num_of_lines; ++i) {
        for (j = 0; j < num_of_culs ; ++j) {


            printf("%.4f,",matrix[i][j]);

        }
        printf("\n");
    }

}
void print_jacobi(jacobiMatrix* jacobi,int num_of_lines){/* print eigen values and eigen vectors */
    int i;
    for (i = 0; i <num_of_lines ; ++i) {
        if(num_of_lines-1==i){
            printf("%.4f",jacobi->eigenValues[i]);
        }
        else{
            printf("%.4f,",jacobi->eigenValues[i]);

        }
    }
    printf("\n");
    print_matrix(jacobi->eigenVectors,num_of_lines);
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
    printf("count_lines: %d\n",count_lines);
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
    printf("count_commas: %d\n",count_commas);
    return count_commas;
}

int* line_lengths (char *filename , int N){
    int* line_lens;
    int index;
    char chr;
    int count=0;
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
    }
    print_matrix(vector_array,num_of_lines);
    return vector_array;
}

double** matBuild(int num_of_rows,int num_of_culs) {
    double** matrix;
    int i;
    matrix = (double **) safe_malloc(sizeof(double *) * num_of_rows);
    for (i = 0; i < num_of_rows; ++i) {
        matrix[i] = (double *) safe_calloc(num_of_culs);
    }
    return  matrix;
}

double calc_exp_norm(double* vector1,double* vector2,int dim){ /* calc norm with exponent*/
    int i;
    double norm=0;
    printf("dim: %d\n",dim);
    for (i = 0; i <dim ; ++i) {
        norm=norm+pow((vector1[i]-vector2[i]),2);
        printf("vector1[i]: %f vector[2][i]: %f i: %d norm: %f\n", vector1[i],vector2[i],i,norm);
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
    printf("this is matrix: \n");
    print_vector_array(matrix,num_of_rows,num_of_culs);
    printf("end of matrix");
    printf("rows: %d culs: %d\n",num_of_rows,num_of_culs);
    adjMatrix = matBuild(num_of_rows,num_of_rows);
    for(i=0;i<num_of_rows; ++i){
        adjMatrix[i][i]=0;
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
    diagMatrix = matBuild(num_of_rows,num_of_rows);
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
    lnormMatrix=matBuild(num_of_rows,num_of_rows);
    for ( i = 0; i < num_of_rows; ++i) {

        for ( j = 0; j < num_of_rows; ++j) {
            if(i==j){
                lnormMatrix[i][j]=1-diagMatrix[i][i]*adjMatrix[i][j]*diagMatrix[j][j];/* each element multiplied by i diag and j diag*/
            }
            else{
                lnormMatrix[i][j]=-diagMatrix[i][i]*adjMatrix[i][j]*diagMatrix[j][j];
            }
        }
    }
    printf("lnorm matrix: \n");
    print_matrix(lnormMatrix,num_of_rows);
    printf("end of lnorm matrix\n");
    return lnormMatrix;
}

void PMatrix(double**P,double** A,double** V,int num_of_rows,int num_of_culs,pivoter* rotator) {/* create rotation matrix */
        int i;
        int j;
        int k;
        double max;
        double theta;
        double c;
        double t;
        double s;
        int sign;
        double temp;
        double **tempMatrix = NULL;
        int pivotRow=0;
        int pivotCul=0;
        max = 0;

        for (i = 0; i < num_of_rows; ++i) { /* finding max element */
            for (j = 0; j < num_of_culs; ++j) {
                if (i == j) {
                    continue;
                }
                if (fabs(A[i][j]) >= max) {
                    max = fabs(A[i][j]);
                    pivotRow = i;
                    pivotCul = j;
                }
            }
        }
        theta = (A[pivotCul][pivotCul] - A[pivotRow][pivotRow]) / (2 * A[pivotRow][pivotCul]);
        if(theta>=0){
            sign=1;
        }
        else{
            sign=0;
        }
        t = sign / (fabs(theta) + sqrt(pow(theta, 2) + 1));
        c = 1 / (sqrt(pow(t, 2) + 1));
        s = t * c;
        P[pivotRow][pivotRow] = c;
        P[pivotCul][pivotCul] = c;
        P[pivotRow][pivotCul] = s;
        P[pivotCul][pivotRow] = -s;
        printf("A[pivotCul][pivotCul]: %f A[pivotRow][pivotRow] %f A[pivotRow][pivotCul]: %f\n",A[pivotCul][pivotCul],A[pivotRow][pivotRow],A[pivotRow][pivotCul] );
        printf("this is c: %f s: %f theta: %f sign: %d t: %f \n",c,s,theta,sign,t);
        rotator->c = c;
        rotator->s = s;
        rotator->pivotRow = pivotRow;
        rotator->pivotCul = pivotCul;
        tempMatrix=matBuild( num_of_rows, num_of_culs);
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
                AtagMatrix[i][j]=(rotator->c*A[i][j])-rotator->s*A[i][rotator->pivotCul];
            }
            else if(j!=rotator->pivotRow && j!=rotator->pivotCul && i==rotator->pivotRow){/*symetry*/
                AtagMatrix[i][j]=A[j][i];
            }
            else if(i!=rotator->pivotRow && i!=rotator->pivotCul && j==rotator->pivotCul){
                AtagMatrix[i][j]=(rotator->c*A[i][j])+rotator->s*A[i][rotator->pivotRow];
            }
            else if(j!=rotator->pivotRow && j!=rotator->pivotCul && i==rotator->pivotCul){/*symetry*/
                AtagMatrix[i][j]=AtagMatrix[j][i];
            }
            else if(i==rotator->pivotRow && j==rotator->pivotRow){
                AtagMatrix[i][j]=(pow(rotator->c,2)*A[i][i])+ (pow(rotator->s,2)*A[rotator->pivotCul][rotator->pivotCul])
                                 -(2*rotator->s*rotator->c*A[i][rotator->pivotCul]);
            }
            else if(i==rotator->pivotCul && j==rotator->pivotCul){
                AtagMatrix[i][j]=(pow(rotator->s,2)*A[rotator->pivotRow][rotator->pivotRow])+(pow(rotator->c,2)*A[j][j])
                                 +(2*rotator->s*rotator->c*A[rotator->pivotRow][j]);
            }
            else if(i==rotator->pivotRow && j==rotator->pivotCul){
                AtagMatrix[i][j]=((pow(rotator->c,2)-pow(rotator->s,2))
                        *A[i][j])+((rotator->s*rotator->c)*
                        (A[i][i]-A[j][j]));
                printf("AtagMatrix[i][j]: %f\n",AtagMatrix[i][j]);
                AtagMatrix[i][j]=0;
            }
            else{
                AtagMatrix[i][j]=A[i][j];
            }
        }
    }
    printf("AtagMatrix start:\n");
    print_matrix(AtagMatrix,num_of_rows);
    printf("AtagMatrix end\n");
    flag=checkJacobiConverge(AtagMatrix,A,num_of_rows,num_of_culs);
    for ( i = 0; i <num_of_rows; ++i) { /*change A to Atag*/
        for (j = 0; j <num_of_culs ; ++j) {
            A[i][j]=AtagMatrix[i][j];
        }


    }

return flag;
}


jacobiMatrix * jacobi(double** A,int num_of_rows,int num_of_culs){
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
    ATag=matBuild(num_of_rows,num_of_culs);
    for ( i = 0; i <100; ++i) {
        PMatrix(P,A,ret->eigenVectors,num_of_rows,num_of_culs,rotator);
        flag=buildATag(ATag,A,rotator,num_of_rows,num_of_culs);
        printf("this is A:\n");
        print_matrix(A,num_of_rows);
        printf("end of A:\n");
        if(flag==1){
            break;
        }
    }
    for (i = 0; i <num_of_rows; ++i) {
        ret->eigenValues[i]=A[i][i];

    }
    printf("jacobi matrix: \n");
    print_jacobi(ret,num_of_rows);
    printf("end of jacobi matrix: \n");
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
eigenTuple** buildEigenTuple(jacobiMatrix* eigens,int num_of_rows){
    int i;
    eigenTuple** egarr;
    egarr=(eigenTuple**) safe_malloc(sizeof(eigenTuple*)*num_of_rows);
    for ( i = 0; i <num_of_rows; ++i) {
        egarr[i]=(eigenTuple*) safe_malloc(sizeof(eigenTuple));
        egarr[i]->val=eigens->eigenValues[i];
        egarr[i]->index=i;
    }
    return egarr;
}
int getK(int num_of_rows,eigenTuple** egarr,int k) /* get number of clusters */{
    int i;
    double temp;
    double max;
    int indexMax=0;

    qsort(egarr, num_of_rows, sizeof(eigenTuple*), cmpfunc);
    if (k>0){
        return k;
    }
    max=-1;
    for ( i = 0; i <num_of_rows/2; ++i) {
        temp=egarr[i]->val-egarr[i+1]->val;
        if(temp>max){
            max=temp;
            indexMax=i+1;
        }

    }
    return indexMax;

}
double** getT(jacobiMatrix* eigens,int num_of_rows,eigenTuple** arr,int k){/* build T matrix */
    double** T=NULL;
    int i;
    int j;
    int index;
    double norm;
    T=matBuild(num_of_rows,k);
    for (j = 0; j < k; ++j) {
        index=arr[j]->index;
        for (i = 0; i <num_of_rows; ++i) {
            T[i][j]=eigens->eigenVectors[i][index];
        }

    }

    for (j = 0; j <k ; ++j) {
        norm=0;
        for (i = 0; i <num_of_rows; ++i) {
            norm=norm+pow(T[i][j],2);

        }
        norm=sqrt(norm);
        for (i = 0; i <num_of_rows; ++i) {
            T[i][j]=T[i][j]/norm;

        }
    }
    printf("this is T matrix after norm\n");
    print_vector_array(T,num_of_rows,k);
    printf("this is T matrix after norm\n");
    return T;
}

int main(int argc, char** argv){
    jacobiMatrix * matJacobi;
    double** vector_array;
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
    /*print_matrix(vector_array,num_of_lines);*/
    free(lines_lens);
    free_mat(file_lines,num_of_lines);/*free string matrix*/
    if(strcmp(argv[1],"wam")==0){
        adjMatrix=wam(vector_array,num_of_lines,dimention);
        print_matrix(adjMatrix,num_of_lines);
        free_matrix(adjMatrix,num_of_lines);
    }
    else if (strcmp(argv[1],"ddg")==0){
        adjMatrix=wam(vector_array,num_of_lines,dimention);
        diagMatrix=ddg(adjMatrix,num_of_lines,0);
        print_matrix(diagMatrix,num_of_lines);
        free_matrix(adjMatrix,num_of_lines);
        free_matrix(diagMatrix,num_of_lines);

    }
    else if (strcmp(argv[1],"lnorm")==0){
        adjMatrix=wam(vector_array,num_of_lines,dimention);
        diagMatrix=ddg(adjMatrix,num_of_lines,1);
        lnormMatrix = lnorm(diagMatrix,adjMatrix ,num_of_lines);
        print_matrix(lnormMatrix,num_of_lines);
        free_matrix(adjMatrix,num_of_lines);
        free_matrix(diagMatrix,num_of_lines);
        free_matrix(lnormMatrix,num_of_lines);
    }
    else if (strcmp(argv[1],"jacobi")==0){

        matJacobi = jacobi(vector_array,num_of_lines,dimention);
        print_jacobi(matJacobi,num_of_lines);
        free_jacobi(matJacobi,num_of_lines);
    }
    free(vector_array);
    return 0;

}
/* from here its kmeans!!!!!!++*/

void calc_norm(double* vector,int dim, double**centroids, int k){
    double min_dis;
    double temp_distance;
    int i;
    int j;
    temp_distance=0;
    min_dis=-1;
    for (i = 0; i <k; i++) {/* amount of clusters*/
        for ( j = 0; j < dim; j++) {/* dimention*/
            temp_distance=pow((vector[j]-centroids[i][j]),2)+temp_distance;

        }
        if(min_dis==-1||temp_distance<min_dis){
            min_dis=temp_distance;
            vector[dim]=i+1;/* update cluster*/
        }
        temp_distance=0;


    }
}
void update_centroids(double** vector_array,int k,int dim,int num_of_vectors, double** centroids){
    int vector_to_centroid;
    int i;
    int j;
    for (i = 0; i < k; i++) {/* initilize centroids*/
        for (j = 0; j < dim+1; j++) {
            centroids[i][j]=0;

        }

    }
    for ( i = 0; i < num_of_vectors; i++) {
        vector_to_centroid= (int) vector_array[i][dim]-1;/* which centroid*/
        centroids[vector_to_centroid][dim]++;
        for ( j = 0; j < dim; j++) {
            centroids[vector_to_centroid][j]=centroids[vector_to_centroid][j]+vector_array[i][j];

        }

    }
    for ( i = 0; i < k; i++) {
        for ( j = 0; j < dim; j++) {
            centroids[i][j]=centroids[i][j]/centroids[i][dim];

        }

    }
}
int check_convergence(double** vector_array,int dim,int num_of_vectors, double** centroids,double epsilon){
    double temp_distance;
    int cluster;
    int i;
    int j;
    epsilon=pow(epsilon,2);
    temp_distance=0;
    for ( i = 0; i <num_of_vectors; ++i) {

        cluster=(int) vector_array[i][dim]-1;
        for ( j = 0; j < dim; j++) {/* dimention*/
            temp_distance=pow((vector_array[i][j]-centroids[cluster][j]),2)+temp_distance;
        }
        if(temp_distance>=epsilon){
            return 0;/* false*/
        }

        temp_distance=0;
    }
    return 1;
}




