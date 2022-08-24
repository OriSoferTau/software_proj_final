
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spkmeans.h"

int memCheck=0;
int numOfMallocs=0;
int numOfFrees=0;


void exit_func(int idx){
    if(idx == 1){
        printf("Invalid Input!\n");
    }

    if (idx == 0){
        printf("An Error Has Occurred\n");
    }
    exit(1);
}

void safe_free(void* element_to_free){
    free(element_to_free);
    memCheck--;
    numOfFrees--;
}

void free_matrix(void** matrix, int num_of_lines){/*free matrix */
    int i;
    for (i = 0; i < num_of_lines;i++) {
        safe_free(matrix[i]);
    }
    safe_free(matrix);
}

void free_jacobi(jacobiMatrix* jacobi,int num_of_lines){
    safe_free(jacobi->eigenValues);
    free_matrix((void **) jacobi->eigenVectors, num_of_lines);
    safe_free(jacobi);
}

void* safe_malloc(size_t size) {
    void * ptr = malloc(size);
    memCheck++;
    numOfMallocs++;
    if(ptr == NULL){
        exit_func(0);
    }
    return ptr;
}

void* safe_calloc(size_t size) {
    void * ptr = calloc(size,sizeof(double));
    memCheck++;
    numOfMallocs++;
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
        for (j = 0; j < num_of_lines-1 ; ++j) {
                printf("%.4f,",matrix[i][j]);

        }

            printf("%.4f\n",matrix[i][num_of_lines-1]);
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
    fclose(fp);
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

    for (i = 0; i <dim ; ++i) {
        norm=norm+pow((vector1[i]-vector2[i]),2);

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
    adjMatrix=matBuild(num_of_rows,num_of_rows);
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
    return lnormMatrix;
}
int isDiagonal(double**matrix,int num_of_rows){
    int i;
    int j;
    for ( i = 0; i <num_of_rows ; ++i) {
        for ( j = 0; j < num_of_rows ; ++j) {
            if (i!=j && matrix[i][j]!=0){
                return 0;
            }
        }
    }
    return 1;
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
                if (fabs(A[i][j]) > max) {
                    max = fabs(A[i][j]);
                    pivotRow = i;
                    pivotCul = j;
                }
            }
        }

    theta = (A[pivotCul][pivotCul] - A[pivotRow][pivotRow]) / (2.0 * A[pivotRow][pivotCul]);
        sign = theta >= 0 ? 1 : -1;
        t = sign / (fabs(theta) + sqrt(theta*theta + 1));
        c = 1.0 / (sqrt(t*t + 1));
        s = t * c;
        P[pivotRow][pivotRow] = c;
        P[pivotCul][pivotCul] = c;
        P[pivotRow][pivotCul] = s;
        P[pivotCul][pivotRow] = -s;
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
        /*printf("start of P:\n");
        print_matrix(P,num_of_rows);
        printf("end of P:\n");*/
        P[pivotRow][pivotRow] = 1;
        P[pivotCul][pivotCul] = 1;
        P[pivotRow][pivotCul] = 0;
        P[pivotCul][pivotRow] = 0;
        free_matrix((void **) tempMatrix, num_of_rows);


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
    int r;
    int j;
    int isConverged;
    double temp;
    double c=rotator->c;
    double s=rotator->s;
    int row=rotator->pivotRow;
    int cul=rotator->pivotCul;
    for ( i = 0; i <num_of_rows; ++i) {
        for (j = 0; j < num_of_culs; ++j) {
            AtagMatrix[i][j] = A[i][j];
        }
    }
        for (r = 0; r <num_of_rows; ++r) {
            temp=c*A[r][row]-s*A[r][cul];
            AtagMatrix[r][row]=temp;
            AtagMatrix[row][r]=temp;

            temp=c*A[r][cul]+s*A[r][row];
            AtagMatrix[r][cul]=temp;
            AtagMatrix[cul][r]=temp;
        }
    AtagMatrix[row][row]=c*c*A[row][row]+s*s*A[cul][cul]-2*s*c*A[row][cul];
    AtagMatrix[cul][cul]=s*s*A[row][row]+c*c*A[cul][cul]+2*s*c*A[row][cul];
    AtagMatrix[row][cul]=0;
    AtagMatrix[cul][row]=0;
    isConverged=checkJacobiConverge(AtagMatrix, A, num_of_rows, num_of_culs);
    for ( i = 0; i <num_of_rows; ++i) {
        for (j = 0; j <num_of_culs ; ++j) {
            A[i][j]=AtagMatrix[i][j];
        }


    }

    return isConverged;

}



jacobiMatrix * jacobi(double** A,int num_of_rows,int num_of_culs){
    int i;
    double** P;
    int isConverged;
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
        if (isDiagonal(A,num_of_rows)==1){
            break;
        }
        PMatrix(P,A,ret->eigenVectors,num_of_rows,num_of_culs,rotator);
        isConverged=buildATag(ATag, A, rotator, num_of_rows, num_of_culs);

        if(isConverged == 1){
            break;
        }
    }
    for (i = 0; i <num_of_rows; ++i) {
        ret->eigenValues[i]=A[i][i];

    }


    free_matrix((void **) ATag, num_of_rows);
    free_matrix((void **) P, num_of_rows);

    safe_free(rotator);
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
    safe_free(lines_lens);
    free_matrix((void **) file_lines, num_of_lines);/*free string matrix*/
    if(strcmp(argv[1],"wam")==0){
        adjMatrix=wam(vector_array,num_of_lines,dimention);
        print_matrix(adjMatrix,num_of_lines);
        free_matrix((void **) adjMatrix, num_of_lines);
    }
    else if (strcmp(argv[1],"ddg")==0){
        adjMatrix=wam(vector_array,num_of_lines,dimention);
        diagMatrix=ddg(adjMatrix,num_of_lines,0);
        print_matrix(diagMatrix,num_of_lines);
        free_matrix((void **) adjMatrix, num_of_lines);
        free_matrix((void **) diagMatrix, num_of_lines);

    }
    else if (strcmp(argv[1],"lnorm")==0){
        adjMatrix=wam(vector_array,num_of_lines,dimention);
        diagMatrix=ddg(adjMatrix,num_of_lines,1);
        lnormMatrix = lnorm(diagMatrix,adjMatrix ,num_of_lines);
        print_matrix(lnormMatrix,num_of_lines);
        free_matrix((void **) adjMatrix, num_of_lines);
        free_matrix((void **) diagMatrix, num_of_lines);
        free_matrix((void **) lnormMatrix, num_of_lines);
    }
    else if (strcmp(argv[1],"jacobi")==0){
/*        printf("this is input matrix:\n");
        print_matrix(vector_array,num_of_lines);
        printf("end of input matrix:\n");*/
        matJacobi = jacobi(vector_array,num_of_lines,dimention);
        print_jacobi(matJacobi,num_of_lines);
        free_jacobi(matJacobi,num_of_lines);
    }
    free_matrix((void **) vector_array, num_of_lines);
    /*printf("memcheck= %d num of mallocs= %d num of frees= %d",memCheck,numOfMallocs,numOfFrees);*/
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




