#include <stdio.h>
#include "matrix.h"
#include "qr.h"

#include <time.h>
#include <gsl/gsl_matrix.h>



void a_b_(matrix* a){
    printf("The original matrix A:\n");
    printMatrix(a);
    // gram-schimit decomposition
    matrix* gramSchimidtQ =NULL;
    matrix* gramSchimidtR= NULL;

    gram_schmidt(a, &gramSchimidtQ,&gramSchimidtR);
    printf("The gram schimidt decomposition:");
    printf("Q: \n");
    printMatrix(gramSchimidtQ);
    printf("R: \n"); 
    printMatrix(gramSchimidtR);

    houseHolderFactor* hhf = houseHolderQR(a);
    matrix* houseHolderQ = NULL;
    matrix* houseHolderR = NULL;
    getExplicitQRFromHouseholder(hhf, &houseHolderQ, &houseHolderR);
    printf("The Householder decomposition:");
    printf("Q: \n");
    printMatrix(houseHolderQ);
    printf("R: \n"); 
    printMatrix(houseHolderR);
    printf("\n\n\n");
}

void c_(houseHolderFactor* hhf, matrix* x, matrix* y){
    matrix* r = NULL;
    matrix* q = NULL;

    getExplicitQRFromHouseholder(hhf, &q, &r);
    printf("Q: \n");
    printMatrix(q);

    printf("x: \n");
    printMatrix(x);
    printf("explicit computation of QX: \n");
    printMatrix(multiplyMatrix(q,x));
    printf("Implicit computation of Qx: \n");
    implilcitQx(hhf,x);
    printMatrix(x);

    printf("\ny: \n");
    printMatrix(y);
    printf("explicit computation of QTy: \n");
    printMatrix(multiplyMatrix(transposeMatrix(q),y));
    printf("implicit computation of QTy: \n");
    implilcitQTx(hhf, y);
    printMatrix(y);
}


// Function to generate a random upper triangular matrix of size n x n
void generateRandomUpperTriangularSquareMatrix(int n, matrix** r, gsl_matrix ** gsl_r) {
    *r = makeMatrix(n , n);
    *gsl_r = gsl_matrix_alloc(n, n);
    printf("gslr size = %d\n",(*gsl_r)->size2);   
    // Use the current time as the seed for the random number generator
    srand(time(NULL));
    
    // Generate the upper triangular matrix
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            // Assign a random number to the upper triangular part (including the diagonal)
            int randn = rand() % 100;
            (*r)->data[i*n + j] = randn;  // Generate a random integer between 0 and 99
            gsl_matrix_set (*gsl_r, i, j, (double)randn);
        }
    }
}

// Function to generate a random upper triangular matrix of size n x n
void generateRandomSquareMatrix(int n, matrix** r, gsl_matrix ** gsl_r) {
    *r = makeMatrix(n , n);
    *gsl_r = gsl_matrix_alloc(n, n);
    
    // Use the current time as the seed for the random number generator
    srand(time(NULL));
    
    // Generate the upper triangular matrix
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            // Assign a random number to the upper triangular part (including the diagonal)
            int randn = rand() % 100;
            (*r)->data[i*n + j] = randn;  // Generate a random integer between 0 and 99
            gsl_matrix_set (*gsl_r, i, j, (double)randn);
        }
    }
}

// Function to print a gsl_matrix
void printGslMatrix(const gsl_matrix *matrix) {
    size_t rows = matrix->size1;
    size_t cols = matrix->size2;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            // Access the matrix element using gsl_matrix_get
            printf("%6.2f ", gsl_matrix_get(matrix, i, j));
        }
        printf("\n");
    }

    
}

void getNormalMatrixFromGsl(matrix* a, gsl_matrix* b){
    for (size_t i = 0; i < a->height; i++) {
        for (size_t j = 0; j < a->width; j++) {
            // Access the matrix element using gsl_matrix_get
            a->data[i* a->width + j] = gsl_matrix_get(b, i,j);
        }
    }
}

double calculateFrobeniusNorm(matrix* a) {
    double norm = 0.0;
    for (int i = 0; i < a->height; i++) {
        for (int j = 0; j < a->width; j++) {
            norm += a->data[i* a->width + j] *a->data[i* a->width + j];
        }
    }
    return sqrt(norm);
}

double calculateGslFrobeniusNorm(gsl_matrix* a) {
    double norm = 0.0;
    for (int i = 0; i < a->size1; i++) {
        for (int j = 0; j < a->size2; j++) {
            norm += gsl_matrix_get(a, i,j) *  gsl_matrix_get(a, i,j) ;
        }
    }
    return sqrt(norm);
}

void d_(){
    int n = 3;

    // setting up the random sample 
    matrix* r = NULL;
    matrix* b = NULL;

    gsl_matrix *gsl_r = NULL;
    gsl_matrix *gsl_b = NULL;

    generateRandomUpperTriangularSquareMatrix(n, &r, &gsl_r);
    generateRandomSquareMatrix(n, &b, &gsl_b);
    printMatrix(r);
    printf("gsl_r: \n");
    printGslMatrix(gsl_r);

    gsl_matrix *gsl_t = gsl_matrix_alloc(n,n);
    gsl_linalg_QR_decomp_r(gsl_b, gsl_t);

    gsl_matrix* gsl_bq = gsl_matrix_alloc(n,n);
    gsl_matrix* gsl_br = gsl_matrix_alloc(n,n);

    gsl_linalg_QR_unpack_r(gsl_b, gsl_t,gsl_bq,gsl_br);

    printf("gsl_bq: \n");
    printGslMatrix(gsl_bq);

    gsl_matrix* gsl_a = gsl_matrix_alloc(n,n);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, gsl_bq, gsl_r, 0.0 , gsl_a);
    printf("gsl_a: \n");
    printGslMatrix(gsl_a);
    matrix* a = makeMatrix(n,n);
    getNormalMatrixFromGsl(a, gsl_a);
    printMatrix(a);

    // my QR decomposition
    houseHolderFactor* hhf = houseHolderQR(a);
    matrix* q;
    // matrix* r;
    getExplicitQRFromHouseholder(hhf, &q, &r);
    // error R
    // double myErrorR = calculateFrobeniusNorm();


    
    // gsl_linalg_QR_decomp_r()
}


int main(int argc, char** argv) {
    matrix* a = makeMatrix(3,4);
    for (int i = 0; i < a->height; i++)
        for (int j = 0; j < a->width; j++)
            a->data[i*3 + j] = i*3 + j;

    // a_b_(a); 

    matrix* x = makeMatrix(1, 4);
    for (size_t i = 0; i < x->height; i++)
        x->data[i] = i + 1;
    
    matrix* y = copyMatrix(x);

    c_(houseHolderQR(a),x, y);

    // d_();

    return 0;
}
