#include <stdio.h>
#include "matrix.h"
#include "qr.h"

int main(int argc, char** argv) {
    
    // house vector test
    matrix* vecx = makeMatrix(1, 3);
    vecx->data[0] = 1;
    vecx->data[1] = 3;
    vecx->data[2] = 1;

    matrix *v = copyMatrix(vecx);
    v->data[0] = 1;
    double beta = 0;
    house(vecx, v, &beta);

    printMatrix(v);
    printf("beta = %lf", beta);

    // matrix* q = NULL;
    // matrix* r = NULL;
    // matrix* p = NULL;
    // matrix* a = readMatrix(argv[1]);

    // printf("Original -----------------\n");
    // printMatrix(a);

    // gram_schmidt(a, &q, &r);
    // p = multiplyMatrix(q,r);

    // printf("A again ------------------\n");
    // printMatrix(a);
    // printf("Q ------------------------\n");
    // printMatrix(q);
    // printf("R ------------------------\n");
    // printMatrix(r);
    // printf("Product ------------------\n");
    // printMatrix(p);

    // freeMatrix(a);
    // freeMatrix(p);
    // freeMatrix(q);
    // freeMatrix(r);
    return 0;
}
