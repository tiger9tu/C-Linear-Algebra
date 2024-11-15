#include <stdio.h>
#include "matrix.h"
#include "qr.h"

int main(int argc, char** argv) {
    

    matrix* a = makeMatrix(3,3);
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            a->data[i*3 + j] = i*3 + j;
        }
        
    }
    printf("A:\n");
    printMatrix(a);
    // house vector test
    houseHolderFactor* vod = houseHolderQR(a);
    printMatrix(vod->qrT);
    restoreFromHouseholderFactor(vod);

    // matrix* vecx = makeMatrix(1, 3);
    // printf("vecx = %p\n", vecx);
    // vecx->data[0] = 1;
    // vecx->data[1] = 3;
    // vecx->data[2] = 1;

    // matrix* v = copyMatrix(vecx);
    // double beta = 0;

    // house(vecx, &beta);
    // printMatrix(vecx);
    // printf("beta = %lf\n", beta);

    // matrix* I = eyeMatrix(3);
    // matrix* vvt = multiplyMatrix(vecx, transposeMatrix(vecx));

    // matrix* P = addMatrix(I, scaleMatrix(vvt, - beta));
    // matrix* Px = multiplyMatrix(P, v);
    // printMatrix(Px);

    // printMatrix(v);
    // printf("beta = %lf", beta);

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
