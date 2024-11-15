#ifndef __ISOMAP_QR
#define __ISOMAP_QR
#include "matrix.h"
typedef struct {
    double* betas;
    matrix* qrT;
}houseHolderFactor;
void gram_schmidt(matrix* a, matrix** q, matrix** r);
matrix* unitVectorRows(matrix* a);
matrix* unitVectorColumns(matrix* a);
void house( matrix* v, double* beta);
houseHolderFactor* houseHolderQR(matrix* a);
void restoreFromHouseholderFactor(houseHolderFactor* hhf);
#endif
