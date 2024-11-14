#ifndef __ISOMAP_QR
#define __ISOMAP_QR
#include "matrix.h"

void gram_schmidt(matrix* a, matrix** q, matrix** r);
matrix* unitVectorRows(matrix* a);
matrix* unitVectorColumns(matrix* a);
void house(matrix* x /*vector x*/, matrix* v, double* beta);
#endif
