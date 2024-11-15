#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "qr.h"
/*===========================================================================
 * gram_schmidt
 * This algorithm performs the Modified Gram-Schmidt algorithm for computing
 * the orthanormalized vectors of a matrix.
 *
 * This algorithm is explained in "Matrix Computations" (3rd Edition) by
 * Golub and Loan on pages 230-232. It closely follows the explicit
 * algorithm printed at the top of page 232.
 *=========================================================================*/
void gram_schmidt(matrix* a, matrix** q, matrix** r) {
    int i, j, k;
    double l2norm;
    double sum;
    double *vPtr;
    double *qPtr;
    double *rPtr;
    matrix *v;

    assert(a->width == a->height, "Currently, QR decomposition only suports square matricies.");
    assert(*q == NULL && *r == NULL, "The q and r matrices must be initalized to NULL.");

    // Because the Modified Gram-Schmidt algorithm overwrites and damages the original
    // 'A' matrix, we copy all of the values into matrix 'V' and use it instead.
    v = copyMatrix(a);
    *q = makeMatrix(a->width, a->height);
    *r = makeMatrix(a->width, a->height);

    // For each column in A (now called V)
    for (k = 0; k < v->width; k++) {

        vPtr = v->data + k;
        // Step 1: Get the L2-Norm of column k
        l2norm = 0;
        for (i = 0; i < v->height; i++) {
            l2norm += *vPtr * *vPtr;
            vPtr += v->width;
        }
        l2norm = sqrt(l2norm);
        if(isnan(l2norm))
          l2norm=0;

        // Store this value in R(k,k)
        // The nice thing about the rPtr variable is that
        // it only has to be readjusted at the beginning of the
        // first 'for' loop. After each use, we just increment by 1.
        rPtr = (*r)->data + (k * (*r)->width) + k;
        *rPtr = l2norm;
        rPtr++;

        // Step 2: Normalize A's column k and store it in Q.
        vPtr = v->data + k;
        qPtr = (*q)->data + k;
        for (i = 0; i < v->height; i++) {
            *qPtr=(l2norm==0? 0 : *vPtr / l2norm);
            vPtr += v->width;
            qPtr += (*q)->width;
        }

        // Step 3: 2 parts. For each column after K, do the following:
        for (j = k + 1; j < v->width; j++) {
            // Step 3a: Dot Product Q's column K with A's column J,
            // storing the result of the Dot Product at R(k,j)
            qPtr = (*q)->data + k;
            vPtr = v->data + j;
            sum = 0;
            for (i = 0; i < a->height; i++) {
                sum += *vPtr * *qPtr;
                vPtr += v->width;
                qPtr += (*q)->width;
            }
            *rPtr = sum;
            rPtr++;

            // Step 3b: Multiply Q's column K with R(k,j)
            // (which is stored at *rPtr). Then take A's column J
            // and subtract from it Q's K * R(k,j). Take this
            // result and store it back into A's column J.
            // R(k,j) is represented by 'sum' calculated in step 3b.
            vPtr = v->data + j;
            qPtr = (*q)->data + k;
            for (i = 0; i < v->height; i++) {
                *vPtr = *vPtr - (*qPtr * sum);
                vPtr += v->width;
                qPtr += (*q)->width;
            }
        }
    }
    freeMatrix(v);
}

/*===========================================================================
 * unitVectorRows
 * This algorithm normalizes each row of a matrix as if it were a vector.
 *=========================================================================*/
matrix* unitVectorRows(matrix* a) {
    matrix* out = copyMatrix(a);
    int i, j;
    double l2norm;
    double *aPtr = a->data;
    double *outPtr = out->data;

    for (i = 0; i < a->height; i++) {
        l2norm = 0;

        for (j = 0; j < a->width; j++) {
            l2norm += *aPtr * *aPtr;
            aPtr++;
        }

        l2norm = sqrt(l2norm);

        for (j = 0; j < out->width; j++) {
            *outPtr = *outPtr / l2norm;
            outPtr++;
        }
    }

    return out;
}

/*===========================================================================
 * unitVectorColumns
 * This algorithm normalizes each column of a matrix as if it were a vector.
 *=========================================================================*/
matrix* unitVectorColumns(matrix* a) {
    matrix* out = copyMatrix(a);
    int i, j;
    double l2norm;
    double *aPtr = a->data;
    double *outPtr = out->data;

    for (i = 0; i < a->width; i++) {
        aPtr = a->data + i;
        outPtr = out->data + i;
        l2norm = 0;

        for (j = 0; j < a->height; j++) {
            l2norm += *aPtr * *aPtr;
            aPtr += a->width;
        }

        l2norm = sqrt(l2norm);

        for (j = 0; j < out->height; j++) {
            *outPtr = *outPtr / l2norm;
            outPtr += out->width;
        }
    }

    return out;
}

void house(matrix* v, double* beta){
    /* v is initialized as x */
    int m = v->height;
    matrix* subx = subVectorRef(v,1,m);
    double sigma = innerProductVector(subx,subx);
    if(sigma == 0 && v->data[0] >= 0) {
        *beta = 0;
    }
    else if(sigma == 0 && v->data[0] < 0) {
        *beta = -2;
    }
    else{
        // calculate x1 - ||x|| while reduce cancelation error
        double normx = sqrt(sigma + (v->data[0])* (v->data[0]));
        if(v->data[0] <= 0) {
            v->data[0] = v->data[0] - normx;
        }
        else {
            v->data[0] = -sigma / (v->data[0] + normx);
        }
        *beta = 2 * v->data[0] * v->data[0] / (sigma + v->data[0] * v->data[0]);
        rescaleMatrix(v, 1 / v->data[0]);
    }

}




houseHolderFactor* makeHouseHolderFactor(int n){
    houseHolderFactor* out = (houseHolderFactor*) malloc(sizeof(houseHolderFactor));
    out->betas = (double*) malloc(sizeof(double) * (n));
    return out;    
}

/* x = xscale*x + yscale* y*/
void rescaleVecAdd(matrix* x, matrix*y, double xScale, double yScale){
    int length = x->height > x->width ? x->height : x->width;
    for (size_t i = 0; i < length; i++)
    {
        x->data[i] = xScale * x->data[i] + yScale*y->data[i]; 
    }
}

/*===========================================================================
 * Implicit Q.c, where Q is the one Householder factor Q, and c is a vector
 *=========================================================================*/
void Qc(matrix* v, double beta, matrix*c){
    double projLength = innerProductVector(v, c);
    rescaleVecAdd(c, v, 1, -beta * projLength);
}

matrix* removeV(houseHolderFactor* hhf, int vrank){
    matrix* R = transposeMatrix(hhf->qrT);
    for (size_t i = 0; i < vrank && i < hhf->qrT->width; i++)
    {
        for (size_t j = i + 1; j < R->height; j++)
        {
            R->data[j * R->width + i] = 0;
        }   
    }
    return R;
}

void printTranspose(matrix* a){
    matrix* tmp = transposeMatrix(a);
    printMatrix(tmp);
    freeMatrix(tmp);
}

/*===========================================================================
 * Implicit 
 *=========================================================================*/
houseHolderFactor* houseHolderQR(matrix* a){
    // allocate factor
    houseHolderFactor* hhf = makeHouseHolderFactor(a->width);
    hhf->qrT = transposeMatrix(copyMatrix(a));

    matrix* tmpx1 = subVectorRef(hhf->qrT, 0, 3);
    matrix* x1 = copyMatrix(tmpx1);
    printf("x1 = \n");
    printMatrix(x1);
matrix* I3 = eyeMatrix(3);
    
    matrix* xbuffer = makeMatrix(1, a->height);
    for(int i = 0; i < 2; i++){
        matrix* I = eyeMatrix(a->height - i);
        matrix* vi = subVectorRef(hhf->qrT, i*hhf->qrT->width + i, (i+1)*hhf->qrT->width);

        for (size_t l = 0; l < a->height - i; l++)
        {
            xbuffer->data[l] = vi->data[l];
        }
        

        house(vi, &(hhf->betas[i]));
              printf("house vi height = %d: \n", vi->height);
        printMatrix(vi);
        // printf("vi[2] = %lf\n", vi->data[2]);
        printf("end house vi\n");
        printf("beta = %lf\n", hhf->betas[i]);
        double scaleback = sqrt(hhf->betas[i] / 2);
        matrix* ve = scaleMatrix(vi,  scaleback);
        printf("ve: \n");
        printMatrix(ve);
        // now check I - beta vvT x = ||x||e1
        matrix* viT = transposeMatrix(vi);
        matrix* vvT = multiplyMatrix(vi, viT);
        printf("outer: \n");
        printMatrix(vvT);
        rescaleMatrix(vvT, -hhf->betas[i]);
        matrix* mirror = addMatrix(I, vvT);
        printf("P: \n");
        printMatrix(mirror);
        printf("x = \n");
        matrix* xsub = subVectorRef(xbuffer, 0, a->height - i);
        printMatrix(xsub);
        matrix* flakx = multiplyMatrix(mirror,xsub );
        printf("flacx: \n");
        printMatrix(flakx);

        // update the submatrix
        for (size_t j = i + 1; j < a->width; j++)
        {
            matrix* vj = subVectorRef(hhf->qrT, j*hhf->qrT->width + i, (j+1)*hhf->qrT->width);
            Qc(vi, hhf->betas[i], vj);
        }


        // update the current column's top entry R(i, i)
        double vTc = 0;
        for (size_t l = 0; l < a->height - i; l++)
        {
            vTc += xbuffer->data[l] * vi->data[l];
        }
        
        vi->data[0] = xbuffer->data[0] - hhf->betas[i] * vTc;

        printf("updated:\n");
        printTranspose(hhf->qrT);
        //debug
        // matrix* debug = copyMatrix(hhf->qrT);
        matrix* Rdebug = removeV(hhf, i+ 1);
        printf("Rdebug:  \n");
        printMatrix(Rdebug);
        matrix* tmpvi = makeMatrix(1, 3);
        for (size_t f = i; f < 3; f++)
        {
            tmpvi->data[f] = vi->data[f-i];
        }
        
        
        printf("tmpvi = \n");
        tmpvi->data[i] = 1;

        printMatrix(tmpvi);
        matrix* tmpviT = transposeMatrix(tmpvi);
        matrix* bvvt = multiplyMatrix(tmpvi, tmpviT);
        rescaleMatrix(bvvt,- hhf->betas[i]);
        matrix* QT = addMatrix(I3, bvvt);
        matrix* QTR = multiplyMatrix(QT, Rdebug);
      printf("QTF = \n");
        printMatrix(QTR);
        freeMatrix(Rdebug);
        freeMatrix(tmpvi);
        freeMatrix(tmpviT);
        freeMatrix(bvvt);
        freeMatrix(QT);
        freeMatrix(QTR);
        freeMatrix(I);
  
        
        printf("Rdebug end\n");
    }
    freeMatrix(xbuffer);
    freeMatrix(I3);
    return hhf;
}

void restoreFromHouseholderFactor(houseHolderFactor* hhf){
    matrix* a = transposeMatrix(hhf->qrT);
    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = i + 1; j < a->height; j++)
        {
            a->data[j * a->width + i] = 0;
        }
    }
    printf("a = \n");
    printMatrix(a);
    matrix* vs = copyMatrix(hhf->qrT);
   for (size_t i = 0; i < vs->height; i++)
    {
        vs->data[i * vs->width + i] = 1;
        for (size_t j = 0; j < i; j++)
        {
            vs->data[i * vs->width + j] = 0;
        }
    }

    matrix* I = eyeMatrix(hhf->qrT->width);

    for (int i =  1; i >=0; i--)
    {
        matrix* aref = a;
        matrix* vi = subVectorRef(vs, i*vs->width, (i + 1)*vs->width);

        printf("voi");
        printMatrix(vi);
        matrix* vvT = multiplyMatrix(vi, transposeMatrix(vi));
        rescaleMatrix(vvT, - hhf->betas[i]); // T -> +
        matrix* Q = addMatrix(I, vvT);
        a = multiplyMatrix(Q, aref);
        freeMatrix(vvT);
        freeMatrix(Q);
        printf("a in i = %d\n", i);
printMatrix(a);

    }
   freeMatrix(I);
   freeMatrix(vs); 

   printMatrix(a);

}