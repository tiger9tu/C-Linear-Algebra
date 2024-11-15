#include <stdio.h>
#include "matrix.h"
#include "qr.h"

int main(int argc, char** argv) {
    

    matrix* a = makeMatrix(3,4);
    for (int i = 0; i < a->height; i++)
    {
        for (int j = 0; j < a->width; j++)
        {
            a->data[i*3 + j] = i*3 + j;
        }
        
    }
    printf("A:\n");
    printMatrix(a);
    // house vector test
    houseHolderFactor* vod = houseHolderQR(a);
    printMatrix(vod->qrT);
    
    matrix* Q = makeMatrix(4,4);
    matrix* R = makeMatrix(3,4);

    getExplicitQRFromHouseholder(vod,Q,R);
    printf("matrix Q = \n");
    printMatrix(Q);
    printf("R= \n");
    printMatrix(R);
    printf("QR= \n");
    printMatrix(multiplyMatrix(Q,R));    
    freeMatrix(vod->qrT);
    freeMatrix(a);

    return 0;
}
