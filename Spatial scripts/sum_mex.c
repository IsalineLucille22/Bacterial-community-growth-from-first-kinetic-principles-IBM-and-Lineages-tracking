//
//  sum_mex.c
//  hellomex
//
//  Created by Isaline Lucille Guex on 14.10.23.
//

#include "sum_mex.h"
#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get the size of the input matrices
    mwSize n = mxGetN(prhs[0]); // Number of columns of S
    mwSize m = mxGetM(prhs[0]); // Number of rows of S

    // Get pointers to the input matrices
    double *matrix_S = mxGetPr(prhs[0]);
    double *sum_type = mxGetPr(prhs[1]);
    
    int dim;
    double *sum;
    if (sum_type[0] == 1) { // Sum on column
        dim = n;
        sum = (double*)malloc(dim * sizeof(double));
        for (mwIndex j = 0; j < n; j++) {
            double sum_temp = 0.0;
            for (mwIndex i = 0; i < m; i++) {
                sum_temp += matrix_S[i + j * m];
            }
            sum[j] = sum_temp;
        }
    }
    else { // Sum on row
        dim = m;
        sum = (double*)malloc(dim * sizeof(double));
        for (mwIndex i = 0; i < m; i++) {
            double sum_temp = 0.0;
            for (mwIndex j = 0; j < n; j++) {
                sum_temp += matrix_S[i + j * m];
            }
            sum[i] = sum_temp;
        }
        
    }

    // Create the output matrix
    plhs[0] = mxCreateDoubleMatrix(dim, 1, mxREAL);
    double *outputData = mxGetPr(plhs[0]);

    for (int i = 0; i < dim; i++) {
        outputData[i] = sum[i];
    }
    
    //free(sum); // Don't forget to free the dynamically allocated memory.
}
/*for (mwIndex col = 0; col < n; col++) {
    double sum = 0.0;
    for (mwIndex row = 0; row < m; row++) {
        sum += matrix[row + col * m];
    }
    outputData[col] = sum;
}

for (mwIndex row = 0; row < m; row++) {
    double sum = 0.0;
    for (mwIndex col = 0; col < n; col++) {
        sum += matrix[row + col * m];
    }
    outputData[row] = sum;
}*/
