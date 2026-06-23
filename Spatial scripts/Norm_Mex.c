//
//  Norm_Mex.c
//  hellomex
//
//  Created by Isaline Lucille Guex on 11.10.23.
//

#include "Norm_Mex.h"
#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get the size of the input matrices
    mwSize n = mxGetN(prhs[0]); //Nb of columns of S


    // Create the output matrix
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    double *outputData = mxGetPr(plhs[0]);

    // Get pointers to the input matrices
    double *matrix_S = mxGetPr(prhs[0]);
    double *matrix_R = mxGetPr(prhs[1]);
    double sum = 0;

    for (mwIndex i = 0; i < n; i++) {
        double diff = matrix_S[i] - matrix_R[i];
        sum += diff*diff;
    }
    sum = sqrt(sum);
    for (mwIndex i = 0; i < n; i++) {
        double diff = matrix_S[i] - matrix_R[i];
        outputData[i]= diff/sum;
    }
}
