//
//  Comp_Wise_Sub.c
//  hellomex
//
//  Created by Isaline Lucille Guex on 11.10.23.
//

#include "Comp_Wise_Sub.h"
#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get the size of the input matrices
    mwSize n = mxGetN(prhs[0]); //Nb of columns of S
    mwSize n1 = mxGetN(prhs[1]); //Nb of columns of R


    // Create the output matrix
    plhs[0] = mxCreateDoubleMatrix(n, n1, mxREAL);
    double *outputData = mxGetPr(plhs[0]);

    // Get pointers to the input matrices
    double *matrix_S = mxGetPr(prhs[0]);
    double *matrix_R = mxGetPr(prhs[1]);

    for (mwIndex i = 0; i < n; i++) {
        for (mwIndex j = 0; j < n1; j++) {
            double diff = (matrix_S[i] - matrix_R[j]);
            outputData[i + j * n] = diff * diff; //sqrt
        }
    }
}
