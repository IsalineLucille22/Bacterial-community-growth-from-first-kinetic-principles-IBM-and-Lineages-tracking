//
//  distEuclid_Mex.c
//  hellomex
//
//  Created by Isaline Lucille Guex on 11.10.23.
//

#include "distEuclid_Mex.h"
#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get the size of the input matrices
    mwSize n = mxGetN(prhs[0]); //Nb of columns of S
    mwSize n1 = mxGetN(prhs[2]); //Nb of columns of R


    // Create the output matrix
    plhs[0] = mxCreateDoubleMatrix(n, n1, mxREAL);
    double *outputData = mxGetPr(plhs[0]);

    // Get pointers to the input matrices
    double *matrix_Sx = mxGetPr(prhs[0]);
    double *matrix_Sy = mxGetPr(prhs[1]);
    double *matrix_Rx = mxGetPr(prhs[2]);
    double *matrix_Ry = mxGetPr(prhs[3]);

    for (mwIndex i = 0; i < n; i++) {
        for (mwIndex j = 0; j < n1; j++) {
            double diff_x = (matrix_Sx[i] - matrix_Rx[j]);
            double diff_y = (matrix_Sy[i] - matrix_Ry[j]);
            outputData[i + j * n] = diff_x * diff_x + diff_y * diff_y; //sqrt
        }
    }
}
