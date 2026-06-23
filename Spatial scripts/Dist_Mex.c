//
//  Dist_Mex.c
//  hellomex
//
//  Created by Isaline Lucille Guex on 11.10.23.
//

#include "Dist_Mex.h"
#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check the number of input arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:myFunction:nrhs",
                          "Two input arguments required.");
    }

    // Get the size of the input matrices
    mwSize m = mxGetM(prhs[0]);
    mwSize n = mxGetN(prhs[0]);

    // Check that the input matrices have the same dimensions
    if (m != mxGetM(prhs[1]) || n != mxGetN(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:myFunction:inputDimensions",
                          "Input matrices must have the same dimensions.");
    }

    // Create the output matrix
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    double *outputData = mxGetPr(plhs[0]);

    // Get pointers to the input matrices
    double *matrix1 = mxGetPr(prhs[0]);
    double *matrix2 = mxGetPr(prhs[1]);

    // Calculate the component-wise Euclidean distance and store it in the output matrix
    for (mwIndex i = 0; i < m; i++) {
        for (mwIndex j = 0; j < n; j++) {
            double diff = sqrt(matrix1[i + j * m]*matrix1[i + j * m] + matrix2[i + j * m]*matrix2[i + j * m]);
            outputData[i + j * m] = diff;
        }
    }
}
