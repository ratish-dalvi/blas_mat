#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mkl_blas.h>
#include <blas.h>

void lumatvec_(double *B, int *pldb, int *pm, int *pn, int *ipiv, int *perrflag)
{

//===================================================================================
//
// LU factorization with partial pivoting on an m x n matrix B. This version should
// be based on matrix-vector products, and calls BLAS-2 (BLAS level 2) routines.
// The parameter small below determines when partial pivoting has failed from the
// diagonal element B(k,k) being too small to safely divide by. When that happens
// errflag is set to k and the function bails out. The array B can be rectangular
// (viz., m /= n). This assumes that the indexing of the matrix B starts from
// 1, not 0. Otherwise a failure on the first step is indistinguishable from a
// no-error condition.
//
// B  is an array of doubles containing the matrix B in column-major order. It
//    is (part of) an array of declared leading dimension ldb in the calling
//    function. 
//   m is the number of rows in B 
//   n is the number of columns in B 
// ipiv is an array of length at least m which contains the pivot indices
// errflag is what would normally be returned by a function in C. Values:
// errflag = 0: success
//         < 0: if errflag = -k, the k-th argument had an illegal value
//         > 0: if errflag =  k, B(k,k) is too small to rely upon and the 
//              factorization probably failed.
//

// The BLAS requires each vector argument also has an "increment" giving how far
// apart consecutive entries of the vector are in memory. E.g., if a matrix H is
// stored in the upper 128 x 128 part of an array declared as G[128][256] in
// row-major order, consecutive entries in a row of G have an increment of 1.
// Accessing consecutive entries in a column of G will have an increment of 256.
//
// Similarly with column major ordering, if a matrix H stored in the first 128x128
// entries of an array declared as G(256,128), then accessing elements along a
// row of G have increment = 256, and elements along a column have increment = 1.
//
//===================================================================================

    // Value to use to bail out because pivot is too small
    double small = ((double) 1000.0)*2.216e-16;

    // To avoid messing with stars and ampersands, make local variables
    int ldb     = *pldb;
    int n       = *pn;
    int m       = *pm;

    double ONE   =  1.0;
    double ZERO  =  0.0;
    double MONE  = -1.0;

    // Because quantities like m-k cannot be passed by address to the BLAS, use
    // some temporary variables, e.g., nrows = ldb-j and then pass in &nrows
    int nrows, ncols, whatever;
    *perrflag = 0;

}
