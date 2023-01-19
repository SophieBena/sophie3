//
//  MathSupport.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/25/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>


@interface MathSupport : NSObject {

}

#define HOUSEHOLDER_CALC_MODE_H1     1
#define HOUSEHOLDER_CALC_MODE_H2     2

/**
 Solution of a least-squares problem with possible rank deficiency.

 Algorithm HFTI (Lawson and Hanson, page 81-82).  See source documentation.

 @param a
   i/o: matrix of coefficients of the least squares problem

 @param m
   input row dimension of a, i.e. a[m][n]

 @param n
   input column dimension of a, i.e. a[m][n]

 @param b
   i/o: [m] right-hand-side vectors of the least-squares prob

 @param nb
   inp: number of vectors (columns) in the matrix b

 @param tau
   inp: absolute tolerance param for determ of pseudorank

 @param k
   out: pseudorank

 @param rnorm
   out: [nb] residual norms for the nb solution vectors

 @param h
   out: [n] pivot scalars for Householder transformations Qi

 @param g
   out: [n] pivot scalars for Householder transformations Ki

 @param p
   out: [n] interchange record
*/

+(void)hfti:(double **)a m:(NSInteger)m n:(NSInteger)n b:(double **)b nb:(NSInteger)nb tau:(double)tau k:(NSInteger *)k
	  rnorm:(double *)rnorm h:(double *)h g:(double *)g p:(NSInteger *)p;
/**
 Householder decomposition of a vector and application of the resulting transformation to column or row vectors of a matrix.

 Algorithm H12 (Lawson and Hanson, page 57). See source documentation

 @param mode
   Calculation mode, i.e. algorithm H1 or H2

 @param p
   index of the pivot element of the vector stored in v

 @param l
   range of indices to zero elements of vector stored in v

 @param m
   TBD

 @param v
   matrix containing the vector v stored as v[][vCol]

 @param vCol
   index number of column of v that contains vector

 @param h
   extra storage for the transformed pth element of v (u[p])

 @param c
  matrix containing the set of column vectors to apply the transformation, i.e. c[][cColStart] -> c[][cColEnd]

 @param cColStart
   starting column

 @param cColEnd
   ending column
 */
+(void)householderColCol:(NSInteger)mode p:(NSInteger)p l:(NSInteger)l m:(NSInteger)m v:(double **)v vCol:(NSInteger)vCol
					   h:(double *)h c:(double **)c cColStart:(NSInteger)cColStart cColEnd:(NSInteger)cColEnd;
/**
 see documentation for householderColCol
 */
+(void)householderColRow:(NSInteger)mode p:(NSInteger)p l:(NSInteger)l m:(NSInteger)m v:(double **)v vCol:(NSInteger)vCol
					   h:(double *)h c:(double **)c cRowStart:(NSInteger)cRowStart cRowEnd:(NSInteger)cRowEnd;

/**
 see documentation for householderColCol
 */
+(void)householderRowCol:(NSInteger)mode p:(NSInteger)p l:(NSInteger)l m:(NSInteger)m v:(double **)v vRow:(NSInteger)vRow
					   h:(double *)h c:(double **)c cColStart:(NSInteger)cColStart cColEnd:(NSInteger)cColEnd;

/**
 see documentation for householderColCol
 */
+(void)householderRowRow:(NSInteger)mode p:(NSInteger)p l:(NSInteger)l m:(NSInteger)m v:(double **)v vRow:(NSInteger)vRow
					   h:(double *)h c:(double **)c cRowStart:(NSInteger)cRowStart cRowEnd:(NSInteger)cRowEnd;

#define MIN1D_SUCCESS          0
#define MIN1D_BAD_INITIAL      1
#define MIN1D_ITERS_EXCEEDED   2

#define MODMRT_SUCCESS         0
#define MODMRT_BAD_INITIAL     1
#define MODMRT_ITERS_EXCEEDED  2

#define RQMCG_SUCCESS          0
#define RQMCG_ITERS_EXCEEDED   1
#define RQMCG_SINGULAR_B       2
#define RQMCG_ZERO_DETERMINANT 3

#define VMMIN_SUCCESS          0
#define VMMIN_BAD_INITIAL      1

/**
 One-dimensional minimisation of a function using success-failure search and parabolic inverse interpolation.

 Algorithm 17 of Nash (1979)

 @param bb
   initial guess to minimum, resulting minimum position

 @param st
   initial and final step-length

 @param reltest
   maximal |error| allowed in computing value of bb

 @param ifn
   input: maximum No. of func eval; output: actual No.

 @param fnminval
   minimum function value on return

 @param fn1d
  code block that returns the function value at bb.  Arguments:
    - bb: argument
    - notcomp: flag to indicate the function is not computable

 @return
   mode flag (as defined in MathSupport.h):
     - MIN1D_SUCCESS: Successful termination.
     - MIN1D_BAD_INITIAL: The initial value of the function is invalid.
     - MIN1D_ITERS_EXCEEDED: Number of function iterations exceeds limit.
 */
+(NSInteger)min1d:(double *)bb st:(double *)st reltest:(double)reltest ifn:(NSInteger *)ifn fnminval:(double *)fnminval fn1d:(double (^)(double bb, NSInteger *notcomp))fn1d;

/**
 Modified Marquardt method for minimising a nonlinear sum-of-squares function

 Algorithm 23 of Nash (1979)

 @param n
   number of parameters in nonlinear least-squares problem

 @param m
   number of residuals in nonlinear least-squares problem

 @param Bvec
   vector of parameters in nonlinear least-squares problem

 @param Fmin
   returned value, minimum value of sum-of-squares function

 @param reltest
   maximal |error| allowed in computing elements of Bvec

 @param maxIter
   maximum number of times the func (*nlres) may be called

 @param nlres
   code block for residual function: (*nlres)(NSInteger i, NSInteger n, double *Bvec, NSInteger *notcomp)

 @param nljac
   code block for gradient function: (*nljac)(NSInteger i, NSInteger n, double *Bvec, double *X))

 @return
   mode flag (as defined in MathSupport.h):
     - MODMRT_SUCCESS: Successful termination.
     - MODMRT_BAD_INITIAL: The initial value(s) for the solution is(are) invalid.
     - MODMRT_ITERS_EXCEEDED: Number of function iterations exceeds limit.
 */
+(NSInteger)modmrt:(NSInteger)n m:(NSInteger)m Bvec:(double *)Bvec Fmin:(double *)Fmin reltest:(double)reltest maxIter:(NSInteger)maxIter
	   nlres:(double (^)(NSInteger i, NSInteger n, double *Bvec, NSInteger *notcomp))nlres nljac:(void (^)(NSInteger i, NSInteger n, double *Bvec, double *X))nljac;

/**
 Rayleigh quotient minimisation by conjugate gradients

 Algorithm 25 of Nash (1979). The quotient is defined by x(T)Ax / x(T)Bx, where A and B are square matrices and B is positive definite.

 @param n
   input:  order of the square matrices A and B

 @param A
   input: square matrix in numerator of the quotient

 @param B
   input: square matrix in denominator of the quotient

 @param X
   i/o: approximation to minimal eigenvector

 @param ipr
   i/o: limit on number of iterations/ actual number

 @param rq
   output: Rayleigh quotient at minimum (eigenvalue)

 @return
   mode flag (as defined in MathSupport.h):
     - RQMCG_SUCCESS: Successful termination.
     - RQMCG_ITERS_EXCEEDED: Number of iterations exceeded.
     - RQMCG_SINGULAR_B: The denominator matrix in the Rayleigh quotient is singular.
     - RQMCG_ZERO_DETERMINANT: The determinant of the Rayleigh quotient is zero.
 */
+(NSInteger)rqmcg:(NSInteger)n A:(double **)A B:(double **)B X:(double *)X ipr:(NSInteger *)ipr rq:(double *)rq;

/**
 Variable metric minimiser

 Algorithm 21 of Nash (1979)

 @param n
   number of parameters in function to be minimized

 @param Bvec
   input: initial guess; output: best estimate of solution

 @param Fmin
   function value at the minimum

 @param reltest
   user-initialized convergence tolerance for linear search

 @param fminfn
   code block for function (NSInteger n, double *Bvec, NSInteger *notcomp)

 @param fmingr
   code block for gradient (NSInteger n, double *Bvec, double *g))

 @return
   - VMMIN_SUCCESS: Successful termination.
   - VMMIN_BAD_INITIAL: The initial solution guess supplied to the method is invalid.

*/
+(NSInteger)vmmin:(NSInteger)n Bvec:(double *)Bvec Fmin:(double *)Fmin reltest:(double)reltest
	 fminfn:(double (^)(NSInteger n, double *Bvec, NSInteger *notcomp))fminfn fmingr:(void (^)(NSInteger n, double *Bvec, double *g))fmingr;

/**
 SUBROUTINE G1 (A,B,CTERM,STERM,SIG)

 COMPUTE ORTHOGONAL ROTATION MATRIX..

 The original version of this code was developed by
 Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratoryc  1973 JUN 12, and published in the book
 "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.

 Revised FEB 1995 to accompany reprinting of the book by SIAM.

 COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))
 (-S,C)         (-S,C)(B)   (   0          )
 COMPUTE SIG = SQRT(A**2+B**2)
 SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT
 SIG MAY BE IN THE SAME LOCATION AS A OR B .
 */
+(void)G1withA:(double)a andB:(double)b andCterm:(double *)cterm andSterm:(double *)sterm andSig:(double *)sig;

/**
 SUBROUTINE G2    (CTERM,STERM,X,Y)

 APPLY THE ROTATION COMPUTED BY G1 TO (X,Y).

 The original version of this code was developed by
 Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
 1972 DEC 15, and published in the book
 "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
 Revised FEB 1995 to accompany reprinting of the book by SIAM.
 */
+(void)G2withCterm:(double)cterm andSterm:(double)sterm andX:(double *)x andY:(double *)y;

/**
 Algorithm NNLS (Non-negative least-squares)
 Given an m by n matrix A, and an m-vector B, computes an n-vector X, that solves the least squares problem A * X = B, subject to X>=0

 @param a
 On entry, a[ 0... N ][ 0 ... M ] contains the M by N matrix A.
 On exit, a[][] contains the product matrix Q*A, where Q is an m by n orthogonal matrix generated implicitly by this function.

 @param m
 Matrix dimension

 @param n
 Matrix dimension

 @param b
 On entry, b[] must contain the m-vector B.
 On exit, b[] contains Q*B

 @param x
 On exit, x[] will contain the solution vector

 @param rnorm
 On exit, rnorm contains the Euclidean norm of the residual vector.

 @param w
 An n-array of working space, wp[].
 On exit, wp[] will contain the dual solution vector. wp[i]=0.0 for all i in set p and wp[i]<=0.0 for all i in set z.

 @param zz
 An m-array of working space, zz[].

 @param index
 An n-array of working space, index[].

 @return
 True if succesful. False if iteration count exceeded 3*N or initial conditions are invalid.
 */
+(BOOL)nnlsWithConstraintMatrix:(double **)a
				andRowDimension:(NSInteger)m
			 andColumnDimension:(NSInteger)n
			andConstraintVector:(double *)b
			  andSolutionVector:(double *)x
				andResidualNorm:(double *)rnorm
		  andDualSolutionVector:(double *)w
		 andDoubleWorkingVector:(double *)zz
	   andUIntegerWorkingVector:(NSInteger *)index;

+(void)simplx:(double **)a m:(NSInteger)m n:(NSInteger)n m1:(NSInteger)m1 m2:(NSInteger)m2 m3:(NSInteger)m3 icase:(NSInteger *)icase izrov:(NSInteger *)izrov iposv:(NSInteger *)iposv;

/**
 Cholesky decomposition of a positive definite matrix.

 @param a
 On entry, a[0..n][0..n] contains a symmetric square matrix
 On exit, the Cholesky decomposition of a if the input matrix is positive definite, else output is undefined.

 @param n
 atrix dimension

 @return
 YES if successful (matrix is positive definite), NO if the decomposition cannot be computed.
 */
+(BOOL)choldc:(double **)a dimension:(NSUInteger)n;

+(void)polint:(double *)xa ya:(double *)ya n:(NSInteger)n x:(double)x y:(double *)y dy:(double *)dy;

+(void)polin2:(double *)x1a x2a:(double *)x2a ya:(double **)ya m:(NSInteger)m n:(NSInteger)n x1:(double)x1 x2:(double)x2 y:(double *)y dy:(double *)dy;

/**
 Algorithm LDP (Least distance programming)
 Given an m by n matrix G, and an m-vector H, computes an n-vector X that solves the least squares problem minimize ||x||, subject to G X >= H

 @param g
 On entry g[][] contains the m by n matrix G.  Contents are unmodified on exit.

 @param m
 Matrix dimension.

 @param n
 Matrix dimension.

 @param h
 On entry, h[] must contain the m-vector H

 @param x
 On entry, x[] is an n-vector that need not be initialized.
 On exit, x[] contains the solution vector, if the computation is successful.

 @param xnorm
 On exit, the norm of the solution vector, if the computation is successful.

 @param w
 On entry, a vector of length (m+2)*(n+1) + 2*m of working storage.

 @param index
 On entry, an m-vector of working storage.

 @return
 YES if successful, NO if NNLS iteration limit exceeded or if inequality constraints are incompatible
 */
+(BOOL)ldpWithConstraintMatrix:(double **)g
               andRowDimension:(NSInteger)m
            andColumnDimension:(NSInteger)n
           andConstraintVector:(double *)h
             andSolutionVector:(double *)x
         andSolutionVectorNorm:(double *)xnorm
        andDoubleWorkingVector:(double *)w
        andDoubleWorkingMatrix:(double **)e
       andIntegerWorkingVector:(NSInteger *)index;

@end
