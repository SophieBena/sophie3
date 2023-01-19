//
//  MathSupport.m
//  PhasePlot
//
//  Source file for functions translated from algorithms published in:
//
//  Lawson, Charles L, and Hanson, Richard J. (1974)
//  Solving Least Squares Problems
//  Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 340 pp
//
//  Created by Mark Ghiorso on 6/25/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "MathSupport.h"
#import "DoubleVector.h"
#import "IntegerVector.h"

@implementation MathSupport

#undef MAX
#undef MIN

#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define MIN(a,b)  ((a) < (b) ? (a) : (b))

#define SQUARE(x) ((x)*(x))

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/* ==========================================================================
 Algorithm HFTI (page 81-82)
 Solution of a least-squares problem with possible rank deficiency.

 The LS problem is specified by designating the matrices a and b
 which are stored as pointers to pointers to double. Note that the
 LS solution can be provided for multiple right-hand-side vectors.

 -----------n-----------    -----nb------      -----nb------
 |                     |    |           |      |           |
 |                     |    |           |      |           |
 |                     |    |           |      |           |
 |                     |    |           |      |           |
 |                     |    |           |      |           |
 |                     |    |           |      |           |
 m         a           |    n     x     |  =   m     b     |
 |                     |    |           |      |           |
 |                     |    |           |      |           |
 |                     |    |           |      |           |
 |                     |    |           |      |           |
 |                     |    |           |      |           |
 |                     |    |           |      |           |
 -----------------------    -------------      -------------

 The Method involves Householder decomposition of the columns of a,
 accompanied by suitable permutations (P) of the elements to render the
 resulting matrix (R) upper triangular with strictly decreasing diagonal
 elements. The resulting transformations are in turn applied to the columns
 of b, i.e.

 Q a P = R   Q b  = c.

 The diagonal elements of R are then evaluated against an absolute
 tolerance parameter (tau) to determine the pseudorank (k) of the
 resulting solution. Thus:

 ------------             ------
 | R11  R12 |             | c1 |
 Q a P = R = |          |   Q b = c = |    |
 | 0    R22 |             | c2 |
 ------------             ------

 where R11 is k by k; R12 k by n-k, R22 m-k by n-k, c1 is k by nb
 and c2 is m-k by nb. If the problem is rank deficient, Householder
 tranformations (K) are applied to the rows of [ R11 : R12 ] to yield
 a matrix W:

 [ R11 : R12 ] K = [ W : 0 ] .

 If the solution is full rank (k = n), W is simply R and K is the
 identity matrix. A solution (y1) is then computed (W y1 = c1) to the stable
 subproblem, where if the system is full rank, c1 is taken to be c. The
 final (possibly rank deficient) solution is computed by setting y2 = 0
 (if k < n, otherwise y2 has zero length), and solving the system:

 ------
 | y1 |
 x = P K |    | = P K y .
 | y2 |
 ------

 Note that the orthogonal transformation (K) is applied to the columns of
 the solution vector even though it was derived from the rows of R.

 In the following algorithm, we let the solution vector(s) and the
 right-hand-side vector(s) occupy the same storage, denoted b.

 rnrom is the residual norm (= || c2 ||), p is a vector representation
 of the permutation matrix P, and h and g are vectors storing quantities
 related to the Q and K orthogonal decompositions, respectively. The
 remainder of the information related to the orthogonal transformations
 Q and K, is returned in the storage for a.

 ========================================================================== */

+(void)hfti:(double **)a m:(NSInteger)m n:(NSInteger)n b:(double **)b nb:(NSInteger)nb tau:(double)tau k:(NSInteger *)k rnorm:(double *)rnorm h:(double *)h g:(double *)g p:(NSInteger *)p {
	double hBar = 0.0, tmp;
	NSInteger i, j, l, lambda = 0;

	NSInteger mu = MIN(m, n);

	for (j=0; j<mu; j++) {
		if (j > 0) {
			for (l=j; l<n; l++) h[l] -= SQUARE(a[j-1][l]);
			for (l=j, lambda=j; l<n; l++) if (h[l] > h[lambda]) lambda = l;
		}
		if (j == 0 || (hBar + 0.001*h[lambda]) <= hBar) {
			for (l=j; l<n; l++) {
				for (i=j, h[l]=0.0; i<m; i++) h[l] += SQUARE(a[i][l]);
			}
			for (l=j, lambda=j; l<n; l++) if (h[l] > h[lambda]) lambda = l;
			hBar = h[lambda];
		}
		p[j] = lambda;
		if (p[j] != j) {
			for (i=0; i<m; i++) {
				tmp = a[i][j]; a[i][j] = a[i][lambda]; a[i][lambda] = tmp;
			}
			h[lambda] = h[j];
		}
		[MathSupport householderColCol:HOUSEHOLDER_CALC_MODE_H1 p:j l:j+1 m:m-1 v:a vCol:j h:&h[j] c:a cColStart:j+1 cColEnd:n-1];
		//householderColCol(           HOUSEHOLDER_CALC_MODE_H1,  j,  j+1,  m-1,  a,     j,  &h[j],  a,          j+1,        n-1);
		[MathSupport householderColCol:HOUSEHOLDER_CALC_MODE_H2 p:j l:j+1 m:m-1 v:a vCol:j h:&h[j] c:b cColStart:0 cColEnd:nb-1];
		//householderColCol(           HOUSEHOLDER_CALC_MODE_H2,  j,  j+1,  m-1,  a,     j,  &h[j],  b,          0,        nb-1);
	}

	if (tau < 0.0) (*k) = (NSInteger) fabs(tau); /* special case -tau = known pseudorank */
	else {
		for (i=0, *k=0; i<n; i++) if (fabs(a[i][i]) > tau) (*k)++;
	}

	for (l=0; l<nb; l++) {
		for (i=(*k), rnorm[l]=0.0; i<m; i++) rnorm[l] += SQUARE(b[i][l]);
		rnorm[l] = sqrt(rnorm[l]);
	}

	if (*k == 0) {
		for (i=0; i<n; i++) { for (j=0; j<nb; j++) b[i][j] = 0.0; } return; }

	if (*k < n) for (i=(*k)-1; i>=0; i--)
		[MathSupport householderRowRow:HOUSEHOLDER_CALC_MODE_H1 p:i l:*k m:n-1 v:a vRow:i h:&g[i] c:a cRowStart:0 cRowEnd:i-1];
		//householderRowRow(           HOUSEHOLDER_CALC_MODE_H1,  i,  *k,  n-1,  a,     i,  &g[i],  a,          0,        i-1);

	for (l=0; l<nb; l++) {
		b[*k-1][l] /= a[*k-1][*k-1];
		for (i=(*k)-2; i>=0; i--) {
			for (j=i+1; j<(*k); j++) b[i][l] -= a[i][j]*b[j][l];
			b[i][l] /= a[i][i];
		}
		if (*k < n) for (i=(*k); i<n; i++) b[i][l] = 0.0;
	}
	if (*k < n) for (i=0; i<(*k); i++)
		[MathSupport householderRowCol:HOUSEHOLDER_CALC_MODE_H2 p:i l:*k m:n-1 v:a vRow:i h:&g[i] c:b cColStart:0 cColEnd:nb-1];
		//householderRowCol(           HOUSEHOLDER_CALC_MODE_H2,  i,  *k,  n-1,  a,     i,  &g[i],  b,          0,        nb-1);

	for (j=mu-1; j>=0; j--) if (p[j] != j) {
		for (l=0; l<nb; l++) {
			tmp = b[j][l]; b[j][l] = b[p[j]][l]; b[p[j]][l] = tmp;
		}
	}

}

/* ==========================================================================
 Algorithm H12 (page 57)
 Householder decomposition of a vector and application of the resulting
 transformation to column or row vectors of a matrix.

 If mode = HOUSEHOLDER_CALC_MODE_H1

 Householder decomposition of a vector v into:

 -        -
 |  v[0]  |
 |    .   |
 | v[p-1] |
 |  y[p]  |
 Q v  =  | v[p+1] |  =  y (stored in v)
 |    .   |
 | v[l-1] |
 |    0   |
 |    .   |
 |    0   |
 -        -

 with y[p] = -s*(v[p]^2 + sum(l<=i<=m) v[i]^2)^1/2. On return, the pth
 element of v is used to store s, and the storage location h is used to
 contain the quantity v[p] - s. The vector v may be stored in a column
 of the input matrix (in which case one uses the householderCol* routines
 and indicates the designated column, vCol), or as a row of the input matrix
 (in which case one uses the householderRow* routines and indicates the row,
 vRow). The elements of v are in either case stored in a pointer to pointer
 to double matrix.

 If mode = HOUSEHOLDER_CALC_MODE_H1 or HOUSEHOLDER_CALC_MODE_H2

 The transformation is applied to the set of vectors stored in the matrix c.
 These may be column vectors beginning with cColStart and ending with
 cColEnd, for which the householder*Col routines are suitable. Alternatively,
 the transformation may be applied to row vectors beginning with cRowStart
 and ending with cRowEnd, in which case the householder*Row routines are
 appropriate. Note that if end < start, c is ignored (in effect the
 identity transformation is performed). The information used to tranform
 c is provided in the elements of v and the scaler h.

 For both modes the user must input the pivot element p, and the indices
 of the elements of the vector v that are to be zeroed (H1 mode) or
 transformed in c (H1 mode, if end >= start, or H2 mode).

 Note that all vectors and matrices have zero-based indices and that
 on input 0 <= p < l and that l <= m. The elemnets of v indexed l through
 m are zero in H1 or transformed in H1/2.

 ========================================================================== */

+(void)householderColCol:(NSInteger)mode p:(NSInteger)p l:(NSInteger)l m:(NSInteger)m v:(double **)v vCol:(NSInteger)vCol h:(double *)h c:(double **)c cColStart:(NSInteger)cColStart cColEnd:(NSInteger)cColEnd {
	double b, s;
	NSInteger i, j;

	if (0 > p || p >= l || l > m) return;
	switch (mode) {

		case HOUSEHOLDER_CALC_MODE_H1:
			for (i=l, s=SQUARE(v[p][vCol]); i<=m; i++) s += SQUARE(v[i][vCol]);
			s = sqrt(s);
			if (v[p][vCol] > 0.0) s *= -1.0;
			*h = v[p][vCol] - s; v[p][vCol] = s;

		case HOUSEHOLDER_CALC_MODE_H2:
			b = v[p][vCol]*(*h);
			if (b == 0.0) return;

			for (j=cColStart; j<=cColEnd; j++) {
				for (i=l, s=c[p][j]*(*h); i<=m; i++) s += c[i][j]*v[i][vCol];
				s /= b;
				c[p][j] += s*(*h);
				for (i=l; i<=m; i++) c[i][j] += s*v[i][vCol];
			}

	} /* end switch */
}

+(void)householderColRow:(NSInteger)mode p:(NSInteger)p l:(NSInteger)l m:(NSInteger)m v:(double **)v vCol:(NSInteger)vCol h:(double *)h c:(double **)c cRowStart:(NSInteger)cRowStart cRowEnd:(NSInteger)cRowEnd {
	double b, s;
	NSInteger i, j;

	if (0 > p || p >= l || l > m) return;
	switch (mode) {

		case HOUSEHOLDER_CALC_MODE_H1:
			for (i=l, s=SQUARE(v[p][vCol]); i<=m; i++) s += SQUARE(v[i][vCol]);
			s = sqrt(s);
			if (v[p][vCol] > 0.0) s *= -1.0;
			*h = v[p][vCol] - s; v[p][vCol] = s;

		case HOUSEHOLDER_CALC_MODE_H2:
			b = v[p][vCol]*(*h);
			if (b == 0.0) return;

			for (j=cRowStart; j<=cRowEnd; j++) {
				for (i=l, s=c[j][p]*(*h); i<=m; i++) s += c[j][i]*v[i][vCol];
				s /= b;
				c[j][p] += s*(*h);
				for (i=l; i<=m; i++) c[j][i] += s*v[i][vCol];
			}

	} /* end switch */
}

+(void)householderRowCol:(NSInteger)mode p:(NSInteger)p l:(NSInteger)l m:(NSInteger)m v:(double **)v vRow:(NSInteger)vRow h:(double *)h c:(double **)c cColStart:(NSInteger)cColStart cColEnd:(NSInteger)cColEnd {
	double b, s;
	NSInteger i, j;

	if (0 > p || p >= l || l > m) return;
	switch (mode) {

		case HOUSEHOLDER_CALC_MODE_H1:
			for (i=l, s=SQUARE(v[vRow][p]); i<=m; i++) s += SQUARE(v[vRow][i]);
			s = sqrt(s);
			if (v[vRow][p] > 0.0) s *= -1.0;
			*h = v[vRow][p] - s; v[vRow][p] = s;

		case HOUSEHOLDER_CALC_MODE_H2:
			b = v[vRow][p]*(*h);
			if (b == 0.0) return;

			for (j=cColStart; j<=cColEnd; j++) {
				for (i=l, s=c[p][j]*(*h); i<=m; i++) s += c[i][j]*v[vRow][i];
				s /= b;
				c[p][j] += s*(*h);
				for (i=l; i<=m; i++) c[i][j] += s*v[vRow][i];
			}

	} /* end switch */
}

+(void)householderRowRow:(NSInteger)mode p:(NSInteger)p l:(NSInteger)l m:(NSInteger)m v:(double **)v vRow:(NSInteger)vRow h:(double *)h c:(double **)c cRowStart:(NSInteger)cRowStart cRowEnd:(NSInteger)cRowEnd {
	double b, s;
	NSInteger i, j;

	if (0 > p || p >= l || l > m) return;
	switch (mode) {

		case HOUSEHOLDER_CALC_MODE_H1:
			for (i=l, s=SQUARE(v[vRow][p]); i<=m; i++) s += SQUARE(v[vRow][i]);
			s = sqrt(s);
			if (v[vRow][p] > 0.0) s *= -1.0;
			*h = v[vRow][p] - s; v[vRow][p] = s;

		case HOUSEHOLDER_CALC_MODE_H2:
			b = v[vRow][p]*(*h);
			if (b == 0.0) return;

			for (j=cRowStart; j<=cRowEnd; j++) {
				for (i=l, s=c[j][p]*(*h); i<=m; i++) s += c[j][i]*v[vRow][i];
				s /= b;
				c[j][p] += s*(*h);
				for (i=l; i<=m; i++) c[j][i] += s*v[vRow][i];
			}

	} /* end switch */
}

#define A1  1.5
#define A2 -0.25

+(NSInteger)min1d:(double *)bb st:(double *)st reltest:(double)reltest ifn:(NSInteger *)ifn fnminval:(double *)fnminval fn1d:(double (^)(double, NSInteger *))fn1d {
	double fii, s0, s1, s2, tt0, tt1, tt2, x0, x1, x2, xii;
	NSInteger notcomp, tripleok;
	NSInteger ifnMax = *ifn;
	double BIG = sqrt(DBL_MAX);

	*ifn = 0;
	x1   = *bb;
	s0 = fn1d(x1, &notcomp); if (notcomp) return MIN1D_BAD_INITIAL; (*ifn)++;

	do {
		x0 = x1;
		*bb = x0;
		x1 = x0 + (*st);
		s1 = fn1d(x1, &notcomp); if (notcomp) s1 = BIG; (*ifn)++;
		if (*ifn > ifnMax) return MIN1D_ITERS_EXCEEDED;
		tripleok = FALSE;
		if (s1 < s0) {
			do {
				*st *= A1;
				x2 = x1 + (*st);
				s2 = fn1d(x2, &notcomp); if (notcomp) s2 = BIG; (*ifn)++;
				if (*ifn > ifnMax) return MIN1D_ITERS_EXCEEDED;
				if (s2 < s1) {
					s0 = s1; s1 = s2;
					x0 = x1; x1 = x2;
				} else {
					tripleok = TRUE;
				}
			} while (!tripleok);
		} else {
			*st *= A2;
			tt2 = s0; s0 = s1; s1 = tt2;
			tt2 = x0; x0 = x1; x1 = tt2;
			do {
				x2 = x1 + (*st);
				s2 = fn1d(x2, &notcomp); if (notcomp) s2 = BIG; (*ifn)++;
				if (*ifn > ifnMax) return MIN1D_ITERS_EXCEEDED;
				if (s2 < s1) {
					s0 = s1; s1 = s2; x0 = x1; x1 = x2;
					*st *= A1;
				} else {
					tripleok = TRUE;
				}
			} while (!tripleok);
		}
		tt0 = x0 - x1;
		tt1 = (s0 - s1)*(*st); tt2 = (s2 - s1)*tt0;
		if (tt1 != tt2) {
			*st = 0.5*( (tt2/(tt2 - tt1))*tt0 - (tt1/(tt2 - tt1))*(*st) );
			xii = x1 + (*st);
			if (fabs(xii-x1) > reltest) {
				fii = fn1d(xii, &notcomp); if (notcomp) fii = BIG; (*ifn)++;
				if (*ifn > ifnMax) return MIN1D_ITERS_EXCEEDED;
				if (fii < s1) {
					s1 = fii; x1 = xii;
				}
			}
		}
		s0 = s1;
	} while (fabs(*bb - x1) > reltest);

	*fnminval = s1;
	// (void) fn1d(*bb, &notcomp);
	BOOL debug = [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.MAIN.VERBOSE"];
	if (debug) NSLog(@"min1d: s0-s1 %13.6g s2-s1 %13.6g x0 %13.6g x1 %13.6g x2 %13.6g", s0-s1, s2-s1, x0, x1, x2);
	return MIN1D_SUCCESS;
}

#undef A1
#undef A2

#define STEPREDN  0.2
#define ACCTOL    0.0001

+(NSInteger)vmmin:(NSInteger)n Bvec:(double *)Bvec Fmin:(double *)Fmin reltest:(double)reltest
	 fminfn:(double (^)(NSInteger, double *, NSInteger *))fminfn fmingr:(void (^)(NSInteger, double *, double *))fmingr {
	double **B, *c, D1, D2, f, *g, gradproj, s, steplength, *t, *X;
	NSInteger accpoint, count, funcount, gradcount, i, ilast, j, notcomp, iterations;

	f = fminfn(n, Bvec, &notcomp); if (notcomp) return VMMIN_BAD_INITIAL;

	B = (double **) malloc((unsigned) n*sizeof(double *));
	for (i=0; i<n; i++) B[i] = (double *) malloc((unsigned) n*sizeof(double));
	c = (double *) malloc((unsigned) n*sizeof(double));
	g = (double *) malloc((unsigned) n*sizeof(double));
	t = (double *) malloc((unsigned) n*sizeof(double));
	X = (double *) malloc((unsigned) n*sizeof(double));

	*Fmin = f; funcount = 1; gradcount = 1;
	fmingr(n, Bvec, g); ilast = gradcount;
	iterations = 0;

	do {
		iterations++; if (iterations > 10000.0) {
            for (i=0; i<n; i++) free(B[i]); free(B); free(c); free(g); free(t); free(X);
            return VMMIN_BAD_INITIAL;
        }
		if (ilast == gradcount)
			for (i=0; i<n; i++) { for (j=0; j<n; j++) B[i][j] = 0.0; B[i][i] = 1.0; }
		for (i=0; i<n; i++) { X[i] = Bvec[i]; c[i] = g[i]; }
		for (i=0, gradproj=0.0; i<n; i++) {
			for (j=0, s=0.0; j<n; j++) s -= B[i][j]*g[j];
			t[i] = s; gradproj += s*g[i];
		}
		if (gradproj < 0.0) {
			steplength = 1.0; accpoint = FALSE;
			do {
				for (i=0, count=0; i<n; i++) {
					Bvec[i] = X[i] + steplength*t[i];
					if (fabs(X[i]-Bvec[i]) < reltest) count++;
				}
				if (count < n) {
					f = fminfn(n, Bvec, &notcomp); funcount++;
					accpoint = (!notcomp) && (f <= *Fmin+gradproj*steplength*ACCTOL);
					if (!accpoint) steplength *= STEPREDN;
				}
			} while (count != n && !accpoint);
			if (count < n) {
				*Fmin = f;
				fmingr(n, Bvec, g); gradcount++;
				for (i=0, D1=0.0; i<n; i++)
				{ t[i] *= steplength; c[i] = g[i] - c[i]; D1 += t[i]*c[i]; }
				if (D1 > 0.0) {
					for (i=0, D2=0.0; i<n; i++) {
						for (j=0, s=0.0; j<n; j++) s += B[i][j]*c[j];
						X[i] = s; D2 += s*c[i];
					}
					D2 = 1.0 + D2/D1;
					for (i=0; i<n; i++) for (j=0; j<n; j++)
						B[i][j] -= (t[i]*X[j] + X[i]*t[j] - D2*t[i]*t[j])/D1;
				} else ilast = gradcount;
			} else {
				if (ilast < gradcount) { count = 0; ilast = gradcount; }
			}
		} else { count = 0; ilast = gradcount; }
	} while (count != n || ilast != gradcount);

	for (i=0; i<n; i++) free(B[i]); free(B); free(c); free(g); free(t); free(X);
	return VMMIN_SUCCESS;
}

#undef STEPREDN
#undef ACCTOL

#define DEC  0.4
#define INC 10.0
#define PHI  1.0

/******************************************************************************
 * Algorithm 7.  Choleski decomposition in compact storage
 ******************************************************************************/

+(void)Choldcmp:(NSInteger)n a:(double *)a singmat:(NSInteger *)singmat
{
	NSInteger i, j, k, m, q;
	double s;

	*singmat = FALSE;
	for (j=1; j<=n; j++) {                                         /* STEP  1 */
		q = j*(j+1)/2;                                             /* STEP  2 */
		if (j > 1) {                                               /* STEP  3 */
			for (i=j; i<=n; i++) {                                 /* STEP  4 */
				m = i*(i-1)/2 + j; s = a[m-1];
				for (k=1; k<=(j-1); k++) s -= a[m-k-1]*a[q-k-1];
				a[m-1] = s;
			}
		}
		if (a[q-1] <= 0.0) {                                       /* STEP  5 */
			*singmat = TRUE;
			a[q-1] = 0.0;
		}
		s = sqrt(a[q-1]);                                          /* STEP 7 */
		for (i=j; i<=n; i++) {                                     /* STEP 8 */
			m = i*(i-1)/2 + j;
			a[m-1] = (s == 0.0) ? 0.0 : a[m-1]/s;
		}
	}
}

/******************************************************************************
 * Algorithm 8.
 * Choleski back-substitution
 ******************************************************************************/

+(void)Cholback:(NSInteger)n a:(double *)a x:(double *)x
{
	NSInteger i, j, q;

	x[0] = (a[0] == 0.0) ? 0.0 : x[0]/a[0];                          /* STEP  1 */
	if (n > 1) {                                                     /* STEP  2 */
		q = 1;                                                       /* STEP  3 */
		for (i=2; i<=n; i++) {                                       /* STEP  4 */
			for (j=1; j<=(i-1); j++) {q++; x[i-1] -= a[q-1]*x[j-1];} /* STEP  5 */
			q++;                                                     /* STEP  6 */
			x[i-1] = (a[q-1] == 0.0) ? 0.0 : x[i-1]/a[q-1];          /* STEP  7 */
		}                                                            /* STEP  8 */
	}
	                                                                 /* STEP  9 */
	x[n-1] = (a[n*(n+1)/2-1] == 0.0) ? 0.0 : x[n-1]/a[n*(n+1)/2-1];
	if (n > 1) {                                                     /* STEP 10 */
		for (i=n; i>=2; i--) {                                       /* STEPS 11/12 */
			q = i*(i-1)/2;
			for (j=1; j<=(i-1); j++) x[j-1] -= x[i-1]*a[q+j-1];      /* STEP 13 */
			x[i-2] = (a[q-1] == 0.0) ? 0.0 : x[i-2]/a[q-1];          /* STEP 14 */
		}                                                            /* STEP 15 */
	}                                                                /* STEP 16 */
}

+(NSInteger)modmrt:(NSInteger)n m:(NSInteger)m Bvec:(double *)Bvec Fmin:(double *)Fmin reltest:(double)reltest maxIter:(NSInteger)maxIter
	   nlres:(double (^)(NSInteger, NSInteger, double *, NSInteger *))nlres nljac:(void (^)(NSInteger, NSInteger, double *, double *))nljac {

	/* The residual function returns the residual of the ith equation and a
	 flag (notcomp) which is set to true if the residual cannot be computed.
     The gradient function returns a vector corresponding to the ith row
	 of the Jacobian of the system of nonlinear equations                   */

	double *a, *c, *delta, *res, *v, *X;
	double lambda, p;
	NSInteger count, i, ifn, igrad, j, k, nn2, q;
	NSInteger notcomp, singmat, calcmat;

	*Fmin   = DBL_MAX;                                               /* STEP  0 */
	lambda  = 0.0001;
	ifn     = 0;
	igrad   = 0;
	calcmat = TRUE;
	nn2     = n*(n+1)/2;

	if (n == 0) return MODMRT_SUCCESS;

	a     = (double *) calloc((unsigned) nn2, sizeof(double));
	c     = (double *) calloc((unsigned) nn2, sizeof(double));
	delta = (double *) calloc((unsigned) n,   sizeof(double));
	res   = (double *) calloc((unsigned) m,   sizeof(double));
	v     = (double *) calloc((unsigned) n,   sizeof(double));
	X     = (double *) calloc((unsigned) n,   sizeof(double));

	for (i=0, p=0.0; i<m; i++) {                                        /* STEP  1 */
		res[i] = nlres(i, n, Bvec, &notcomp);
		if (notcomp) {
			free(a); free(c); free(delta); free(res); free(v); free(X);
			return MODMRT_BAD_INITIAL;
		}
		p += SQUARE(res[i]);
	}
	ifn++; *Fmin = p; count = 0;

	while (count < n) {                                                 /* STEP  2 */
		if (calcmat) {
			igrad++;                                                    /* STEP  3 */
			for (j=0; j<nn2; j++) a[j] = 0.0;
			for (j=0; j<n;   j++) v[j] = 0.0;
			for (i=0; i<m; i++) {                                       /* STEP  4 */
				nljac(i, n, Bvec, X);
				for (j=1; j<=n; j++) {
					v[j-1] = v[j-1] + X[j-1]*res[i];
					q = j*(j-1)/2;
					for (k=1; k<=j; k++) a[q+k-1] += X[j-1]*X[k-1];
				}
			}
			for (j=0; j<nn2; j++) c[j] = a[j];                          /* STEP  5 */
			for (j=0; j<n;   j++) X[j] = Bvec[j];
		}
		for (j=1; j<=n; j++) {                                          /* STEP  6 */
			q = j*(j+1)/2;
			a[q-1] = c[q-1]*(1.0+lambda) + PHI*lambda;
			delta[j-1] = -v[j-1];
			if (j > 1) for (i=1; i<=(j-1); i++) a[q-i-1] = c[q-i-1];
		}
		notcomp = FALSE;
		[self Choldcmp:n a:a singmat:&singmat];                         /* STEP 7 */
		if (!singmat) {
			[self Cholback:n a:a x:delta];                              /* STEP  8 */
			for (i=0, count = 0; i<n; i++) {                            /* STEP  9 */
				Bvec[i] = X[i] + delta[i];
				if (fabs(Bvec[i]-X[i]) <= reltest) count++;
			}
			if (count < n) {                                            /* STEP 10 */
				p = 0.0; i = 0;
				do {
					res[i] = nlres(i, n, Bvec, &notcomp);
					if (!notcomp) p += SQUARE(res[i]);
					i++;
				} while (!notcomp && i < m); ifn++;
				if (ifn > maxIter) {
					free(a); free(c); free(delta); free(res); free(v); free(X);
					return MODMRT_ITERS_EXCEEDED;
				}
			}
		}
		if (count < n) {
			if (!singmat && !notcomp && p < *Fmin) {
				lambda *= DEC;
				*Fmin = p;
				calcmat = TRUE;
			} else {
				lambda *= INC;
				if (lambda < SQUARE(DBL_EPSILON)) lambda = DBL_EPSILON;
				calcmat = FALSE;
			}
		}
	}

	free(a); free(c); free(delta); free(res); free(v); free(X);
	return MODMRT_SUCCESS;
}

#undef DEC
#undef INC
#undef PHI

+(NSInteger)rqmcg:(NSInteger)n A:(double **)A B:(double **)B X:(double *)X ipr:(NSInteger *)ipr rq:(double *)rq {
	NSInteger conv, count, i, itn, itlimit, j;
	double *avec, *bvec, *yvec, *zvec, *g, *t;
	double beta, d, gg, pa, pn, step, ta, tabt, tat,
	tbt, tol, u, v, w, xat, xax, xbt, xbx;

	avec = (double *) malloc((unsigned) n*sizeof(double));
	bvec = (double *) malloc((unsigned) n*sizeof(double));
	yvec = (double *) malloc((unsigned) n*sizeof(double));
	zvec = (double *) malloc((unsigned) n*sizeof(double));
	g    = (double *) malloc((unsigned) n*sizeof(double));
	t    = (double *) malloc((unsigned) n*sizeof(double));

	itlimit = *ipr; conv = FALSE;
	*ipr = 0; tol = SQUARE(n)*DBL_EPSILON; pa = DBL_MAX;

	while (*ipr <= itlimit && !conv) {                                           /* STEP  1 */
		for (i=0; i<n; i++) for (j=0, avec[i]=0.0, bvec[i]=0.0; j<n; j++)
		{ avec[i] += A[i][j]*X[j]; bvec[i] += B[i][j]*X[j]; } (*ipr)++;

		for (i=0, xax=0.0, xbx=0.0; i<n; i++)                                    /* STEP  2 */
		{ xax += X[i]*avec[i]; xbx += X[i]*bvec[i]; }
		if (xbx <= tol) {                                                        /* STEP  3 */
			free(avec); free(bvec); free(yvec); free(zvec); free(g); free(t);
			return RQMCG_SINGULAR_B;
		}
		*rq = xax/xbx;                                                           /* STEP  4 */

		if (*rq < pa) {                                                          /* STEP  5 */
			pa = *rq;
			for (i=0, gg=0.0; i<n; i++)
			{ g[i] = 2.0*(avec[i] - (*rq)*bvec[i])/xbx; gg += g[i]*g[i]; }
			if (gg > tol) {                                                      /* STEP  7 */
				for (i=0; i<n; i++) t[i] = -g[i];                                /* STEP  8 */
				itn = 0;                                                         /* STEP  9 */
				do {
					itn++;                                                       /* STEP 10 */
					for (i=0; i<n; i++) for (j=0, yvec[i]=0.0, zvec[i]=0.0; j<n; j++)
					{ yvec[i] += A[i][j]*t[j]; zvec[i] += B[i][j]*t[j]; } (*ipr)++;
					for (i=0, tat = 0.0, tbt = 0.0, xat = 0.0, xbt = 0.0; i<n; i++) {
						xat += X[i]*yvec[i]; tat += t[i]*yvec[i];
						xbt += X[i]*zvec[i]; tbt += t[i]*zvec[i];
					}
					u = tat*xbt - xat*tbt; v = tat*xbx - xax*tbt;                /* STEP 12 */
					w = xat*xbx - xax*xbt; d = v*v - 4.0*u*w;
					if (d < 0.0) {                                               /* STEP 13 */
						free(avec); free(bvec); free(yvec); free(zvec); free(g); free(t);
						return RQMCG_ZERO_DETERMINANT;
					}
					d = sqrt(d);                                                 /* STEP 14 */
					step = (v > 0.0) ? -2.0*w/(v+d) : 0.5*(d-v)/u;
					for (i=0, count=0, xax=0.0, xbx=0.0; i<n; i++) {             /* STEP 15 */
						avec[i] += step*yvec[i]; bvec[i] += step*zvec[i];
						w = X[i]; X[i] = w + step*t[i];
						if (fabs(X[i]-w) < 10.0*DBL_EPSILON) count++;
						xax += X[i]*avec[i];  xbx += X[i]*bvec[i];
					}
					if (xbx <= tol) {                                            /* STEP 16 */
						free(avec); free(bvec); free(yvec); free(zvec); free(g); free(t);
						return RQMCG_SINGULAR_B;
					} else pn = xax/xbx;
					if (count < n && pn < *rq) {                         /* STEPS 17 and 18 */
						*rq = pn;                                                /* STEP 19 */
						for (i=0, gg=0.0; i<n; i++)
						{ g[i] = 2.0*(avec[i]-pn*bvec[i])/xbx; gg += g[i]*g[i]; }
						if (gg > tol) {                                          /* STEP 20 */
							for (i=0, xbt=0.0; i<n; i++) xbt += X[i]*zvec[i];    /* STEP 21 */
							for (i=0, tabt=0.0, beta=0.0; i<n; i++) {            /* STEP 22 */
								w = yvec[i] - pn*zvec[i]; tabt += t[i]*w;
								beta += g[i]*(w - g[i]*xbt);
							}
							beta /= tabt;                                        /* STEP 23 */
							for (i=0; i<n; i++) t[i] = beta*t[i] - g[i];
						}
					} else {
						if (itn == 1) conv = TRUE;
						itn = n + 1;
					}
				} while (itn < n && count != n && gg > tol && !conv);
			} else conv = TRUE;
		} else conv = TRUE;
		for (i=0, ta=0.0; i<n; i++) ta += SQUARE(X[i]);
		ta = 1.0/sqrt(ta);
		for (i=0; i<n; i++) X[i] *= ta;
	}

	free(avec); free(bvec); free(yvec); free(zvec); free(g); free(t);
	if (*ipr > itlimit) return RQMCG_ITERS_EXCEEDED;
	else                return RQMCG_SUCCESS;

}

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
+(void)G1withA:(double)a andB:(double)b andCterm:(double *)cterm andSterm:(double *)sterm andSig:(double *)sig {
	if (fabs(a) > fabs(b)) {
		double xr = b/a;
		double yr = sqrt(1.0 + xr*xr);
		*cterm = (a > 0.0) ? 1.0/yr : -1.0/yr;
		*sterm = (*cterm)*xr;
		*sig   = fabs(a)*yr;
	} else if (b != 0.0) {
		double xr = a/b;
		double yr = sqrt(1.0 + xr*xr);
		*sterm = (b > 0.0) ? 1.0/yr : -1.0/yr;
		*cterm = (*sterm)*xr;
		*sig   = fabs(b)*yr;
	} else {
		*sig   = 0.0;
		*cterm = 0.0;
		*sterm = 1.0;
	}
}

/**
 SUBROUTINE G2    (CTERM,STERM,X,Y)

 APPLY THE ROTATION COMPUTED BY G1 TO (X,Y).

 The original version of this code was developed by
 Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
 1972 DEC 15, and published in the book
 "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
 Revised FEB 1995 to accompany reprinting of the book by SIAM.
*/
+(void)G2withCterm:(double)cterm andSterm:(double)sterm andX:(double *)x andY:(double *)y {
	double xr = cterm*(*x) + sterm*(*y);
	*y = -sterm*(*x) + cterm*(*y);
    *x = xr;
}

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

@param wp
 An n-array of working space, wp[].
 On exit, wp[] will contain the dual solution vector. wp[i]=0.0 for all i in set p and wp[i]<=0.0 for all i in set z.

@param zzp
 An m-array of working space, zz[].

@param indexp
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
	   andUIntegerWorkingVector:(NSInteger *)index {

  	if (m <= 0 || n <= 0) return false;

	for (NSInteger i=0; i<n; i++) { x[i] = 0.0; index[i] = i; }
	NSInteger itmax = 3*n;
	NSInteger iter  = 0;
	NSInteger iz2   = n - 1;
	NSInteger iz1   = 0;
	NSInteger nsetp = 0;
	NSInteger npp1  = 0;

	BOOL (^finish)(BOOL,NSInteger) = ^(BOOL state, NSInteger npivot) {
		double sm = 0.0;
		if (npivot < m) for (NSInteger i=npivot; i<m; i++) sm += b[i]*b[i];
		else            for (NSInteger i=0; i<n; i++) w[i] = 0.0;
		*rnorm = sqrt(sm);
		return state;
	};

	BOOL mainLoop = true;
	while (mainLoop) {
		BOOL debug = [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.MAIN.VERBOSE"];
		if (debug) NSLog(@"npp1 %ld, nsetp %ld, iter %ld, iz2 %ld, iz1 %ld", npp1, nsetp, iter, iz2, iz1);

		// quit if all coefficients are already in the solution or if m cols of a have been triangularized.
		if (iz1 > iz2 || nsetp >= m) return finish(true, npp1);

		// compute components of the dual (negative gradient) vector w().
		for (NSInteger iz=iz1; iz<=iz2; iz++) {
			double sm = 0.0;
			for (NSInteger l=npp1; l<m; l++) sm += a[l][index[iz]]*b[l];
			w[index[iz]] = sm;
		}

		// find largest positive w(j).
		BOOL wLoop = true;
		NSInteger izmax = iz1;
		double up;
		while (wLoop) {
			double wmax = 0.0;
			for (NSInteger iz=iz1; iz<=iz2; iz++) {
				if (w[index[iz]] > wmax) {
					wmax  = w[index[iz]];
					izmax = iz;
				}
			}

			// if wmax <= 0.0 then terminate. This condition indicates satisfaction of the kuhn-tucker conditions.
			if (wmax <= 0.0)  return finish(true, npp1);

			// the sign of w(index[izmax]) is ok for index[izmax] to be moved to set p.
			// begin the transformation and check new diagonal element to avoid near linear dependence.
			// original FORTRAN code: call h12 (1,npp1,npp1+1,m,a(1,index[izmax]),1,up,dummy,1,1,0)
			double asave = a[npp1][index[izmax]];
			[MathSupport householderColCol:HOUSEHOLDER_CALC_MODE_H1 p:npp1 l:npp1+1 m:m-1 v:a vCol:index[izmax] h:&up c:NULL cColStart:0 cColEnd:-1];

			double unorm = 0.0;
			if (nsetp != 0) for (NSInteger l=0; l<nsetp; l++) unorm += a[l][index[izmax]]*a[l][index[izmax]];
			unorm = sqrt(unorm);

			if (fabs(a[npp1][index[izmax]]*0.01) > DBL_EPSILON*fabs(unorm)) {
				// col index[izmax] is sufficiently independent.  copy b into zz, update zz and solve for ztest ( = proposed new value for x(index[izmax]) ).
				// original FORTRAN code: call h12 (2,npp1,npp1+1,m,a(1,index[izmax]),1,up,zz,1,1,1)
				for (NSInteger l=0; l<m; l++) zz[l] = b[l];
				[MathSupport householderColRow:HOUSEHOLDER_CALC_MODE_H2 p:npp1 l:npp1+1 m:m-1 v:a vCol:index[izmax] h:&up c:&zz cRowStart:0 cRowEnd:0];

	            if (zz[npp1]/a[npp1][index[izmax]] > 0.0) break; // out of wLoop
			}

			// reject index[izmax] as a candidate to be moved from set z to set p.
			// restore a(npp1,index[izmax]), set w(index[izmax])=0., and loop back to test dual coeffs again.
			a[npp1][index[izmax]] = asave;
			w[index[izmax]]       = 0.0;
		} // end while on wloop

		// The index index(izmax) has been selected to be moved from set z to set p.
		// Update b, update indices, apply householder transformations to cols in new set z,  zero subdiagonal elts in col index(izmax),  set w(index(izmax))=0.
		for (NSInteger l=0; l<m; l++) b[l] = zz[l];
		NSInteger indexTemp = index[izmax];
		index[izmax] = index[iz1];
		index[iz1]   = indexTemp;
		nsetp        = npp1+1;
	    iz1++;
		npp1++;

		// call h12 (2,nsetp,npp1,m,a(1,indexTemp),1,up,a(1,index[jz]),1,mda,1)
		if (iz1 <= iz2) for (NSInteger jz=iz1; jz<=iz2; jz++)
			[MathSupport householderColCol:HOUSEHOLDER_CALC_MODE_H2 p:nsetp-1 l:npp1 m:m-1 v:a vCol:indexTemp h:&up c:a cColStart:index[jz] cColEnd:index[jz]];

		if (nsetp != m) for (NSInteger l=npp1; l<m; l++) a[l][indexTemp] = 0.0;
		w[indexTemp] = 0.0;

		//solve the triangular system. store the solution temporarily in zz().
		for (NSInteger l=0; l<nsetp; l++) {
			if (l != 0) for (NSInteger ii=0; ii<=(nsetp-(l+1)); ii++) zz[ii] -= a[ii][index[nsetp-(l+1)+1]]*zz[nsetp-(l+1)+1];
			zz[nsetp-(l+1)] /= a[nsetp-(l+1)][index[nsetp-(l+1)]];
		}

		BOOL secondaryLoop = true;
		while (secondaryLoop) {
			iter++;
			if (iter > itmax) return finish(false, npp1);

			// see if all new constrained coeffs are feasible. if not compute alpha.
			double alpha = 2.0;
			NSInteger alphaIndex = 0;
			for (NSInteger ip=0; ip<nsetp; ip++) {
				if (zz[ip] <= 0.0) {
					double t = -x[index[ip]]/(zz[ip] - x[index[ip]]);
					if (alpha > t) {
						alpha      = t;
						alphaIndex = ip-1;
					}
				}
			}

			// if all new constrained coeffs are feasible then alpha will still = 2. if so exit from secondary loop to main loop.
			if (alpha == 2.0) break; // out of the while loop

			// otherwise use alpha which will be between 0. and 1. to interpolate between the old x and the new zz.
			for (NSInteger ip=0; ip<nsetp; ip++) x[index[ip]] += alpha*(zz[ip]-x[index[ip]]);

			// modify a and b and the index arrays to move coefficient i from set p to set z.
			BOOL loop = true;
			NSInteger superAlphaIndex = index[alphaIndex+1];
			while (loop) {
				x[superAlphaIndex] = 0.0;

				if (alphaIndex != (nsetp-1)) {
					alphaIndex++;
					for (NSInteger j=(alphaIndex+1); j<nsetp; j++) {
						index[j-1] = index[j];
						// call g1 (a(j-1,index[j]),a(j,index[j]),cc,ss,a(j-1,index[j]))
						double cc, ss;
						[MathSupport G1withA:a[j-1][index[j]] andB:a[j][index[j]] andCterm:&cc andSterm:&ss andSig:&a[j-1][index[j]]];
						a[j][index[j]] = 0.0;
						for (NSInteger l=0; l<n; l++) {
							if (l != index[j]) {
								// apply procedure g2 (cc,ss,a(j-1,l),a(j,l))
								double temp = a[j-1][l];
								a[j-1][l] =  cc*temp + ss*a[j][l];
								a[j][l]   = -ss*temp + cc*a[j][l];
							}
						}

						// apply procedure g2 (cc,ss,b(j-1),b(j))
						double temp = b[j-1];
						b[j-1] =  cc*temp + ss*b[j];
						b[j]   = -ss*temp + cc*b[j];
					}
				}

				npp1 = nsetp-1;
				nsetp--;
				iz1--;
				index[iz1] = superAlphaIndex;

				// see if the remaining coeffs in set p are feasible.  they should be because of the way alpha was determined.
				// if any are infeasible it is due to round-off error.  any that are nonpositive will be set to zero and moved from set p to set z.
				loop = false;
				for (NSInteger jj=0; jj<nsetp; jj++) {
					superAlphaIndex = index[jj];
					if (x[superAlphaIndex] <= 0.0) { loop = true; break; }
				}
			}

			// copy b( ) into zz( ).  then solve again and loop back.
			for (NSInteger i=0; i<m; i++) zz[i] = b[i];
			for (NSInteger l=0; l<nsetp; l++) {
				if (l != 0) for (NSInteger ii=0; ii<=(nsetp-(l+1)); ii++) zz[ii] -= a[ii][index[nsetp-(l+1)+1]]*zz[nsetp-(l+1)+1];
				zz[nsetp-(l+1)] /= a[nsetp-(l+1)][index[nsetp-(l+1)]];
			}

		} // end secondary while loop

		for (NSInteger ip=0; ip<nsetp; ip++) x[index[ip]] = zz[ip];
		//all new coeffs are positive.  loop back to beginning.
	} // end primary while loop

	return finish(false, npp1);

}

+(void)simp1:(double **)a mm:(NSInteger)mm ll:(NSInteger *)ll nll:(NSInteger)nll iabf:(NSInteger)iabf kp:(NSInteger *)kp bmax:(double *)bmax {
	NSInteger k;
	double test;

	*kp = ll[1];
	*bmax = a[mm+1][*kp+1];
	for (k=2; k<=nll; k++) {
		if (iabf == 0)
			test = a[mm+1][ll[k]+1] - (*bmax);
		else
			test = fabs(a[mm+1][ll[k]+1]) - fabs(*bmax);
		if (test > 0.0) {
			*bmax = a[mm+1][ll[k]+1];
			*kp = ll[k];
		}
	}
}

#define EPS 1.0e-6

+(void)simp2:(double **)a n:(NSInteger)n l2:(NSInteger *)l2 nl2:(NSInteger)nl2 ip:(NSInteger *)ip kp:(NSInteger)kp q1:(double *)q1 {
	NSInteger k,ii,i;
	double qp = 0.0, q0 = 0.0, q;

	*ip = 0;
	for (i=1; i<=nl2; i++) {
		if (a[l2[i]+1][kp+1] < -EPS) {
			*q1 = -a[l2[i]+1][1]/a[l2[i]+1][kp+1];
			*ip = l2[i];
			for (i=i+1; i<=nl2; i++) {
				ii = l2[i];
				if (a[ii+1][kp+1] < -EPS) {
					q = -a[ii+1][1]/a[ii+1][kp+1];
					if (q < *q1) {
						*ip = ii;
						*q1 = q;
					} else if (q == *q1) {
						for (k=1; k<=n; k++) {
							qp = -a[*ip+1][k+1]/a[*ip+1][kp+1];
							q0 = -a[ii+1][k+1]/a[ii+1][kp+1];
							if (q0 != qp) break;
						}
						if (q0 < qp) *ip=ii;
					}
				}
			}
		}
	}
}

#undef EPS

+(void)simp3:(double **)a i1:(NSInteger)i1 k1:(NSInteger)k1 ip:(NSInteger)ip kp:(NSInteger)kp {
	NSInteger kk, ii;
	double piv;

	piv = 1.0/a[ip+1][kp+1];
	for (ii=1; ii<=i1+1; ii++)
		if (ii-1 != ip) {
			a[ii][kp+1] *= piv;
			for (kk=1; kk<=k1+1; kk++)
				if (kk-1 != kp)
					a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1];
		}
	for (kk=1; kk<=k1+1; kk++)
		if (kk-1 != kp) a[ip+1][kk] *= -piv;
	a[ip+1][kp+1] = piv;
}

#define EPS 1.0e-6

+(void)simplx:(double **)a m:(NSInteger)m n:(NSInteger)n m1:(NSInteger)m1 m2:(NSInteger)m2 m3:(NSInteger)m3 icase:(NSInteger *)icase izrov:(NSInteger *)izrov iposv:(NSInteger *)iposv {
	NSInteger i, ip, ir, is, k, kh, kp, m12, nl1, nl2;
	NSInteger *l1, *l2, *l3;
	double q1, bmax;

	if (m != (m1+m2+m3)) NSLog(@"Bad input constraint counts in simplx.");
    IntegerVector *l1Wrapper = [[IntegerVector alloc] initWithSize:n+2 andInitialValue:0];
    IntegerVector *l2Wrapper = [[IntegerVector alloc] initWithSize:m+1 andInitialValue:0];
    IntegerVector *l3Wrapper = [[IntegerVector alloc] initWithSize:m+1 andInitialValue:0];
	l1 = [l1Wrapper pointerToInteger];
	l2 = [l2Wrapper pointerToInteger];
	l3 = [l3Wrapper pointerToInteger];
	nl1 = n;
	for (k=1; k<=n; k++) l1[k] = izrov[k] = k;
	nl2 = m;
	for (i=1; i<=m; i++) {
		if (a[i+1][1] < 0.0) NSLog(@"Bad input tableau in simplx.");
		l2[i] = i;
		iposv[i] = n+i;
	}
	for (i=1; i<=m2; i++) l3[i] = 1;
	if (m2+m3) {
		ir = 1;
		for (k=1; k<=(n+1); k++) {
			q1 = 0.0;
			for (i=m1+1; i<=m; i++) q1 += a[i+1][k];
			a[m+2][k] = -q1;
		}
		do {
			// simp1(a, m+1, l1, nl1, 0, &kp, &bmax);
			[MathSupport simp1:a mm:m+1 ll:l1 nll:nl1 iabf:0 kp:&kp bmax:&bmax];
			if (bmax <= EPS && a[m+2][1] < -EPS) {
				*icase = -1;
				return;
			} else if (bmax <= EPS && a[m+2][1] <= EPS) {
				m12 = m1 + m2 + 1;
				if (m12 <= m) {
					for (ip=m12; ip<=m; ip++) {
						if (iposv[ip] == (ip+n)) {
							[MathSupport simp1:a mm:ip ll:l1 nll:nl1 iabf:1 kp:&kp bmax:&bmax];
							// simp1(a, ip, l1, nl1, 1, &kp, &bmax);
							if (bmax > 0.0)
								goto one;
						}
					}
				}
				--m12;
				if (m1+1 <= m12)
					for (i=m1+1; i<=m12; i++)
						if (l3[i-m1] == 1)
							for (k=1; k<=n+1; k++)
								a[i+1][k] = -a[i+1][k];
				break;
			}
			[MathSupport simp2:a n:n l2:l2 nl2:nl2 ip:&ip kp:kp q1:&q1];
			// simp2(a, n, l2, nl2, &ip, kp, &q1);
			if (ip == 0) {
				*icase = -1;
				return;
			}
		one: [MathSupport simp3:a i1:m+1 k1:n ip:ip kp:kp];
			// simp3(a, m+1, n, ip, kp);
			if (iposv[ip] >= (n+m1+m2+1)) {
				for (k=1; k<=nl1; k++)
					if (l1[k] == kp) break;
				--nl1;
				for (is=k; is<=nl1; is++) l1[is] = l1[is+1];
				++a[m+2][kp+1];
				for (i=1; i<=m+2; i++) a[i][kp+1] = -a[i][kp+1];
			} else {
				if (iposv[ip] >= (n+m1+1)) {
					kh = iposv[ip] - m1 - n;
					if (l3[kh]) {
						l3[kh] = 0;
						++a[m+2][kp+1];
						for (i=1; i<=m+2; i++)
							a[i][kp+1] = -a[i][kp+1];
					}
				}
			}
			is = izrov[kp];
			izrov[kp] = iposv[ip];
			iposv[ip] = is;
		} while (ir);
	}
	for (;;) {
		[MathSupport simp1:a mm:0 ll:l1 nll:nl1 iabf:0 kp:&kp bmax:&bmax];
		// simp1(a, 0, l1, nl1, 0, &kp, &bmax);
		if (bmax <= 0.0) {
			*icase = 0;
			return;
		}
		[MathSupport simp2:a n:n l2:l2 nl2:nl2 ip:&ip kp:kp q1:&q1];
		// simp2(a, n, l2, nl2, &ip, kp, &q1);
		if (ip == 0) {
			*icase = 1;
			return;
		}
		[MathSupport simp3:a i1:m k1:n ip:ip kp:kp];
		// simp3(a, m, n, ip, kp);
		is = izrov[kp];
		izrov[kp] = iposv[ip];
		iposv[ip] = is;
	}
}
#undef EPS
#undef NRANSI

+(BOOL)choldc:(double **)a dimension:(NSUInteger)n {
	for (NSUInteger i=0; i<n; i++) {
		double p = 0.0;
		for (NSUInteger j=i; j<n; j++) {
			double sum = a[i][j];
			if (i > 0) for (NSInteger k=i-1; k>=0; k--) sum -= a[i][k]*a[j][k];
			if (i == j) {
				if (sum <= 0.0) return NO;
				p = sqrt(sum);
			} else a[j][i] = sum/p;
		}
	}
	return YES;
}

+(void)polint:(double *)xa ya:(double *)ya n:(NSInteger)n x:(double)x y:(double *)y dy:(double *)dy {
	NSInteger ns=1;
	double dif = fabs(x - xa[0]);
    DoubleVector *cWrapper = [[DoubleVector alloc] initWithSize:n andInitialValue:0.0];
    DoubleVector *dWrapper = [[DoubleVector alloc] initWithSize:n andInitialValue:0.0];
    double *c = [cWrapper pointerToDouble];
	double *d = [dWrapper pointerToDouble];
	for (NSInteger i=1; i<=n; i++) {
		double dift = fabs(x - xa[i-1]);
		if (dift < dif) {
			ns = i;
			dif = dift;
		}
		c[i-1] = ya[i-1];
		d[i-1] = ya[i-1];
	}
	*y = ya[ns-- - 1];
	for (NSInteger m=1; m<n; m++) {
		for (NSInteger i=1; i<=(n-m); i++) {
			double ho = xa[i-1] - x;
			double hp = xa[i+m-1] - x;
			double w = c[i] - d[i-1];
			double den = ho - hp;
			if (den == 0.0) NSLog(@"Error in routine POLINT");
			den = w/den;
			d[i-1] = hp*den;
			c[i-1] = ho*den;
		}
		*y += (*dy = (2*ns < (n-m) ? c[ns] : d[ns-- -1]) );
	}
}

+(void)polin2:(double *)x1a x2a:(double *)x2a ya:(double **)ya m:(NSInteger)m n:(NSInteger)n x1:(double)x1 x2:(double)x2 y:(double *)y dy:(double *)dy {
    DoubleVector *ymtmpWrapper = [[DoubleVector alloc] initWithSize:m andInitialValue:0.0];
    double *ymtmp = [ymtmpWrapper pointerToDouble];
	for (NSInteger j=0; j<m; j++) [MathSupport polint:x2a ya:ya[j] n:n x:x2 y:&ymtmp[j] dy:dy];
	[MathSupport polint:x1a ya:ymtmp n:m x:x1 y:y dy:dy];
}

/**
 Algorithm LDP (Least distance programming)
 Given an m by n matrix G, and an m-vector H, computes an n-vector X that solves the least squares problem minimize ||x||, subject to G X >= H

 @param g
 On entry g[m][n+1] contains the m by n matrix G.  Contents are unmodified on exit.

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
 On entry, a vector of length 2*n + 2*m + 2 of working storage.

 @param e
 On entry a n+1 by m matrix of working storage.

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
       andIntegerWorkingVector:(NSInteger *)index {

    if (n <= 0) return NO;
    for (NSUInteger i=0; i<n; i++) x[i] = 0.0;
    *xnorm = 0.0;
    if (m <= 0) return YES;

    for (NSUInteger j=0; j<m; j++) {
        for (NSUInteger i=0; i<n; i++) e[i][j] = g[j][i];
        e[n][j] = h[j];
    }
    for (NSUInteger i=0; i<n; i++) w[i] = 0.0;
    w[n] = 1.0;

    double rNorm = 0.0;
    if (![MathSupport nnlsWithConstraintMatrix:e
                               andRowDimension:n+1
                            andColumnDimension:m
                           andConstraintVector:w
                             andSolutionVector:&w[2*n+2]
                               andResidualNorm:&rNorm
                         andDualSolutionVector:&w[2*n+m+2]
                        andDoubleWorkingVector:&w[n+1]
                      andUIntegerWorkingVector:index]) return NO;
    if (rNorm <= 0.0) return NO;

    double factor = 1.0;
    for (NSUInteger i=0; i<m; i++) factor -= h[i]*w[2*n+2+i];
    if (factor == 0.0) return NO;

    for (NSUInteger j=0; j<n; j++) {
        for (NSUInteger i=0; i<m; i++) x[j] += g[i][j]*w[2*n+2+i];
        x[j] /= factor;
    }

    for (NSUInteger j=0; j<n; j++) *xnorm += x[j]*x[j];
    *xnorm = sqrt(*xnorm);
    return YES;
}



@end
