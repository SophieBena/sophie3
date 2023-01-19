//
//  HornblendeBerman.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/10/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import "HornblendeBerman.h"
#import "BermanProperties.h"
#import "BermanStoichiometricPhases.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"

@implementation HornblendeBerman

static NSArray *endmembers;

#define NR         2    // Independent composition variables
#define NA         3    // Endmember compositions
#define NATOMS  42.0    // Number of atoms in the formula unit

#pragma mark -
#pragma mark class methods

+(void)initialize {
	if (self == [HornblendeBerman class]) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"Initialize(HornblendeBerman) - entry ...");
		NSMutableArray *mutableEndmembers = [NSMutableArray arrayWithCapacity:NA];

		BermanProperties *pargasite = [[BermanProperties alloc] initWithH:-3016624.0*4.184
																		S:160.0*4.184
																	   k0:1267.25
																	   k1:-66.5434e2
																	   k2:-303.787e5
																	   k3:391.353e7
																	   v0:27.35
																	   v1:-1.392e-6
																	   v2:3.481e-12
																	   v3:24.374e-6
																	   v4:98.338e-10];
		[pargasite setPhaseFormula:@"NaCa2Mg4AlAl2Si6O22(OH)2"];
		[pargasite setPhaseName:@"pargasite"];
		[mutableEndmembers addObject:pargasite];
		if (debug) NSLog(@"... allocated pargasite ...");

#define DHGrnCum (-9623300.0-(-12067517.0))

		BermanProperties *ferropargasite = [[BermanProperties alloc] initWithH:-3016624.0*4.184+4.0*DHGrnCum/7.0
																			 S:185.5*4.184
																			k0:1342.61
																			k1:-83.4862e2
																			k2:-247.604e5
																			k3:348.507e7
																			v0:27.989
																			v1:-1.392e-6
																			v2:3.481e-12
																			v3:24.374e-6
																			v4:98.338e-10];
		[ferropargasite setPhaseFormula:@"NaCa2Fe4AlAl2Si6O22(OH)2"];
		[ferropargasite setPhaseName:@"ferropargasite"];
		[mutableEndmembers addObject:ferropargasite];
		if (debug) NSLog(@"... allocated ferropargasite ...");

#undef DHGrnCum
#define DHBfABf ((-2836709.0+7932.05)-(-3275265.0-8565.18))

		BermanProperties *magnesiohastingsite = [[BermanProperties alloc] initWithH:-3016624.0*4.184+DHBfABf
																				  S:163.8*4.184
																				 k0:1273.66
																				 k1:-67.1606e2
																				 k2:-280.331e5
																				 k3:350.697e7
																				 v0:27.38
																				 v1:-1.392e-6
																				 v2:3.481e-12
																				 v3:24.374e-6
																				 v4:98.338e-10];
		[magnesiohastingsite setPhaseFormula:@"NaCa2Mg4FeAl2Si6O22(OH)2"];
		[magnesiohastingsite setPhaseName:@"magnesiohastingsite"];
		[mutableEndmembers addObject:magnesiohastingsite];
		if (debug) NSLog(@"... allocated magnesiohastingsite ...");

#undef DHBfABf

		endmembers = [NSArray arrayWithArray:mutableEndmembers];
	}
}

#pragma mark -
#pragma mark instance methods

-(id)init {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(HornblendeBerman) ...");

		[self setPhaseName:@"Hornblende"];
		computeMixingQuantities = NO;

		if (debug) NSLog(@"... exiting.");
	}
	return self;
}

/**
 "C" code from MELTS
 */

#define FIRST       00000001 /* octal for binary 00000000000000000001 */
#define SECOND      00000002 /* octal for binary 00000000000000000010 */
#define THIRD       00000004 /* octal for binary 00000000000000000100 */
#define FOURTH      00000010 /* octal for binary 00000000000000001000 */
#define FIFTH       00000020 /* octal for binary 00000000000000010000 */
#define SIXTH       00000040 /* octal for binary 00000000000000100000 */
#define SEVENTH     00000100 /* octal for binary 00000000000001000000 */
#define EIGHTH      00000200 /* octal for binary 00000000000010000000 */
#define NINTH       00000400 /* octal for binary 00000000000100000000 */
#define TENTH       00001000 /* octal for binary 00000000001000000000 */
#define ELEVENTH    00002000 /* octal for binary 00000000010000000000 */
#define TWELFTH     00004000 /* octal for binary 00000000100000000000 */
#define THIRTEENTH  00010000 /* octal for binary 00000001000000000000 */
#define FOURTEENTH  00020000 /* octal for binary 00000010000000000000 */
#define FIFTEENTH   00040000 /* octal for binary 00000100000000000000 */
#define SIXTEENTH   00100000 /* octal for binary 00001000000000000000 */
#define SEVENTEENTH 00200000 /* octal for binary 00010000000000000000 */
#define EIGHTEENTH  00400000 /* octal for binary 00100000000000000000 */
#define NINETEENTH  01000000 /* octal for binary 01000000000000000000 */
#define TWENTIETH   02000000 /* octal for binary 10000000000000000000 */

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * hornblende (pargasite-hastingite) reciprocal solution
 */

#define WFEMG   4.0*1.68*4.184*1000.0  /* 4*W12 pyroxenes joules */
#define WFEAL  16.78*1000.0            /* W34   pyroxenes joules */
#define DGR     0.0                    /* joules */

/* Change nothing below this line */

#define G0   0.0
#define GX   (WFEMG)
#define GY   (WFEAL)
#define GXX -(WFEMG)
#define GYY -(WFEAL)
#define GXY  (DGR)

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conHrn defines the conversion from m[i], to r[j]
 */
#define FR0(i)     (i == 1) ? 1.0 - r[0] : - r[0]
#define FR1(i)     (i == 2) ? 1.0 - r[1] : - r[1]
#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */
#define R          8.3143
#define S          -R*(4.0*xMgM12*log(xMgM12) + 4.0*xFe2M12*log(xFe2M12) \
+ xFe3M3*log(xFe3M3) + xAlM3*log(xAlM3))
#define H          (G0)+(GX)*r[0]+(GY)*r[1]+(GXX)*r[0]*r[0] \
+(GYY)*r[1]*r[1]+(GXY)*r[0]*r[1]
#define V          0.0
#define G          (H) - t*(S) + (p-1.0)*(V)

#define DGDR0      R*t*(4.0*log(xFe2M12) - 4.0*log(xMgM12)) \
+ (GX) + 2.0*(GXX)*r[0] + (GXY)*r[1]
#define DGDR1      R*t*(log(xFe3M3) - log(xAlM3)) \
+ (GY) + 2.0*(GYY)*r[1] + (GXY)*r[0]
#define DGDT       -(S)
#define DGDP       (V)

#define D2GDR0R0   R*t*(4.0/xFe2M12 + 4.0/xMgM12) + 2.0*(GXX)
#define D2GDR0R1   (GXY)
#define D2GDR1R1   R*t*(1.0/xFe3M3 + 1.0/xAlM3) + 2.0*(GYY)
#define D2GDR0DT   R*(4.0*log(xFe2M12) - 4.0*log(xMgM12))
#define D2GDR1DT   R*(log(xFe3M3) - log(xAlM3))
#define D2GDR0DP   0.0
#define D2GDR1DP   0.0
#define D2GDT2     0.0
#define D2GDTDP    0.0
#define D2GDP2     0.0

#define D3GDR0R0R0 R*t*(-4.0/(xFe2M12*xFe2M12) + 4.0/(xMgM12*xMgM12))
#define D3GDR0R0R1 0.0
#define D3GDR0R0DT R*(4.0/xFe2M12 + 4.0/xMgM12)
#define D3GDR0R0DP 0.0
#define D3GDR0R1R1 0.0
#define D3GDR0R1DT 0.0
#define D3GDR0R1DP 0.0
#define D3GDR1R1R1 R*t*(-1.0/(xFe3M3*xFe3M3) + 1.0/(xAlM3*xAlM3))
#define D3GDR1R1DT R*(1.0/xFe3M3 + 1.0/xAlM3)
#define D3GDR1R1DP 0.0
#define D3GDT3     0.0
#define D3GDT2DP   0.0
#define D3GDTDP2   0.0
#define D3GDP3     0.0

#define D3GDR0DT2  0.0
#define D3GDR1DT2  0.0
#define D3GDR0DTDP 0.0
#define D3GDR1DTDP 0.0
#define D3GDR0DP2  0.0
#define D3GDR1DP2  0.0

/*
 *=============================================================================
 * Public functions:
 *    mask  -  bitwise mask for selecting output
 *    t     -  Temperature (K)
 *    p     -  Pressure (bars)
 *    *x    -  (pointer to x[]) Array of independent compositional variables
 */

-(BOOL)test:(int)mask
		  t:(double)t
		  p:(double)p
		 na:(int)na           // Expected number of endmember components
		 nr:(int)nr           // Expected number of independent compositional variables
	  names:(char **)names    // array of strings of names of endmember components
   formulas:(char **)formulas // array of strings of formulas of endmember components
		  r:(double [NR])r    // array of indepependent compos variables, check bounds
		  m:(double [NA])m    // array of moles of endmember components, check bounds
{
	const char *NAMES[NA]    = { "pargasite",  "ferropargasite",       "magnesiohastingsite" };
	const char *FORMULAS[NA] = { "NaCa2Mg4AlAl2Si6O22(OH)2", "NaCa2Fe4AlAl2Si6O22(OH)2",
		"NaCa2Mg4FeAl2Si6O22(OH)2"  };
	BOOL result = YES;
	int i;

	if (mask & FIRST) {
		result = result && (na == NA);
	}
	if (mask & SECOND) {
		result = result && (nr == NR);
	}
	if (mask & THIRD) {
		for (i=0; i<NA; i++) {
			result = result && (strcmp(names[i],NAMES[i]) == 0);
		}
	}
	if (mask & FOURTH) {
		for (i=0; i<NA; i++) {
			result = result && (strcmp(formulas[i],FORMULAS[i]) == 0);
		}
	}
	/* Check bounds on the independent compositional variables */
	if (mask & FIFTH) {
		for (i=0; i<NR; i++) {
			result = result && (r[i] >= 0.0) && (r[i] <= 1.0);
		}
	}
	/* Check bounds on moles of endmember components */
	if (mask & SIXTH) {
		for (i=1; i<NA; i++) result = result && (m[i] >= 0.0);
		result = result && (m[0] >= -m[2]);
	}

	return result;
}

-(void)convert:(int)inpMask
	   outMask:(int)outMask
			 t:(double)t
			 p:(double)p
			 e:(double [107])e              // comp of olivine in moles of elements
			 m:(double [NA])m               // comp of olivine in moles of endmember components
			 r:(double [NR])r               // comp of olivine in terms of the independent comp var
			 x:(double [NA])x               // comp of olivine in mole fractions of endmember comp
			dm:(double [NR][NA])dm          // Jacobian matrix: dm[i][j] = dr[i]/dm[j]
		   d2m:(double [NR][NA][NA])d2m     // vector of matrices: d2m[i][j][k] = d2r[i]/dm[j]dm[k]
			dr:(double [NA][NR])dr          // Jacobian matrix: dr[i][j] = dx[i]/dr[j]
		   d3m:(double [NR][NA][NA][NA])d3m // 3rd deriv matrix: d3m[i][j][k][l]=d3r[i]/dm[j]dm[k]dm[l]
{
	/*---------------------------------------------------------------------------
	 Not all combinations of inpMask and outMask are feasible. Valid
	 combinations are:

	 inpMask          outMask
	 (1)  FIRST            SECOND
	 (2)  SECOND           THIRD  | FOURTH  | FIFTH | SIXTH | EIGHTH
	 (3)  THIRD            FOURTH | SEVENTH

	 (1) converts a vector of moles of elements into a vector of moles of
	 endmember leucite components.
	 (2) calculates from a vector of moles of endmember components, one or
	 all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
	 (3) calculates from a vector of independent compositional variables
	 mole fractions of endmember components and/or the Jacobian matrix
	 dx[]/dr[]

	 In this routine it is assumed that the elements are in the order of atomic
	 numbers and that the order of leucite components has been verified as:
	 m[0] = pargasite           (NaCa2Mg4AlAl2Si6O22(OH)2),
	 m[1] = ferropargasite      (NaCa2Fe4AlAl2Si6O22(OH)2),
	 m[2] = magnesiohastingsite (NaCa2Mg4FeAl2Si6O22(OH)2).

	 ----------------------------------------------------------------------------*/

	int i, j, k;

	if (inpMask == FIRST && outMask == SECOND) {
		double sumchg, fe2, fe3;
		static const int Hy =  1;
		static const int O  =  8;
		static const int Na = 11;
		static const int Mg = 12;
		static const int Al = 13;
		static const int Si = 14;
		static const int Ca = 20;
		static const int Fe = 26;

		sumchg  = e[Hy] + e[Na] + 2.0*e[Mg] + 3.0*e[Al] + 4.0*e[Si] + 2.0*e[Ca];
		fe3 = 2.0*e[O] - sumchg - 2.0*e[Fe];
		fe2 = e[Fe] - fe3;

		if (fe3 < 0.0) { fe3 = 0.0; fe2 = e[Fe]; }
		if (fe2 < 0.0) { fe2 = 0.0; fe3 = e[Fe]; }
		if (fe3 > e[Fe]) { fe3 = e[Fe]; fe2 = 0.0; }
		if (fe2 > e[Fe]) { fe2 = e[Fe]; fe3 = 0.0; }

		/* Projection into the  ternary */
		m[0] = (e[Mg]-fe3*4.0)/4.0; /* moles of NaCa2Mg4AlAl2Si6O22(OH)2  */
		m[1] = fe2/4.0;             /* moles of NaCa2Fe4AlAl2Si6O22(OH)2  */
		m[2] = fe3;                 /* moles of NaCa2Mg4FeAl2Si6O22(OH)2  */

	} else if (inpMask == SECOND) {
		double sum;

		if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
			NSLog(@"Illegal call to conFld with inpMask = %o and outMask = %o\n", inpMask, outMask);

		for (i=0, sum=0.0; i<NA; i++) sum += m[i];

		if (outMask & THIRD) {
			for (i=0; i<NR; i++) r[i] = (sum != 0.0) ? m[i+1]/sum : 0.0;
		}

		if (outMask & FOURTH) {
			/* Converts a vector of moles of end-member components (m) into a vector
			 of mole fractions of endmember components                            */
			for (i=0; i<NA; i++) x[i] = (sum != 0.0) ? m[i]/sum : 0.0;
		}

		if (outMask & FIFTH) {
			/* Calculates the matrix dr[i]/dm[j] using m[] as input                 */
			if (sum == 0.0) {
				for (i=0; i<NR; i++) { for (j=0; j<NA; j++) dm[i][j] = 0.0; }
			} else {
				for (i=0; i<NR; i++) {
					for (j=0; j<NA; j++)
						dm[i][j] = (i+1 == j) ? (1.0-m[i+1]/sum)/sum : - m[i+1]/SQUARE(sum);
				}
			}

		}

		if (outMask & SIXTH) {
			/* Calculates the matrix d2r[i]/dm[j]dm[k] using m[] as input           */

			if (sum == 0.0) {
				for (i=0; i<NR; i++) {
					for (j=0; j<NA; j++)  {
						for (k=0; k<NA; k++) d2m[i][j][k] = 0.0;
					}
				}
			} else {
				for (i=0; i<NR; i++) {
					for (j=0; j<NA; j++)  {
						for (k=0; k<NA; k++) {
							d2m[i][j][k]  = 2.0*m[i+1]/CUBE(sum);
							d2m[i][j][k] -= (i+1 == j) ? 1.0/SQUARE(sum) : 0.0;
							d2m[i][j][k] -= (i+1 == k) ? 1.0/SQUARE(sum) : 0.0;
						}
					}
				}
			}

		}

		if (outMask & EIGHTH) {
			/* calculates the matrix d3r[i]/dm[j]dm[k]dm[l] using m[] as input */
			int l;

			if (sum == 0.0) {
				for (i=0; i<NR; i++) {
					for (j=0; j<NA; j++)  {
						for (k=0; k<NA; k++)  {
							for (l=0; l<NA; l++) d3m[i][j][k][l] = 0.0;
						}
					}
				}
			} else {
				for (i=0; i<NR; i++) {
					for (j=0; j<NA; j++)  {
						for (k=0; k<NA; k++)  {
							for (l=0; l<NA; l++)  {
								d3m[i][j][k][l]  = -6.0*m[i+1]/QUARTIC(sum);
								d3m[i][j][k][l] += (i+1 == j) ? 2.0/CUBE(sum) : 0.0;
								d3m[i][j][k][l] += (i+1 == k) ? 2.0/CUBE(sum) : 0.0;
								d3m[i][j][k][l] += (i+1 == l) ? 2.0/CUBE(sum) : 0.0;
							}
						}
					}
				}
			}

		}

	} else if (inpMask == THIRD) {

		if (outMask & ~(FOURTH | SEVENTH))
			NSLog(@"Illegal call to conFld with inpMask = %o and outMask = %o\n", inpMask, outMask);

		if (outMask & FOURTH) {
			/* Converts a vector of independent compositional variables (r)
			 into a vector of mole fractions of endmember components (x).         */

			for (i=0, x[0]=1.0; i<NR; i++) { x[i+1] = r[i]; x[0] -= r[i]; }
		}

		if (outMask & SEVENTH) {
			/* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
			for (i=0; i<NR; i++) for (j=0; j<NR; j++) dr[i+1][j] = (i == j) ? 1.0 : 0.0;
			for (j=0; j<NR; j++) dr[0][j] = -1.0;
		}

	} else  {
		NSLog(@"Illegal call to conFld with inpMask = %o and outMask = %o\n", inpMask, outMask);
	}

}

-(NSString *)displayFormula:(double)t
						  p:(double)p
						  r:(double [NA])r
{
	double totMg, totFe2, totAl, totFe3;

	totMg  = (1.0-r[0])*4.0;
	totFe2 = r[0]*4.0;
	totAl  = 1.0-r[1];
	totFe3 = r[1];

	return [NSString stringWithFormat:@"NaCa2Mg%4.2fFe2+%4.2fAl%4.2fFe3+%4.2fAl2Si6O22(OH)2", totMg, totFe2, totAl, totFe3];
}

-(void)activity:(int)mask
			  t:(double)t
			  p:(double)p
			  r:(double [NA])r
			  a:(double [NA])a      // (pointer to a[]) activities              BINARY MASK: 0001
			 mu:(double [NA])mu     // (pointer to mu[]) chemical potentials    BINARY MASK: 0010
			 dx:(double [NA][NR])dx // (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100
{
	double xMgM12  = 1.0-r[0];
	double xFe2M12 = r[0];
	double xAlM3   = 1.0-r[1];
	double xFe3M3  = r[1];
	if (xMgM12  <= 0.0) xMgM12  = DBL_EPSILON;
	if (xFe2M12 <= 0.0) xFe2M12 = DBL_EPSILON;
	if (xAlM3   <= 0.0) xAlM3   = DBL_EPSILON;
	if (xFe3M3  <= 0.0) xFe3M3  = DBL_EPSILON;

	double g, dgdr[NR], fr[NA][NR];
	int i, j;

	for(i=0; i<NA; i++) fr[i][0] = FR0(i);
	for(i=0; i<NA; i++) fr[i][1] = FR1(i);

	g       = G;
	dgdr[0] = DGDR0;
	dgdr[1] = DGDR1;

	if (mask & FIRST) {
		for(i=0; i<NA; i++) {
			for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
			a[i] = exp(a[i]/(R*t));
		}
	}

	if (mask & SECOND) {
		for(i=0; i<NA; i++) {
			for (mu[i]=g, j=0; j<NR; j++) mu[i] += fr[i][j]*dgdr[j];
		}
	}

	if (mask & THIRD) {
		double d2gdr2[NR][NR], dfrdr[NA][NR], sum;
		int k;

		d2gdr2[0][0] = D2GDR0R0;
		d2gdr2[0][1] = D2GDR0R1; d2gdr2[1][0] = d2gdr2[0][1];
		d2gdr2[1][1] = D2GDR1R1;

		for(i=0; i<NA; i++) dfrdr[i][0] = DFR0DR0(i);
		for(i=0; i<NA; i++) dfrdr[i][1] = DFR1DR1(i);

		for (i=0; i<NA; i++) {
			for (k=0; k<NR; k++) {
				for (dx[i][k]=g, j=0; j<NR; j++) dx[i][k] += fr[i][j]*dgdr[j];
				dx[i][k] = exp(dx[i][k]/(R*t));
				sum = (1.0+dfrdr[i][k])*dgdr[k];
				for (j=0; j<NR; j++) sum += fr[i][j]*d2gdr2[j][k];
				dx[i][k] *= sum/(R*t);
			}
		}
	}

}

-(void)gmix:(int)mask
		  t:(double)t
		  p:(double)p
		  r:(double [NR])r
	   gmix:(double *)gmix           // Gibbs energy of mixing             BINARY MASK: 0001
		 dx:(double [NR])dx          // (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010
		dx2:(double [NR][NR])dx2     // (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100
		dx3:(double [NR][NR][NR])dx3 // (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000
{
	double xMgM12  = 1.0-r[0];
	double xFe2M12 = r[0];
	double xAlM3   = 1.0-r[1];
	double xFe3M3  = r[1];
	if (xMgM12  <= 0.0) xMgM12  = DBL_EPSILON;
	if (xFe2M12 <= 0.0) xFe2M12 = DBL_EPSILON;
	if (xAlM3   <= 0.0) xAlM3   = DBL_EPSILON;
	if (xFe3M3  <= 0.0) xFe3M3  = DBL_EPSILON;

	if (mask & FIRST) {
		*gmix = G;
	}

	if(mask & SECOND) {
		dx[0] = DGDR0;
		dx[1] = DGDR1;
	}

	if(mask & THIRD) {
		double d2gdr2[NR][NR];
		int i, j;

		d2gdr2[0][0] = D2GDR0R0;
		d2gdr2[0][1] = D2GDR0R1; d2gdr2[1][0] = d2gdr2[0][1];
		d2gdr2[1][1] = D2GDR1R1;

		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) dx2[i][j] = d2gdr2[i][j];
		}
	}
	if (mask & FOURTH) {
		double d3gdr3[NR][NR][NR];
		int i, j, k;

		d3gdr3[0][0][0] = D3GDR0R0R0;
		d3gdr3[0][0][1] = D3GDR0R0R1; d3gdr3[0][1][0] = d3gdr3[0][0][1]; d3gdr3[1][0][0] = d3gdr3[0][0][1];
		d3gdr3[0][1][1] = D3GDR0R1R1; d3gdr3[1][0][1] = d3gdr3[0][1][1]; d3gdr3[1][1][0] = d3gdr3[0][1][1];
		d3gdr3[1][1][1] = D3GDR1R1R1;

		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) {
				for (k=0; k<NR; k++) dx3[i][j][k] = d3gdr3[i][j][k];
			}
		}
	}

}

-(void)hmix:(int)mask
		  t:(double)t
		  p:(double)p
		  r:(double [NR])r
	   hmix:(double *)hmix // Enthalpy of mixing BINARY MASK: 1
/*  This function calculates enthalpy of mixing corrected to 1 bar. */
{
	double xMgM12  = 1.0-r[0];
	double xFe2M12 = r[0];
	double xAlM3   = 1.0-r[1];
	double xFe3M3  = r[1];
	if (xMgM12  <= 0.0) xMgM12  = DBL_EPSILON;
	if (xFe2M12 <= 0.0) xFe2M12 = DBL_EPSILON;
	if (xAlM3   <= 0.0) xAlM3   = DBL_EPSILON;
	if (xFe3M3  <= 0.0) xFe3M3  = DBL_EPSILON;

	*hmix = (G) + t*(S);
}

-(void)smix:(int)mask
		  t:(double)t
		  p:(double)p
		  r:(double [NR])r
	   smix:(double *)smix        // Entropy of mixing                  BINARY MASK: 001
		 dx:(double [NR])dx       // (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010
		dx2:(double [NR][NR])dx2  // (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100
{
	double xMgM12  = 1.0-r[0];
	double xFe2M12 = r[0];
	double xAlM3   = 1.0-r[1];
	double xFe3M3  = r[1];
	if (xMgM12  <= 0.0) xMgM12  = DBL_EPSILON;
	if (xFe2M12 <= 0.0) xFe2M12 = DBL_EPSILON;
	if (xAlM3   <= 0.0) xAlM3   = DBL_EPSILON;
	if (xFe3M3  <= 0.0) xFe3M3  = DBL_EPSILON;

	if (mask & FIRST) {
		*smix = S;
	}

	if(mask & SECOND) {
		double d2gdrdt[NR];
		int i;

		d2gdrdt[0] = D2GDR0DT;
		d2gdrdt[1] = D2GDR1DT;

		for (i=0; i<NR; i++) dx[i] = - d2gdrdt[i];
	}

	if(mask & THIRD) {
		double d3gdr2dt[NR][NR];
		int i, j;

		d3gdr2dt[0][0] = D3GDR0R0DT;
		d3gdr2dt[0][1] = D3GDR0R1DT; d3gdr2dt[1][0] = d3gdr2dt[0][1];
		d3gdr2dt[1][1] = D3GDR1R1DT;

		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) dx2[i][j] = - d3gdr2dt[i][j];
		}
	}

}

-(void)cpmix:(int)mask
		   t:(double)t
		   p:(double)p
		   r:(double [NR])r
	   cpmix:(double *)cpmix // Heat capacity of mixing         BINARY MASK: 001
		  dt:(double *)dt    // d(cp)/d(t)                      BINARY MASK: 010
		  dx:(double [NR])dx // d(cp)/d(x[])d(t)                BINARY MASK: 100
{
	double d2gdt2 = D2GDT2;

	if (mask & FIRST) {
		*cpmix = - t*d2gdt2;
	}

	if(mask & SECOND) {
		double d3gdt3   = D3GDT3;

		*dt = -t*d3gdt3 - d2gdt2;
	}

	if(mask & THIRD) {
		double d3gdrdt2[NR];
		int i;

		d3gdrdt2[0] = D3GDR0DT2;
		d3gdrdt2[1] = D3GDR1DT2;

		for (i=0; i<NR; i++) dx[i] = -t*d3gdrdt2[i];
	}

}

-(void)vmix:(int)mask
		  t:(double)t
		  p:(double)p
		  r:(double [NR])r
	   vmix:(double *)vmix       // Volume of mixing                BINARY MASK: 0000000001
		 dx:(double [NR])dx      // pointer to dx[]) d(v)/d(x[])    BINARY MASK: 0000000010
		dx2:(double [NR][NR])dx2 // pointer to dx2[][]) d(v)/d(x[])2BINARY MASK: 0000000100
		 dt:(double *)dt         // d(v)/d(t)                       BINARY MASK: 0000001000
		 dp:(double *)dp         // d(v)/d(p)                       BINARY MASK: 0000010000
		dt2:(double *)dt2        // d2(v)/d(t)2                     BINARY MASK: 0000100000
	   dtdp:(double *)dtdp       // d2(v)/d(t)d(p)                  BINARY MASK: 0001000000
		dp2:(double *)dp2        // d2(v)/d(p)2                     BINARY MASK: 0010000000
	   dxdt:(double [NR])dxdt    // d2(v)/d(x[])d(t)                BINARY MASK: 0100000000
	   dxdp:(double [NR])dxdp    // d2(v)/d(x[])d(p)                BINARY MASK: 1000000000
{
	if (mask & FIRST) {
		*vmix = DGDP;
	}

	if(mask & SECOND) {
		double d2gdrdp[NR];
		int i;

		d2gdrdp[0] = D2GDR0DP;
		d2gdrdp[1] = D2GDR1DP;

		for (i=0; i<NR; i++) dx[i] = d2gdrdp[i];
	}

	if(mask & THIRD) {
		double d3gdr2dp[NR][NR];
		int i, j;

		d3gdr2dp[0][0] = D3GDR0R0DP;
		d3gdr2dp[0][1] = D3GDR0R1DP; d3gdr2dp[1][0] = d3gdr2dp[0][1];
		d3gdr2dp[1][1] = D3GDR1R1DP;

		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) dx2[i][j] = d3gdr2dp[i][j];
		}
	}

	if(mask & FOURTH) {
		*dt = D2GDTDP;
	}

	if(mask & FIFTH) {
		*dp = D2GDP2;
	}

	if(mask & SIXTH) {
		*dt2 = D3GDT2DP;
	}

	if(mask & SEVENTH) {
		*dtdp = D3GDTDP2;
	}

	if(mask & EIGHTH) {
		*dp2 = D3GDP3;
	}

	if(mask & NINTH) {
		double d3gdrdtdp[NR];
		int i;

		d3gdrdtdp[0] = D3GDR0DTDP;
		d3gdrdtdp[1] = D3GDR1DTDP;

		for (i=0; i<NR; i++) dxdt[i] = d3gdrdtdp[i];
	}

	if(mask & TENTH) {
		double d3gdrdp2[NR];
		int i;

		d3gdrdp2[0] = D3GDR0DP2;
		d3gdrdp2[1] = D3GDR1DP2;

		for (i=0; i<NR; i++) dxdp[i] = d3gdrdp2[i];
	}

}

#define NAS 4

-(NSUInteger)numberOfSolutionSpecies {
	return NAS;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
	switch (index) {
		case 0:
			return @"pargasite";           // NaCa2Mg4AlAl2Si6O22(OH)2
			break;
		case 1:
			return @"ferropargasite";      // NaCa2Fe4AlAl2Si6O22(OH)2
			break;
		case 2:
			return @"magnesiohastingsite"; // NaCa2Mg4FeAl2Si6O22(OH)2
			break;
		case 3:
			return @"ferrohastingsite";    // NaCa2Fe4FeAl2Si6O22(OH)2
			break;
		default:
			return @"";
			break;
	}
}

-(DoubleVector *)convertMolesOfSpeciesToMolesOfComponents:(double *)mSpecies {
	DoubleVector *mComponentWrapper =[[DoubleVector alloc] initWithSize:NA];
	double *mComponents = [mComponentWrapper pointerToDouble];
	mComponents[0] = mSpecies[0] - mSpecies[3];
	mComponents[1] = mSpecies[1] + mSpecies[3];
	mComponents[2] = mSpecies[2] + mSpecies[3];
	return mComponentWrapper;
}

-(DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	DoubleVector *muSpeciesWrapper = [[DoubleVector alloc] initWithSize:NAS];
	double *muSpecies = [muSpeciesWrapper pointerToDouble];
	double *muComponents = [[self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p] pointerToDouble];
	for (NSUInteger i=0; i<NA; i++) muSpecies[i] = muComponents[i];
	// NaCa2Fe4FeAl2Si6O22(OH)2
	muSpecies[3] = muSpecies[2] + muSpecies[1] - muSpecies[0];
	return muSpeciesWrapper;
}

-(DoubleVector *)elementalCompositionOfSpeciesAtIndex:(NSUInteger)index {
	DoubleVector *speciesElementArrayWrapper = nil;

	if (index < NA) speciesElementArrayWrapper = [[endmembers objectAtIndex:index] formulaAsElementArray];
	else if (index == 3) {
		speciesElementArrayWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
		double *speciesElementArray = [speciesElementArrayWrapper pointerToDouble];
        for (NSUInteger i=1; i<107; i++) speciesElementArray[i] = 0.0;
		double *component = [[[endmembers objectAtIndex:2] formulaAsElementArray] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] += component[i];
		component = [[[endmembers objectAtIndex:1] formulaAsElementArray] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] += component[i];
		component = [[[endmembers objectAtIndex:0] formulaAsElementArray] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] -= component[i];
	}
	return speciesElementArrayWrapper;
}

-(void)correctActivityCoefficients:(double [NAS])gamma forComposition:(double [NAS])x { }

// --> SolutionPhaseProtocol public function
-(NSArray *)affinityAndCompositionFromLiquidChemicalPotentialSum:(double *)chemicalPotentials andT:(double)t andP:(double)p {
	NSMutableArray *results = [NSMutableArray arrayWithCapacity:NA+1];
	double mu0[NAS], deltaMu[NAS], xNz[NAS], x[NAS], gamma[NAS], xLast[NAS], gammaLast[NAS], affinity = 0.0;
	NSUInteger i, j, nz = 0, index[NAS];

	BOOL debugS = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
	BOOL debugV = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
	if (debugV) NSLog(@"Entering [... affinityAndCompositionFromLiquidChemicalPotentialSum] ...");
    for (i=0; i<NAS; i++) gamma[i] = 0.0;

	// Compute solid -> liquid delta mus and deflate composition space
	for (i=0; i<NA; i++) {
		if (chemicalPotentials[i] != 0.0) {
			mu0[i] = [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
			deltaMu[nz] = chemicalPotentials[i] - mu0[i];
			index[nz] = i;
			gamma[nz] = 1.0;
			nz++;
		} else mu0[i] = 0.0;
		x[i] = 0.0;
		xLast[i] = 0.0;
	}

	// NaCa2Fe4FeAl2Si6O22(OH)2
	if ((chemicalPotentials[0] != 0.0) && (chemicalPotentials[1] != 0.0) && (chemicalPotentials[2] != 0.0)) {
		mu0[3] = mu0[2] + mu0[1] - mu0[0] + (DGR);
		deltaMu[nz] = chemicalPotentials[2] + chemicalPotentials[1] - chemicalPotentials[0] - mu0[3];
		index[nz] = 3;
		gamma[nz] = 1.0;
		nz++;
	} else mu0[3] = 0.0;
	x[3] = 0.0;
	xLast[3] = 0.0;

	// There are no non-zero chemical potentials, so the phase can never form.  Return an affinity of zero and zero the mole fractions
	if (nz == 0) {
		[results addObject:[NSNumber numberWithDouble:affinity]];
		for (i=0; i<NA; i++) [results addObject:[NSNumber numberWithDouble:x[i]]];
		[results addObject:[NSNumber numberWithBool:YES]];      // convergence flag
		[results addObject:[NSNumber numberWithUnsignedInteger:0]]; // iteration count
		[results addObject:[NSNumber numberWithDouble:NATOMS]]; // number of atoms used to scale affinity
		[results addObject:[NSNumber numberWithDouble:0.0]];    // likely error in affinity

		if (debugV) NSLog(@"... Terminated. Trivial case.");
		return [NSArray arrayWithArray:results];
	}

	NSUInteger count = 0;
	BOOL converged = NO;
	double affinityLast, xReduced[NA];
	do {
		affinityLast = affinity;
		// Solve for mole ractions in the deflated composition space
		if (nz == 1) {
			xNz[0] = 1.0;
			affinity = -(deltaMu[0]-R*t*log(gamma[0]));
		} else {
			double sum = 1.0;
			if (nz > 2) for (i=0; i<(nz-2); i++) {
				xNz[i] = exp(((deltaMu[i]-R*t*log(gamma[i]))-(deltaMu[nz-2]-R*t*log(gamma[nz-2])))/(R*t));
				sum += xNz[i];
			}
			xNz[nz-2] = exp(((deltaMu[nz-2]-R*t*log(gamma[nz-2]))-(deltaMu[nz-1]-R*t*log(gamma[nz-1])))/(R*t));

			xNz[nz-2] /= 1.0 + xNz[nz-2]*sum;
			xNz[nz-1] = 1.0 - xNz[nz-2];
			if (nz > 2) for (i=0; i<(nz-2); i++) {
				xNz[i] *= xNz[nz-2];
				xNz[nz-1] -= xNz[i];
			}
			for (i=0; i<nz; i++) if (xNz[i] <= DBL_EPSILON) xNz[i] = DBL_EPSILON;

			/* compute the chemical affinity (choice of mu[] is arbitrary) */
			affinity = -(deltaMu[0]-R*t*log(gamma[0])) + R*t*log(xNz[0]);
		}

		// Reinflate the solution
		for (i=0; i<nz; i++) x[index[i]] = xNz[i];

		// Determine activity coefficients
		double a[NAS], mu[NA], r[NR];
		xReduced[0] = x[0] - x[3];
		xReduced[1] = x[1] + x[3];
		xReduced[2] = x[2] + x[3];
		[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:xReduced r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
		[self activity:FIRST | SECOND t:t p:p r:r a:a mu:mu dx:NULL];
		if (![self testPermissibleValuesOfComponents:xReduced]) {
			NSLog(@"Composition estimate is infeasible. T = %f, P = %f, count = %lu", t-273.15, p/10.0, count);
			NSLog(@"species X0 %@ = %g", [[self nameOfSolutionSpeciesAtIndex:0] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[0]);
			NSLog(@"species X1 %@ = %g", [[self nameOfSolutionSpeciesAtIndex:1] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[1]);
			NSLog(@"species X2 %@ = %g", [[self nameOfSolutionSpeciesAtIndex:2] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[2]);
			NSLog(@"species X3 %@ = %g", [[self nameOfSolutionSpeciesAtIndex:3] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[3]);
			NSLog(@"%@ X0 = %13.6g, a = %13.6g, mu = %13.6g", [[self nameOfSolutionSpeciesAtIndex:0] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], xReduced[0], a[0], mu[0]);
			NSLog(@"%@ X1 = %13.6g, a = %13.6g, mu = %13.6g", [[self nameOfSolutionSpeciesAtIndex:1] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], xReduced[1], a[1], mu[1]);
			NSLog(@"%@ X2 = %13.6g, a = %13.6g, mu = %13.6g", [[self nameOfSolutionSpeciesAtIndex:2] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], xReduced[2], a[2], mu[2]);
		}

		for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) gamma[j++] = a[i]/x[i];

		if ((chemicalPotentials[0] != 0.0) && (chemicalPotentials[1] != 0.0) && (chemicalPotentials[2] != 0.0)) {
			a[3] = exp((mu[2] + mu[1] - mu[0] + mu0[2] + mu0[1] - mu0[0] - mu0[3])/(R*t));
			gamma[j] = a[3]/x[3];
		} else gamma[j] = 0.0;

		if (debugV) {
			NSLog(@"Iteration %lu", count);
			NSLog(@"%13.6g %13.6g %13.6g %13.6g %13.6g", a[0], a[1], a[2], a[3], affinity);
			NSLog(@"%13.6g %13.6g %13.6g %13.6g", x[0], x[1], x[2], x[3]);
            double g[NAS];
			for (i=0, j=0; i<NAS; i++) g[i] = (x[i] != 0.0) ? gamma[j++] : 0.0;
            NSLog(@"%13.6g %13.6g %13.6g %13.6g", g[0], g[1], g[2], g[3]);
		}

		if (count > 25) { // pure empiricism
			for (i=0; i<nz; i++) gamma[i] = (gamma[i]+gammaLast[i])/2.0;
		}
		for (i=0; i<nz; i++) gammaLast[i] = gamma[i];

		[self correctActivityCoefficients:gamma forComposition:x];
		converged = (fabs(affinity-affinityLast) < 0.1);
		count++;

	} while (count < 50 && !converged);

	if (debugS) {
		NSLog(@"... Terminated (converged %@) for phase %@ in %lu iterations with affinity %f J (delta %f) for %f atoms.",
			  converged ? @"YES" : @"NO", [self phaseName], count, affinity, fabs(affinity-affinityLast), NATOMS);
		for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) NSLog(@"... ... Activity coefficient of %@ is %13.6g with mole fraction %13.6g (reduced %13.6g)",
														 [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
														 gamma[j++], x[i], xReduced[i]);
		NSLog(@"... ... Activity coefficient of %@ is %13.6g with mole fraction %13.6g", @"ferrohastingite", gamma[j], x[3]);
	}
	if (debugV) NSLog(@"Exiting [... affinityAndCompositionFromLiquidChemicalPotentialSum].");


	[results addObject:[NSNumber numberWithDouble:affinity]];                         // affinity in J
	for (i=0; i<NA; i++) [results addObject:[NSNumber numberWithDouble:xReduced[i]]]; // composition in mole fraction of endmembers
	[results addObject:[NSNumber numberWithBool:converged]];                          // convergence flag
	[results addObject:[NSNumber numberWithUnsignedInteger:count]];                       // iteration count
	[results addObject:[NSNumber numberWithDouble:NATOMS]];                           // number of atoms used to scale affinity
	[results addObject:[NSNumber numberWithDouble:fabs(affinity-affinityLast)]];      // likely error in affinity

	return [NSArray arrayWithArray:results];
}

#import "SolutionPhase.h"

@end
