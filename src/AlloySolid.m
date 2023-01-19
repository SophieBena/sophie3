//
//  AlloySolid.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/11/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import "AlloySolid.h"
#import "BermanProperties.h"
#import "Metals.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

@implementation AlloySolid

static NSArray *endmembers;

#define NR        1    // Independent composition variables
#define NA        2    // Endmember compositions
#define NATOMS  1.0    // Number of atoms in the formula unit

#pragma mark -
#pragma mark class methods

+(void)initialize {
	if (self == [AlloySolid class]) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"Initialize(AlloySolid) - entry ...");
		NSMutableArray *mutableEndmembers = [NSMutableArray arrayWithCapacity:NA];

		FeSolid *feSolid = [[FeSolid alloc] init];
		[feSolid setPhaseFormula:@"Fe"];
		[feSolid setPhaseName:@"Fe solid"];
		[mutableEndmembers addObject:feSolid];
		if (debug) NSLog(@"... allocated Fe solid ...");

		NiSolid *niSolid = [[NiSolid alloc] init];
		[niSolid setPhaseFormula:@"Ni"];
		[niSolid setPhaseName:@"Ni solid"];
		[mutableEndmembers addObject:niSolid];
		if (debug) NSLog(@"... allocated Ni solid ...");

		endmembers = [NSArray arrayWithArray:mutableEndmembers];
	}
}


#pragma mark -
#pragma mark instance methods

-(id)init {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(AlloySolid) ...");

		[self setPhaseName:@"Solid Alloy"];
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
 * Fe-Ni Alloy solution models from Tomiska and Tomiska & Neckel
 */
#define BH1 - 14441.0
#define BH2 -  9182.0
#define BH3 -  1520.0

#define BS1 - 3.2620
#define BS2 - 3.6564
#define BS3 - 1.4098

/*
 #define BH1Liq - 17769.0
 #define BH2Liq -  8942.0
 #define BH3Liq -  1460.0

 #define BS1Liq - 4.1938
 #define BS1Liq - 3.3928
 #define BS1Liq - 0.7806
 */

#define BG1 (BH1) - t*(BS1)
#define BG2 (BH2) - t*(BS2)
#define BG3 (BH3) - t*(BS3)

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conAlS defines the conversion from m[i], to r[j]
 */
#define NR         1
#define NS         0
#define NA         2
#define FR0(i)     (i == 1) ? 1.0 - xni : - xni
#define DFR0DR0(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */
#define R          8.3143
#define S          - R*(xni*log(xni)+xfe*log(xfe)) \
+ xni*xfe*((BS1) + (BS2)*(xni-xfe) + (BS3)*(xni-xfe)*(xni-xfe))
#define H          xni*xfe*((BH1) + (BH2)*(xni-xfe) + (BH3)*(xni-xfe)*(xni-xfe))
#define V          0.0
#define G          H - t*(S) + (p-1.0)*(V)

#define DGDR0      R*t*(log(xni)-log(xfe)) \
+ (xfe-xni)*((BG1) + (BG2)*(xni-xfe) + (BG3)*(xni-xfe)*(xni-xfe)) \
+ xfe*xni*(2.0*(BG2) + 4.0*(BG3)*(xni-xfe))
#define DGDT       -(S)
#define DGDP       (V)

#define D2GDR0R0   R*t*(1.0/xni+1.0/xfe) \
- 2.0*((BG1) + (BG2)*(xni-xfe) + (BG3)*(xni-xfe)*(xni-xfe)) \
+ 2.0*(xfe-xni)*(2.0*(BG2) + 4.0*(BG3)*(xni-xfe)) \
+ xfe*xni*(8.0*(BG3))
#define D2GDR0DT   R*(log(xni)-log(xfe)) \
- (xfe-xni)*((BS1) + (BS2)*(xni-xfe) + (BS3)*(xni-xfe)*(xni-xfe)) \
- xfe*xni*(2.0*(BS2) + 4.0*(BS3)*(xni-xfe))
#define D2GDR0DP   0.0
#define D2GDT2     0.0
#define D2GDTDP    0.0
#define D2GDP2     0.0

#define D3GDR0R0R0 R*t*(1.0/SQUARE(xfe) - 1.0/SQUARE(xni)) \
- 6.0*(2.0*(BG2) + 4.0*(BG3)*(xni-xfe)) \
+ 2.0*(xfe-xni)*(8.0*(BG3)) \
+ (xfe-xni)*(8.0*(BG3))
#define D3GDR0R0DT R*(1.0/xni+1.0/xfe) \
+ 2.0*((BS1) + (BS2)*(xni-xfe) + (BS3)*(xni-xfe)*(xni-xfe)) \
- 2.0*(xfe-xni)*(2.0*(BS2) + 4.0*(BS3)*(xni-xfe)) \
- xfe*xni*(8.0*(BS3))
#define D3GDR0R0DP 0.0
#define D3GDT3     0.0
#define D3GDT2DP   0.0
#define D3GDTDP2   0.0
#define D3GDP3     0.0

#define D3GDR0DT2  0.0
#define D3GDR0DTDP 0.0
#define D3GDR0DP2  0.0

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
	const char *NAMES[NA]    = { "Fe-metal", "Ni-metal" };
	const char *FORMULAS[NA] = { "Fe", "Ni" };
    NSInteger i;
	BOOL result = TRUE;
	double sum;

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
		for (i=0, sum=0.0; i<NR; i++) {
			result = result && (r[i] >= 0.0) && (r[i] <= 1.0);
			sum += r[i];
		}
		result = result && (sum <= 1.0);
	}
	/* Check bounds on moles of endmember components */
	if (mask & SIXTH) {
		for (i=0; i<NA; i++) result = result && (m[i] >= 0.0);
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
	 endmember alloy components.
	 (2) calculates from a vector of moles of endmember components, one or
	 all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
	 (3) calculates from a vector of independent compositional variables
	 mole fractions of endmember components and/or the Jacobian matrix
	 dx[]/dr[]

	 In this routine it is assumed that the elements are in the order of atomic
	 numbers and that the order of biotite components has been verified as:
	 m[0] = Fe-metal (Fe),
	 m[1] = Ni-metal (Ni) .

	 ----------------------------------------------------------------------------*/

	int i, j, k;

	if (inpMask == FIRST && outMask == SECOND) {
		static const int Fe = 26;
		static const int Ni = 28;

		/* Projection into the Fe-Ni binary */
		m[0] = e[Fe]; /* moles of Fe                      */
		m[1] = e[Ni]; /* moles of Ni                      */

	} else if (inpMask == SECOND) {
		double sum;

		if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
			NSLog(@"Illegal call to conFld with inpMask = %o and outMask = %o\n", inpMask, outMask);

		for (i=0, sum=0.0; i<NA; i++) sum += m[i];

		if (outMask & THIRD) {
			r[0] = (sum != 0.0) ? m[1]/sum : 0.0;
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
				for (j=0; j<NA; j++) dm[0][j] = (1 == j) ? (1.0-m[1]/sum)/sum : - m[1]/SQUARE(sum);
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
			/* calculates the matrix d3r[i]/dm[j]dm[k]dm[l] using m[] as input	*/
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

			x[1] = r[0];
            x[0] = 1.0 - r[0];
		}

		if (outMask & SEVENTH) {
			/* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
			dr[1][0] =  1.0;
			dr[0][0] = -1.0;
		}

	} else  {
		NSLog(@"Illegal call to conFld with inpMask = %o and outMask = %o\n", inpMask, outMask);
	}

}

-(NSString *)displayFormula:(double)t
						  p:(double)p
						  r:(double [NA])r
{
	double totFe, totNi;
	totFe = 1.0 - r[0];
	totNi = r[0];

	return [NSString stringWithFormat:@"Fe%4.2fNi%4.2f", totFe, totNi];
}

-(void)activity:(int)mask
			  t:(double)t
			  p:(double)p
			  r:(double [NA])r
			  a:(double [NA])a      // (pointer to a[]) activities              BINARY MASK: 0001
			 mu:(double [NA])mu     // (pointer to mu[]) chemical potentials    BINARY MASK: 0010
			 dx:(double [NA][NR])dx // (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100
{
	double xni = (r[0]     > DBL_EPSILON) ? r[0]     : DBL_EPSILON;
	double xfe = (1.0-r[0] > DBL_EPSILON) ? 1.0-r[0] : DBL_EPSILON;

	double g, dgdr[NR], fr[NA][NR];
	int i, j;

	for(i=0; i<NA; i++) fr[i][0] = FR0(i);

	g       = G;
	dgdr[0] = DGDR0;

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

		for(i=0; i<NA; i++) dfrdr[i][0] = DFR0DR0(i);

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
	double xni = (r[0]     > DBL_EPSILON) ? r[0]     : DBL_EPSILON;
	double xfe = (1.0-r[0] > DBL_EPSILON) ? 1.0-r[0] : DBL_EPSILON;

	if (mask & FIRST) {
		*gmix = G;
	}

	if(mask & SECOND) {
		dx[0] = DGDR0;
	}

	if(mask & THIRD) {
		double d2gdr2[NR][NR];
		int i, j;

		d2gdr2[0][0] = D2GDR0R0;

		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) dx2[i][j] = d2gdr2[i][j];
		}
	}

	if (mask & FOURTH) {
		double d3gdr3[NR][NR][NR];
		int i, j, k;

		d3gdr3[0][0][0] = D3GDR0R0R0;

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
	double xni = (r[0]     > DBL_EPSILON) ? r[0]     : DBL_EPSILON;
	double xfe = (1.0-r[0] > DBL_EPSILON) ? 1.0-r[0] : DBL_EPSILON;

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
	double xni = (r[0]     > DBL_EPSILON) ? r[0]     : DBL_EPSILON;
	double xfe = (1.0-r[0] > DBL_EPSILON) ? 1.0-r[0] : DBL_EPSILON;

	if (mask & FIRST) {
		*smix = S;
	}

	if(mask & SECOND) {
		double d2gdrdt[NR];
		int i;

		d2gdrdt[0] = D2GDR0DT;

		for (i=0; i<NR; i++) dx[i] = - d2gdrdt[i];
	}

	if(mask & THIRD) {
		double d3gdr2dt[NR][NR];
		int i, j;

		d3gdr2dt[0][0] = D3GDR0R0DT;

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

		for (i=0; i<NR; i++) dx[i] = d2gdrdp[i];
	}

	if(mask & THIRD) {
		double d3gdr2dp[NR][NR];
		int i, j;

		d3gdr2dp[0][0] = D3GDR0R0DP;

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

		for (i=0; i<NR; i++) dxdt[i] = d3gdrdtdp[i];
	}

	if(mask & TENTH) {
		double d3gdrdp2[NR];
		int i;

		d3gdrdp2[0] = D3GDR0DP2;

		for (i=0; i<NR; i++) dxdp[i] = d3gdrdp2[i];
	}

}

-(NSUInteger)numberOfSolutionSpecies {
	return NA;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
	switch (index) {
		case 0:
			return @"Fe solid";
			break;
		case 1:
			return @"Ni solid";
			break;
		default:
			return @"";
			break;
	}
}

-(DoubleVector *)convertMolesOfSpeciesToMolesOfComponents:(double *)mSpecies {
    DoubleVector *mComponentWrapper =[[DoubleVector alloc] initWithSize:NA];
	double *mComponents = [mComponentWrapper pointerToDouble];
	for (NSUInteger i=0; i<NA; i++) mComponents[i] = mSpecies[i];
	return mComponentWrapper;
}

-(DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	return [self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p];
}

-(DoubleVector *)elementalCompositionOfSpeciesAtIndex:(NSUInteger)index {
	return [[endmembers objectAtIndex:index] formulaAsElementArray];
}

-(void)correctActivityCoefficients:(double [NA])gamma forComposition:(double [NA])x { }

// --> SolutionPhaseProtocol public function
-(NSArray *)affinityAndCompositionFromLiquidChemicalPotentialSum:(double *)chemicalPotentials andT:(double)t andP:(double)p {
	NSMutableArray *results = [NSMutableArray arrayWithCapacity:NA+1];
	double deltaMu[NA], xNz[NA], x[NA], gamma[NA], xLast[NA], affinity = 0.0;
	NSUInteger i, j, nz = 0, index[NA];
    for (i=0; i<NA; i++) gamma[i] = 0.0;

	BOOL debugS = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
	BOOL debugV = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
	if (debugV) NSLog(@"Entering [... affinityAndCompositionFromLiquidChemicalPotentialSum] ...");

	// Compute solid -> liquid delta mus and deflate composition space
	for (i=0; i<NA; i++) {
		if (chemicalPotentials[i] != 0.0) {
			deltaMu[nz] = chemicalPotentials[i] - [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
			index[nz] = i;
			gamma[nz] = 1.0;
			nz++;
		}
		x[i] = 0.0;
		xLast[i] = 0.0;
	}

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
	double affinityLast;
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

			/* compute the chemical affinity (choice of mu[] is arbitrary) */
			affinity = -(deltaMu[0]-R*t*log(gamma[0])) + R*t*log(xNz[0]);
		}

		// Reinflate the solution
		for (i=0; i<nz; i++) x[index[i]] = xNz[i];

		// Determine activity coefficients
		double a[NA], r[NR];
		[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:x r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
		[self activity:FIRST t:t p:p r:r a:a mu:NULL dx:NULL];

		for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) gamma[j++] = a[i]/x[i];
		[self correctActivityCoefficients:gamma forComposition:x];
		converged = (fabs(affinity-affinityLast) < 0.1);
		count++;

	} while (count < 50 && !converged);

	if (debugS) {
		NSLog(@"... Terminated (converged %@) for phase %@ in %lu iterations with affinity %f (delta %f) J for %f atoms.",
			  converged ? @"YES" : @"NO", [self phaseName], count, affinity, fabs(affinity-affinityLast), NATOMS);
		for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) NSLog(@"... ... Activity coefficient of component %@ is %f with mole fraction %f",
														 [[endmembers objectAtIndex:i] phaseName], gamma[j++], x[i]);
	}

	[results addObject:[NSNumber numberWithDouble:affinity]];                    // affinity in J
	for (i=0; i<NA; i++) [results addObject:[NSNumber numberWithDouble:x[i]]];   // composition in mole fraction of endmembers
	[results addObject:[NSNumber numberWithBool:converged]];                     // convergence flag
	[results addObject:[NSNumber numberWithUnsignedInteger:count]];                  // iteration count
	[results addObject:[NSNumber numberWithDouble:NATOMS]];                      // number of atoms used to scale affinity
	[results addObject:[NSNumber numberWithDouble:fabs(affinity-affinityLast)]]; // likely error in affinity

	if (debugV) NSLog(@"Exiting [... affinityAndCompositionFromLiquidChemicalPotentialSum].");
	return [NSArray arrayWithArray:results];
}

#import "SolutionPhase.h"

@end
