//
//  GarnetBerman.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/3/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import "GarnetBerman.h"
#import "BermanProperties.h"
#import "BermanStoichiometricPhases.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"
#import "MathSupport.h"

@implementation GarnetBerman

static NSArray *endmembers;

#define NR        2 // Two independent composition variables
#define NA        3 // Three endmember compositions
#define NATOMS 20.0 // Number of atoms in the formula unit

#pragma mark -
#pragma mark class methods

+(void)initialize {
	if (self == [GarnetBerman class]) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"Initialize(GarnetBerman) - entry ...");
		NSMutableArray *mutableEndmembers = [NSMutableArray arrayWithCapacity:NA];

		BermanProperties *almandine = [[BermanProperties alloc] initWithH:-5267216.0
																		S:340.007
																	   k0:573.96
																	   k1:-14.831E2
																	   k2:-292.920E5
																	   k3:502.208E7
																	   v0:11.511
																	   v1:-0.558E-6
																	   v2:0.321E-12
																	   v3:18.613E-6
																	   v4:74.539E-10];
		[almandine setPhaseFormula:@"Fe3Al2Si3O12"];
		[almandine setPhaseName:@"almandine"];
		[mutableEndmembers addObject:almandine];
		if (debug) NSLog(@"... allocated almandine ...");

		BermanProperties *grossular = [[BermanProperties alloc] initWithH:-6632859.0
																		S:255.150
																	   k0:573.43
																	   k1:-20.394E2
																	   k2:-188.872E5
																	   k3:231.931E7
																	   v0:12.538
																	   v1:-0.654E-6
																	   v2:1.635E-12
																	   v3:18.994E-6
																	   v4:79.756E-10];
		[grossular setPhaseFormula:@"Ca3Al2Si3O12"];
		[grossular setPhaseName:@"grossular"];
		[mutableEndmembers addObject:grossular];
		if (debug) NSLog(@"... allocated grossular ...");

		BermanProperties *pyrope = [[BermanProperties alloc] initWithH:-6286548.0
																	 S:266.359
																	k0:640.72
																	k1:-45.421E2
																	k2:-47.019E5
																	k3:0.0E7
																	v0:11.316
																	v1:-0.576E-6
																	v2:0.442E-12
																	v3:22.519E-6
																	v4:37.044E-10];
		[pyrope setPhaseFormula:@"Mg3Al2Si3O12"];
		[pyrope setPhaseName:@"pyrope"];
		[mutableEndmembers addObject:pyrope];
		if (debug) NSLog(@"... allocated pyrope ...");

		endmembers = [NSArray arrayWithArray:mutableEndmembers];
	}
}

#pragma mark -
#pragma mark instance methods

-(id)init {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(GarnetBerman) ... entry ...");

		[self setPhaseName:@"Garnet"];
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
 * Garnet solution parameters:
 *
 * Berman, R.G.
 * Mixing properties of Ca-Mg-Fe-Mn garnets
 * American Mineralogist 75, 328-344
 *
 * Berman, R.G., Koziol, A.M.
 * Ternary excess properties of grossular-pyrope-almandine garnet and
 * their influence in geothermobarometry
 * American Mineralogist 76, 1223-1231
 *
 * 1 == Grossular, 2 == Pyrope, 3 == Almandine
 */
#define WH112 21560.0   /* joules     */
#define WS112    18.79  /* joules/K   */
#define WV112     0.10  /* joules/bar */
#define WH122 69200.0   /* joules     */
#define WS122    18.79  /* joules/K   */
#define WV122     0.10  /* joules/bar */
#define WH113 20320.0   /* joules     */
#define WS113     5.08  /* joules/K   */
#define WV113     0.17  /* joules/bar */
#define WH133  2620.0   /* joules     */
#define WS133     5.08  /* joules/K   */
#define WV133     0.09  /* joules/bar */
#define WH223   230.0   /* joules     */
#define WS223     0.0   /* joules/K   */
#define WV223     0.01  /* joules/bar */
#define WH233  3720.0   /* joules     */
#define WS233     0.0   /* joules/K   */
#define WV233     0.06  /* joules/bar */
#define WH123     0.0   /* joules     */
#define WS123     0.0   /* joules/K   */
#define WV123     0.00  /* joules/bar */

#define WG112  (WH112-t*WS112+p*WV112)
#define WG122  (WH122-t*WS122+p*WV122)
#define WG113  (WH113-t*WS113+p*WV113)
#define WG133  (WH133-t*WS133+p*WV133)
#define WG223  (WH223-t*WS223+p*WV223)
#define WG233  (WH233-t*WS233+p*WV233)
#define WG123  (WH123-t*WS123+p*WV123)

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conFld defines the conversion from m[i], to r[j]
 */
#define FR0(i)     (i == 0) ? 1.0 - xal : - xal
#define FR1(i)     (i == 1) ? 1.0 - xgr : - xgr
#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */
#define R         8.3143
#define S         - 3.0*R*(xgr*log(xgr) + xpy*log(xpy) + xal*log(xal)) + \
WS112*xgr*xpy*(xgr+xal/2.0) + WS122*xgr*xpy*(xpy+xal/2.0) + \
WS113*xgr*xal*(xgr+xpy/2.0) + WS133*xgr*xal*(xal+xpy/2.0) + \
WS223*xpy*xal*(xpy+xgr/2.0) + WS233*xpy*xal*(xal+xgr/2.0) + \
WS123*xgr*xpy*xal
#define H         WH112*xgr*xpy*(xgr+xal/2.0) + WH122*xgr*xpy*(xpy+xal/2.0) + \
WH113*xgr*xal*(xgr+xpy/2.0) + WH133*xgr*xal*(xal+xpy/2.0) + \
WH223*xpy*xal*(xpy+xgr/2.0) + WH233*xpy*xal*(xal+xgr/2.0) + \
WH123*xgr*xpy*xal
#define V         WV112*xgr*xpy*(xgr+xal/2.0) + WV122*xgr*xpy*(xpy+xal/2.0) + \
WV113*xgr*xal*(xgr+xpy/2.0) + WV133*xgr*xal*(xal+xpy/2.0) + \
WV223*xpy*xal*(xpy+xgr/2.0) + WV233*xpy*xal*(xal+xgr/2.0) + \
WV123*xgr*xpy*xal
#define G         H - t*(S) + p*(V)

#define DGDR0     3.0*R*t*(log(xal) - log(xpy)) + \
WG112*( xgr*xpy/2.0 - xgr*(xgr+xal/2.0)) + \
WG122*(-xgr*xpy/2.0 - xgr*(xpy+xal/2.0)) + \
WG113*(-xgr*xal/2.0 + xgr*(xgr+xpy/2.0)) + \
WG133*( xgr*xal/2.0 + xgr*(xal+xpy/2.0)) + \
WG223*(-xpy*xal     + (xpy-xal)*(xpy+xgr/2.0)) + \
WG233*( xpy*xal     + (xpy-xal)*(xal+xgr/2.0)) + \
WG123*xgr*(xpy-xal)
#define DGDR1     3.0*R*t*(log(xgr) - log(xpy)) + \
WG112*( xgr*xpy     + (xpy-xgr)*(xgr+xal/2.0)) + \
WG122*(-xgr*xpy     + (xpy-xgr)*(xpy+xal/2.0)) + \
WG113*( xgr*xal/2.0 + xal*(xgr+xpy/2.0)) + \
WG133*(-xgr*xal/2.0 + xal*(xal+xpy/2.0)) + \
WG223*(-xpy*xal/2.0 - xal*(xpy+xgr/2.0)) + \
WG233*( xpy*xal/2.0 - xal*(xal+xgr/2.0)) + \
WG123*(xpy-xgr)*xal
#define DGDP      (V)

#define D2GDR0R0  3.0*R*t*(1.0/xal + 1.0/xpy) + \
xgr*(WG122-WG112+WG133-WG113) + 2.0*(xpy-xal)*(WG233-WG223) \
- 2.0*WG223*(xpy+xgr/2.0) - 2.0*WG233*(xal+xgr/2.0) \
- 2.0*WG123*xgr
#define D2GDR0R1  3.0*R*t*(1.0/xpy) + \
WG112*( (xpy-xgr)/2.0 - (xgr+xal/2.0) - xgr) + \
WG122*(-(xpy-xgr)/2.0 - (xpy+xal/2.0) + xgr) + \
WG113*(-xal/2.0 + (xgr+xpy/2.0) + xgr/2.0) + \
WG133*( xal/2.0 + (xal+xpy/2.0) - xgr/2.0) + \
WG223*( xal     - (xpy+xgr/2.0) - (xpy-xal)/2.0) + \
WG233*(-xal     - (xal+xgr/2.0) + (xpy-xal)/2.0) + \
WG123*(xpy-xal-xgr)
#define D2GDR0DT  3.0*R*(log(xal) - log(xpy)) - (\
WS112*( xgr*xpy/2.0 - xgr*(xgr+xal/2.0)) + \
WS122*(-xgr*xpy/2.0 - xgr*(xpy+xal/2.0)) + \
WS113*(-xgr*xal/2.0 + xgr*(xgr+xpy/2.0)) + \
WS133*( xgr*xal/2.0 + xgr*(xal+xpy/2.0)) + \
WS223*(-xpy*xal     + (xpy-xal)*(xpy+xgr/2.0)) + \
WS233*( xpy*xal     + (xpy-xal)*(xal+xgr/2.0)) + \
WS123*xgr*(xpy-xal) )
#define D2GDR0DP  WV112*( xgr*xpy/2.0 - xgr*(xgr+xal/2.0)) + \
WV122*(-xgr*xpy/2.0 - xgr*(xpy+xal/2.0)) + \
WV113*(-xgr*xal/2.0 + xgr*(xgr+xpy/2.0)) + \
WV133*( xgr*xal/2.0 + xgr*(xal+xpy/2.0)) + \
WV223*(-xpy*xal     + (xpy-xal)*(xpy+xgr/2.0)) + \
WV233*( xpy*xal     + (xpy-xal)*(xal+xgr/2.0)) + \
WV123*xgr*(xpy-xal)
#define D2GDR1R1  3.0*R*t*(1.0/xgr + 1.0/xpy) + \
2.0*(xpy-xgr)*(WG112-WG122) - 2.0*WG112*(xgr+xal/2.0) \
- 2.0*WG122*(xpy+xal/2.0) + xal*(WG113-WG133+WG223-WG233) \
- 2.0*WG123*xal
#define D2GDR1DT  3.0*R*(log(xgr) - log(xpy)) - ( \
WS112*( xgr*xpy     + (xpy-xgr)*(xgr+xal/2.0)) + \
WS122*(-xgr*xpy     + (xpy-xgr)*(xpy+xal/2.0)) + \
WS113*( xgr*xal/2.0 + xal*(xgr+xpy/2.0)) + \
WS133*(-xgr*xal/2.0 + xal*(xal+xpy/2.0)) + \
WS223*(-xpy*xal/2.0 - xal*(xpy+xgr/2.0)) + \
WS233*( xpy*xal/2.0 - xal*(xal+xgr/2.0)) + \
WS123*(xpy-xgr)*xal )
#define D2GDR1DP  WV112*( xgr*xpy     + (xpy-xgr)*(xgr+xal/2.0)) + \
WV122*(-xgr*xpy     + (xpy-xgr)*(xpy+xal/2.0)) + \
WV113*( xgr*xal/2.0 + xal*(xgr+xpy/2.0)) + \
WV133*(-xgr*xal/2.0 + xal*(xal+xpy/2.0)) + \
WV223*(-xpy*xal/2.0 - xal*(xpy+xgr/2.0)) + \
WV233*( xpy*xal/2.0 - xal*(xal+xgr/2.0)) + \
WV123*(xpy-xgr)*xal
#define D2GDT2    0.0
#define D2GDTDP   0.0
#define D2GDP2    0.0

#define D3GDR0R0R0 3.0*R*t*(1.0/SQUARE(xpy)-1.0/SQUARE(xal)) + \
6.0*WG223 - 6.0*WG233
#define D3GDR0R0R1 3.0*R*t*(1.0/SQUARE(xpy)) + WG122 - WG112 \
+ WG133 - WG113 - 3.0*WG233 + 3.0*WG223 - 2.0*WG123
#define D3GDR0R0DT 3.0*R*(1.0/xal + 1.0/xpy) - ( \
xgr*(WS122-WS112+WS133-WS113) \
+ 2.0*(xpy-xal)*(WS233-WS223) - 2.0*WS223*(xpy+xgr/2.0) \
- 2.0*WS233*(xal+xgr/2.0) - 2.0*WS123*xgr )
#define D3GDR0R0DP xgr*(WV122-WV112+WV133-WV113) \
+ 2.0*(xpy-xal)*(WV233-WV223) - 2.0*WV223*(xpy+xgr/2.0) \
- 2.0*WV233*(xal+xgr/2.0) - 2.0*WV123*xgr
#define D3GDR0R1R1 3.0*R*t*(1.0/SQUARE(xpy)) + WG113 - WG133 \
+ WG223 - WG233 - 3.0*WG112 + 3.0*WG122 - 2.0*WG123
#define D3GDR0R1DT 3.0*R*(1.0/xpy) - ( \
WS112*( (xpy-xgr)/2.0 - (xgr+xal/2.0) - xgr) + \
WS122*(-(xpy-xgr)/2.0 - (xpy+xal/2.0) + xgr) + \
WS113*(-xal/2.0 + (xgr+xpy/2.0) + xgr/2.0) + \
WS133*( xal/2.0 + (xal+xpy/2.0) - xgr/2.0) + \
WS223*( xal     - (xpy+xgr/2.0) - (xpy-xal)/2.0) + \
WS233*(-xal     - (xal+xgr/2.0) + (xpy-xal)/2.0) + \
WS123*(xpy-xal-xgr) )
#define D3GDR0R1DP WV112*( (xpy-xgr)/2.0 - (xgr+xal/2.0) - xgr) + \
WV122*(-(xpy-xgr)/2.0 - (xpy+xal/2.0) + xgr) + \
WV113*(-xal/2.0 + (xgr+xpy/2.0) + xgr/2.0) + \
WV133*( xal/2.0 + (xal+xpy/2.0) - xgr/2.0) + \
WV223*( xal     - (xpy+xgr/2.0) - (xpy-xal)/2.0) + \
WV233*(-xal     - (xal+xgr/2.0) + (xpy-xal)/2.0) + \
WV123*(xpy-xal-xgr)
#define D3GDR1R1R1 3.0*R*t*(1.0/SQUARE(xpy)-1.0/SQUARE(xgr)) + \
6.0*WG122 - 6.0*WG112
#define D3GDR1R1DT 3.0*R*(1.0/xgr + 1.0/xpy) - ( \
2.0*(xpy-xgr)*(WS112-WS122) - 2.0*WS112*(xgr+xal/2.0) \
- 2.0*WS122*(xpy+xal/2.0) + xal*(WS113-WS133+WS223-WS233) \
- 2.0*WS123*xal )
#define D3GDR1R1DP 2.0*(xpy-xgr)*(WV112-WV122) - 2.0*WV112*(xgr+xal/2.0) \
- 2.0*WV122*(xpy+xal/2.0) + xal*(WV113-WV133+WV223-WV233) \
- 2.0*WV123*xal
#define D3GDT3     0.0
#define D3GDT2DP   0.0
#define D3GDTDP2   0.0
#define D3GDP3     0.0

#define D3GDR0DT2  0.0
#define D3GDR0DTDP 0.0
#define D3GDR0DP2  0.0
#define D3GDR1DT2  0.0
#define D3GDR1DTDP 0.0
#define D3GDR1DP2  0.0

/*
 *=============================================================================
 *    mask  -  bitwise mask for selecting output
 *    t     -  Temperature (K)
 *    p     -  Pressure (bars)
 *    r[NR] -  Array of independent compositional variables
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
	const char *NAMES[NA]    = { "almandine", "grossular", "pyrope" };
	const char *FORMULAS[NA] = { "Fe3Al2Si3O12", "Ca3Al2Si3O12", "Mg3Al2Si3O12" };
	int i;
	BOOL result = YES;
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
	 endmember garnet components.
	 (2) calculates from a vector of moles of endmember components, one or
	 all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
	 (3) calculates from a vector of independent compositional variables
	 mole fractions of endmember components and/or the Jacobian matrix
	 dx[]/dr[]

	 In this routine it is assumed that the elements are in the order of atomic
	 numbers and that the order of garnet components has been verified as:
	 m[0] = almandine (Fe3Al2Si3O12) ,
	 m[1] = grossular (Ca3Al2Si3O12) and
	 m[2] = pyrope    (Mg3Al2Si3O12)

	 ----------------------------------------------------------------------------*/

	int i, j, k;

	if (inpMask == FIRST && outMask == SECOND) {
		static const int Mg = 12;
		static const int Ca = 20;
		static const int Fe = 26;

		/* Projection into the Fe, Ca, Mg triangle */
		m[0] = e[Fe]/3.0; /* moles of Fe3Al2Si3O12                   */
		m[1] = e[Ca]/3.0; /* Moles of Ca3Al2Si3O12                   */
		m[2] = e[Mg]/3.0; /* Moles of Mg3Al2Si3O12                   */

	} else if (inpMask == SECOND) {
		double sum;

		if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
			NSLog(@"Illegal call to conGrn with inpMask = %o and outMask = %o", inpMask, outMask);

		for (i=0, sum=0.0; i<NA; i++) sum += m[i];

		if (outMask & THIRD) {
			for (i=0; i<NR; i++) r[i] = (sum != 0.0) ? m[i]/sum : 0.0;
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
						dm[i][j] = (i == j) ? (1.0-m[i]/sum)/sum : - m[i]/SQUARE(sum);
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
							d2m[i][j][k]  = 2.0*m[i]/CUBE(sum);
							d2m[i][j][k] -= (i == j) ? 1.0/SQUARE(sum) : 0.0;
							d2m[i][j][k] -= (i == k) ? 1.0/SQUARE(sum) : 0.0;
						}
					}
				}
			}

		}

		if (outMask & EIGHTH) {
			/* Calculates the matrix d3r[i]/dm[j]dm[k]dm[l] using m[] as input      */
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
					for (j=0; j<NA; j++) {
						for (k=0; k<NA; k++) {
							for (l=0; l<NA; l++) {
								d3m[i][j][k][l]  = -6.0*m[i]/QUARTIC(sum);
								d3m[i][j][k][l] += (i == j) ? 2.0/CUBE(sum) : 0.0;
								d3m[i][j][k][l] += (i == k) ? 2.0/CUBE(sum) : 0.0;
								d3m[i][j][k][l] += (i == l) ? 2.0/CUBE(sum) : 0.0;
							}
						}
					}
				}
			}

		}
	} else if (inpMask == THIRD) {

		if (outMask & ~(FOURTH | SEVENTH))
			NSLog(@"Illegal call to conGrn with inpMask = %o and outMask = %o", inpMask, outMask);

		if (outMask & FOURTH) {
			/* Converts a vector of independent compositional variables (r)
			 into a vector of mole fractions of endmember components (x).         */

			for (i=0, x[2]=1.0; i<NR; i++) { x[i] = r[i]; x[2] -= r[i]; }
		}

		if (outMask & SEVENTH) {
			/* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
			for (i=0; i<NR; i++) for (j=0; j<NR; j++) dr[i][j] = (i == j) ? 1.0 : 0.0;
			for (j=0; j<NR; j++) dr[2][j] = -1.0;
		}

	} else  {
		NSLog(@"Illegal call to conGrn with inpMask = %o and outMask = %o", inpMask, outMask);
	}

}

-(NSString *)displayFormula:(double)t
						  p:(double)p
						  r:(double [NA])r
{
	double totCa, totFe2, totMg;

	totCa  = r[1];
	totFe2 = r[0];
	totMg  = 1.0 - r[0] - r[1];

	return [NSString stringWithFormat:@"(Ca%4.2fFe''%4.2fMg%4.2f)3Al2Si3O12", totCa, totFe2, totMg];
}

-(void)activity:(int)mask
			  t:(double)t
			  p:(double)p
			  r:(double [NA])r
			  a:(double [NA])a      // (pointer to a[]) activities              BINARY MASK: 0001
			 mu:(double [NA])mu     // (pointer to mu[]) chemical potentials    BINARY MASK: 0010
			 dx:(double [NA][NR])dx // (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100
{
	double xal = (r[0]          > DBL_EPSILON) ? r[0]          : DBL_EPSILON;
	double xgr = (r[1]          > DBL_EPSILON) ? r[1]          : DBL_EPSILON;
	double xpy = (1.0-r[0]-r[1] > DBL_EPSILON) ? 1.0-r[0]-r[1] : DBL_EPSILON;

	double g, dgdr[NR], fr[NA][NR];
	int i, j;

	for(i=0; i<NA; i++) {
		fr[i][0] = FR0(i);
		fr[i][1] = FR1(i);
	}

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
		d2gdr2[0][1] = D2GDR0R1;
		d2gdr2[1][0] = d2gdr2[0][1];
		d2gdr2[1][1] = D2GDR1R1;

		for(i=0; i<NA; i++) {
			dfrdr[i][0] = DFR0DR0(i);
			dfrdr[i][1] = DFR1DR1(i);
		}

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
	double xal = (r[0]          > DBL_EPSILON) ? r[0]          : DBL_EPSILON;
	double xgr = (r[1]          > DBL_EPSILON) ? r[1]          : DBL_EPSILON;
	double xpy = (1.0-r[0]-r[1] > DBL_EPSILON) ? 1.0-r[0]-r[1] : DBL_EPSILON;

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
		d2gdr2[0][1] = D2GDR0R1;
		d2gdr2[1][0] = d2gdr2[0][1];
		d2gdr2[1][1] = D2GDR1R1;

		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) dx2[i][j] = d2gdr2[i][j];
		}
	}

	if(mask & FOURTH) {
		double d3gdr3[NR][NR][NR];
		int i, j, k;

		d3gdr3[0][0][0] = D3GDR0R0R0;
		d3gdr3[0][0][1] = D3GDR0R0R1;
		d3gdr3[0][1][0] = d3gdr3[0][0][1];
		d3gdr3[0][1][1] = D3GDR0R1R1;
		d3gdr3[1][0][0] = d3gdr3[0][0][1];
		d3gdr3[1][0][1] = d3gdr3[0][1][1];
		d3gdr3[1][1][0] = d3gdr3[0][1][1];
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
	double xal = (r[0]          > DBL_EPSILON) ? r[0]          : DBL_EPSILON;
	double xgr = (r[1]          > DBL_EPSILON) ? r[1]          : DBL_EPSILON;
	double xpy = (1.0-r[0]-r[1] > DBL_EPSILON) ? 1.0-r[0]-r[1] : DBL_EPSILON;

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
	double xal = (r[0]          > DBL_EPSILON) ? r[0]          : DBL_EPSILON;
	double xgr = (r[1]          > DBL_EPSILON) ? r[1]          : DBL_EPSILON;
	double xpy = (1.0-r[0]-r[1] > DBL_EPSILON) ? 1.0-r[0]-r[1] : DBL_EPSILON;

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
		d3gdr2dt[0][1] = D3GDR0R1DT;
		d3gdr2dt[1][0] = d3gdr2dt[0][1];
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
	double xal = (r[0]          > DBL_EPSILON) ? r[0]          : DBL_EPSILON;
	double xgr = (r[1]          > DBL_EPSILON) ? r[1]          : DBL_EPSILON;
	double xpy = (1.0-r[0]-r[1] > DBL_EPSILON) ? 1.0-r[0]-r[1] : DBL_EPSILON;

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
		d3gdr2dp[0][1] = D3GDR0R1DP;
		d3gdr2dp[1][0] = d3gdr2dp[0][1];
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

-(NSUInteger)numberOfSolutionSpecies {
	return NA;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
	switch (index) {
		case 0:
			return @"almandine";
			break;
		case 1:
			return @"grossular";
			break;
		case 2:
			return @"pyrope";
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
				xNz[i] = exp(((deltaMu[i]-R*t*3.0*log(gamma[i]))-(deltaMu[nz-2]-R*t*3.0*log(gamma[nz-2])))/(3.0*R*t));
				sum += xNz[i];
			}
			xNz[nz-2] = exp(((deltaMu[nz-2]-R*t*3.0*log(gamma[nz-2]))-(deltaMu[nz-1]-R*t*3.0*log(gamma[nz-1])))/(3.0*R*t));

			xNz[nz-2] /= 1.0 + xNz[nz-2]*sum;
			xNz[nz-1] = 1.0 - xNz[nz-2];
			if (nz > 2) for (i=0; i<(nz-2); i++) {
				xNz[i] *= xNz[nz-2];
				xNz[nz-1] -= xNz[i];
			}

			for (i=0; i<nz; i++) if (xNz[i] <= DBL_EPSILON) xNz[i] = DBL_EPSILON;

			/* compute the chemical affinity (choice of mu[] is arbitrary) */
			affinity = -(deltaMu[0]-R*t*3.0*log(gamma[0])) + R*t*3.0*log(xNz[0]);
		}

		// Reinflate the solution
		for (i=0; i<nz; i++) x[index[i]] = xNz[i];

		// Determine activity coefficients
		double a[NA], r[NR];
		[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:x r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
		[self activity:FIRST t:t p:p r:r a:a mu:NULL dx:NULL];

		for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) gamma[j++] = pow(a[i], 1.0/3.0)/x[i];
		[self correctActivityCoefficients:gamma forComposition:x];

		if (debugV) {
			NSLog(@"Iteration %lu", count);
			NSLog(@"%10.3g %10.3g %10.3g %10.3g", a[0], a[1], a[2], affinity);
			NSLog(@"%10.3g %10.3g %10.3g", x[0], x[1], x[2]);
            double g[NA];
			for (i=0, j=0; i<NA; i++) g[i] = (x[i] != 0.0) ? gamma[j++] : 0.0;
            NSLog(@"%10.3g %10.3g %10.3g", g[0], g[1], g[2]);
		}

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
