//
//  FeldsparBerman.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/30/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "FeldsparBerman.h"
#import "BermanProperties.h"
#import "BermanStoichiometricPhases.h"
#import "BermanAlbite.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"
#import "IntegerVector.h"
#import "MathSupport.h"


@implementation FeldsparBerman

static NSCountedSet *instanceSet;

#ifdef DEBUG
#undef DEBUG
#endif

#define NR         2    // Two independent composition variables
#define NA         3    // Three endmember compositions
#define NATOMS  13.0    // Number of atoms in the formula unit

@synthesize whabor;
@synthesize wsabor;
@synthesize wvabor;
@synthesize whorab;
@synthesize wsorab;
@synthesize wvorab;
@synthesize whaban;
@synthesize whanab;
@synthesize whoran;
@synthesize whanor;
@synthesize wvanor;
@synthesize whabanor;
@synthesize wvabanor;

@synthesize endmembers;

#pragma mark -
#pragma mark class methods

/*
 Note that in this class the endmembers array is not a global static but is defined as an instance variable.
 This is because calculations for albite have to be made with instance variables (ordering and teh like), and
 therefore the BermanAlbite class must be instantiated for each class instance.  The endmember instance
 variable is created in the init() method.
 */
+(void)initialize {
	if (self == [FeldsparBerman class]) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"Initialize(FeldsparBerman) - entry ...");
		instanceSet = [[NSCountedSet alloc] initWithCapacity:1];
	}
}

#pragma mark -
#pragma mark instance methods

@synthesize operationParent;

-(id)initWithCompositionConstraint:(NSString *)name {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(FeldsparBerman) with %@ composition ... entry ...", name);

		[self setPhaseName:[NSString stringWithString:name]];
		computeMixingQuantities = NO;
		operationParent = @"";

		isMELTS = YES;
		if ([[NSUserDefaults standardUserDefaults] integerForKey:@"PPCalculationDatabase"] != 1) isMELTS = NO;

		NSMutableArray *mutableEndmembers = [NSMutableArray arrayWithCapacity:NA];
		[mutableEndmembers addObject:[[BermanAlbite alloc] init]];
		if (debug) NSLog(@"... allocated albite ...");

		BermanProperties *anorthite = nil;
		if (isMELTS) {
			if (debug) NSLog(@"FeldsparBerman is configured for MELTS.");
			anorthite = [[BermanProperties alloc] initWithH:-4228730.0+3.7*4184.0
														  S:200.186+3.7*4184.0/2200.0
														 k0:439.37
														 k1:-37.341E2
														 k2:0.0
														 k3:-31.702E7
														 v0:10.075
														 v1:-1.272E-6
														 v2:3.176E-12
														 v3:10.918E-6
														 v4:41.985E-10];
		} else {
			if (debug) NSLog(@"FeldsparBerman is configured for pMELTS.");
			anorthite = [[BermanProperties alloc] initWithH:-4228730.0
														  S:200.186
														 k0:439.37
														 k1:-37.341E2
														 k2:0.0
														 k3:-31.702E7
														 v0:10.075
														 v1:-1.272E-6
														 v2:3.176E-12
														 v3:10.918E-6
														 v4:41.985E-10];
		}
		[anorthite setPhaseFormula:@"CaAl2Si2O8"];
		[anorthite setPhaseName:@"anorthite"];
		[mutableEndmembers addObject:anorthite];
		if (debug) NSLog(@"... allocated anorthite ...");

		[mutableEndmembers addObject:[[SanidineBerman alloc] init]];
		if (debug) NSLog(@"... allocated sanidine ...");

		endmembers = [NSArray arrayWithArray:mutableEndmembers];

        /*
         *=============================================================================
         * Feldspar solution parameters:
         * Elkins, Linda T., Grove, Timothy L.
         * Ternary feldspar experiments and thermodynamic models
         * American Mineralogist 75, 544-559
         */
        whabor   = 18810.0;  /* joules     */
        wsabor   = 10.3;     /* joules/K   */
        wvabor   = 0.4602;   /* joules/bar */
        whorab   = 27320.0;  /* joules     */
        wsorab   = 10.3;     /* joules/K   */
        wvorab   = 0.3264;   /* joules/bar */
        whaban   = 7924.0;   /* joules     */
        whanab   = 0.0;      /* joules     */
        whoran   = 40317.0;  /* joules     */
        whanor   = 38974.0;  /* joules     */
        wvanor   = -0.1037;  /* joules/bar */
        whabanor = 12545.0;  /* joules     */
        wvabanor = -1.095;   /* joules/bar */

		if (debug) NSLog(@"... exiting.");
	}
	return self;
}

- (id)init {
    return [self initWithCompositionConstraint:@"Feldspar"];
}

static NSString *kwhabor = @"whabor";
static NSString *kwsabor = @"wsabor";
static NSString *kwvabor = @"wvabor";
static NSString *kwhorab = @"whorab";
static NSString *kwsorab = @"wsorab";
static NSString *kwvorab = @"wvorab";
static NSString *kwhaban = @"whaban";
static NSString *kwhanab = @"whanab";
static NSString *kwhoran = @"whoran";
static NSString *kwhanor = @"whanor";
static NSString *kwvanor = @"wvanor";
static NSString *kwhabanor = @"whabanor";
static NSString *kwvabanor = @"wvabanor";
static NSString *kendmembers = @"endmembers";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        whabor = [aDecoder decodeDoubleForKey:kwhabor];
        wsabor = [aDecoder decodeDoubleForKey:kwsabor];
        wvabor = [aDecoder decodeDoubleForKey:kwvabor];
        whorab = [aDecoder decodeDoubleForKey:kwhorab];
        wsorab = [aDecoder decodeDoubleForKey:kwsorab];
        wvorab = [aDecoder decodeDoubleForKey:kwvorab];
        whaban = [aDecoder decodeDoubleForKey:kwhaban];
        whanab = [aDecoder decodeDoubleForKey:kwhanab];
        whoran = [aDecoder decodeDoubleForKey:kwhoran];
        whanor = [aDecoder decodeDoubleForKey:kwhanor];
        wvanor = [aDecoder decodeDoubleForKey:kwvanor];
        whabanor = [aDecoder decodeDoubleForKey:kwhabanor];
        wvabanor = [aDecoder decodeDoubleForKey:kwvabanor];
#ifdef __APPLE__
        endmembers = (NSArray *) [aDecoder decodeObjectOfClass:[NSArray class] forKey:kendmembers];
#else
        endmembers = (NSArray *) [aDecoder decodeObjectForKey:kendmembers];
#endif

        computeMixingQuantities = NO;
        operationParent = @"";

        isMELTS = YES;
        if ([[NSUserDefaults standardUserDefaults] integerForKey:@"PPCalculationDatabase"] != 1) isMELTS = NO;
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:whabor forKey:kwhabor];
        [aCoder encodeDouble:wsabor forKey:kwsabor];
        [aCoder encodeDouble:wvabor forKey:kwvabor];
        [aCoder encodeDouble:whorab forKey:kwhorab];
        [aCoder encodeDouble:wsorab forKey:kwsorab];
        [aCoder encodeDouble:wvorab forKey:kwvorab];
        [aCoder encodeDouble:whaban forKey:kwhaban];
        [aCoder encodeDouble:whanab forKey:kwhanab];
        [aCoder encodeDouble:whoran forKey:kwhoran];
        [aCoder encodeDouble:whanor forKey:kwhanor];
        [aCoder encodeDouble:wvanor forKey:kwvanor];
        [aCoder encodeDouble:whabanor forKey:kwhabanor];
        [aCoder encodeDouble:wvabanor forKey:kwvabanor];
        [aCoder encodeObject:endmembers forKey:kendmembers];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", [self className]];
}

-(void)incrementInstanceCountOfPhase {
	[instanceSet addObject:[[self phaseName] stringByAppendingString:[self operationParent]]];
}

-(void)decrementInstanceCountOfPhase {
	[instanceSet removeObject:[[self phaseName] stringByAppendingString:[self operationParent]]];
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

#define FR0(i)     (i == 0) ? 1.0 - xab : - xab
#define FR1(i)     (i == 1) ? 1.0 - xan : - xan
#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0

#define LOG(x) ((x > 0.0) ? log(x) : log(DBL_EPSILON))

/*
 * Global (to this file): derivative definitions
 */
#define R       8.3143
#define S       - R*(xab*LOG(xab) + xan*LOG(xan) + xor*LOG(xor)) + wsabor*xab*xor*(xor+xan/2.0) + wsorab*xab*xor*(xab+xan/2.0)
#define H       whaban*xab*xan*(xan+xor/2.0) + whanab*xab*xan*(xab+xor/2.0) + whabor*xab*xor*(xor+xan/2.0) + whorab*xab*xor*(xab+xan/2.0) + \
                whanor*xan*xor*(xor+xab/2.0) + whoran*xan*xor*(xan+xab/2.0) + whabanor*xab*xan*xor
#define V       wvabor*xab*xor*(xor+xan/2.0) + wvorab*xab*xor*(xab+xan/2.0) + wvanor*xan*xor*(xor+xab/2.0) + wvabanor*xab*xan*xor
#define G       H - t*(S) + (p-1.0)*(V)

#define DGDR0   R*t*(LOG(xab) - LOG(xor)) + whaban*(xan*(xan+xor/2.0) - 0.5*xab*xan) + whanab*(xan*(xab+xor/2.0) + 0.5*xab*xan) + \
                (whabor-t*wsabor+(p-1.0)*wvabor)* ((xor-xab)*(xor+xan/2.0) - xab*xor) + (whorab-t*wsorab+(p-1.0)*wvorab)* \
                ((xor-xab)*(xab+xan/2.0) + xab*xor) + (whanor+(p-1.0)*wvanor)*(- xan*(xor+xab/2.0) - 0.5*xan*xor) + \
                whoran*(- xan*(xan+xab/2.0) + 0.5*xan*xor) + (whabanor+(p-1.0)*wvabanor)*xan*(xor - xab)
#define DGDR1   R*t*(LOG(xan) - LOG(xor)) + whaban*(xab*(xan+xor/2.0) + 0.5*xab*xan) + whanab*(xab*(xab+xor/2.0) - 0.5*xab*xan) + \
                (whabor-t*wsabor+(p-1.0)*wvabor)* (- xab*(xor+xan/2.0) - 0.5*xab*xor) + (whorab-t*wsorab+(p-1.0)*wvorab)* \
                (- xab*(xab+xan/2.0) + 0.5*xab*xor) + (whanor+(p-1.0)*wvanor)*((xor-xan)*(xor+xab/2.0) - xan*xor) + \
                whoran*((xor-xan)*(xan+xab/2.0) + xan*xor) + (whabanor+(p-1.0)*wvabanor)*xab*(xor - xan)
#define DGDP    (V)

#define D2GDR0R0  R*t*(1.0/xab + 1.0/xor) + (whanab - whaban)*xan - 2.0*(whabor-t*wsabor+(p-1.0)*wvabor)* ((xor+xan/2.0) + (xor-xab)) \
                  - 2.0*(whorab-t*wsorab+(p-1.0)*wvorab)* ((xab+xan/2.0) - (xor-xab)) + (whanor+(p-1.0)*wvanor - whoran)*xan + \
                  -2.0*(whabanor+(p-1.0)*wvabanor)*xan
#define D2GDR0R1  R*t/xor + whaban*((xan+xor/2.0) + 0.5*xan - 0.5*xab) + whanab*((xab+xor/2.0) - 0.5*xan + 0.5*xab) + \
                  0.5*(whabor-t*wsabor+(p-1.0)*wvabor)*(3.0*xab-xan-3.0*xor) + 0.5*(whorab-t*wsorab+(p-1.0)*wvorab)*(xor - 5.0*xab - xan) + \
                  0.5*(whanor+(p-1.0)*wvanor)*(3.0*xan - 3.0*xor - xab) + 0.5*whoran*(xor - 5.0*xan - xab) + (whabanor+(p-1.0)*wvabanor)*(xor - xab - xan)
#define D2GDR0DT  R*(LOG(xab) - LOG(xor)) - wsabor*((xor-xab)*(xor+xan/2.0) - xab*xor) - wsorab*((xor-xab)*(xab+xan/2.0) + xab*xor)
#define D2GDR0DP  wvabor*((xor-xab)*(xor+xan/2.0) - xab*xor) + wvorab*((xor-xab)*(xab+xan/2.0) + xab*xor) + wvanor*(- xan*(xor+xab/2.0) - 0.5*xan*xor) + \
                  wvabanor*xan*(xor - xab)

#define D2GDR1R1  R*t*(1.0/xan + 1.0/xor) + xab*(whaban - whanab) + (whabor-t*wsabor+(p-1.0)*wvabor)*xab - (whorab-t*wsorab+(p-1.0)*wvorab)*xab + \
                  (whanor+(p-1.0)*wvanor)*(2.0*xan - 4.0*xor - xab) + whoran*(2.0*xor - 4.0*xan - xab) - 2.0*(whabanor+(p-1.0)*wvabanor)*xab
#define D2GDR1DT  R*(LOG(xan) - LOG(xor)) - wsabor*(- xab*(xor+xan/2.0) - 0.5*xab*xor) - wsorab*(- xab*(xab+xan/2.0) + 0.5*xab*xor)
#define D2GDR1DP  wvabor*(- xab*(xor+xan/2.0) - 0.5*xab*xor) + wvorab*(- xab*(xab+xan/2.0) + 0.5*xab*xor) + wvanor*((xor-xan)*(xor+xab/2.0) - xan*xor) + \
                  wvabanor*xab*(xor - xan)
#define D2GDT2    0.0
#define D2GDTDP   0.0
#define D2GDP2    0.0

#define D3GDR0R0R0 R*t*(1.0/SQUARE(xor)-1.0/SQUARE(xab)) + 6.0* (whabor-t*wsabor+(p-1)*wvabor) - 6.0*(whorab - t*wsorab+(p-1)*wvorab)
#define D3GDR0R0R1 R*t*(1.0/SQUARE(xor)) + (whanab-whaban) + 3.0*(whabor-t*wsabor+(p-1)*wvabor) - 3.0*(whorab - t*wsorab+(p-1)*wvorab) + \
                   (whanor+(p-1)*wvanor - whoran) - 2*(whabanor + (p-1)*wvabanor)
#define D3GDR0R0DT R*(1.0/xab + 1.0/xor) + 2.0*wsabor*((xor+xan/2.0) + (xor-xab)) + 2.0*wsorab*((xab+xan/2.0) - (xor-xab))
#define D3GDR0R0DP - 2.0*wvabor*((xor+xan/2.0) + (xor-xab)) - 2.0*wvorab*((xab+xan/2.0) - (xor-xab)) + wvanor*xan - 2.0*wvabanor*xan
#define D3GDR0R1DT R/xor - 0.5*wsabor*(3.0*xab - xan - 3.0*xor) - 0.5*wsorab*(xor - 5.0*xab - xan)
#define D3GDR0R1DP 0.5*wvabor*(3.0*xab - xan - 3.0*xor) + 0.5*wvorab*(xor - 5.0*xab - xan) + 0.5*wvanor*(3.0*xan - 3.0*xor - xab) + wvabanor*(xor - xab - xan)
#define D3GDR1R1R1 R*t*(1.0/SQUARE(xor)-1.0/SQUARE(xan)) + 6.0* (whanor+(p-1.0)*wvanor) - 6.0*whoran
#define D3GDR1R1R0 R*t*(1.0/SQUARE(xor)) + (whaban-whanab) + (whabor-t*wsabor+(p-1)*wvabor) - (whorab - t*wsorab+(p-1)*wvorab) + \
                   3.0*(whanor+(p-1.0)*wvanor) - 3.0*whoran - 2*(whabanor + (p-1)*wvabanor)
#define D3GDR1R1DT R*(1.0/xan + 1.0/xor) + (wsorab-wsabor)*xab
#define D3GDR1R1DP (wvabor - wvorab)*xab + wvanor*(2.0*xan - 4.0*xor - xab) - 2.0*wvabanor*xab
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
	const char *NAMES[NA]    = { "albite", "anorthite", "sanidine" };
	const char *FORMULAS[NA] = { "NaAlSi3O8", "CaAl2Si2O8", "KAlSi3O8" };
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
	 endmember feldspar components.
	 (2) calculates from a vector of moles of endmember components, one or
	 all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
	 (3) calculates from a vector of independent compositional variables
	 mole fractions of endmember components and/or the Jacobian matrix
	 dx[]/dr[]

	 In this routine it is assumed that the elements are in the order of atomic
	 numbers and that the order of feldspar components has been verified as:
	 m[0] = albite    (NaAlSi3O8) ,
	 m[1] = anorthite (CaAl2Si2O8) and
	 m[2] = sanidine  (KAlSi3O8)

	 ----------------------------------------------------------------------------*/

	int i, j, k;

	if (inpMask == FIRST && outMask == SECOND) {
		static const int Na = 11;
		static const int K  = 19;
		static const int Ca = 20;

		/* Projection into the Na, Ca, K  triangle */
		m[0] = e[Na]; /* moles of NaAlSi3O8                      */
		m[1] = e[Ca]; /* Moles of CaAl2Si2O8                     */
		m[2] = e[K];  /* Moles of KAlSi3O8                       */

	} else if (inpMask == SECOND) {
		double sum;

		if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
			NSLog(@"Illegal call to conFld with inpMask = %o and outMask = %o\n", inpMask, outMask);

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
			/* calculates the matrix d3r[i]/dm[j]dm[k]dm[l] using m[] as input        */
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
			NSLog(@"Illegal call to conFld with inpMask = %o and outMask = %o\n", inpMask, outMask);

		if (outMask & FOURTH) {
			/* Converts a vector of independent compositional variables (r)
			 into a vector of mole fractions of endmember components (x).         */

			for (i=0, x[2]=1.0; i<NR; i++) { x[i] = r[i]; x[2] -= r[i]; }
			if (fabs(x[2]) < sqrt(DBL_EPSILON)) x[2] = 0.0;
		}

		if (outMask & SEVENTH) {
			/* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
			for (i=0; i<NR; i++) for (j=0; j<NR; j++) dr[i][j] = (i == j) ? 1.0 : 0.0;
			for (j=0; j<NR; j++) dr[2][j] = -1.0;
		}

	} else  {
		NSLog(@"Illegal call to conFld with inpMask = %o and outMask = %o\n", inpMask, outMask);
	}

}

-(NSString *)displayFormula:(double)t
						  p:(double)p
						  r:(double [NA])r
{
	double totAl, totCa, totNa, totK, totSi;

	totK   = 1.0 - r[0] - r[1];
    totNa  = r[0];
    totCa  = r[1];
    totAl  = 1.0 + r[1];
    totSi  = 3.0 - r[1];

	return [NSString stringWithFormat:@"K%4.2fNa%4.2fCa%4.2fAl%4.2fSi%4.2fO8", totK, totNa, totCa, totAl, totSi];
}

-(void)activity:(int)mask
			  t:(double)t
			  p:(double)p
			  r:(double [NA])r
			  a:(double [NA])a      // (pointer to a[]) activities              BINARY MASK: 0001
			 mu:(double [NA])mu     // (pointer to mu[]) chemical potentials    BINARY MASK: 0010
			 dx:(double [NA][NR])dx // (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100
{
	double xab = (r[0]          >  DBL_EPSILON) ? r[0]          : DBL_EPSILON;
	double xan = (r[1]          >  DBL_EPSILON) ? r[1]          : DBL_EPSILON;
	double xor = (1.0-r[0]-r[1] >  DBL_EPSILON) ? 1.0-r[0]-r[1] : DBL_EPSILON;

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
	double xab = (r[0]          >  DBL_EPSILON) ? r[0]          : DBL_EPSILON;
	double xan = (r[1]          >  DBL_EPSILON) ? r[1]          : DBL_EPSILON;
	double xor = (1.0-r[0]-r[1] >  DBL_EPSILON) ? 1.0-r[0]-r[1] : DBL_EPSILON;

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

	if (mask & FOURTH) {
		double d3gdr3[NR][NR][NR];
		int i, j, k;

		d3gdr3[0][0][0] = D3GDR0R0R0;
		d3gdr3[0][0][1] = D3GDR0R0R1;
		d3gdr3[0][1][0] = d3gdr3[0][0][1];
		d3gdr3[0][1][1] = D3GDR1R1R0;
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
	double xab = (r[0]          >  DBL_EPSILON) ? r[0]          : DBL_EPSILON;
	double xan = (r[1]          >  DBL_EPSILON) ? r[1]          : DBL_EPSILON;
	double xor = (1.0-r[0]-r[1] >  DBL_EPSILON) ? 1.0-r[0]-r[1] : DBL_EPSILON;

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
	double xab = (r[0]          >  DBL_EPSILON) ? r[0]          : DBL_EPSILON;
	double xan = (r[1]          >  DBL_EPSILON) ? r[1]          : DBL_EPSILON;
	double xor = (1.0-r[0]-r[1] >  DBL_EPSILON) ? 1.0-r[0]-r[1] : DBL_EPSILON;

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
	double xab = (r[0]          >  DBL_EPSILON) ? r[0]          : DBL_EPSILON;
	double xan = (r[1]          >  DBL_EPSILON) ? r[1]          : DBL_EPSILON;
	double xor = (1.0-r[0]-r[1] >  DBL_EPSILON) ? 1.0-r[0]-r[1] : DBL_EPSILON;

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
			return @"albite";
			break;
		case 1:
			return @"anorthite";
			break;
		case 2:
			return @"sanidine";
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
			for (i=0; i<nz; i++) if (xNz[i] <= DBL_EPSILON) xNz[i] = DBL_EPSILON;

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

	double KOverTotalMoles = DBL_MAX;
	if ((x[0]+x[1]+x[2]) != 0.0) KOverTotalMoles =  x[2]/(x[0]+x[1]+x[2]);

	if      ((KOverTotalMoles > 0.15) && [[self phaseName] isEqualToString:@"Plagioclase"]) affinity = 0.0;
	else if ((KOverTotalMoles < 0.15) && [[self phaseName] isEqualToString:@"Sanidine"])    affinity = 0.0;

	[results addObject:[NSNumber numberWithDouble:affinity]];                    // affinity in J
	for (i=0; i<NA; i++) [results addObject:[NSNumber numberWithDouble:x[i]]];   // composition in mole fraction of endmembers
	[results addObject:[NSNumber numberWithBool:converged]];                     // convergence flag
	[results addObject:[NSNumber numberWithUnsignedInteger:count]];                  // iteration count
	[results addObject:[NSNumber numberWithDouble:NATOMS]];                      // number of atoms used to scale affinity
	[results addObject:[NSNumber numberWithDouble:fabs(affinity-affinityLast)]]; // likely error in affinity

	if (debugV) NSLog(@"Exiting [... affinityAndCompositionFromLiquidChemicalPotentialSum].");
	return [NSArray arrayWithArray:results];
}

// --> SolutionPhaseProtocol public function

-(NSDictionary *)checkForAndDetermineCompositionOfCoexistingImmisciblePhase:(double *)refMoles andT:(double)t andP:(double)p {
	BOOL debugS = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
	BOOL debugV = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
	if (debugV) NSLog(@"Entering [... checkForAndDetermineCompositionOfCoexistingImmisciblePhase] ...");

	NSMutableDictionary *results = [NSMutableDictionary dictionaryWithCapacity:6];

	if ((refMoles[1] == 0.0) || (refMoles[2] == 0.0) ||
		([[self phaseName] isEqualToString:@"Feldspar"]) ||
		(([[self phaseName] isEqualToString:@"Sanidine"]) &&
		 ([instanceSet countForObject:[@"Plagioclase" stringByAppendingString:[self operationParent]]] > 0)) ||
		(([[self phaseName] isEqualToString:@"Plagioclase"]) &&
		 ([instanceSet countForObject:[@"Sanidine" stringByAppendingString:[self operationParent]]] > 0))
		) {

		NSMutableArray *composition = [NSMutableArray arrayWithCapacity:NA];
		for (NSUInteger i=0; i<NA; i++) [composition addObject:[NSNumber numberWithDouble:refMoles[i]]];
		[results setObject:[NSNumber numberWithBool:NO] forKey:@"additionalPhaseDetected"];
		[results setObject:@"" forKey:@"nameOfCoexistingPhase"];
		[results setObject:[NSNumber numberWithBool:NO] forKey:@"converged"];
		[results setObject:[NSNumber numberWithUnsignedInteger:0] forKey:@"iterations"];
		[results setObject:[NSNumber numberWithDouble:0.0] forKey:@"affinity"];
		[results setObject:[NSArray arrayWithArray:composition] forKey:@"composition"];

		return [NSDictionary dictionaryWithDictionary:results];
	}

	double *referenceActivities = [[self getActivityFromMolesOfComponents:refMoles andT:t andP:p] pointerToDouble];

	double moles[NA];
	double refTotalMoles = 0.0;
	for (NSUInteger i=0; i<NA; i++) refTotalMoles += refMoles[i];
	if (refTotalMoles <= 0.0) refTotalMoles = 1.0;

	// Specific to feldspar ...
	moles[0] = refMoles[0]/refTotalMoles; // Ab
	moles[1] = refMoles[2]/refTotalMoles; // switch in Sn for An
	moles[2] = refMoles[1]/refTotalMoles; // switch in An for Sn

	if (debugS) {
		NSLog(@"... d = %d, r = %d, rNorm = %10.3e Ab = %10.3e An = %10.3e Sn = %10.3e Affinity/RT = %10.3e Ab = %10.3e An = %10.3e Sn = %10.3e",
			  NA, NA, 0.0, 0.0, 0.0, 0.0, 0.0, refMoles[0]/refTotalMoles, refMoles[1]/refTotalMoles, refMoles[2]/refTotalMoles);
		NSLog(@"... d = %d, r = %d, rNorm = %10.3e Ab = %10.3e An = %10.3e Sn = %10.3e Affinity/RT = %10.3e Ab = %10.3e An = %10.3e Sn = %10.3e",
			  NA, NA, 0.0, 0.0, 0.0, 0.0, 0.0, moles[0], moles[1], moles[2]);
	}
	// ... end

	NSUInteger iters = 0;
	BOOL converged = NO;
	double tolerance = sqrt(DBL_EPSILON);
	NSUInteger maxIters = 50;
	double affinityScaledByRT = 0.0;
	double totalMoles = 1.0;

	DoubleMatrix *dMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:NA+1 andWithColumnSize:1 andInitialValue:0.0];
	double **dMatrix = [dMatrixWrapper pointerToPointerToDouble];
    DoubleMatrix *cMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:NA+1 andWithColumnSize:NA+1 andInitialValue:0.0];
	double **cMatrix = [cMatrixWrapper pointerToPointerToDouble];
	DoubleVector *hVectorWrapper = [[DoubleVector alloc] initWithSize:NA+1];
	double *hVector  = [hVectorWrapper pointerToDouble];
    DoubleVector *gVectorWrapper = [[DoubleVector alloc] initWithSize:NA+1];
	double *gVector  = [gVectorWrapper pointerToDouble];
    IntegerVector *pVectorWrapper = [[IntegerVector alloc] initWithSize:NA+1 andInitialValue:0];
	NSInteger *pVector     = [pVectorWrapper pointerToInteger];

	while (!converged && (iters < maxIters)) {
        DoubleVector *activitiesWrapper = [self getActivityFromMolesOfComponents:moles andT:t andP:p];
		double *activities = [activitiesWrapper pointerToDouble];
        DoubleMatrix *dadmWrapper = [self getDaDmFromMolesOfComponents:moles andT:t andP:p];
		double **dadm = [dadmWrapper pointerToPointerToDouble];
		double rNorm = 0.0;
		double expAffinityScaledByRT = exp(affinityScaledByRT/totalMoles);

		for (NSUInteger i=0; i<NA; i++) {
			if ((refMoles[i] != 0.0) && (moles[i] != 0.0)) {
				dMatrix[i][0] = 1.0 - activities[i]*expAffinityScaledByRT/referenceActivities[i];
				for (NSUInteger j=0; j<NA; j++) cMatrix[i][j] = dadm[i][j]*expAffinityScaledByRT/referenceActivities[i]
					- (activities[i]/referenceActivities[i])*expAffinityScaledByRT*affinityScaledByRT/pow(totalMoles, 2.0);
				cMatrix[i][NA] = (activities[i]/referenceActivities[i])*expAffinityScaledByRT/totalMoles;
			} else {
				dMatrix[i][0] = 0.0;
				for (NSUInteger j=0; j<NA; j++) cMatrix[i][j] = 0.0;
				cMatrix[i][NA] = 0.0;
			}
			rNorm += pow(dMatrix[i][0], 2.0);
			cMatrix[NA][i] = 1.0;
		}
		rNorm = sqrt(rNorm);
		cMatrix[NA][NA] = 0.0;
		dMatrix[NA][0] = 0.0;

		if (rNorm < tolerance) converged = YES;
		if (converged || (iters > maxIters)) break;

		double HFTItolerance = 10.0*DBL_EPSILON, HFTIrNorm;
		NSInteger pseudorank;
		[MathSupport hfti:cMatrix
						m:NA+1
						n:NA+1
						b:dMatrix
					   nb:1
					  tau:HFTItolerance
						k:&pseudorank
					rnorm:&HFTIrNorm
						h:hVector
						g:gVector
						p:pVector];

		BOOL compositionOK = NO;
		double tempMoles[NA], stepLength = 1.0;
		while (!compositionOK) {
			for (NSUInteger i=0; i<NA; i++) tempMoles[i] = moles[i] + stepLength*dMatrix[i][0];
			compositionOK = [self testPermissibleValuesOfComponents:tempMoles];
			if (!compositionOK) stepLength /= 2.0;
            if (stepLength < DBL_EPSILON) {
                converged = NO;
                compositionOK = YES;
                iters = maxIters;
                if (debugS) {
                    NSLog(@"... iter = %lu, computation aborted because steplength (%10.3e) is < DBL_EPSILON", iters, stepLength);
                }
            }
		}
		totalMoles = 0.0;
		for (NSUInteger i=0; i<NA; i++) {
			moles[i] = tempMoles[i];
			totalMoles += moles[i];
		}
		affinityScaledByRT += dMatrix[NA][0];

		// Specific to feldspar ...
		if (debugS) {
			NSLog(@"... d = %d, r = %ld, rNorm = %10.3e Ab = %10.3e An = %10.3e Sn = %10.3e Affinity/RT = %10.3e Ab = %10.3e An = %10.3e Sn = %10.3e sL = %10.3e",
				  NA, pseudorank, rNorm, dMatrix[0][0], dMatrix[1][0], dMatrix[2][0], affinityScaledByRT, moles[0], moles[1], moles[2], stepLength);
		}
		// ... end

		iters++;
	}

	if (debugS) {
		NSLog(@"... Soln: Reference composition is a %@", [self phaseName]);
		NSLog(@"... There are %lu instances of Plagioclase instantiated for this class",
			  [instanceSet countForObject:[@"Plagioclase" stringByAppendingString:[self operationParent]]]);
		NSLog(@"... There are %lu instances of Sanidine    instantiated for this class",
			  [instanceSet countForObject:[@"Sanidine" stringByAppendingString:[
                                                                                                            self operationParent]]]);
		NSLog(@"... There are %lu instances of Feldspar    instantiated for this class",
			  [instanceSet countForObject:[@"Feldspar" stringByAppendingString:[self operationParent]]]);
		NSLog(@"... Soln: affinity = %g", affinityScaledByRT*8.314472*t);
		for (NSUInteger i=0; i<NA; i++) NSLog(@"... Soln: %6.3f Ref: %6.3f delta: %13.6e", moles[i], refMoles[i]/refTotalMoles, moles[i]-refMoles[i]/refTotalMoles);
	}

	NSMutableArray *composition = [NSMutableArray arrayWithCapacity:NA];
	for (NSUInteger i=0; i<NA; i++) [composition addObject:[NSNumber numberWithDouble:moles[i]]];
	BOOL additionalPhaseDetected = (affinityScaledByRT > tolerance) ? YES : NO;
	if (!converged) additionalPhaseDetected = NO;

	NSString *nameOfCoexistingPhase = @"";
	if (additionalPhaseDetected) {
		if      ([[self phaseName] isEqualToString:@"Sanidine"])    nameOfCoexistingPhase = @"Plagioclase";
		else if ([[self phaseName] isEqualToString:@"Plagioclase"]) nameOfCoexistingPhase = @"Sanidine";
	}

	[results setObject:[NSNumber numberWithBool:additionalPhaseDetected] forKey:@"additionalPhaseDetected"];
	[results setObject:nameOfCoexistingPhase forKey:@"nameOfCoexistingPhase"];
	[results setObject:[NSNumber numberWithBool:converged] forKey:@"converged"];
	[results setObject:[NSNumber numberWithUnsignedInteger:iters] forKey:@"iterations"];
	[results setObject:[NSNumber numberWithDouble:affinityScaledByRT*8.314472*t] forKey:@"affinity"];
	[results setObject:[NSArray arrayWithArray:composition] forKey:@"composition"];

	return [NSDictionary dictionaryWithDictionary:results];
}

-(NSString *)nameOfPhaseWithComposition:(double *)refMoles {
	// order NaAlSi3O8, CaAl2Si2O8, KAlSi3O8
	double totNa = refMoles[0];
	double totCa = refMoles[1];
	double totK  = refMoles[2];

	if ((totNa+totK+totCa) == 0.0) return [NSString stringWithString:[self phaseName]];
	if ((totK/(totNa+totK+totCa)) < 0.15) return @"Plagioclase";
	else                                  return @"Sanidine";
}

#import "SolutionPhase.h"

@end
