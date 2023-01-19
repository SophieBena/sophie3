//
//  ClinoamphiboleBerman.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/8/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import "ClinoamphiboleBerman.h"
#import "BermanProperties.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "IntegerVector.h"
#import "DoubleTensor.h"
#import "MathSupport.h"

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

@implementation ClinoamphiboleBerman

static const int clino = TRUE;
static NSArray *endmembers;
static NSCountedSet *instanceSet;

#define NR        2  /* Six independent composition variables       */
#define NS        2  /* Two ordering parameters                     */
#define NA        3  /* Seven endmember compositions                */
#define NATOMS 41.0  /* Average number of atoms in the formula unit */

#define MAX_ITER 100

#pragma mark -
#pragma mark class methods

+(void)initialize {
	if (self == [ClinoamphiboleBerman class]) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"Initialize(ClinoamphiboleBerman) - entry ...");
		NSMutableArray *mutableEndmembers = [NSMutableArray arrayWithCapacity:NA];

		BermanProperties *cummingtonite = [[BermanProperties alloc] initWithH:-12067920.38
																			S:540.2587
																		   k0:1233.79200
																		   k1:-71.33980e2
																		   k2:-221.63800e5
																		   k3:233.393749e7
																		   v0:26.34
																		   v1:-1.1393847e-6
																		   v2:0.0
																		   v3:28.1048234e-6
																		   v4:62.894e-10];
		[cummingtonite setPhaseFormula:@"Mg7Si8O22(OH)2"];
		[cummingtonite setPhaseName:@"cummingtonite"];
		[mutableEndmembers addObject:cummingtonite];
		if (debug) NSLog(@"... allocated cummingtonite ...");

		BermanProperties *grunerite = [[BermanProperties alloc] initWithH:-9623550.0
																		S:725.0
																	   k0:1347.83
																	   k1:-93.5691e2
																	   k2:-202.2848e5
																	   k3:303.919e7
																	   v0:27.840
																	   v1:-1.670299e-6
																	   v2:8.68919e-12
																	   v3:28.400e-6
																	   v4:0.0];
		[grunerite setPhaseFormula:@"Fe7Si8O22(OH)2"];
		[grunerite setPhaseName:@"grunerite"];
		[mutableEndmembers addObject:grunerite];
		if (debug) NSLog(@"... allocated grunerite ...");

		BermanProperties *tremolite = [[BermanProperties alloc] initWithH:-12307863.0
																		S:549.500
																	   k0:1229.36
																	   k1:-64.019e2
																	   k2:-320.899e5
																	   k3:420.881e7
																	   v0:27.312
																	   v1:-1.4763e-6
																	   v2:8.9462e-12
																	   v3:24.374e-6
																	   v4:98.338e-10];
		[tremolite setPhaseFormula:@"Ca2Mg5Si8O22(OH)2"];
		[tremolite setPhaseName:@"tremolite"];
		[mutableEndmembers addObject:tremolite];
		if (debug) NSLog(@"... allocated tremolite ...");

		endmembers = [NSArray arrayWithArray:mutableEndmembers];

		instanceSet = [[NSCountedSet alloc] initWithCapacity:1];
	}
}


#pragma mark -
#pragma mark instance methods

@synthesize operationParent;

-(id)initWithCompositionConstraint:(NSString *)name {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(CpxBerman) with %@ composition ... entry ...", name);

		[self setPhaseName:[NSString stringWithString:name]];
		computeMixingQuantities = NO;
		operationParent = @"";

		tOld = -9999.0;
		pOld = -9999.0;
		for (NSUInteger i=0; i<NR; i++) rOld[i] = -9999.0;
		for (NSUInteger i=0; i<NS; i++) sOld[i] = 2.0;
		xfe2m13 = 0.0;
		xmg2m13 = 0.0;
		xfe2m2  = 0.0;
		xmg2m2  = 0.0;
		xfe2m4  = 0.0;
		xmg2m4  = 0.0;
		xca2m4  = 0.0;

		if (debug) NSLog(@"... exiting.");
	}
	return self;
}

-(id)init {
	return [self initWithCompositionConstraint:@"Clinoamphibole"];
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

/*
 *=============================================================================
 * Amphibole solution parameters:
 * Ghiorso, M.S. and B.W. Evans (1996)
 *   Work in progress.
 */

#define  R       8.3143

#define  DH1      (-12073132.08) - (-12067920.38) /* Ortho-Clino: Mg7Si8O22(OH)2     (joules) */
#define  DS1      (535.2587) - (540.2587)         /*                               (joules/K) */
#define  DV1      (26.310) - (26.340)             /*                             (joules/bar) */
#define  DH2      (-9627014.85) - (-9623550.0)    /* Ortho-Clino: Fe7Si8O22(OH)2     (joules) */
#define  DS2      (720.0) - (725.0)               /*                               (joules/K) */
#define  DV2      (27.810) - (27.840)             /*                             (joules/bar) */
#define  DH3      (-12307863.0) - (-12307863.0)   /* Ortho-Clino: Ca2Mg5Si8O22(OH)2  (joules) */
#define  DS3      (544.500) - (549.500)           /*                               (joules/K) */
#define  DV3      (27.312) - (27.312)             /*                             (joules/bar) */

#define  GEX        -25.8556    * 1000.0 /* joules */
#define  GDEX         5.90595   * 1000.0 /* joules */

#define  cGX123_4    12.1027    * 1000.0 /* joules */
#define  cGX13_24     3.29327   * 1000.0 /* joules */
#define  cGX134_2    15.1135    * 1000.0 /* joules */
#define  cW123        5.11564   * 1000.0 /* joules */
#define  cW4          3.55954   * 1000.0 /* joules */
#define  cDZ         -0.781732  * 1000.0 /* joules */
#define  cGXCA_13    17.801     * 1000.0 /* joules */
#define  cGXCA_2    -11.395     * 1000.0 /* joules */
#define  cWCAMG      28.100     * 1000.0 /* joules */
#define  cWCAFE      19.500     * 1000.0 /* joules */
#define  cDWCAMG      0.0000    * 1000.0 /* joules */
#define  cDWCAFE      0.0000    * 1000.0 /* joules */

#define  oGX123_4    13.7354    * 1000.0 /* joules */
#define  oGX13_24     5.25246   * 1000.0 /* joules */
#define  oGX134_2    16.7461    * 1000.0 /* joules */
#define  oW123        5.52380   * 1000.0 /* joules */
#define  oW4          3.72280   * 1000.0 /* joules */
#define  oDZ         -0.210304  * 1000.0 /* joules */
#define  oGXCA_13    17.801     * 1000.0 /* joules */
#define  oGXCA_2    -11.395     * 1000.0 /* joules */
#define  oWCAMG      30.100     * 1000.0 /* joules */
#define  oWCAFE      19.500     * 1000.0 /* joules */
#define  oDWCAMG      0.0000    * 1000.0 /* joules */
#define  oDWCAFE      0.0000    * 1000.0 /* joules */

/*
 *=============================================================================
 * Dependent parameters (DO NOT change definitions below this line:
 */

#define  DG1     (DH1)-t*(DS1)+(p-1.0)*DV1
#define  DG2     (DH2)-t*(DS2)+(p-1.0)*DV2
#define  DG3     (DH3)-t*(DS3)+(p-1.0)*DV3

#define  GX123_4  (clino) ? (cGX123_4) : (oGX123_4)
#define  GX13_24  (clino) ? (cGX13_24) : (oGX13_24)
#define  GX134_2  (clino) ? (cGX134_2) : (oGX134_2)
#define  W123     (clino) ? (cW123)    : (oW123)
#define  W4       (clino) ? (cW4)      : (oW4)
#define  DZ       (clino) ? (cDZ)      : (oDZ)
#define  GXCA_13  (clino) ? (cGXCA_13) : (oGXCA_13)
#define  GXCA_2   (clino) ? (cGXCA_2)  : (oGXCA_2)
#define  WCAMG    (clino) ? (cWCAMG)   : (oWCAMG)
#define  WCAFE    (clino) ? (cWCAFE)   : (oWCAFE)
#define  DWCAMG   (clino) ? (cDWCAMG)  : (oDWCAMG)
#define  DWCAFE   (clino) ? (cDWCAFE)  : (oDWCAFE)

/*
 * Vertices of composition space
 */

/* Mg7Si8O22(OH)2 */
#define  H1      (clino) ? 0.0 : (DH1)
#define  S1      (clino) ? 0.0 : (DS1)
#define  V1      (clino) ? 0.0 : (DV1)
#define  G1      (clino) ? 0.0 : (DH1)-t*(DS1)+(p-1.0)*(DV1)
/* Fe7Si8O22(OH)2 */
#define  H2      (clino) ? 0.0 : (DH2)
#define  S2      (clino) ? 0.0 : (DS2)
#define  V2      (clino) ? 0.0 : (DV2)
#define  G2      (clino) ? 0.0 : (DH2)-t*(DS2)+(p-1.0)*(DV2)
/* Ca2Mg5Si8O22(OH)2 */
#define  H3      (clino) ? 0.0 : (DH3)
#define  S3      (clino) ? 0.0 : (DS3)
#define  V3      (clino) ? 0.0 : (DV3)
#define  G3      (clino) ? 0.0 : (DH3)-t*(DS3)+(p-1.0)*(DV3)

/*
 * Definitions of Taylor expansion coefficients in terms of solution
 * parameters. Independent variables are r, p, s13, s2
 */

#define  H0         (H1)/2.0 + (H2)/2.0 + 5.0*(W123)/4.0 + (W4)/2.0 + (GX123_4)/4.0
#define  HR        -(H1)/2.0 + (H2)/2.0
#define  HS13       3.0*(GEX)/5.0 + (GDEX)
#define  HS2        2.0*(GEX)/5.0 - (GDEX)
#define  HRR       -5.0*(W123)/4.0 - (W4)/2.0 - (GX123_4)/4.0
#define  HRS13      5.0*(W123)/14.0 - 6.0*(W4)/7.0 - 5.0*(GX123_4)/28.0 + (DZ)/2.0
#define  HRS2       15.0*(W123)/14.0 - 4.0*(W4)/7.0 - (GX123_4)/28.0 - (DZ)/2.0
#define  HS13S13   -125.0*(W123)/98.0 - 18.0*(W4)/49.0 - 43.0*(GX123_4)/196.0 \
+ (GX13_24)/2.0 - (DZ)/14.0
#define  HS13S2     115.0*(W123)/49.0 - 24.0*(W4)/49.0 + 30.0*(GX123_4)/49.0 \
- (GX13_24)/2.0 - (GX134_2)/2.0 - (DZ)/7.0
#define  HS2S2     -145.0*(W123)/98.0 - 8.0*(W4)/49.0 - 37.0*(GX123_4)/196.0 \
+ (GX134_2)/2.0 + 3.0*(DZ)/14.0

#define  HP        -(H1) + (H3) - (W4) + (WCAFE) + (WCAMG) + (GXCA_13)/2.0 \
+ (GXCA_2) - (DWCAFE)/2.0 - (DWCAMG)/2.0
#define  HRP       -(W4) + (WCAFE) - (WCAMG) + (GXCA_13)/2.0 + (GXCA_2) \
- (DWCAFE) + (DWCAMG)
#define  HPS13     -6.0*(W4)/7.0 + 6.0*(WCAFE)/7.0 - 6.0*(WCAMG)/7.0 \
- 4.0*(GXCA_13)/7.0 - (GXCA_2)/7.0 - 6.0*(DWCAFE)/7.0 + 6.0*(DWCAMG)/7.0
#define  HPS2      -4.0*(W4)/7.0 + 4.0*(WCAFE)/7.0 - 4.0*(WCAMG)/7.0 \
+ 2.0*(GXCA_13)/7.0 - 3.0*(GXCA_2)/7.0 - 4.0*(DWCAFE)/7.0 + 4.0*(DWCAMG)/7.0
#define  HPP       -2.0*(WCAMG) + (DWCAFE) + 3.0*(DWCAMG)

#define  HPS13S13  -18.0*(DWCAFE)/49.0 - 18.0*(DWCAMG)/49.0
#define  HPPS13     6.0*(DWCAFE)/7.0 - 18.0*(DWCAMG)/7.0
#define  HPPS2      4.0*(DWCAFE)/7.0 - 12.0*(DWCAMG)/7.0
#define  HRPS13    -6.0*(DWCAFE)/7.0 - 6.0*(DWCAMG)/7.0
#define  HRPS2     -4.0*(DWCAFE)/7.0 - 4.0*(DWCAMG)/7.0
#define  HPS2S2    -8.0*(DWCAFE)/49.0 - 8.0*(DWCAMG)/49.0
#define  HPS13S2   -24.0*(DWCAFE)/49.0 - 24.0*(DWCAMG)/49.0
#define  HRRP      -(DWCAFE)/2.0 - (DWCAMG)/2.0
#define  HRPP       (DWCAFE) - 3.0*(DWCAMG)
#define  HPPP      -4.0*(DWCAMG)

#define  S0         (S1)/2.0 + (S2)/2.0
#define  SR        -(S1)/2.0 + (S2)/2.0
#define  SS13       0.0
#define  SS2        0.0
#define  SRR        0.0
#define  SRS13      0.0
#define  SRS2       0.0
#define  SS13S13    0.0
#define  SS13S2     0.0
#define  SS2S2      0.0

#define  SP        -(S1) + (S3)
#define  SRP        0.0
#define  SPS13      0.0
#define  SPS2       0.0
#define  SPP        0.0

#define  SPS13S13   0.0
#define  SPPS13     0.0
#define  SPPS2      0.0
#define  SRPS13     0.0
#define  SRPS2      0.0
#define  SPS2S2     0.0
#define  SPS13S2    0.0
#define  SRRP       0.0
#define  SRPP       0.0
#define  SPPP       0.0

#define  V0         (V1)/2.0 + (V2)/2.0
#define  VR        -(V1)/2.0 + (V2)/2.0
#define  VS13       0.0
#define  VS2        0.0
#define  VRR        0.0
#define  VRS13      0.0
#define  VRS2       0.0
#define  VS13S13    0.0
#define  VS13S2     0.0
#define  VS2S2      0.0

#define  VP        -(V1) + (V3)
#define  VRP        0.0
#define  VPS13      0.0
#define  VPS2       0.0
#define  VPP        0.0

#define  VPS13S13   0.0
#define  VPPS13     0.0
#define  VPPS2      0.0
#define  VRPS13     0.0
#define  VRPS2      0.0
#define  VPS2S2     0.0
#define  VPS13S2    0.0
#define  VRRP       0.0
#define  VRPP       0.0
#define  VPPP       0.0

#define  G0         (G1)/2.0 + (G2)/2.0 + 5.0*(W123)/4.0 + (W4)/2.0 + (GX123_4)/4.0
#define  GR        -(G1)/2.0 + (G2)/2.0
#define  GS13       3.0*(GEX)/5.0 + (GDEX)
#define  GS2        2.0*(GEX)/5.0 - (GDEX)
#define  GRR       -5.0*(W123)/4.0 - (W4)/2.0 - (GX123_4)/4.0
#define  GRS13      5.0*(W123)/14.0 - 6.0*(W4)/7.0 - 5.0*(GX123_4)/28.0 + (DZ)/2.0
#define  GRS2       15.0*(W123)/14.0 - 4.0*(W4)/7.0 - (GX123_4)/28.0 - (DZ)/2.0
#define  GS13S13   -125.0*(W123)/98.0 - 18.0*(W4)/49.0 - 43.0*(GX123_4)/196.0 \
+ (GX13_24)/2.0 - (DZ)/14.0
#define  GS13S2     115.0*(W123)/49.0 - 24.0*(W4)/49.0 + 30.0*(GX123_4)/49.0 \
- (GX13_24)/2.0 - (GX134_2)/2.0 - (DZ)/7.0
#define  GS2S2     -145.0*(W123)/98.0 - 8.0*(W4)/49.0 - 37.0*(GX123_4)/196.0 \
+ (GX134_2)/2.0 + 3.0*(DZ)/14.0

#define  GP        -(G1) + (G3) - (W4) + (WCAFE) + (WCAMG) + (GXCA_13)/2.0 \
+ (GXCA_2) - (DWCAFE)/2.0 - (DWCAMG)/2.0
#define  GRP       -(W4) + (WCAFE) - (WCAMG) + (GXCA_13)/2.0 + (GXCA_2) \
- (DWCAFE) + (DWCAMG)
#define  GPS13     -6.0*(W4)/7.0 + 6.0*(WCAFE)/7.0 - 6.0*(WCAMG)/7.0 \
- 4.0*(GXCA_13)/7.0 - (GXCA_2)/7.0 - 6.0*(DWCAFE)/7.0 + 6.0*(DWCAMG)/7.0
#define  GPS2      -4.0*(W4)/7.0 + 4.0*(WCAFE)/7.0 - 4.0*(WCAMG)/7.0 \
+ 2.0*(GXCA_13)/7.0 - 3.0*(GXCA_2)/7.0 - 4.0*(DWCAFE)/7.0 + 4.0*(DWCAMG)/7.0
#define  GPP       -2.0*(WCAMG) + (DWCAFE) + 3.0*(DWCAMG)

#define  GPS13S13  -18.0*(DWCAFE)/49.0 - 18.0*(DWCAMG)/49.0
#define  GPPS13     6.0*(DWCAFE)/7.0 - 18.0*(DWCAMG)/7.0
#define  GPPS2      4.0*(DWCAFE)/7.0 - 12.0*(DWCAMG)/7.0
#define  GRPS13    -6.0*(DWCAFE)/7.0 - 6.0*(DWCAMG)/7.0
#define  GRPS2     -4.0*(DWCAFE)/7.0 - 4.0*(DWCAMG)/7.0
#define  GPS2S2    -8.0*(DWCAFE)/49.0 - 8.0*(DWCAMG)/49.0
#define  GPS13S2   -24.0*(DWCAFE)/49.0 - 24.0*(DWCAMG)/49.0
#define  GRRP      -(DWCAFE)/2.0 - (DWCAMG)/2.0
#define  GRPP       (DWCAFE) - 3.0*(DWCAMG)
#define  GPPP      -4.0*(DWCAMG)

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conAph defines the conversion from m[i], to r[j]
 */
/* Order: r, p */
#define FR(i)     (i == 0 || i == 2) ? -(1.0+r[0]) : (1.0-r[0])
#define FP(i)     (i == 2)           ?   1.0-r[1]  :     -r[1]

/* Order: S13, S2 */
#define FS13(i)    - s[0]
#define FS2(i)     - s[1]

#define DFRDR(i)     - 1.0
#define DFPDP(i)     - 1.0

#define DFS13DS13(i) - 1.0
#define DFS2DS2(i)   - 1.0

/*
 * Global (to this file): derivative definitions
 */

#define SIC -R*(2.0*xfe2m4*log(xfe2m4) + 2.0*xmg2m4*log(xmg2m4) \
+ 2.0*xca2m4*log(xca2m4) + 3.0*xfe2m13*log(xfe2m13) \
+ 3.0*xmg2m13*log(xmg2m13) + 2.0*xfe2m2*log(xfe2m2) \
+ 2.0*xmg2m2*log(xmg2m2) )
#define S   (SIC) + \
(S0) + (SR)*r[0] + (SS13)*s[0] + (SS2)*s[1] + (SRR)*r[0]*r[0] \
+ (SRS13)*r[0]*s[0] + (SRS2)*r[0]*s[1] + (SS13S13)*s[0]*s[0] \
+ (SS13S2)*s[0]*s[1] + (SS2S2)*s[1]*s[1] + (SP)*r[1] + (SRP)*r[0]*r[1] \
+ (SPS13)*r[1]*s[0] + (SPS2)*r[1]*s[1] + (SPP)*r[1]*r[1] \
+ (SPS13S13)*r[1]*s[0]*s[0] + (SPPS13)*r[1]*r[1]*s[0] \
+ (SPPS2)*r[1]*r[1]*s[1] + (SRPS13)*r[0]*r[1]*s[0] + (SRPS2)*r[0]*r[1]*s[1] \
+ (SPS2S2)*r[1]*s[1]*s[1] + (SPS13S2)*r[1]*s[0]*s[1] + (SRRP)*r[0]*r[0]*r[1] \
+ (SRPP)*r[0]*r[1]*r[1] + (SPPP)*r[1]*r[1]*r[1]
#define H   (H0) + (HR)*r[0] + (HS13)*s[0] + (HS2)*s[1] + (HRR)*r[0]*r[0] \
+ (HRS13)*r[0]*s[0] + (HRS2)*r[0]*s[1] + (HS13S13)*s[0]*s[0] \
+ (HS13S2)*s[0]*s[1] + (HS2S2)*s[1]*s[1] + (HP)*r[1] + (HRP)*r[0]*r[1] \
+ (HPS13)*r[1]*s[0] + (HPS2)*r[1]*s[1] + (HPP)*r[1]*r[1] \
+ (HPS13S13)*r[1]*s[0]*s[0] + (HPPS13)*r[1]*r[1]*s[0] \
+ (HPPS2)*r[1]*r[1]*s[1] + (HRPS13)*r[0]*r[1]*s[0] + (HRPS2)*r[0]*r[1]*s[1] \
+ (HPS2S2)*r[1]*s[1]*s[1] + (HPS13S2)*r[1]*s[0]*s[1] + (HRRP)*r[0]*r[0]*r[1] \
+ (HRPP)*r[0]*r[1]*r[1] + (HPPP)*r[1]*r[1]*r[1]
#define V   (V0) + (VR)*r[0] + (VS13)*s[0] + (VS2)*s[1] + (VRR)*r[0]*r[0] \
+ (VRS13)*r[0]*s[0] + (VRS2)*r[0]*s[1] + (VS13S13)*s[0]*s[0] \
+ (VS13S2)*s[0]*s[1] + (VS2S2)*s[1]*s[1] + (VP)*r[1] + (VRP)*r[0]*r[1] \
+ (VPS13)*r[1]*s[0] + (VPS2)*r[1]*s[1] + (VPP)*r[1]*r[1] \
+ (VPS13S13)*r[1]*s[0]*s[0] + (VPPS13)*r[1]*r[1]*s[0] \
+ (VPPS2)*r[1]*r[1]*s[1] + (VRPS13)*r[0]*r[1]*s[0] + (VRPS2)*r[0]*r[1]*s[1] \
+ (VPS2S2)*r[1]*s[1]*s[1] + (VPS13S2)*r[1]*s[0]*s[1] + (VRRP)*r[0]*r[0]*r[1] \
+ (VRPP)*r[0]*r[1]*r[1] + (VPPP)*r[1]*r[1]*r[1]
#define G   -t*(SIC) + \
(G0) + (GR)*r[0] + (GS13)*s[0] + (GS2)*s[1] + (GRR)*r[0]*r[0] \
+ (GRS13)*r[0]*s[0] + (GRS2)*r[0]*s[1] + (GS13S13)*s[0]*s[0] \
+ (GS13S2)*s[0]*s[1] + (GS2S2)*s[1]*s[1] + (GP)*r[1] + (GRP)*r[0]*r[1] \
+ (GPS13)*r[1]*s[0] + (GPS2)*r[1]*s[1] + (GPP)*r[1]*r[1] \
+ (GPS13S13)*r[1]*s[0]*s[0] + (GPPS13)*r[1]*r[1]*s[0] \
+ (GPPS2)*r[1]*r[1]*s[1] + (GRPS13)*r[0]*r[1]*s[0] + (GRPS2)*r[0]*r[1]*s[1] \
+ (GPS2S2)*r[1]*s[1]*s[1] + (GPS13S2)*r[1]*s[0]*s[1] + (GRRP)*r[0]*r[0]*r[1] \
+ (GRPP)*r[0]*r[1]*r[1] + (GPPP)*r[1]*r[1]*r[1]

/*----------------------------------------------------------------------------*/

#define DGDR0 R*t*(log(xfe2m4) - log(xmg2m4) + 3.0*log(xfe2m13)/2.0 \
- 3.0*log(xmg2m13)/2.0 + log(xfe2m2) - log(xmg2m2) ) \
+ (GR) + 2.0*(GRR)*r[0] + (GRS13)*s[0] + (GRS2)*s[1] \
+ (GRP)*r[1] + (GRPS13)*r[1]*s[0] + (GRPS2)*r[1]*s[1] \
+ 2.0*(GRRP)*r[0]*r[1] + (GRPP)*r[1]*r[1]
#define DGDR1 R*t*(-2.0*log(xmg2m4) + 2.0*log(xca2m4) ) \
+ (GP) + (GRP)*r[0] + (GPS13)*s[0] + (GPS2)*s[1] + 2.0*(GPP)*r[1] \
+ (GPS13S13)*s[0]*s[0] + 2.0*(GPPS13)*r[1]*s[0] \
+ 2.0*(GPPS2)*r[1]*s[1] + (GRPS13)*r[0]*s[0] + (GRPS2)*r[0]*s[1] \
+ (GPS2S2)*s[1]*s[1] + (GPS13S2)*s[0]*s[1] + (GRRP)*r[0]*r[0] \
+ 2.0*(GRPP)*r[0]*r[1] + 3.0*(GPPP)*r[1]*r[1]
#define DGDS0 R*t*(6.0*log(xfe2m4)/7.0 - 6.0*log(xmg2m4)/7.0 \
- 12.0*log(xfe2m13)/7.0 + 12.0*log(xmg2m13)/7.0 \
+ 6.0*log(xfe2m2)/7.0 - 6.0*log(xmg2m2)/7.0 ) \
+ (GS13) + (GRS13)*r[0] + 2.0*(GS13S13)*s[0] + (GS13S2)*s[1] \
+ (GPS13)*r[1] + 2.0*(GPS13S13)*r[1]*s[0] + (GPPS13)*r[1]*r[1] \
+ (GRPS13)*r[0]*r[1] + (GPS13S2)*r[1]*s[1]
#define DGDS1 R*t*(4.0*log(xfe2m4)/7.0 - 4.0*log(xmg2m4)/7.0 \
+ 6.0*log(xfe2m13)/7.0 - 6.0*log(xmg2m13)/7.0 \
- 10.0*log(xfe2m2)/7.0 + 10.0*log(xmg2m2)/7.0 ) \
+ (GS2) + (GRS2)*r[0] + (GS13S2)*s[0] + 2.0*(GS2S2)*s[1] \
+ (GPS2)*r[1] + (GPPS2)*r[1]*r[1] + (GRPS2)*r[0]*r[1] \
+ 2.0*(GPS2S2)*r[1]*s[1] + (GPS13S2)*r[1]*s[0]
#define DGDT  -(S)
#define DGDP   (V)

/*----------------------------------------------------------------------------*/

#define D2GDR0R0 R*t*(1.0/(2.0*xfe2m4) + 1.0/(2.0*xmg2m4) + 3.0/(4.0*xfe2m13) \
+ 3.0/(4.0*xmg2m13) + 1.0/(2.0*xfe2m2) + 1.0/(2.0*xmg2m2) ) \
+ 2.0*(GRR) + 2.0*(GRRP)*r[1]
#define D2GDR0R1 R*t*(1.0/xmg2m4) + (GRP) + (GRPS13)*s[0] + (GRPS2)*s[1] \
+ 2.0*(GRRP)*r[0] + 2.0*(GRPP)*r[1]
#define D2GDR0S0 R*t*(3.0/(7.0*xfe2m4) + 3.0/(7.0*xmg2m4) - 6.0/(7.0*xfe2m13) \
- 6.0/(7.0*xmg2m13) + 3.0/(7.0*xfe2m2) + 3.0/(7.0*xmg2m2) ) \
+ (GRS13) + (GRPS13)*r[1]
#define D2GDR0S1 R*t*(2.0/(7.0*xfe2m4) + 2.0/(7.0*xmg2m4) + 3.0/(7.0*xfe2m13) \
+ 3.0/(7.0*xmg2m13) - 5.0/(7.0*xfe2m2) - 5.0/(7.0*xmg2m2) ) \
+ (GRS2) + (GRPS2)*r[1]
#define D2GDR0DT R*(log(xfe2m4) - log(xmg2m4) + 3.0*log(xfe2m13)/2.0 \
- 3.0*log(xmg2m13)/2.0 + log(xfe2m2) - log(xmg2m2) ) \
- (SR) - 2.0*(SRR)*r[0] - (SRS13)*s[0] - (SRS2)*s[1] \
- (SRP)*r[1] - (SRPS13)*r[1]*s[0] - (SRPS2)*r[1]*s[1] \
- 2.0*(SRRP)*r[0]*r[1] - (SRPP)*r[1]*r[1]
#define D2GDR0DP (VR) + 2.0*(VRR)*r[0] + (VRS13)*s[0] + (VRS2)*s[1] \
+ (VRP)*r[1] + (VRPS13)*r[1]*s[0] + (VRPS2)*r[1]*s[1] \
+ 2.0*(VRRP)*r[0]*r[1] + (VRPP)*r[1]*r[1]

#define D2GDR1R1 R*t*(2.0/(xmg2m4) + 2.0/(xca2m4) ) \
+ 2.0*(GPP) + 2.0*(GPPS13)*s[0] + 2.0*(GPPS2)*s[1] \
+ 2.0*(GRPP)*r[0] + 6.0*(GPPP)*r[1]
#define D2GDR1S0 R*t*(6.0/(7.0*xmg2m4) ) \
+ (GPS13) + 2.0*(GPS13S13)*s[0] + 2.0*(GPPS13)*r[1] \
+ (GRPS13)*r[0] + (GPS13S2)*s[1]
#define D2GDR1S1 R*t*(4.0/(7.0*xmg2m4) ) \
+ (GPS2) + 2.0*(GPPS2)*r[1] + (GRPS2)*r[0] \
+ 2.0*(GPS2S2)*s[1] + (GPS13S2)*s[0]
#define D2GDR1DT R*(-2.0*log(xmg2m4) + 2.0*log(xca2m4) ) \
- (SP) - (SRP)*r[0] - (SPS13)*s[0] - (SPS2)*s[1] - 2.0*(SPP)*r[1] \
- (SPS13S13)*s[0]*s[0] - 2.0*(SPPS13)*r[1]*s[0] \
- 2.0*(SPPS2)*r[1]*s[1] - (SRPS13)*r[0]*s[0] - (SRPS2)*r[0]*s[1] \
- (SPS2S2)*s[1]*s[1] - (SPS13S2)*s[0]*s[1] - (SRRP)*r[0]*r[0] \
- 2.0*(SRPP)*r[0]*r[1] - 3.0*(SPPP)*r[1]*r[1]
#define D2GDR1DP (VP) + (VRP)*r[0] + (VPS13)*s[0] + (VPS2)*s[1] + 2.0*(VPP)*r[1] \
+ (VPS13S13)*s[0]*s[0] + 2.0*(VPPS13)*r[1]*s[0] \
+ 2.0*(VPPS2)*r[1]*s[1] + (VRPS13)*r[0]*s[0] + (VRPS2)*r[0]*s[1] \
+ (VPS2S2)*s[1]*s[1] + (VPS13S2)*s[0]*s[1] + (VRRP)*r[0]*r[0] \
+ 2.0*(VRPP)*r[0]*r[1] + 3.0*(VPPP)*r[1]*r[1]

#define D2GDS0S0 R*t*(6.0*3.0/(7.0*7.0*xfe2m4) + 6.0*3.0/(7.0*7.0*xmg2m4) \
+ 12.0*4.0/(7.0*7.0*xfe2m13) + 12.0*4.0/(7.0*7.0*xmg2m13) \
+ 6.0*3.0/(7.0*7.0*xfe2m2) + 6.0*3.0/(7.0*7.0*xmg2m2) ) \
+ 2.0*(GS13S13) + 2.0*(GPS13S13)*r[1]
#define D2GDS0S1 R*t*(6.0*2.0/(7.0*7.0*xfe2m4) + 6.0*2.0/(7.0*7.0*xmg2m4) \
- 12.0*2.0/(7.0*7.0*xfe2m13) - 12.0*2.0/(7.0*7.0*xmg2m13) \
- 6.0*5.0/(7.0*7.0*xfe2m2) - 6.0*5.0/(7.0*7.0*xmg2m2) ) \
+ (GS13S2) + (GPS13S2)*r[1]
#define D2GDS0DT R*(6.0*log(xfe2m4)/7.0 - 6.0*log(xmg2m4)/7.0 \
- 12.0*log(xfe2m13)/7.0 + 12.0*log(xmg2m13)/7.0 \
+ 6.0*log(xfe2m2)/7.0 - 6.0*log(xmg2m2)/7.0 ) \
- (SS13) - (SRS13)*r[0] - 2.0*(SS13S13)*s[0] - (SS13S2)*s[1] \
- (SPS13)*r[1] - 2.0*(SPS13S13)*r[1]*s[0] - (SPPS13)*r[1]*r[1] \
- (SRPS13)*r[0]*r[1] - (SPS13S2)*r[1]*s[1]
#define D2GDS0DP (VS13) + (VRS13)*r[0] + 2.0*(VS13S13)*s[0] + (VS13S2)*s[1] \
+ (VPS13)*r[1] + 2.0*(VPS13S13)*r[1]*s[0] + (VPPS13)*r[1]*r[1] \
+ (VRPS13)*r[0]*r[1] + (VPS13S2)*r[1]*s[1]

#define D2GDS1S1 R*t*(4.0*2.0/(7.0*7.0*xfe2m4) + 4.0*2.0/(7.0*7.0*xmg2m4) \
+ 6.0*2.0/(7.0*7.0*xfe2m13) + 6.0*2.0/(7.0*7.0*xmg2m13) \
+ 10.0*5.0/(7.0*7.0*xfe2m2) + 10.0*5.0/(7.0*7.0*xmg2m2) ) \
+ 2.0*(GS2S2) + 2.0*(GPS2S2)*r[1]
#define D2GDS1DT R*(4.0*log(xfe2m4)/7.0 - 4.0*log(xmg2m4)/7.0 \
+ 6.0*log(xfe2m13)/7.0 - 6.0*log(xmg2m13)/7.0 \
- 10.0*log(xfe2m2)/7.0 + 10.0*log(xmg2m2)/7.0 ) \
- (SS2) - (SRS2)*r[0] - (SS13S2)*s[0] - 2.0*(SS2S2)*s[1] \
- (SPS2)*r[1] - (SPPS2)*r[1]*r[1] - (SRPS2)*r[0]*r[1] \
- 2.0*(SPS2S2)*r[1]*s[1] - (SPS13S2)*r[1]*s[0]
#define D2GDS1DP (VS2) + (VRS2)*r[0] + (VS13S2)*s[0] + 2.0*(VS2S2)*s[1] \
+ (VPS2)*r[1] + (VPPS2)*r[1]*r[1] + (VRPS2)*r[0]*r[1] \
+ 2.0*(VPS2S2)*r[1]*s[1] + (VPS13S2)*r[1]*s[0]

#define D2GDT2   0.0
#define D2GDTDP  0.0
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/

#define D3GDR0R0R0 R*t*(-1.0/(4.0*xfe2m4*xfe2m4) + 1.0/(4.0*xmg2m4*xmg2m4) \
- 3.0/(8.0*xfe2m13*xfe2m13) + 3.0/(8.0*xmg2m13*xmg2m13) \
- 1.0/(4.0*xfe2m2*xfe2m2) + 1.0/(4.0*xmg2m2*xmg2m2) )
#define D3GDR0R0R1 R*t*(1.0/(2.0*xmg2m4*xmg2m4) ) + 2.0*(GRRP)
#define D3GDR0R0S0 R*t*(-3.0/(2.0*7.0*xfe2m4*xfe2m4) + 3.0/(2.0*7.0*xmg2m4*xmg2m4) \
+ 3.0*4.0/(4.0*7.0*xfe2m13*xfe2m13) - 3.0*4.0/(4.0*7.0*xmg2m13*xmg2m13) \
- 3.0/(2.0*7.0*xfe2m2*xfe2m2) + 3.0/(2.0*7.0*xmg2m2*xmg2m2) )
#define D3GDR0R0S1 R*t*(-2.0/(2.0*7.0*xfe2m4*xfe2m4) + 2.0/(2.0*7.0*xmg2m4*xmg2m4) \
- 3.0*2.0/(4.0*7.0*xfe2m13*xfe2m13) + 3.0*2.0/(4.0*7.0*xmg2m13*xmg2m13) \
+ 5.0/(2.0*7.0*xfe2m2*xfe2m2) - 5.0/(2.0*7.0*xmg2m2*xmg2m2) )
#define D3GDR0R0DT R*(1.0/(2.0*xfe2m4) + 1.0/(2.0*xmg2m4) + 3.0/(4.0*xfe2m13) \
+ 3.0/(4.0*xmg2m13) + 1.0/(2.0*xfe2m2) + 1.0/(2.0*xmg2m2) ) \
- 2.0*(SRR) - 2.0*(SRRP)*r[1]
#define D3GDR0R0DP 2.0*(VRR) + 2.0*(VRRP)*r[1]

#define D3GDR0R1R1 R*t*(1.0/(xmg2m4*xmg2m4)) + 2.0*(GRPP)
#define D3GDR0R1S0 R*t*(3.0/(7.0*xmg2m4*xmg2m4)) + (GRPS13)
#define D3GDR0R1S1 R*t*(2.0/(7.0*xmg2m4*xmg2m4)) + (GRPS2)
#define D3GDR0R1DT R*(1.0/xmg2m4) - (SRP) - (SRPS13)*s[0] - (SRPS2)*s[1] \
- 2.0*(SRRP)*r[0] - 2.0*(SRPP)*r[1]
#define D3GDR0R1DP (VRP) + (VRPS13)*s[0] + (VRPS2)*s[1] + 2.0*(VRRP)*r[0] \
+ 2.0*(VRPP)*r[1]

#define D3GDR0S0S0 R*t*(-3.0*3.0/(7.0*7.0*xfe2m4*xfe2m4) + 3.0*3.0/(7.0*7.0*xmg2m4*xmg2m4) \
- 6.0*4.0/(7.0*7.0*xfe2m13*xfe2m13) + 6.0*4.0/(7.0*7.0*xmg2m13*xmg2m13) \
- 3.0*3.0/(7.0*7.0*xfe2m2*xfe2m2) + 3.0*3.0/(7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDR0S0S1 R*t*(-3.0*2.0/(7.0*7.0*xfe2m4*xfe2m4) + 3.0*2.0/(7.0*7.0*xmg2m4*xmg2m4) \
+ 6.0*2.0/(7.0*7.0*xfe2m13*xfe2m13) - 6.0*2.0/(7.0*7.0*xmg2m13*xmg2m13) \
+ 3.0*5.0/(7.0*7.0*xfe2m2*xfe2m2) - 3.0*5.0/(7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDR0S0DT R*(3.0/(7.0*xfe2m4) + 3.0/(7.0*xmg2m4) - 6.0/(7.0*xfe2m13) \
- 6.0/(7.0*xmg2m13) + 3.0/(7.0*xfe2m2) + 3.0/(7.0*xmg2m2) ) \
- (SRS13) - (SRPS13)*r[1]
#define D3GDR0S0DP (VRS13) + (VRPS13)*r[1]

#define D3GDR0S1S1 R*t*(-2.0*2.0/(7.0*7.0*xfe2m4*xfe2m4) + 2.0*2.0/(7.0*7.0*xmg2m4*xmg2m4) \
- 3.0*2.0/(7.0*7.0*xfe2m13*xfe2m13) + 3.0*2.0/(7.0*7.0*xmg2m13*xmg2m13) \
- 5.0*5.0/(7.0*7.0*xfe2m2*xfe2m2) + 5.0*5.0/(7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDR0S1DT R*(2.0/(7.0*xfe2m4) + 2.0/(7.0*xmg2m4) + 3.0/(7.0*xfe2m13) \
+ 3.0/(7.0*xmg2m13) - 5.0/(7.0*xfe2m2) - 5.0/(7.0*xmg2m2) ) \
- (SRS2) - (SRPS2)*r[1]
#define D3GDR0S1DP (VRS2) + (VRPS2)*r[1]

#define D3GDR1R1R1 R*t*(2.0/(xmg2m4*xmg2m4) - 2.0/(xca2m4*xca2m4) ) + 6.0*(GPPP)
#define D3GDR1R1S0 R*t*(2.0*3.0/(7.0*xmg2m4*xmg2m4) ) + 2.0*(GPPS13)
#define D3GDR1R1S1 R*t*(2.0*2.0/(7.0*xmg2m4*xmg2m4) ) + 2.0*(GPPS2)
#define D3GDR1R1DT R*(2.0/(xmg2m4) + 2.0/(xca2m4) ) - 2.0*(SPP) \
- 2.0*(SPPS13)*s[0] - 2.0*(SPPS2)*s[1] - 2.0*(SRPP)*r[0] - 6.0*(SPPP)*r[1]
#define D3GDR1R1DP 2.0*(VPP) + 2.0*(VPPS13)*s[0] + 2.0*(VPPS2)*s[1] \
+ 2.0*(VRPP)*r[0] + 6.0*(VPPP)*r[1]

#define D3GDR1S0S0 R*t*(6.0*3.0/(7.0*7.0*xmg2m4*xmg2m4) ) + 2.0*(GPS13S13)
#define D3GDR1S0S1 R*t*(6.0*2.0/(7.0*7.0*xmg2m4*xmg2m4) ) + (GPS13S2)
#define D3GDR1S0DT R*(6.0/(7.0*xmg2m4) ) \
- (SPS13) - 2.0*(SPS13S13)*s[0] - 2.0*(SPPS13)*r[1] \
- (SRPS13)*r[0] - (SPS13S2)*s[1]
#define D3GDR1S0DP (VPS13) + 2.0*(VPS13S13)*s[0] + 2.0*(VPPS13)*r[1] \
+ (VRPS13)*r[0] + (VPS13S2)*s[1]

#define D3GDR1S1S1 R*t*(4.0*2.0/(7.0*7.0*xmg2m4*xmg2m4) ) + 2.0*(GPS2S2)
#define D3GDR1S1DT R*(4.0/(7.0*xmg2m4) ) \
- (SPS2) - 2.0*(SPPS2)*r[1] - (SRPS2)*r[0] \
- 2.0*(SPS2S2)*s[1] - (SPS13S2)*s[0]
#define D3GDR1S1DP (VPS2) + 2.0*(VPPS2)*r[1] + (VRPS2)*r[0] \
+ 2.0*(VPS2S2)*s[1] + (VPS13S2)*s[0]

#define D3GDS0S0S0 R*t*(-6.0*3.0*3.0/(7.0*7.0*7.0*xfe2m4*xfe2m4) + 6.0*3.0*3.0/(7.0*7.0*7.0*xmg2m4*xmg2m4) \
+ 12.0*4.0*4.0/(7.0*7.0*7.0*xfe2m13*xfe2m13) - 12.0*4.0*4.0/(7.0*7.0*7.0*xmg2m13*xmg2m13) \
- 6.0*3.0*3.0/(7.0*7.0*7.0*xfe2m2*xfe2m2) + 6.0*3.0*3.0/(7.0*7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDS0S0S1 R*t*(-6.0*3.0*2.0/(7.0*7.0*7.0*xfe2m4*xfe2m4) + 6.0*3.0*2.0/(7.0*7.0*7.0*xmg2m4*xmg2m4) \
- 12.0*4.0*2.0/(7.0*7.0*7.0*xfe2m13*xfe2m13) + 12.0*4.0*2.0/(7.0*7.0*7.0*xmg2m13*xmg2m13) \
+ 6.0*3.0*5.0/(7.0*7.0*7.0*xfe2m2*xfe2m2) - 6.0*3.0*5.0/(7.0*7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDS0S0DT R*(6.0*3.0/(7.0*7.0*xfe2m4) + 6.0*3.0/(7.0*7.0*xmg2m4) \
+ 12.0*4.0/(7.0*7.0*xfe2m13) + 12.0*4.0/(7.0*7.0*xmg2m13) \
+ 6.0*3.0/(7.0*7.0*xfe2m2) + 6.0*3.0/(7.0*7.0*xmg2m2) ) \
- 2.0*(SS13S13) - 2.0*(SPS13S13)*r[1]
#define D3GDS0S0DP 2.0*(VS13S13) + 2.0*(VPS13S13)*r[1]

#define D3GDS0S1S1 R*t*(-6.0*2.0*2.0/(7.0*7.0*7.0*xfe2m4*xfe2m4) + 6.0*2.0*2.0/(7.0*7.0*7.0*xmg2m4*xmg2m4) \
+ 12.0*2.0*2.0/(7.0*7.0*7.0*xfe2m13*xfe2m13) - 12.0*2.0*2.0/(7.0*7.0*7.0*xmg2m13*xmg2m13) \
- 6.0*5.0*5.0/(7.0*7.0*7.0*xfe2m2*xfe2m2) + 6.0*5.0*5.0/(7.0*7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDS0S1DT R*(6.0*2.0/(7.0*7.0*xfe2m4) + 6.0*2.0/(7.0*7.0*xmg2m4) \
- 12.0*2.0/(7.0*7.0*xfe2m13) - 12.0*2.0/(7.0*7.0*xmg2m13) \
- 6.0*5.0/(7.0*7.0*xfe2m2) - 6.0*5.0/(7.0*7.0*xmg2m2) ) \
- (SS13S2) - (SPS13S2)*r[1]
#define D3GDS0S1DP (VS13S2) + (VPS13S2)*r[1]

#define D3GDS0DT2  0.0
#define D3GDS0DTDP 0.0
#define D3GDS0DP2  0.0

#define D3GDS1S1S1 R*t*(-4.0*2.0*2.0/(7.0*7.0*7.0*xfe2m4*xfe2m4) + 4.0*2.0*2.0/(7.0*7.0*7.0*xmg2m4*xmg2m4) \
- 6.0*2.0*2.0/(7.0*7.0*7.0*xfe2m13*xfe2m13) + 6.0*2.0*2.0/(7.0*7.0*7.0*xmg2m13*xmg2m13) \
+ 10.0*5.0*5.0/(7.0*7.0*7.0*xfe2m2*xfe2m2) - 10.0*5.0*5.0/(7.0*7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDS1S1DT R*(4.0*2.0/(7.0*7.0*xfe2m4) + 4.0*2.0/(7.0*7.0*xmg2m4) \
+ 6.0*2.0/(7.0*7.0*xfe2m13) + 6.0*2.0/(7.0*7.0*xmg2m13) \
+ 10.0*5.0/(7.0*7.0*xfe2m2) + 10.0*5.0/(7.0*7.0*xmg2m2) ) \
- 2.0*(SS2S2) - 2.0*(SPS2S2)*r[1]
#define D3GDS1S1DP 2.0*(VS2S2) + 2.0*(VPS2S2)*r[1]

#define D3GDS1DT2  0.0
#define D3GDS1DTDP 0.0
#define D3GDS1DP2  0.0

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
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1; \
d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;

#define fillD2GDRDS \
d2gdrds[0][0] = D2GDR0S0; d2gdrds[0][1] = D2GDR0S1; \
d2gdrds[1][0] = D2GDR1S0; d2gdrds[1][1] = D2GDR1S1;

#define fillD2GDRDT \
d2gdrdt[0] = D2GDR0DT; d2gdrdt[1] = D2GDR1DT;

#define fillD2GDRDP \
d2gdrdp[0] = D2GDR0DP; d2gdrdp[1] = D2GDR1DP;

#define fillD2GDS2 \
d2gds2[0][0] = D2GDS0S0;     d2gds2[0][1] = D2GDS0S1; \
d2gds2[1][0] = d2gds2[0][1]; d2gds2[1][1] = D2GDS1S1;

#define fillD2GDSDT \
d2gdsdt[0] = D2GDS0DT;  d2gdsdt[1] = D2GDS1DT;

#define fillD2GDSDP \
d2gdsdp[0] = D2GDS0DP;  d2gdsdp[1] = D2GDS1DP;

#define fillD3GDR3 \
d3gdr3[0][0][0] = D3GDR0R0R0;		d3gdr3[0][0][1] = D3GDR0R0R1; \
d3gdr3[0][1][0] = d3gdr3[0][0][1];	d3gdr3[0][1][1] = D3GDR0R1R1; \
d3gdr3[1][0][0] = d3gdr3[0][0][1];	d3gdr3[1][0][1] = d3gdr3[0][1][1]; \
d3gdr3[1][1][0] = d3gdr3[0][1][1];	d3gdr3[1][1][1] = D3GDR1R1R1;

#define fillD3GDR2DS \
d3gdr2ds[0][0][0] = D3GDR0R0S0;        d3gdr2ds[0][0][1] = D3GDR0R0S1; \
d3gdr2ds[0][1][0] = D3GDR0R1S0;        d3gdr2ds[0][1][1] = D3GDR0R1S1; \
d3gdr2ds[1][0][0] = d3gdr2ds[0][1][0]; d3gdr2ds[1][0][1] = d3gdr2ds[0][1][1]; \
d3gdr2ds[1][1][0] = D3GDR1R1S0;        d3gdr2ds[1][1][1] = D3GDR1R1S1;

#define fillD3GDR2DT \
d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT; \
d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT;

#define fillD3GDR2DP \
d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP; \
d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP;

#define fillD3GDRDS2 \
d3gdrds2[0][0][0] = D3GDR0S0S0;        d3gdrds2[0][0][1] = D3GDR0S0S1; \
d3gdrds2[0][1][0] = d3gdrds2[0][0][1]; d3gdrds2[0][1][1] = D3GDR0S1S1; \
d3gdrds2[1][0][0] = D3GDR1S0S0;        d3gdrds2[1][0][1] = D3GDR1S0S1; \
d3gdrds2[1][1][0] = d3gdrds2[1][0][1]; d3gdrds2[1][1][1] = D3GDR1S1S1;

#define fillD3GDRDSDT \
d3gdrdsdt[0][0] = D3GDR0S0DT;  d3gdrdsdt[0][1] = D3GDR0S1DT; \
d3gdrdsdt[1][0] = D3GDR1S0DT;  d3gdrdsdt[1][1] = D3GDR1S1DT;

#define fillD3GDRDSDP \
d3gdrdsdp[0][0] = D3GDR0S0DP; d3gdrdsdp[0][1] = D3GDR0S1DP; \
d3gdrdsdp[1][0] = D3GDR1S0DP; d3gdrdsdp[1][1] = D3GDR1S1DP;

#define fillD3GDS3 \
d3gds3[0][0][0] = D3GDS0S0S0;      d3gds3[0][0][1] = D3GDS0S0S1; \
d3gds3[0][1][0] = d3gds3[0][0][1]; d3gds3[0][1][1] = D3GDS0S1S1; \
d3gds3[1][0][0] = d3gds3[0][0][1]; d3gds3[1][0][1] = d3gds3[0][1][1]; \
d3gds3[1][1][0] = d3gds3[0][1][1]; d3gds3[1][1][1] = D3GDS1S1S1;

#define fillD3GDS2DT \
d3gds2dt[0][0] = D3GDS0S0DT;     d3gds2dt[0][1] = D3GDS0S1DT; \
d3gds2dt[1][0] = d3gds2dt[0][1]; d3gds2dt[1][1] = D3GDS1S1DT;

#define fillD3GDS2DP \
d3gds2dp[0][0] = D3GDS0S0DP;     d3gds2dp[0][1] = D3GDS0S1DP; \
d3gds2dp[1][0] = d3gds2dp[0][1]; d3gds2dp[1][1] = D3GDS1S1DP;

#define fillD3GDSDT2 \
d3gdsdt2[0] = D3GDS0DT2; d3gdsdt2[1] = D3GDS1DT2;

#define fillD3GDSDTDP \
d3gdsdtdp[0] = D3GDS0DTDP; d3gdsdtdp[1] = D3GDS1DTDP;

#define fillD3GDSDP2 \
d3gdsdp2[0] = D3GDS0DP2; d3gdsdp2[1] = D3GDS1DP2;

#define fillD3GDRDT2 \
d3gdrdt2[0] = D3GDR0DT2; d3gdrdt2[1] = D3GDR1DT2;

#define fillD3GDRDTDP \
d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP;

#define fillD3GDRDP2 \
d3gdrdp2[0] = D3GDR0DP2; d3gdrdp2[1] = D3GDR1DP2;

/**
 Matrix inversion routine converted from Numerical recipies to have zero array indexing
 */

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

-(void)gaussj:(double [NS][NS])a {
	int indxc[NS], indxr[NS], ipiv[NS];
	int i, icol = -1, irow = -1, j, k, l,ll;
	double big, dum, pivinv, temp;

	for (j=0; j<NS; j++) ipiv[j]=0;
	for (i=0; i<NS; i++) {
		big=0.0;
		for (j=0; j<NS; j++)
			if (ipiv[j] != 1)
				for (k=0; k<NS; k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big = fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					} else if (ipiv[k] > 1) return;
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0; l<NS; l++) SWAP(a[irow][l],a[icol][l])
				}
		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol][icol] == 0.0) return;
		pivinv = 1.0/a[icol][icol];
		a[icol][icol] = 1.0;
		for (l=0; l<NS; l++) a[icol][l] *= pivinv;
		for (ll=0; ll<NS; ll++)
			if (ll != icol) {
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l=0; l<NS; l++) a[ll][l] -= a[icol][l]*dum;
			}
	}
	for (l=(NS-1); l>=0; l--) {
		if (indxr[l] != indxc[l])
			for (k=0; k<NS; k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
}

/*
 *=============================================================================
 * Local function to compute ordering state and associated derivatives
 */

-(void)order:(int)mask
		   t:(double)t
		   p:(double)p
		   r:(double [NR])r
		   s:(double [NS])s           // s[NS]                BINARY MASK: 0000000001
		  dr:(double [NS][NR])dr      // ds[NS]/dr[NR]        BINARY MASK: 0000000010
		  dt:(double [NS])dt          // ds[NS]/dt            BINARY MASK: 0000000100
		  dp:(double [NS])dp          // ds[NS]/dp            BINARY MASK: 0000001000
		 dr2:(double [NS][NR][NR])dr2 // d2s[NS]/dr[NR]dr[NR] BINARY MASK: 0000010000
		 drt:(double [NS][NR])drt     // d2s[NS]/dr[NR]dt     BINARY MASK: 0000100000
		 drp:(double [NS][NR])drp     // d2s[NS]/dr[NR]dp     BINARY MASK: 0001000000
		 dt2:(double [NS])dt2         // d2s[NS]/dt2          BINARY MASK: 0010000000
		 dtp:(double [NS])dtp         // d2s[NS]/dtp          BINARY MASK: 0100000000
		 dp2:(double [NS])dp2         // d2s[NS]/dp2          BINARY MASK: 1000000000
{

	int i, j, iter=0;

	/* look-up or compute the current ordering state */
	if ( (t != tOld)       || (p != pOld)  || (r[0] != rOld[0]) || (r[1] != rOld[1]) ) {
		double dgds[NS], sNew[NS], xfe, xmg;

		for (i=0; i<NS; i++) { sOld[i] = 2.0; dgds[i] = 0.0; }

		xca2m4  = r[1];
		xfe     = (1.0+r[0])/2.0;
		xmg     = (7.0-2.0*xca2m4-7.0*xfe)/7.0;
		sNew[0] = ((1.0-xca2m4)*xfe - xfe)/(xfe+xmg);
		sNew[1] = ((1.0-xca2m4)*xfe - xfe)/(xfe+xmg);

		if (xfe < 10.0*DBL_EPSILON) {
			xfe2m4  = DBL_EPSILON;
			xfe2m13 = DBL_EPSILON;
			xfe2m2  = DBL_EPSILON;
			xmg2m4  = 1.0 - xca2m4;
			xmg2m13 = 1.0;
			xmg2m2  = 1.0;
			sOld[0] = sNew[0];
			sOld[1] = sNew[1];
		}

		while ( ((fabs(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) ||
				 (fabs(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON)    ) &&
			   (iter < MAX_ITER)) {
			double s[NS], sCorr[NS], lambda = 1.0;

			for (i=0; i<NS; i++) s[i] = sNew[i];

			xfe2m4  = (1.0 + r[0] + 6.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
			xfe2m13 = (1.0 + r[0] - 8.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
			xfe2m2  = (1.0 + r[0] + 6.0*s[0]/7.0 - 10.0*s[1]/7.0)/2.0;
			xmg2m4  = 1.0 - xfe2m4  - xca2m4;
			xmg2m13 = 1.0 - xfe2m13;
			xmg2m2  = 1.0 - xfe2m2;

			if (xfe2m4  <= DBL_EPSILON) xfe2m4  = DBL_EPSILON;
			if (xfe2m13 <= DBL_EPSILON) xfe2m13 = DBL_EPSILON;
			if (xfe2m2  <= DBL_EPSILON) xfe2m2  = DBL_EPSILON;
			if (xmg2m4  <= DBL_EPSILON) xmg2m4  = DBL_EPSILON;
			if (xmg2m13 <= DBL_EPSILON) xmg2m13 = DBL_EPSILON;
			if (xmg2m2  <= DBL_EPSILON) xmg2m2  = DBL_EPSILON;

			dgds[0] = DGDS0;
			dgds[1] = DGDS1;

			invd2gds2[0][0] = D2GDS0S0;
			invd2gds2[0][1] = D2GDS0S1;
			invd2gds2[1][0] = invd2gds2[0][1];
			invd2gds2[1][1] = D2GDS1S1;


			for (i=0; i<NS; i++) sOld[i] = s[i];

			[self gaussj:invd2gds2];

			for (i=0; i<NS; i++) {
				for(j=0, sCorr[i]=0.0; j<NS; j++) sCorr[i] += - invd2gds2[i][j]*dgds[j];
				s[i] = sOld[i] + lambda*sCorr[i];
			}

			xfe2m4  = (1.0 + r[0] + 6.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
			xfe2m13 = (1.0 + r[0] - 8.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
			xfe2m2  = (1.0 + r[0] + 6.0*s[0]/7.0 - 10.0*s[1]/7.0)/2.0;
			xmg2m4  = 1.0 - xfe2m4  - xca2m4;
			xmg2m13 = 1.0 - xfe2m13;
			xmg2m2  = 1.0 - xfe2m2;

			while ((   xfe2m4  < 0.0 || xfe2m4  > 1.0 || xmg2m4  < 0.0 || xmg2m4  > 1.0
					|| xfe2m13 < 0.0 || xfe2m13 > 1.0 || xmg2m13 < 0.0 || xmg2m13 > 1.0
					|| xfe2m2  < 0.0 || xfe2m2  > 1.0 || xmg2m2  < 0.0 || xmg2m2  > 1.0)
				   && lambda > DBL_EPSILON ) {
				lambda /= 2.0;
				for (i=0; i<NS; i++) s[i] = sOld[i] + lambda*sCorr[i];
				xfe2m4  = (1.0 + r[0] + 6.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
				xfe2m13 = (1.0 + r[0] - 8.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
				xfe2m2  = (1.0 + r[0] + 6.0*s[0]/7.0 - 10.0*s[1]/7.0)/2.0;
				xmg2m4  = 1.0 - xfe2m4  - xca2m4;
				xmg2m13 = 1.0 - xfe2m13;
				xmg2m2  = 1.0 - xfe2m2;
			}

			for (i=0; i<NS; i++) sNew[i] = s[i];
			iter++;
		}
		tOld = t;
		pOld = p;
		for (i=0; i<NR; i++) rOld[i] = r[i];

		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
		if (debug) {
			for (i=0; i<NS; i++) {
				if (dgds[i] > sqrt(DBL_EPSILON) && fabs(sOld[i]) > DBL_EPSILON) {
					NSLog(@"ERROR in AMPHIBOLE.C (function ORDER). Failed to converge!");
					NSLog(@"  X2    = %13.6g, X3    = %13.6g", r[0], r[1]);
					NSLog(@"  s1    = %13.6g, s2    = %13.6g", sOld[0], sOld[1]);
					NSLog(@"  dgds1 = %13.6g, dgds2 = %13.6g", dgds[0], dgds[1]);
					NSLog(@"  X Ca2+ m4:  %13.6g", xca2m4);
					NSLog(@"  X Mg2+ m4:  %13.6g  X Mg2+ m13: %13.6g  X Mg2+ m2: %13.6g", xmg2m4, xmg2m13, xmg2m2);
					NSLog(@"  X Fe2+ m4:  %13.6g  X Fe2+ m13: %13.6g  X Fe2+ m2: %13.6g", xfe2m4, xfe2m13, xfe2m2);
					break;
				}
			}
		}

	}

	if (mask & FIRST  ) {   /* return s        */
		for (i=0; i<NS; i++) s[i] = sOld[i];
	}
	if (mask & SECOND ) {   /* compute ds/dr:  */
		double *s = sOld;
		double d2gdrds[NR][NS];
		int k;

		fillD2GDRDS

		for (i=0; i<NS; i++) {
			for (j=0; j<NR; j++) {
				dr[i][j] = 0.0;
				for (k=0; k<NS; k++) dr[i][j] += - invd2gds2[i][k]*d2gdrds[j][k];
			}
		}
	}
	if (mask & THIRD  ) {   /* compute ds/dt:  */
		double *s = sOld;
		double d2gdsdt[NS];

		fillD2GDSDT

		for (i=0; i<NS; i++) {
			dt[i] = 0.0;
			for (j=0; j<NS; j++) dt[i] += - invd2gds2[i][j]*d2gdsdt[j];
		}
	}
	if (mask & FOURTH ) {   /* compute ds/dp:  */
		double *s = sOld;
		double d2gdsdp[NS];

		fillD2GDSDP

		for (i=0; i<NS; i++) {
			dp[i] = 0.0;
			for (j=0; j<NS; j++) dp[i] += - invd2gds2[i][j]*d2gdsdp[j];
		}
	}
	if (mask & FIFTH  ) {   /* compute d2s/dr2 */
		double *s = sOld;
		double d2gdrds[NR][NS], d3gdr2ds[NR][NR][NS], d3gdrds2[NR][NS][NS],
		d3gds3[NS][NS][NS], dsdr[NS][NR], temp[NS];
		int k, l, m, n;

		fillD2GDRDS
		fillD3GDR2DS
		fillD3GDRDS2
		fillD3GDS3

		/* compute dsdr matrix */
		for (i=0; i<NS; i++) {
			for (j=0; j<NR; j++) {
				dsdr[i][j] = 0.0;
				for (k=0; k<NS; k++) dsdr[i][j] += - invd2gds2[i][k]*d2gdrds[j][k];
			}
		}

		/* compute dsdr2 cube */
		for (i=0; i<NS; i++) {
			for (j=0; j<NR; j++) {
				for (k=0; k<NR; k++) {
					for (l=0; l<NS; l++) {
						temp[l] = d3gdr2ds[j][k][l];
						for (m=0; m<NS; m++) {
							temp[l] += d3gdrds2[j][l][m]*dsdr[m][k]
							+ d3gdrds2[k][l][m]*dsdr[m][j];
							for (n=0; n<NS; n++)
								temp[l] += d3gds3[l][m][n]*dsdr[m][j]*dsdr[n][k];
						}
					}
					dr2[i][j][k] = 0.0;
					for (l=0; l<NS; l++) dr2[i][j][k] += - invd2gds2[i][l]*temp[l];
				}
			}
		}

	}
	if (mask & SIXTH  ) {   /* compute d2s/drt */
		double *s = sOld;
		double d2gdrds[NR][NS], d2gdsdt[NS], d3gdrds2[NR][NS][NS],
		d3gdrdsdt[NR][NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS], dsdr[NS][NR],
		dsdt[NS], temp[NS];
		int k, l, m;

		fillD2GDRDS
		fillD2GDSDT
		fillD3GDRDS2
		fillD3GDRDSDT
		fillD3GDS3
		fillD3GDS2DT

		/* compute dsdr matrix */
		for (i=0; i<NS; i++) {
			for (j=0; j<NR; j++) {
				dsdr[i][j] = 0.0;
				for (k=0; k<NS; k++) dsdr[i][j] += - invd2gds2[i][k]*d2gdrds[j][k];
			}
		}

		/* compute dsdt vector */
		for (i=0; i<NS; i++) {
			dsdt[i] = 0.0;
			for (j=0; j<NS; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
		}

		/* compute dsdrdt matrix */
		for (i=0; i<NS; i++) {
			for (j=0; j<NR; j++) {
				for (k=0; k<NS; k++) {
					temp[k] = d3gdrdsdt[j][k];
					for (l=0; l<NS; l++) {
						temp[k] += d3gdrds2[j][k][l]*dsdt[l] + d3gds2dt[k][l]*dsdr[l][j];
						for (m=0; m<NS; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdt[m];
					}
				}
				drt[i][j] = 0.0;
				for (k=0; k<NS; k++) drt[i][j] += - invd2gds2[i][k]*temp[k];
			}
		}

	}
	if (mask & SEVENTH) {   /* compute d2s/drp */
		double *s = sOld;
		double d2gdrds[NR][NS], d2gdsdp[NS], d3gdrds2[NR][NS][NS],
		d3gdrdsdp[NR][NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS], dsdr[NS][NR],
		dsdp[NS], temp[NS];
		int k, l, m;

		fillD2GDRDS
		fillD2GDSDP
		fillD3GDRDS2
		fillD3GDRDSDP
		fillD3GDS3
		fillD3GDS2DP

		/* compute dsdr matrix */
		for (i=0; i<NS; i++) {
			for (j=0; j<NR; j++) {
				dsdr[i][j] = 0.0;
				for (k=0; k<NS; k++) dsdr[i][j] += - invd2gds2[i][k]*d2gdrds[j][k];
			}
		}

		/* compute dsdp vector */
		for (i=0; i<NS; i++) {
			dsdp[i] = 0.0;
			for (j=0; j<NS; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
		}

		/* compute dsdrdp matrix */
		for (i=0; i<NS; i++) {
			for (j=0; j<NR; j++) {
				for (k=0; k<NS; k++) {
					temp[k] = d3gdrdsdp[j][k];
					for (l=0; l<NS; l++) {
						temp[k] += d3gdrds2[j][k][l]*dsdp[l] + d3gds2dp[k][l]*dsdr[l][j];
						for (m=0; m<NS; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdp[m];
					}
				}
				drp[i][j] = 0.0;
				for (k=0; k<NS; k++) drp[i][j] += - invd2gds2[i][k]*temp[k];
			}
		}

	}
	if (mask & EIGHTH ) {   /* compute d2s/dt2 */
		double *s = sOld;
		double d2gdsdt[NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gdsdt2[NS],
		dsdt[NS], temp[NS];
		int k, l;

		fillD2GDSDT
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDSDT2

		/* compute dsdt vector */
		for (i=0; i<NS; i++) {
			dsdt[i] = 0.0;
			for (j=0; j<NS; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
		}

		/* compute dsdt2 vector */
		for (i=0; i<NS; i++) {
			for (j=0; j<NS; j++) {
				temp[j] = d3gdsdt2[j];
				for (k=0; k<NS; k++) {
					temp[j] +=  2.0*d3gds2dt[j][k]*dsdt[k];
					for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdt[l];
				}
			}
			dt2[i] = 0.0;
			for (j=0; j<NS; j++) dt2[i] += - invd2gds2[i][j]*temp[j];
		}

	}
	if (mask & NINTH  ) {   /* compute d2s/dtp */
		double *s = sOld;
		double d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS],
		d3gds2dp[NS][NS], d3gdsdtdp[NS], dsdt[NS], dsdp[NS], temp[NS];
		int k, l;

		fillD2GDSDT
		fillD2GDSDP
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDS2DP
		fillD3GDSDTDP

		/* compute dsdt vector */
		for (i=0; i<NS; i++) {
			dsdt[i] = 0.0;
			for (j=0; j<NS; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
		}

		/* compute dsdp vector */
		for (i=0; i<NS; i++) {
			dsdp[i] = 0.0;
			for (j=0; j<NS; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
		}

		/* compute dsdtp vector */
		for (i=0; i<NS; i++) {
			for (j=0; j<NS; j++) {
				temp[j] = d3gdsdtdp[j];
				for (k=0; k<NS; k++) {
					temp[j] += d3gds2dt[j][k]*dsdp[k] + d3gds2dp[j][k]*dsdt[k];
					for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdp[l];
				}
			}
			dtp[i] = 0.0;
			for (j=0; j<NS; j++) dtp[i] += - invd2gds2[i][j]*temp[j];
		}

	}
	if (mask & TENTH  ) {   /* compute d2s/dp2 */
		double *s = sOld;
		double d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS], d3gdsdp2[NS],
		dsdp[NS], temp[NS];
		int k, l;

		fillD2GDSDP
		fillD3GDS3
		fillD3GDS2DP
		fillD3GDSDP2

		/* compute dsdp vector */
		for (i=0; i<NS; i++) {
			dsdp[i] = 0.0;
			for (j=0; j<NS; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
		}

		/* compute dsdp2 vector */
		for (i=0; i<NS; i++) {
			for (j=0; j<NS; j++) {
				temp[j] = d3gdsdp2[j];
				for (k=0; k<NS; k++) {
					temp[j] +=  2.0*d3gds2dp[j][k]*dsdp[k];
					for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdp[k]*dsdp[l];
				}
			}
			dp2[i] = 0.0;
			for (j=0; j<NS; j++) dp2[i] += - invd2gds2[i][j]*temp[j];
		}

	}

}

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
	const char *NAMES[NA]    = { "cummingtonite", "grunerite", "tremolite" };
	const char *FORMULAS[NA] = { "Mg7Si8O22(OH)2", "Fe7Si8O22(OH)2", "Ca2Mg5Si8O22(OH)2" };
	BOOL result = YES;
	NSUInteger i;
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
		result = result && ((r[0]+1.0)/2.0 >= 0.0)
		&& ((r[0]+1.0)/2.0 <= 1.0-2.0*r[1]/7.0 );         /* Fe2+ */
		result = result && (r[1] >= 0.0) && (r[1] <= 1.0);                /* Ca2+ */
		result = result && (1.0-(r[0]+1.0)/2.0-r[1]+5.0*r[1]/7.0 >= 0.0)  /* Mg2+ */
		&& (1.0-(r[0]+1.0)/2.0-r[1]+5.0*r[1]/7.0 <= 1.0-2.0*r[1]/7.0);
	}
	/* Check bounds on moles of endmember components */
	if (mask & SIXTH) {
		for (i=0, sum=0.0; i<NA; i++) sum += m[i];
		result = result && (sum >= 0.0);
		if (sum > 0.0) {
			result = result && (m[1]/sum >= 0.0) && (m[1]/sum <= 1.0-2.0*m[2]/(7.0*sum)); /* Fe2+ */
			result = result && (m[2]/sum >= 0.0) && (m[2]/sum <= 1.0);                    /* Ca2+ */
			result = result && (m[0]/sum+5.0*m[2]/(7.0*sum) >= 0.0)                       /* Mg2+ */
			                && (m[0]/sum+5.0*m[2]/(7.0*sum) <= 1.0-2.0*m[2]/(7.0*sum));
		}
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
	 endmember amphibole components.
	 (2) calculates from a vector of moles of endmember components, one or
	 all of: r[], x[], dr[]/dm[] d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
	 (3) calculates from a vector of independent compositional variables
	 mole fractions of endmember components and/or the Jacobian matrix
	 dx[]/dr[]

	 In this routine it is assumed that the elements are in the order of atomic
	 numbers and that the order of amphibole components has been verified as:
	 m[0] = cummingtonite  (Mg7Si8O22(OH)2),
	 m[1] = grunerite      (Fe7Si8O22(OH)2),
	 m[2] = tremolite      (Ca2Mg5Si8O22(OH)2),

	 ----------------------------------------------------------------------------*/

	int i, j, k;

	if (inpMask == FIRST && outMask == SECOND) {
		/* Converts a vector of moles of elements into a vector of moles of
		 end-member components.                                                 */
		static const int Mg = 12;
		static const int Ca = 20;
		static const int Fe = 26;

		/* Assign moles of endmembers */
		m[0] = (e[Mg] - 5.0*e[Ca]/2.0)/7.0;
		m[1] =  e[Fe]/7.0;
		m[2] =  e[Ca]/2.0;

	} else if (inpMask == SECOND) {
		double sum;

		if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
			NSLog(@"Illegal call to conAph with inpMask = %o and outMask = %o", inpMask, outMask);

		for (i=0, sum=0.0; i<NA; i++) sum += m[i];

		if (outMask & THIRD) {
			/* Converts a vector of moles of end-member components (m) into a vector
			 of independent compositional variables (r) required as input for the
			 remaining public functions.                                          */
			r[0] = (sum != 0.0) ? 2.0*m[1]/sum - 1.0 : 0.0;
			r[1] = (sum != 0.0) ?     m[2]/sum : 0.0;
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
				for (j=0; j<NA; j++) {
					dm[0][j]  = -2.0*m[1]/SQUARE(sum);
					dm[0][j] += (j == 1) ? 2.0/sum : 0.0;
					dm[1][j]  = -(m[2])/SQUARE(sum);
					dm[1][j] += (j == 2) ? 1.0/sum : 0.0;
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
				for (j=0; j<NA; j++) {
					for (k=0; k<NA; k++) {
						d2m[0][j][k]  = 4.0*m[1]/CUBE(sum);
						d2m[0][j][k] -= (j == 1) ? 2.0/SQUARE(sum) : 0.0;
						d2m[0][j][k] -= (k == 1) ? 2.0/SQUARE(sum) : 0.0;

						d2m[1][j][k]  = 2.0*m[2]/CUBE(sum);
						d2m[1][j][k] -= (j == 2) ? 1.0/SQUARE(sum) : 0.0;
						d2m[1][j][k] -= (k == 2) ? 1.0/SQUARE(sum) : 0.0;
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
				for (j=0; j<NA; j++)  {
					for (k=0; k<NA; k++)  {
						for (l=0; l<NA; l++)  {
							d3m[0][j][k][l]  = -12.0*m[1]/QUARTIC(sum);
							d3m[0][j][k][l] += (j == 1) ? 4.0/CUBE(sum) : 0.0;
							d3m[0][j][k][l] += (k == 1) ? 4.0/CUBE(sum) : 0.0;
							d3m[0][j][k][l] += (l == 1) ? 4.0/CUBE(sum) : 0.0;

							d3m[1][j][k][l]  = -6.0*m[2]/QUARTIC(sum);
							d3m[1][j][k][l] += (j == 2) ? 2.0/CUBE(sum) : 0.0;
							d3m[1][j][k][l] += (k == 2) ? 2.0/CUBE(sum) : 0.0;
							d3m[1][j][k][l] += (l == 2) ? 2.0/CUBE(sum) : 0.0;
						}
					}
				}
			}
		}

	} else if (inpMask == THIRD) {

		if (outMask & ~(FOURTH | SEVENTH))
			NSLog(@"Illegal call to conAph with inpMask = %o and outMask = %o", inpMask, outMask);

		if (outMask & FOURTH) {
			/* Converts a vector of independent compositional variables (r) into a
			 vector of mole fractions of endmember components (x).                */
			x[0] = 1.0 - (1.0+r[0])/2.0 - r[1];
			x[1] = (1.0+r[0])/2.0;
			x[2] = r[1];
		}

		if (outMask & SEVENTH) {
			/* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
			for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
			dr[0][0] = -1.0/2.0; dr[0][1] = -1.0;
			dr[1][0] =  1.0/2.0;
			dr[2][1] =  1.0;
		}

	} else  {
		NSLog(@"Illegal call to conAph with inpMask = %o and outMask = %o", inpMask, outMask);
	}

}

-(NSString *)displayFormula:(double)t
						  p:(double)p
						  r:(double [NA])r
{
	double totCa, totFe2, totMg;

	totCa  = 2.0*r[1];
	totFe2 = 7.0*(r[0]+1.0)/2.0;
	totMg  = (1.0 - (r[0]+1.0)/2.0 - r[1])*7.0 +  5.0*r[1];

	return [NSString stringWithFormat:@"Ca%4.2fFe%4.2fMg%4.2fSi8O22(OH)2", totCa, totFe2, totMg];
}

-(void)activity:(int)mask
			  t:(double)t
			  p:(double)p
			  r:(double [NA])r
			  a:(double [NA])a      // (pointer to a[]) activities              BINARY MASK: 0001
			 mu:(double [NA])mu     // (pointer to mu[]) chemical potentials    BINARY MASK: 0010
			 dx:(double [NA][NR])dx // (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100
{
	double s[NS], g, dgdr[NR], fr[NA][NR];
	int i, j;

	for(i=0; i<NA; i++) {
		fr[i][0] = FR(i); /* R */
		fr[i][1] = FP(i); /* P */
	}

	[self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

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
		double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR],
		dfrdr[NA][NR], gs[NA][NS], dgsds[NA][NS], sum;
		int k, l;

		fillD2GDR2
		fillD2GDRDS
		fillD2GDS2

		for(i=0; i<NA; i++) {
			gs[i][0] = FS13(i);         /* s13 */
			gs[i][1] = FS2(i);          /* s2  */
			dfrdr[i][0] = DFRDR(i);     /* R   */
			dfrdr[i][1] = DFPDP(i);     /* P   */
			dgsds[i][0] = DFS13DS13(i); /* s13 */
			dgsds[i][1] = DFS2DS2(i);   /* s2  */
		}

		[self order:SECOND t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

		for (i=0; i<NA; i++) {
			for (k=0; k<NR; k++) {
				/* compute activity of the i-th component */
				for (dx[i][k]=g, j=0; j<NR; j++) dx[i][k] += fr[i][j]*dgdr[j];
				dx[i][k] = exp(dx[i][k]/(R*t));

				/* compute derivative of i-th activity with respect to r(k) */
				sum = (1.0+dfrdr[i][k])*dgdr[k];
				for (j=0; j<NR; j++) {
					sum += fr[i][j]*d2gdr2[j][k];
					for (l=0; l<NS; l++) sum += fr[i][j]*d2gdrds[j][l]*dsdr[l][k];
				}
				for (j=0; j<NS; j++) {
					sum += gs[i][j]*d2gdrds[k][j];
					for (l=0; l<NS; l++) sum += gs[i][j]*d2gds2[j][l]*dsdr[l][k];
				}
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
	double s[NS];

	[self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

	if (mask & FIRST) {
		*gmix = G;
	}

	if(mask & SECOND) {
		dx[0] = DGDR0;
		dx[1] = DGDR1;
	}

	if(mask & THIRD) {
		double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR];
		int i, j, k, l;

		fillD2GDR2
		fillD2GDRDS
		fillD2GDS2

		[self order:SECOND t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) {
				dx2[i][j] = d2gdr2[i][j];
				for (k=0; k<NS; k++) {
					dx2[i][j] += d2gdrds[i][k]*dsdr[k][j] + d2gdrds[j][k]*dsdr[k][i];
					for (l=0; l<NS; l++) dx2[i][j] += d2gds2[k][l]*dsdr[k][i]*dsdr[l][j];
				}
			}
		}

	}

	if(mask & FOURTH) {
		double d3gdr3[NR][NR][NR], d3gdr2ds[NR][NR][NS], d3gdrds2[NR][NS][NS];
		double d3gds3[NS][NS][NS], dsdr[NS][NR];
		int i, j, k, l, m, n;

		fillD3GDR3
		fillD3GDR2DS
		fillD3GDRDS2
		fillD3GDS3

		[self order:SECOND t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) {
				for (k=0; k<NR; k++) {
					dx3[i][j][k] = d3gdr3[i][j][k];
					for (l=0; l<NS; l++) {
						dx3[i][j][k] += d3gdr2ds[i][j][l]*dsdr[l][k] +
						d3gdr2ds[j][k][l]*dsdr[l][i] + d3gdr2ds[k][i][l]*dsdr[l][j];
						for (m=0; m<NS; m++) {
							dx3[i][j][k] +=
							d3gdrds2[i][l][m]*dsdr[l][j]*dsdr[m][k] +
							d3gdrds2[j][l][m]*dsdr[l][k]*dsdr[m][i] +
							d3gdrds2[k][l][m]*dsdr[l][i]*dsdr[m][j];
							for (n=0; n<NS; n++)
								dx3[i][j][k] +=
								d3gds3[l][m][n]*dsdr[l][i]*dsdr[m][j]*dsdr[n][k];
						}
					}
				}
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
	double s[NS];

	[self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

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
	double s[NS];

	[self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

	if (mask & FIRST) {
		*smix = S;
	}

	if(mask & SECOND) {
		double d2gdrds[NR][NS], d2gdrdt[NR], d2gds2[NS][NS], d2gdsdt[NS],
		dsdr[NS][NR], dsdt[NS];
		int i, k, l;

		fillD2GDRDS
		fillD2GDRDT
		fillD2GDS2
		fillD2GDSDT

		[self order:SECOND | THIRD t:t p:p r:r s:NULL dr:dsdr dt:dsdt dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

		for (i=0; i<NR; i++) {
			dx[i] = d2gdrdt[i];
			for (k=0; k<NS; k++) {
				dx[i] += d2gdrds[i][k]*dsdt[k] + d2gdsdt[k]*dsdr[k][i];
				for (l=0; l<NS; l++) dx[i] += d2gds2[k][l]*dsdt[k]*dsdr[l][i] ;
			}
			dx[i] *= -1.0;
		}
	}

	if(mask & THIRD) {
		double d2gdrds[NR][NS], d2gds2[NS][NS], d2gdsdt[NS], d3gdr2ds[NR][NR][NS],
		d3gdr2dt[NR][NR], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
		d3gds3[NS][NS][NS], d3gds2dt[NS][NS], dsdr[NS][NR], dsdt[NS],
		d2sdr2[NS][NR][NR], d2sdrdt[NS][NR];
		int i, j, k, l, m;

		fillD2GDRDS
		fillD2GDS2
		fillD2GDSDT
		fillD3GDR2DS
		fillD3GDR2DT
		fillD3GDRDS2
		fillD3GDRDSDT
		fillD3GDS3
		fillD3GDS2DT

		[self order:SECOND | THIRD | FIFTH | SIXTH t:t p:p r:r s:NULL dr:dsdr dt:dsdt dp:NULL dr2:d2sdr2 drt:d2sdrdt drp:NULL dt2:NULL dtp:NULL dp2:NULL];

		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) {
				dx2[i][j] = d3gdr2dt[i][j];
				for (k=0; k<NS; k++) {
					dx2[i][j] += d3gdr2ds[i][j][k]*dsdt[k]
					+ d3gdrdsdt[i][k]*dsdr[k][j]
					+ d3gdrdsdt[j][k]*dsdr[k][i]
					+ d2gdsdt[k]*d2sdr2[k][i][j]
					+ d2gdrds[i][k]*d2sdrdt[k][j]
					+ d2gdrds[j][k]*d2sdrdt[k][i];
					for (l=0; l<NS; l++) {
						dx2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdt[l]
						+ d3gdrds2[j][k][l]*dsdr[k][i]*dsdt[l]
						+ d2gds2[k][l]*d2sdr2[k][i][j]*dsdt[l]
						+ d3gds2dt[k][l]*dsdr[k][i]*dsdr[l][j]
						+ d2gds2[k][l]*dsdr[k][i]*d2sdrdt[l][j]
						+ d2gds2[k][l]*dsdr[l][j]*d2sdrdt[k][i];
						for (m=0; m<NS; m++)
							dx2[i][j] += d3gds3[k][l][m]*dsdr[k][i]*dsdr[l][j]*dsdt[m];
					}
				}
				dx2[i][j] *= -1.0;
			}
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
	double s[NS], dsdt[NS], d2gdsdt[NS], d2gds2[NS][NS], d2gdt2;
	int i, j;

	[self order:FIRST | THIRD t:t p:p r:r s:s dr:NULL dt:dsdt dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

	fillD2GDS2
	fillD2GDSDT
	d2gdt2 = D2GDT2;

	if (mask & FIRST) {
		*cpmix = d2gdt2;
		for (i=0; i<NS; i++) {
			*cpmix += 2.0*d2gdsdt[i]*dsdt[i];
			for (j=0; j<NS; j++) *cpmix += d2gds2[i][j]*dsdt[i]*dsdt[j];
		}
		*cpmix *= -t;
	}

	if(mask & SECOND) {
		double d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gdsdt2[NS], d2sdt2[NS],
		temp;
		double d3gdt3 = D3GDT3;
		int k;

		fillD3GDS3
		fillD3GDS2DT
		fillD3GDSDT2

		[self order:EIGHTH t:t p:p r:r s:NULL dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:d2sdt2 dtp:NULL dp2:NULL];

		/* compute d2gdt2 */
		temp = d2gdt2;
		for (i=0; i<NS; i++) {
			temp += 2.0*d2gdsdt[i]*dsdt[i];
			for (j=0; j<NS; j++) temp += d2gds2[i][j]*dsdt[i]*dsdt[j];
		}

		*dt = d3gdt3;
		for (i=0; i<NS; i++) {
			*dt += 3.0*d3gdsdt2[i]*dsdt[i] + 3.0*d2gdsdt[i]*d2sdt2[i];
			for (j=0; j<NS; j++) {
				*dt += 3.0*d2gds2[i][j]*dsdt[i]*d2sdt2[j]
				+ 3.0*d3gds2dt[i][j]*dsdt[i]*dsdt[j];
				for (k=0; k<NS; k++) *dt += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdt[k];
			}
		}
		*dt = -t*(*dt) - temp;
	}

	if(mask & THIRD) {
		double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
		d3gds2dt[NS][NS], d2gdrds[NR][NS], d3gdrdt2[NR], d3gdsdt2[NS],
		dsdr[NS][NR], d2sdrdt[NS][NR], d2sdt2[NS];
		int k, l;

		fillD2GDRDS
		fillD3GDRDS2
		fillD3GDRDSDT
		fillD3GDRDT2
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDSDT2

		[self order:SECOND | SIXTH | EIGHTH t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:NULL dr2:NULL drt:d2sdrdt drp:NULL dt2:d2sdt2 dtp:NULL dp2:NULL];

		for (i=0; i<NR; i++) {
			for (j=0,dx[i]=d3gdrdt2[i]; j<NS; j++) {
				dx[i] += d3gdsdt2[j]*dsdr[j][i] + 2.0*d2gdsdt[j]*d2sdrdt[j][i] +
				2.0*d3gdrdsdt[i][j]*dsdt[j] + d2gdrds[i][j]*d2sdt2[j];
				for (k=0; k<NS; k++) {
					dx[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdt[k] +
					2.0*d2gds2[j][k]*dsdt[j]*d2sdrdt[k][i] +
					2.0*d3gds2dt[j][k]*dsdr[j][i]*dsdt[k] +
					d2gds2[j][k]*dsdr[j][i]*d2sdt2[k];
					for (l=0; l<NS; l++)
						dx[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdt[k]*dsdt[l];
				}
			}
			dx[i] *= -t;
		}
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
	double s[NS];

	[self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

	if (mask & FIRST) {
		*vmix = DGDP;
	}

	if(mask & SECOND) {
		double d2gdrds[NR][NS], d2gdrdp[NR], d2gds2[NS][NS], d2gdsdp[NS],
		dsdr[NS][NR], dsdp[NS];
		int i, j, k;

		fillD2GDRDS
		fillD2GDRDP
		fillD2GDS2
		fillD2GDSDP

		[self order:SECOND | FOURTH t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

		for (i=0; i<NR; i++) {
			dx[i] = d2gdrdp[i];
			for (j=0; j<NS; j++) {
				dx[i] += d2gdrds[i][j]*dsdp[j] + d2gdsdp[j]*dsdr[j][i];
				for (k=0; k<NS; k++) dx[i] += d2gds2[j][k]*dsdp[j]*dsdr[k][i];
			}
		}
	}

	if(mask & THIRD) {
		double d2gdrds[NR][NS], d2gds2[NS][NS], d2gdsdp[NS], d3gdr2ds[NR][NR][NS],
		d3gdr2dp[NR][NR], d3gdrds2[NR][NS][NS], d3gdrdsdp[NR][NS],
		d3gds3[NS][NS][NS], d3gds2dp[NS][NS], dsdr[NS][NR], dsdp[NS],
		d2sdr2[NS][NR][NR], d2sdrdp[NS][NR];
		int i, j, k, l, m;

		fillD2GDRDS
		fillD2GDS2
		fillD2GDSDP
		fillD3GDR2DS
		fillD3GDR2DP
		fillD3GDRDS2
		fillD3GDRDSDP
		fillD3GDS3
		fillD3GDS2DP

		[self order:SECOND | FOURTH | FIFTH | SEVENTH t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:dsdp dr2:d2sdr2 drt:NULL drp:d2sdrdp dt2:NULL dtp:NULL dp2:NULL];

		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) {
				dx2[i][j] = d3gdr2dp[i][j];
				for (k=0; k<NS; k++) {
					dx2[i][j] += d3gdr2ds[i][j][k]*dsdp[k]
					+ d3gdrdsdp[i][k]*dsdr[k][j]
					+ d3gdrdsdp[j][k]*dsdr[k][i]
					+ d2gdsdp[k]*d2sdr2[k][i][j]
					+ d2gdrds[i][k]*d2sdrdp[k][j]
					+ d2gdrds[j][k]*d2sdrdp[k][i];
					for (l=0; l<NS; l++) {
						dx2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdp[l]
						+ d3gdrds2[j][k][l]*dsdr[k][i]*dsdp[l]
						+ d2gds2[k][l]*d2sdr2[k][i][j]*dsdp[l]
						+ d3gds2dp[k][l]*dsdr[k][i]*dsdr[l][j]
						+ d2gds2[k][l]*dsdr[k][i]*d2sdrdp[l][j]
						+ d2gds2[k][l]*dsdr[k][j]*d2sdrdp[l][i];
						for (m=0; m<NS; m++)
							dx2[i][j] += d3gds3[k][l][m]*dsdr[k][i]*dsdr[l][j]*dsdp[m];
					}
				}
			}
		}

	}

	if(mask & FOURTH) {
		double d2gdtdp = D2GDTDP;
		double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], dsdt[NS], dsdp[NS];
		int i, j;

		fillD2GDS2
		fillD2GDSDT
		fillD2GDSDP

		[self order:THIRD | FOURTH t:t p:p r:r s:NULL dr:NULL dt:dsdt dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

		*dt = d2gdtdp;
		for (i=0; i<NS; i++) {
			*dt += d2gdsdt[i]*dsdp[i] + d2gdsdp[i]*dsdt[i];
			for (j=0; j<NS; j++) *dt += d2gds2[i][j]*dsdt[i]*dsdp[j];
		}
	}

	if(mask & FIFTH) {
		double d2gdp2 = D2GDP2;
		double d2gds2[NS][NS], d2gdsdp[NS], dsdp[NS];
		int i,j;

		fillD2GDS2
		fillD2GDSDP

		[self order:FOURTH t:t p:p r:r s:NULL dr:NULL dt:NULL dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

		*dp = d2gdp2;
		for (i=0; i<NS; i++) {
			*dp += 2.0*d2gdsdp[i]*dsdp[i];
			for (j=0; j<NS; j++) *dp += d2gds2[i][j]*dsdp[i]*dsdp[j];
		}
	}

	if(mask & SIXTH) {
		double d3gdt2dp = D3GDT2DP;
		double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS],
		d3gds2dp[NS][NS], d3gds2dt[NS][NS], d3gdsdtdp[NS], d3gdsdt2[NS],
		dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS];
		int i, j, k;

		fillD2GDS2
		fillD2GDSDT
		fillD2GDSDP
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDS2DP
		fillD3GDSDT2
		fillD3GDSDTDP

		[self order:THIRD | FOURTH | EIGHTH | NINTH t:t p:p r:r s:NULL dr:NULL dt:dsdt dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:d2sdt2 dtp:d2sdtdp dp2:NULL];

		*dt2 = d3gdt2dp;
		for (i=0; i<NS; i++) {
			*dt2 += d3gdsdt2[i]*dsdp[i] + 2.0*d2gdsdt[i]*d2sdtdp[i]
            + d2gdsdp[i]*d2sdt2[i] + 2.0*d3gdsdtdp[i]*dsdt[i];
			for (j=0; j<NS; j++) {
				*dt2 += 2.0*d3gds2dt[i][j]*dsdt[i]*dsdp[j]
				+ d2gds2[i][j]*d2sdt2[i]*dsdp[j]
				+ 2.0*d2gds2[i][j]*dsdt[i]*d2sdtdp[j]
				+ d3gds2dp[i][j]*dsdt[i]*dsdt[j];
				for (k=0; k<NS; k++) *dt2 += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdp[k];
			}
		}
	}

	if(mask & SEVENTH) {
		double d3gdtdp2 = D3GDTDP2;
		double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS],
		d3gds2dt[NS][NS], d3gds2dp[NS][NS], d3gdsdtdp[NS], d3gdsdp2[NS],
		dsdt[NS], dsdp[NS], d2sdtdp[NS], d2sdp2[NS];
		int i, j, k;

		fillD2GDS2
		fillD2GDSDT
		fillD2GDSDP
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDS2DP
		fillD3GDSDTDP
		fillD3GDSDP2

		[self order:THIRD | FOURTH | NINTH | TENTH t:t p:p r:r s:NULL dr:NULL dt:dsdt dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:d2sdtdp dp2:d2sdp2];

		*dtdp = d3gdtdp2;
		for (i=0; i<NS; i++) {
			*dtdp += 2.0*d3gdsdtdp[i]*dsdp[i] + d2gdsdt[i]*d2sdp2[i]
			+ 2.0*d2gdsdp[i]*d2sdtdp[i] + d3gdsdp2[i]*dsdt[i];
			for (j=0; j<NS; j++) {
				*dtdp += 2.0*d3gds2dp[i][j]*dsdt[i]*dsdp[j]
				+ d2gds2[i][j]*dsdt[i]*d2sdp2[j]
				+ 2.0*d2gds2[i][j]*d2sdtdp[i]*dsdp[j]
				+ d3gds2dt[i][j]*dsdp[i]*dsdp[j];
				for (k=0; k<NS; k++) *dtdp += d3gds3[i][j][k]*dsdt[i]*dsdp[j]*dsdp[k];
			}
		}
	}

	if(mask & EIGHTH) {
		double d3gdp3 = D3GDP3;
		double d2gds2[NS][NS], d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS],
		d3gdsdp2[NS], dsdp[NS], d2sdp2[NS];
		int i, j, k;

		fillD2GDS2
		fillD2GDSDP
		fillD3GDS3
		fillD3GDS2DP
		fillD3GDSDP2

		[self order:FOURTH | TENTH t:t p:p r:r s:NULL dr:NULL dt:NULL dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:d2sdp2];

		*dp2 = d3gdp3;
		for (i=0; i<NS; i++) {
			*dp2 += 3.0*d3gdsdp2[i]*dsdp[i] + 3.0*d2gdsdp[i]*d2sdp2[i];
			for (j=0; j<NS; j++) {
				*dp2 += 3.0*d2gds2[i][j]*dsdp[i]*d2sdp2[j]
				+ 3.0*d3gds2dp[i][j]*dsdp[i]*dsdp[j];
				for (k=0; k<NS; k++) *dp2 += d3gds3[i][j][k]*dsdp[i]*dsdp[j]*dsdp[k];
			}
		}
	}

	if(mask & NINTH) {
		double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
		d3gds2dp[NS][NS], d2gdrds[NR][NS], d3gdrdtdp[NR], d3gdsdtdp[NS],
		dsdt[NS], dsdp[NS], dsdr[NS][NR], d2sdrdt[NS][NR], d2sdrdp[NS][NR],
		d2gds2[NS][NS], d2gdsdt[NS], d3gdrdsdp[NR][NS], d2gdsdp[NS],
		d2sdtdp[NS], d3gds2dt[NS][NS];
		int i, j, k, l;

		fillD2GDRDS
		fillD2GDS2
		fillD2GDSDT
		fillD2GDSDP
		fillD3GDRDS2
		fillD3GDRDSDT
		fillD3GDRDSDP
		fillD3GDRDTDP
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDSDTDP
		fillD3GDS2DP

		[self order:SECOND | THIRD | FOURTH | SIXTH | SEVENTH | NINTH t:t p:p r:r s:NULL dr:dsdr dt:dsdt dp:dsdp dr2:NULL drt:d2sdrdt drp:d2sdrdp dt2:NULL dtp:d2sdtdp dp2:NULL];

		for (i=0; i<NR; i++) {
			for (j=0,dxdt[i]=d3gdrdtdp[i]; j<NS; j++) {
				dxdt[i] += d3gdsdtdp[j]*dsdr[j][i] + d2gdsdt[j]*d2sdrdp[j][i] +
				d3gdrdsdt[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdtdp[j] +
				d3gdrdsdp[i][j]*dsdt[j] + d2gdsdp[j]*d2sdrdt[j][i];
				for (k=0; k<NS; k++) {
					dxdt[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdp[k] +
					d2gds2[j][k]*dsdt[j]*d2sdrdp[k][i] +
					d2gds2[j][k]*dsdp[j]*d2sdrdt[k][i] +
					d3gds2dt[j][k]*dsdr[j][i]*dsdp[k] +
					d3gds2dp[j][k]*dsdr[j][i]*dsdt[k] +
					d2gds2[j][k]*dsdr[j][i]*d2sdtdp[k];
					for (l=0; l<NS; l++)
						dxdt[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdt[k]*dsdp[l];
				}
			}
		}
	}

	if(mask & TENTH) {
		double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gds2dp[NS][NS],
		d2gdrds[NR][NS], dsdp[NS], dsdr[NS][NR], d2sdrdp[NS][NR], d2gds2[NS][NS],
		d3gdrdsdp[NR][NS], d3gdrdp2[NR], d3gdsdp2[NS], d2gdsdp[NS], d2sdp2[NS];
		int i, j, k, l;

		fillD2GDRDS
		fillD2GDS2
		fillD2GDSDP
		fillD3GDRDS2
		fillD3GDRDSDP
		fillD3GDRDP2
		fillD3GDS3
		fillD3GDS2DP
		fillD3GDSDP2

		[self order:SECOND | FOURTH | SEVENTH | TENTH t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:dsdp dr2:NULL drt:NULL drp:d2sdrdp dt2:NULL dtp:NULL dp2:d2sdp2];

		for (i=0; i<NR; i++) {
			for (j=0,dxdp[i]=d3gdrdp2[i]; j<NS; j++) {
				dxdp[i] += d3gdsdp2[j]*dsdr[j][i] + d2gdsdp[j]*d2sdrdp[j][i] +
				2.0*d3gdrdsdp[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdp2[j] +
				d2gdsdp[j]*d2sdrdp[j][i];
				for (k=0; k<NS; k++) {
					dxdp[i] += d3gdrds2[i][j][k]*dsdp[j]*dsdp[k] +
					2.0*d2gds2[j][k]*dsdp[j]*d2sdrdp[k][i] +
					2.0*d3gds2dp[j][k]*dsdr[j][i]*dsdp[k] +
					d2gds2[j][k]*dsdr[j][i]*d2sdp2[k];
					for (l=0; l<NS; l++)
						dxdp[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdp[k]*dsdp[l];
				}
			}
		}
	}

}

#define NAS 4

-(NSUInteger)numberOfSolutionSpecies {
	return NAS;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
	switch (index) {
		case 0:
			return @"cummingtonite";
			break;
		case 1:
			return @"grunnerite";
			break;
		case 2:
			return @"tremolite";
			break;
		case 3:
			return @"ferroactinolite";
			break;
		default:
			return @"";
			break;
	}
}

-(DoubleVector *)convertMolesOfSpeciesToMolesOfComponents:(double *)mSpecies {
	DoubleVector *mComponentWrapper =[[DoubleVector alloc] initWithSize:NA];
	double *mComponents = [mComponentWrapper pointerToDouble];
	mComponents[0] = mSpecies[0] - 5.0*mSpecies[3]/7.0; // Mg7Si8O22(OH)2
	mComponents[1] = mSpecies[1] + 5.0*mSpecies[3]/7.0; // Fe7Si8O22(OH)2
	mComponents[2] = mSpecies[2] + mSpecies[3];         // Ca2Mg5Si8O22(OH)2
	return mComponentWrapper;
}

-(DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *muSpeciesWrapper = [[DoubleVector alloc] initWithSize:NAS];
	double *muSpecies = [muSpeciesWrapper pointerToDouble];
	double *muComponents = [[self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p] pointerToDouble];
	for (NSUInteger i=0; i<NA; i++) muSpecies[i] = muComponents[i];
	muSpecies[3] = muSpecies[2] + 5.0*muSpecies[1]/7.0 - 5.0*muSpecies[0]/7.0; // Ca2Fe5Si8O22(OH)2
	return muSpeciesWrapper;
}

-(DoubleVector *)elementalCompositionOfSpeciesAtIndex:(NSUInteger)index {
	DoubleVector *speciesElementArrayWrapper = nil;

	if (index < NA) speciesElementArrayWrapper = [[endmembers objectAtIndex:index] formulaAsElementArray];
	else if (index == 3) {
        speciesElementArrayWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
		double *speciesElementArray = [speciesElementArrayWrapper pointerToDouble];
        double *component = [[[endmembers objectAtIndex:2] formulaAsElementArray] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] += component[i];
		component = [[[endmembers objectAtIndex:1] formulaAsElementArray] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] += 5.0*component[i]/7.0;
		component = [[[endmembers objectAtIndex:0] formulaAsElementArray] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] -= 5.0*component[i]/7.0;
	}
	return speciesElementArrayWrapper;
}

-(void)correctActivityCoefficients:(double [NAS])gamma forComposition:(double [NAS])x { }

// --> SolutionPhaseProtocol public function
-(NSArray *)affinityAndCompositionFromLiquidChemicalPotentialSum:(double *)chemicalPotentials andT:(double)t andP:(double)p {
	NSMutableArray *results = [NSMutableArray arrayWithCapacity:NA+1];
	double mu0[NAS], deltaMu[NAS], xNz[NAS], x[NAS], gamma[NAS], gammaLast[NAS], xLast[NAS], affinity = 0.0;
	// double gammaLast[NAS];
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
			//NSLog(@"sum mu liq %g, mu0 %g, delta %g", chemicalPotentials[i], mu0[i], deltaMu[nz]);
			index[nz] = i;
			gamma[nz] = 1.0;
			nz++;
		} else mu0[i] = 0.0;
		x[i] = 0.0;
		xLast[i] = 0.0;
	}

	// Ca2Fe5Si8O22(OH)2
	if ((chemicalPotentials[0] != 0.0) && (chemicalPotentials[1] != 0.0) && (chemicalPotentials[2] != 0.0)) {
		mu0[3] = mu0[2] + 5.0*mu0[1]/7.0 - 5.0*mu0[0]/7.0 + (GXCA_13) + 2.0*(GXCA_2) + (GX123_4)/2.0 - (GEX);
		deltaMu[nz] = chemicalPotentials[2] + 5.0*chemicalPotentials[1]/7.0 - 5.0*chemicalPotentials[0]/7.0 - mu0[3];
		//NSLog(@"sum mu liq %g, mu0 %g, delta %g", chemicalPotentials[2] + 5.0*chemicalPotentials[1]/7.0 - 5.0*chemicalPotentials[0]/7.0, mu0[3], deltaMu[nz]);
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
				xNz[i] = exp(((deltaMu[i]-R*t*5.0*log(gamma[i]))-(deltaMu[nz-2]-R*t*5.0*log(gamma[nz-2])))/(5.0*R*t));
				sum += xNz[i];
			}
			xNz[nz-2] = exp(((deltaMu[nz-2]-R*t*5.0*log(gamma[nz-2]))-(deltaMu[nz-1]-R*t*5.0*log(gamma[nz-1])))/(5.0*R*t));

			xNz[nz-2] /= 1.0 + xNz[nz-2]*sum;
			xNz[nz-1] = 1.0 - xNz[nz-2];
			if (nz > 2) for (i=0; i<(nz-2); i++) {
				xNz[i] *= xNz[nz-2];
				xNz[nz-1] -= xNz[i];
			}

			for (i=0; i<nz; i++) if (xNz[i] <= DBL_EPSILON) xNz[i] = DBL_EPSILON;

			/* compute the chemical affinity (choice of mu[] is arbitrary) */
			affinity = -(deltaMu[0]-R*t*5.0*log(gamma[0])) + R*t*5.0*log(xNz[0]);
		}

		// Reinflate the solution
		for (i=0; i<nz; i++) x[index[i]] = xNz[i];

		// Determine activity coefficients
		double a[NAS], mu[NA], r[NR];
		xReduced[0] = x[0] - 5.0*x[3]/7.0; // Mg7Si8O22(OH)2
		xReduced[1] = x[1] + 5.0*x[3]/7.0; // Fe7Si8O22(OH)2
		xReduced[2] = x[2] + x[3];         // Ca2Mg5Si8O22(OH)2
		[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:xReduced r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
		[self activity:FIRST | SECOND t:t p:p r:r a:a mu:mu dx:NULL];

		for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) gamma[j++] = pow(a[i], 1.0/5.0)/x[i];

		if ((chemicalPotentials[0] != 0.0) && (chemicalPotentials[1] != 0.0) && (chemicalPotentials[2] != 0.0)) {
			a[3] = exp((mu[2] - 5.0*mu[0]/7.0 + 5.0*mu[1]/7.0)/(R*t));
			gamma[j] = pow(a[3], 1.0/5.0)/x[3];
		}

		if (debugV) {
			NSLog(@"Iteration %lu", count);
			NSLog(@"%10.3g %10.3g %10.3g %10.3g %13.6g", a[0], a[1], a[2], a[3], affinity);
			NSLog(@"%10.3g %10.3g %10.3g %10.3g", x[0], x[1], x[2], x[3]);
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

	} while (count < 75 && !converged);

	if (debugS) {
		NSLog(@"... Terminated (converged %@) for phase %@ in %lu iterations with delta affinity %f J for %f atoms.",
			  converged ? @"YES" : @"NO", [self phaseName], count, fabs(affinity-affinityLast), NATOMS);
		for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) NSLog(@"... ... Activity coefficient of component %@ is %13.6g with mole fraction %13.6g",
														 [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0], gamma[j++], x[i]);
		NSLog(@"... ... Activity coefficient of component %@ is %13.6g with mole fraction %13.6g", @"ferroactinolite     ", gamma[j], x[3]);
	}
	if (debugV) NSLog(@"Exiting [... affinityAndCompositionFromLiquidChemicalPotentialSum].");

	double CaPlusNaOnM4 = xReduced[2];

	if      ((CaPlusNaOnM4 < 0.2) && [[self phaseName] isEqualToString:@"Actinolite"])    affinity = 0.0;
	else if ((CaPlusNaOnM4 > 0.2) && [[self phaseName] isEqualToString:@"Cummingtonite"]) affinity = 0.0;

	[results addObject:[NSNumber numberWithDouble:affinity]];                         // affinity in J
	for (i=0; i<NA; i++) [results addObject:[NSNumber numberWithDouble:xReduced[i]]]; // composition in mole fraction of endmembers
	[results addObject:[NSNumber numberWithBool:converged]];                          // convergence flag
	[results addObject:[NSNumber numberWithUnsignedInteger:count]];                       // iteration count
	[results addObject:[NSNumber numberWithDouble:NATOMS]];                           // number of atoms used to scale affinity
	[results addObject:[NSNumber numberWithDouble:fabs(affinity-affinityLast)]];      // likely error in affinity

	if (debugV) NSLog(@"Exiting [... affinityAndCompositionFromLiquidChemicalPotentialSum].");
	return [NSArray arrayWithArray:results];
}

// --> SolutionPhaseProtocol public function

-(NSDictionary *)checkForAndDetermineCompositionOfCoexistingImmisciblePhase:(double *)refMoles andT:(double)t andP:(double)p {
	BOOL debugS = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
	BOOL debugV = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
	if (debugV) NSLog(@"Entering [... checkForAndDetermineCompositionOfCoexistingImmisciblePhase] ...");

	NSMutableDictionary *results = [NSMutableDictionary dictionaryWithCapacity:6];

	if ((refMoles[1] == 0.0) ||
		([[self phaseName] isEqualToString:@"Clinoamphibole"]) ||
		(([[self phaseName] isEqualToString:@"Actinolite"]) &&
		 ([instanceSet countForObject:[@"Cummingtonite" stringByAppendingString:[self operationParent]]] > 0)) ||
		(([[self phaseName] isEqualToString:@"Cummingtonite"]) &&
		 ([instanceSet countForObject:[@"Actinolite" stringByAppendingString:[self operationParent]]] > 0))
		) {

		NSMutableArray *composition = [NSMutableArray arrayWithCapacity:NA];
		for (NSUInteger i=0; i<NA; i++) [composition addObject:[NSNumber numberWithDouble:refMoles[i]]];
		[results setObject:[NSNumber numberWithBool:NO] forKey:@"additionalPhaseDetected"];
		[results setObject:@""forKey:@"nameOfCoexistingPhase"];
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

	// Specific to clinoamphibole ...
	moles[0] = refMoles[0]/refTotalMoles; // Mg7Si8O22(OH)2
	moles[1] = refMoles[1]/refTotalMoles; // Fe7Si8O22(OH)2
	moles[2] = refMoles[2]/refTotalMoles; // Ca2Mg5Si8O22(OH)2

	double FeOverMg = ((7.0*moles[0]+5.0*moles[2]) != 0.0) ? 7.0*moles[1]/(7.0*moles[0]+5.0*moles[2]) : 10000.0;

	BOOL initialGuessOK = NO;
	if ([[self phaseName] isEqualToString:@"Actinolite"]) moles[2] = 0.95;
	else moles[2] = 0.05;
	while (!initialGuessOK) {
		moles[1] = (1.0 - 2.0*moles[2]/7.0)*FeOverMg/(1.0+FeOverMg);
		moles[0] = 1.0 - moles[1] - moles[2];
		initialGuessOK = [self testPermissibleValuesOfComponents:moles];
		if (!initialGuessOK) {
			if ([[self phaseName] isEqualToString:@"Actinolite"]) moles[2] -= 0.05;
			else moles[2] += 0.05;
		}
	}

	double totalMoles = 1.0;

	if (debugS) {
		NSLog(@"... d = %d, r = %d, rNorm = %10.3e Cm = %10.3e Gn = %10.3e Tr = %10.3e Affinity/RT = %10.3e Cm = %10.3e Gn = %10.3e Tr = %10.3e",
			  NA, NA, 0.0, 0.0, 0.0, 0.0, 0.0, refMoles[0]/refTotalMoles, refMoles[1]/refTotalMoles, refMoles[2]/refTotalMoles);
		NSLog(@"... d = %d, r = %d, rNorm = %10.3e Cm = %10.3e Gn = %10.3e Tr = %10.3e Affinity/RT = %10.3e Cm = %10.3e Gn = %10.3e Tr = %10.3e",
			  NA, NA, 0.0, 0.0, 0.0, 0.0, 0.0, moles[0], moles[1], moles[2]);
	}
	// ... end

	NSUInteger iters = 0;
	BOOL converged = NO;
	double tolerance = sqrt(DBL_EPSILON);
	NSUInteger maxIters = 50;
	double affinityScaledByRT = 0.0;

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

		// Specific to clinoamphibole ...
		if (debugS) {
			NSLog(@"... d = %d, r = %ld, rNorm = %10.3e Cm = %10.3e Gn = %10.3e Tr = %10.3e Affinity/RT = %10.3e Cm = %10.3e Gn = %10.3e Tr = %10.3e sL = %10.3e",
				  NA, pseudorank, rNorm, dMatrix[0][0], dMatrix[1][0], dMatrix[2][0], affinityScaledByRT, moles[0], moles[1], moles[2], stepLength);
		}
		// ... end

		iters++;
	}

	if (debugS) {
		NSLog(@"... Soln: Reference composition is a %@", [self phaseName]);
		NSLog(@"... Reference formula: %@", [self getFormulaFromMolesOfComponents:refMoles andT:t andP:p]);
		NSLog(@"... There are %lu instances of Cummingtonite instantiated for this class",
			  [instanceSet countForObject:[@"Cummingtonite" stringByAppendingString:[self operationParent]]]);
		NSLog(@"... There are %lu instances of Actinolite    instantiated for this class",
			  [instanceSet countForObject:[@"Actinolite" stringByAppendingString:[self operationParent]]]);
		NSLog(@"... There are %lu instances of Clinoamphibole instantiated for this class",
			  [instanceSet countForObject:[@"Clinoamphibole" stringByAppendingString:[self operationParent]]]);
		NSLog(@"... Soln: affinity = %g", affinityScaledByRT*8.314472*t);
		NSLog(@"... Soln: formula = %@", [self getFormulaFromMolesOfComponents:moles andT:t andP:p]);
		for (NSUInteger i=0; i<NA; i++) NSLog(@"... Soln: %6.3f Ref: %6.3f delta: %13.6e", moles[i], refMoles[i]/refTotalMoles, moles[i]-refMoles[i]/refTotalMoles);
	}

	NSMutableArray *composition = [NSMutableArray arrayWithCapacity:NA];
	for (NSUInteger i=0; i<NA; i++) [composition addObject:[NSNumber numberWithDouble:moles[i]]];
	BOOL additionalPhaseDetected = (affinityScaledByRT > tolerance) ? YES : NO;
	if (!converged) additionalPhaseDetected = NO;

	NSString *nameOfCoexistingPhase = @"";
	if (additionalPhaseDetected) {
		if      ([[self phaseName] isEqualToString:@"Actinolite"])    nameOfCoexistingPhase = @"Cummingtonite";
		else if ([[self phaseName] isEqualToString:@"Cummingtonite"]) nameOfCoexistingPhase = @"Actinolite";
	}

	[results setObject:[NSNumber numberWithBool:additionalPhaseDetected] forKey:@"additionalPhaseDetected"];
	[results setObject:nameOfCoexistingPhase forKey:@"nameOfCoexistingPhase"];
	[results setObject:[NSNumber numberWithBool:converged] forKey:@"converged"];
	[results setObject:[NSNumber numberWithUnsignedInteger:iters] forKey:@"iterations"];
	[results setObject:[NSNumber numberWithDouble:affinityScaledByRT*8.314472*t] forKey:@"affinity"];
	[results setObject:[NSArray arrayWithArray:composition] forKey:@"composition"];

	return [NSDictionary dictionaryWithDictionary:results];
}

#import "SolutionPhase.h"

@end
