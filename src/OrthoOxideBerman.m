//
//  OrthoOxideBerman.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/8/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import "OrthoOxideBerman.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"

@implementation OrthoOxideBerman

static NSArray *endmembers;

#define NR        2    /* Two independent composition variables       */
#define NS        3    /* Three ordering parameters s3 = X3           */
#define NA        3    /* Three endmember compositions                */
#define NATOMS  8.0    /* Average number of atoms in the formula unit */

#pragma mark -
#pragma mark class methods

+(void)initialize {
	if (self == [OrthoOxideBerman class]) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"Initialize(MeliliteBerman) - entry ...");
		NSMutableArray *mutableEndmembers = [NSMutableArray arrayWithCapacity:NA];

		BermanProperties *pseudobrookite = [[BermanProperties alloc] initWithH:-1754429.0
																			 S:156.4816
																			k0:261.35
																			k1:-15.307e2
																			k2:0.0
																			k3:-23.466e7
																			v0:5.491
																			v1:0.0
																			v2:0.0
																			v3:0.0
																			v4:0.0];
		[pseudobrookite setPhaseFormula:@"Fe2TiO5"];
		[pseudobrookite setPhaseName:@"pseudobrookite"];
		[mutableEndmembers addObject:pseudobrookite];
		if (debug) NSLog(@"... allocated pseudobrookite ...");

		BermanProperties *ferropseudobrookite = [[BermanProperties alloc] initWithH:-2175332.0
																				  S:127.1936-26.94+59.80
																				 k0:232.58-58.196+77.036
																				 k1:-7.555e2-(-1.6114e2)-5.8471e2
																				 k2:-56.608e5-(-14.0458e5)
																				 k3:58.214e7-11.2673e7+0.5558e7
																				 v0:5.578
																				 v1:0.0
																				 v2:0.0
																				 v3:0.0
																				 v4:0.0];
		[ferropseudobrookite setPhaseFormula:@"FeTi2O5"];
		[ferropseudobrookite setPhaseName:@"ferropseudobrookite"];
		[mutableEndmembers addObject:ferropseudobrookite];
		if (debug) NSLog(@"... allocated ferropseudobrookite ...");

		BermanProperties *karrooite = [[BermanProperties alloc] initWithH:-2507053.0
																		S:127.1936
																	   k0:232.58
																	   k1:-7.555e2
																	   k2:-56.608e5
																	   k3:58.214e7
																	   v0:5.482
																	   v1:0.0
																	   v2:0.0
																	   v3:0.0
																	   v4:0.0];
		[karrooite setPhaseFormula:@"MgTi2O5"];
		[karrooite setPhaseName:@"karrooite"];
		[mutableEndmembers addObject:karrooite];
		if (debug) NSLog(@"... allocated karrooite ...");

		endmembers = [NSArray arrayWithArray:mutableEndmembers];
	}
}

#pragma mark -
#pragma mark instance methods

-(id)init {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(OrthoOxideBerman) ... entry ...");
		[self setPhaseName:@"OrthoOxide"];
		computeMixingQuantities = NO;
		tOld = -9999.0;
		pOld = -9999.0;
		tOldPure = -9999.0;
		pOldPure = -9999.0;
		for (NSUInteger i=0; i<NR; i++) rOld[i] = -9999.0;
		for (NSUInteger i=0; i<NS; i++) sOld[i] = 2.0;
		for (NSUInteger i=0; i<NS; i++) sOldPure[i] = 2.0;
		xmg2M1 = 0.0;
		xfe2M1 = 0.0;
		xti4M1 = 0.0;
		xfe3M1 = 0.0;
		xmg2M2 = 0.0;
		xfe2M2 = 0.0;
		xti4M2 = 0.0;
		xfe3M2 = 0.0;

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

#define MAX_ITER 200    /* Maximum number of iterations allowed in order */

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Guess: July 16, 1997 (see notes)
 */

#define R  8.3143

#define G11P   22000.0 /* joules */
#define G22P   22000.0 /* joules */
#define G33P   22000.0 /* joules */
#define W11P    6000.0 /* joules */
#define W12    26000.0 /* joules */
#define W12P   26000.0 /* joules */
#define W13    45000.0 /* joules */
#define W13P   45000.0 /* joules */
#define W1P2   26000.0 /* joules */
#define W1P2P  26000.0 /* joules */
#define W1P3   45000.0 /* joules */
#define W1P3P  45000.0 /* joules */
#define W22P    6000.0 /* joules */
#define W23     8400.0 /* joules */
#define W23P    8400.0 /* joules */
#define W2P3    8400.0 /* joules */
#define W2P3P   8400.0 /* joules */
#define W33P    6000.0 /* joules */

#define G0     (G11P)/2.0+(W11P)/4.0
#define GX2    ((G22P)-(G11P)-(W11P))/2.0 + ((W12)+(W12P)+(W1P2)+(W1P2P))/4.0
#define GX3    ((G33P)-(G11P)-(W11P))/2.0 + ((W13)+(W13P)+(W1P3)+(W1P3P))/4.0
#define GS1    -(G11P)/2.0
#define GS2    -(G22P)/2.0 + ((W12)-(W12P)+(W1P2)-(W1P2P))/4.0
#define GS3    -(G33P)/2.0 + ((W13)-(W13P)+(W1P3)-(W1P3P))/4.0
#define GX2X2  ((W11P)-(W12)-(W12P)-(W1P2)-(W1P2P)+(W22P))/4.0
#define GX2X3  (W11P)/2.0 + (-(W12)-(W12P)-(W13)-(W13P)-(W1P2)-(W1P2P) \
-(W1P3)-(W1P3P)+(W23)+(W23P)+(W2P3)+(W2P3P))/4.0
#define GX2S1  ((W12)+(W12P)-(W1P2)-(W1P2P))/4.0
#define GX2S2  (-(W12)+(W12P)-(W1P2)+(W1P2P))/4.0
#define GX2S3  (-(W13)+(W13P)-(W1P3)+(W1P3P)+(W23)-(W23P)+(W2P3)-(W2P3P))/4.0
#define GX3X3  ((W11P)-(W13)-(W13P)-(W1P3)-(W1P3P)+(W33P))/4.0
#define GX3S1  ((W13)+(W13P)-(W1P3)-(W1P3P))/4.0
#define GX3S2  (-(W12)+(W12P)-(W1P2)+(W1P2P)+(W23)+(W23P)-(W2P3)-(W2P3P))/4.0
#define GX3S3  (-(W13)+(W13P)-(W1P3)+(W1P3P))/4.0
#define GS1S1  -(W11P)/4.0
#define GS1S2  ((W12)-(W12P)-(W1P2)+(W1P2P))/4.0
#define GS1S3  ((W13)-(W13P)-(W1P3)+(W1P3P))/4.0
#define GS2S2  -(W22P)/4.0
#define GS2S3  ((W23)-(W23P)-(W2P3)+(W2P3P))/4.0
#define GS3S3  -(W33P)/4.0

/*
 *=============================================================================
 * Local function to compute ordering state and associated derivatives of
 * pure component endmembers
 */

#define PB_S           - R*(  ((1.0+s[0])/2.0)*log((1.0+s[0])/2.0) \
+ ((1.0-s[0])/2.0)*log((1.0-s[0])/2.0) \
+ ((3.0-s[0])/2.0)*log((3.0-s[0])/4.0) \
+ ((1.0+s[0])/2.0)*log((1.0+s[0])/4.0) )
#define PB_H           (G0) + (GS1)*s[0] + (GS1S1)*s[0]*s[0]
#define PB_G           PB_H - t*(PB_S)

#define DPB_GDS0       R*t*(  log((1.0+s[0])/2.0)/2.0 - log((1.0-s[0])/2.0)/2.0 \
- log((3.0-s[0])/4.0)/2.0 + log((1.0+s[0])/4.0)/2.0 ) \
+ (GS1) + 2.0*(GS1S1)*s[0]
#define DPB_GDT        -(PB_S)
#define DPB_GDP        0.0

#define D2PB_GDS0S0    R*t*(  0.5/(1.0+s[0]) + 0.5/(1.0-s[0]) \
+ 0.5/(3.0-s[0]) + 0.5/(1.0+s[0]) ) + 2.0*(GS1S1)
#define D2PB_GDS0DT    R*(  log((1.0+s[0])/2.0)/2.0 - log((1.0-s[0])/2.0)/2.0 \
- log((3.0-s[0])/4.0)/2.0 + log((1.0+s[0])/4.0)/2.0 )
#define D2PB_GDS0DP    0.0
#define D2PB_GDT2      0.0
#define D2PB_GDTP      0.0
#define D2PB_GDP2      0.0

#define D3PB_GDS0S0S0  R*t*(- 0.5/SQUARE(1.0+s[0]) + 0.5/SQUARE(1.0-s[0]) \
+ 0.5/SQUARE(3.0-s[0]) - 0.5/SQUARE(1.0+s[0]) )
#define D3PB_GDS0S0DT  R*(  0.5/(1.0+s[0]) + 0.5/(1.0-s[0]) \
+ 0.5/(3.0-s[0]) + 0.5/(1.0+s[0]) )
#define D3PB_GDS0S0DP  0.0
#define D3PB_GDS0DT2   0.0
#define D3PB_GDS0DTDP  0.0
#define D3PB_GDS0DP2   0.0
#define D3PB_GDT3      0.0
#define D3PB_GDT2DP    0.0
#define D3PB_GDTDP2    0.0
#define D3PB_GDP3      0.0

#define FE_S           - R*(  ((1.0+s[1])/2.0)*log((1.0+s[1])/2.0) \
+ ((1.0-s[1])/2.0)*log((1.0-s[1])/4.0) )
#define FE_H           (GX2) + (GS2)*s[1] + (GX2X2) + (GX2S2)*s[1] + (GS2S2)*s[1]*s[1]
#define FE_G           FE_H - t*(FE_S)

#define DFE_GDS1       R*t*( log((1.0+s[1])/2.0)/2.0 - log((1.0-s[1])/4.0)/2.0 ) \
+ (GS2) + (GX2S2) + 2.0*(GS2S2)*s[1]
#define DFE_GDT        -(FE_S)
#define DFE_GDP        0.0

#define D2FE_GDS1S1    R*t*( 0.5/(1.0+s[1]) + 0.5/(1.0-s[1]) ) + 2.0*(GS2S2)
#define D2FE_GDS1DT    R*( log((1.0+s[1])/2.0)/2.0 - log((1.0-s[1])/4.0)/2.0 )
#define D2FE_GDS1DP    0.0
#define D2FE_GDT2      0.0
#define D2FE_GDTP      0.0
#define D2FE_GDP2      0.0

#define D3FE_GDS1S1S1  R*t*( -0.5/SQUARE(1.0+s[1]) + 0.5/SQUARE(1.0-s[1]) )
#define D3FE_GDS1S1DT  R*( 0.5/(1.0+s[1]) + 0.5/(1.0-s[1]) )
#define D3FE_GDS1S1DP  0.0
#define D3FE_GDS1DT2   0.0
#define D3FE_GDS1DTDP  0.0
#define D3FE_GDS1DP2   0.0
#define D3FE_GDT3      0.0
#define D3FE_GDT2DP    0.0
#define D3FE_GDTDP2    0.0
#define D3FE_GDP3      0.0

#define MG_S           - R*(  ((1.0+s[2])/2.0)*log((1.0+s[2])/2.0) \
+ ((1.0-s[2])/2.0)*log((1.0-s[2])/4.0) )
#define MG_H           (GX3) + (GS3)*s[2] + (GX3X3) + (GX3S3)*s[2] + (GS3S3)*s[2]*s[2]
#define MG_G           MG_H - t*(MG_S)

#define DMG_GDS2       R*t*( log((1.0+s[2])/2.0)/2.0 - log((1.0-s[2])/4.0)/2.0 ) \
+ (GS3) + (GX3S3) + 2.0*(GS3S3)*s[2]
#define DMG_GDT        -(MG_S)
#define DMG_GDP        0.0

#define D2MG_GDS2S2    R*t*( 0.5/(1.0+s[2]) + 0.5/(1.0-s[2]) ) + 2.0*(GS3S3)
#define D2MG_GDS2DT    R*( log((1.0+s[2])/2.0)/2.0 - log((1.0-s[2])/4.0)/2.0 )
#define D2MG_GDS2DP    0.0
#define D2MG_GDT2      0.0
#define D2MG_GDTP      0.0
#define D2MG_GDP2      0.0

#define D3MG_GDS2S2S2  R*t*( -0.5/SQUARE(1.0+s[2]) + 0.5/SQUARE(1.0-s[2]) )
#define D3MG_GDS2S2DT  R*( 0.5/(1.0+s[2]) + 0.5/(1.0-s[2]) )
#define D3MG_GDS2S2DP  0.0
#define D3MG_GDS2DT2   0.0
#define D3MG_GDS2DTDP  0.0
#define D3MG_GDS2DP2   0.0
#define D3MG_GDT3      0.0
#define D3MG_GDT2DP    0.0
#define D3MG_GDTDP2    0.0
#define D3MG_GDP3      0.0

#define fillD2GDSDT \
d2gdsdt[0] = D2PB_GDS0DT; d2gdsdt[1] = D2FE_GDS1DT; \
d2gdsdt[2] = D2MG_GDS2DT;

#define fillD2GDSDP \
d2gdsdp[0] = D2PB_GDS0DP; d2gdsdp[1] = D2FE_GDS1DP; \
d2gdsdp[2] = D2MG_GDS2DP;

#define fillD3GDS3 \
d3gds3[0] = D3PB_GDS0S0S0; d3gds3[1] = D3FE_GDS1S1S1; \
d3gds3[2] = D3MG_GDS2S2S2;

#define fillD3GDS2DT \
d3gds2dt[0] = D3PB_GDS0S0DT; d3gds2dt[1] = D3FE_GDS1S1DT; \
d3gds2dt[2] = D3MG_GDS2S2DT;

#define fillD3GDS2DP \
d3gds2dp[0] = D3PB_GDS0S0DP; d3gds2dp[1] = D3FE_GDS1S1DP; \
d3gds2dp[2] = D3MG_GDS2S2DP;

#define fillD3GDSDT2 \
d3gdsdt2[0] = D3PB_GDS0DT2; d3gdsdt2[1] = D3FE_GDS1DT2; \
d3gdsdt2[2] = D3MG_GDS2DT2;

#define fillD3GDSDTDP \
d3gdsdtdp[0] = D3PB_GDS0DTDP; d3gdsdtdp[1] = D3FE_GDS1DTDP; \
d3gdsdtdp[2] = D3MG_GDS2DTDP;

#define fillD3GDSDP2 \
d3gdsdp2[0] = D3PB_GDS0DP2; d3gdsdp2[1] = D3FE_GDS1DP2; \
d3gdsdp2[2] = D3MG_GDS2DP2;

-(void)pureOrder:(int)mask
			   t:(double)t
			   p:(double)p
			   s:(double [NS])s		// s[NS]       BINARY MASK: 000001
			  dt:(double [NS])dt	// ds[NS]/dt   BINARY MASK: 000010
			  dp:(double [NS])dp	// ds[NS]/dp   BINARY MASK: 000100
			 dt2:(double [NS])dt2	// d2s[NS]/dt2 BINARY MASK: 001000
			 dtp:(double [NS])dtp	// d2s[NS]/dtp BINARY MASK: 010000
			 dp2:(double [NS])dp2	// d2s[NS]/dp2 BINARY MASK: 100000
{
	int i;

	if ( (t != tOldPure) || (p != pOldPure) ) {
		double dgds[NS], sNew[NS];
		int iter=0;
		for (i=0; i<NS; i++) { sOldPure[i] = 2.0; sNew[i] = 0.9; }
		while (((fabs(sNew[0]-sOldPure[0]) > 10.0*DBL_EPSILON) ||
				(fabs(sNew[1]-sOldPure[1]) > 10.0*DBL_EPSILON) ||
				(fabs(sNew[2]-sOldPure[2]) > 10.0*DBL_EPSILON) ) && (iter < MAX_ITER) ) {
			double s[NS];

			for (i=0; i<NS; i++) s[i] = sNew[i];

			dgds[0] = DPB_GDS0;
			dgds[1] = DFE_GDS1;
			dgds[2] = DMG_GDS2;

			d2gds2Pure[0] = D2PB_GDS0S0;
			d2gds2Pure[1] = D2FE_GDS1S1;
			d2gds2Pure[2] = D2MG_GDS2S2;

			for (i=0; i<NS; i++) sOldPure[i] = s[i];

			for (i=0; i<NS; i++) {
				s[i] += - dgds[i]/d2gds2Pure[i];
				s[i] = MIN(s[i],  1.0 - DBL_EPSILON);
				s[i] = MAX(s[i], -1.0 + DBL_EPSILON);
			}

			for (i=0; i<NS; i++) sNew[i] = s[i];
			iter++;
		}
		tOldPure = t;
		pOldPure = p;
	}

	if (mask & FIRST  ) {   /* return s        */
		for (i=0; i<NS; i++) s[i] = sOldPure[i];
	}

	if (mask & SECOND ) {   /* compute ds/dt:  */
		double *s = sOldPure;
		double d2gdsdt[NS];

		fillD2GDSDT

		for (i=0; i<NS; i++) dt[i] = - d2gdsdt[i]/d2gds2Pure[i];
	}

	if (mask & THIRD  ) {   /* compute ds/dp:  */
		double d2gdsdp[NS];

		fillD2GDSDP

		for (i=0; i<NS; i++) dp[i] = - d2gdsdp[i]/d2gds2Pure[i];
	}

	if (mask & FOURTH ) {   /* compute d2s/dt2 */
		double *s = sOldPure;
		double d2gdsdt[NS], d3gds3[NS], d3gds2dt[NS], d3gdsdt2[NS], dsdt[NS];

		fillD2GDSDT
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDSDT2

		for (i=0; i<NS; i++) dsdt[i] = - d2gdsdt[i]/d2gds2Pure[i];
		for (i=0; i<NS; i++) dt2[i] = - (d3gdsdt2[i] + 2.0*d3gds2dt[i]*dsdt[i]
										 + d3gds3[i]*dsdt[i]*dsdt[i])/d2gds2Pure[i];
	}

	if (mask & FIFTH  ) {   /* compute d2s/dtp */
		double *s = sOldPure;
		double d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS], d3gds2dt[NS], d3gds2dp[NS],
		d3gdsdtdp[NS], dsdt[NS], dsdp[NS];

		fillD2GDSDT
		fillD2GDSDP
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDS2DP
		fillD3GDSDTDP

		for (i=0; i<NS; i++) dsdt[i] = - d2gdsdt[i]/d2gds2Pure[i];
		for (i=0; i<NS; i++) dsdp[i] = - d2gdsdp[i]/d2gds2Pure[i];

		for (i=0; i<NS; i++)
			dtp[i] = - (d3gdsdtdp[i] + d3gds2dt[i]*dsdp[i] + d3gds2dp[i]*dsdt[i]
						+ d3gds3[i]*dsdt[i]*dsdp[i])/d2gds2Pure[i];

	}

	if (mask & SIXTH  ) {   /* compute d2s/dp2 */
		double *s = sOldPure;
		double d2gdsdp[NS], d3gds3[NS], d3gds2dp[NS], d3gdsdp2[NS], dsdp[NS];

		fillD2GDSDP
		fillD3GDS3
		fillD3GDS2DP
		fillD3GDSDP2

		for (i=0; i<NS; i++) dsdp[i] = - d2gdsdp[i]/d2gds2Pure[i];

		for (i=0; i<NS; i++) dp2[i] = - (d3gdsdp2[i] + 2.0*d3gds2dp[i]*dsdp[i]
										 + d3gds3[i]*dsdp[i]*dsdp[i])/d2gds2Pure[i];
	}
}

#undef fillD2GDSDT
#undef fillD2GDSDP
#undef fillD3GDS3
#undef fillD3GDS2DT
#undef fillD3GDS2DP
#undef fillD3GDSDT2
#undef fillD3GDSDTDP
#undef fillD3GDSDP2

#define fillD2GDS2 \
d2gds2[0][0] = D2PB_GDS0S0; d2gds2[0][1] = 0.0;          d2gds2[0][2] = 0.0;        \
d2gds2[1][0] = 0.0;         d2gds2[1][1] = D2FE_GDS1S1;  d2gds2[1][2] = 0.0;        \
d2gds2[2][0] = 0.0;         d2gds2[2][1] = 0.0;          d2gds2[2][2] = D2MG_GDS2S2;

#define fillD2GDSDT \
d2gdsdt[0][0] = D2PB_GDS0DT; d2gdsdt[0][1] = 0.0;         d2gdsdt[0][2] =  0.0; \
d2gdsdt[1][0] = 0.0;         d2gdsdt[1][1] = D2FE_GDS1DT; d2gdsdt[1][2] =  0.0; \
d2gdsdt[2][0] = 0.0;         d2gdsdt[2][1] = 0.0;         d2gdsdt[2][2] =  D2MG_GDS2DT;

#define fillD2GDSDP \
d2gdsdp[0][0] = D2PB_GDS0DP; d2gdsdp[0][1] = 0.0;         d2gdsdp[0][2] =  0.0; \
d2gdsdp[1][0] = 0.0;         d2gdsdp[1][1] = D2FE_GDS1DP; d2gdsdp[1][2] =  0.0; \
d2gdsdp[2][0] = 0.0;         d2gdsdp[2][1] = 0.0;         d2gdsdp[2][2] =  D2MG_GDS2DP;

#define fillD2GDT2 \
d2gdt2[0] = D2PB_GDT2;       d2gdt2[1] = D2FE_GDT2;       d2gdt2[2] = D2MG_GDT2;

#define fillD2GDTDP \
d2gdtdp[0] = D2PB_GDTP;      d2gdtdp[1] = D2FE_GDTP;      d2gdtdp[2] = D2MG_GDTP;

#define fillD2GDP2 \
d2gdp2[0] = D2PB_GDP2;       d2gdp2[1] = D2FE_GDP2;       d2gdp2[2] = D2MG_GDP2;

#define fillD3GDS3 \
d3gds3[0][0] = D3PB_GDS0S0S0; d3gds3[0][1] = 0.0;            d3gds3[0][2] = 0.0; \
d3gds3[1][0] = 0.0;           d3gds3[1][1] = D3FE_GDS1S1S1;  d3gds3[1][2] = 0.0; \
d3gds3[2][0] = 0.0;           d3gds3[2][1] = 0.0;            d3gds3[2][2] = D3MG_GDS2S2S2;

#define fillD3GDS2DT \
d3gds2dt[0][0] = D3PB_GDS0S0DT; d3gds2dt[0][1] = 0.0;            d3gds2dt[0][2] = 0.0; \
d3gds2dt[1][0] = 0.0;           d3gds2dt[1][1] = D3FE_GDS1S1DT;  d3gds2dt[1][2] = 0.0; \
d3gds2dt[2][0] = 0.0;           d3gds2dt[2][1] = 0.0;            d3gds2dt[2][2] = D3MG_GDS2S2DT;

#define fillD3GDS2DP \
d3gds2dp[0][0] = D3PB_GDS0S0DP; d3gds2dp[0][1] = 0.0;            d3gds2dp[0][2] = 0.0; \
d3gds2dp[1][0] = 0.0;           d3gds2dp[1][1] = D3FE_GDS1S1DP;  d3gds2dp[1][2] = 0.0; \
d3gds2dp[2][0] = 0.0;           d3gds2dp[2][1] = 0.0;            d3gds2dp[2][2] = D3MG_GDS2S2DP;


#define fillD3GDSDT2 \
d3gdsdt2[0][0] = D3PB_GDS0DT2; d3gdsdt2[0][1] = 0.0;           d3gdsdt2[0][2] = 0.0; \
d3gdsdt2[1][0] = 0.0;          d3gdsdt2[1][1] = D3FE_GDS1DT2;  d3gdsdt2[1][2] = 0.0; \
d3gdsdt2[2][0] = 0.0;          d3gdsdt2[2][1] = 0.0;           d3gdsdt2[2][2] = D3MG_GDS2DT2;

#define fillD3GDSDTDP \
d3gdsdtdp[0][0] = D3PB_GDS0DTDP; d3gdsdtdp[0][1] = 0.0;           d3gdsdtdp[0][2] = 0.0; \
d3gdsdtdp[1][0] = 0.0;           d3gdsdtdp[1][1] = D3FE_GDS1DTDP; d3gdsdtdp[1][2] = 0.0; \
d3gdsdtdp[2][0] = 0.0;           d3gdsdtdp[2][1] = 0.0;           d3gdsdtdp[2][2] = D3MG_GDS2DTDP;

#define fillD3GDSDP2 \
d3gdsdp2[0][0] = D3PB_GDS0DP2; d3gdsdp2[0][1] = 0.0;           d3gdsdp2[0][2] = 0.0; \
d3gdsdp2[1][0] = 0.0;          d3gdsdp2[1][1] = D3FE_GDS1DP2;  d3gdsdp2[1][2] = 0.0; \
d3gdsdp2[2][0] = 0.0;          d3gdsdp2[2][1] = 0.0;           d3gdsdp2[2][2] = D3MG_GDS2DP2;

#define fillD3GDT3 \
d3gdt3[0] = D3PB_GDT3;         d3gdt3[1] = D3FE_GDT3;          d3gdt3[2] = D3MG_GDT3;

#define fillD3GDT2DP \
d3gdt2dp[0] = D3PB_GDT2DP;     d3gdt2dp[1] = D3FE_GDT2DP;      d3gdt2dp[2] = D3MG_GDT2DP;

#define fillD3GDTDP2 \
d3gdtdp2[0] = D3PB_GDTDP2;     d3gdtdp2[1] = D3FE_GDTDP2;      d3gdtdp2[2] = D3MG_GDTDP2;

#define fillD3GDP3 \
d3gdp3[0] = D3PB_GDP3;         d3gdp3[1] = D3FE_GDP3;          d3gdp3[2] = D3MG_GDP3;

-(void)pureEndMembers:(int)mask
					t:(double)t
					p:(double)p
					a:(double [NA])a        // activities              BINARY MASK: 0000000000001
				   mu:(double [NA])mu       // chemical potentials     BINARY MASK: 0000000000010
				 gmix:(double [NA])gmix     // Gibbs energy            BINARY MASK: 0000000000100
				 hmix:(double [NA])hmix     // Enthalpy of mixing      BINARY MASK: 0000000001000
				 smix:(double [NA])smix     // Entropy of mixing       BINARY MASK: 0000000010000
				cpmix:(double [NA])cpmix    // Heat capacity of mixing BINARY MASK: 0000000100000
			  cpmixdt:(double [NA])cpmixdt  // d(cp)/d(t)              BINARY MASK: 0000001000000
				 vmix:(double [NA])vmix     // Volume of mixing        BINARY MASK: 0000010000000
			   vmixdt:(double [NA])vmixdt   // d(v)/d(t)               BINARY MASK: 0000100000000
			   vmixdp:(double [NA])vmixdp   // d(v)/d(p)               BINARY MASK: 0001000000000
			  vmixdt2:(double [NA])vmixdt2  // d2(v)/d(t)2             BINARY MASK: 0010000000000
			 vmixdtdp:(double [NA])vmixdtdp // d2(v)/d(t)d(p)          BINARY MASK: 0100000000000
			  vmixdp2:(double [NA])vmixdp2  // d2(v)/d(p)2             BINARY MASK: 1000000000000
{
	double s[NS];
	int i, j;

	[self pureOrder:FIRST t:t p:p s:s dt:NULL dp:NULL dt2:NULL dtp:NULL dp2:NULL];

	if (mask & FIRST) {
		a[0] = PB_G;
		a[0] = exp(a[0]/(R*t));
		a[1] = FE_G;
		a[1] = exp(a[1]/(R*t));
		a[2] = MG_G;
		a[2] = exp(a[2]/(R*t));
	}

	if (mask & SECOND) {
		mu[0] = PB_G;
		mu[1] = FE_G;
		mu[2] = MG_G;
	}

	if (mask & THIRD) {
		gmix[0] = PB_G;
		gmix[1] = FE_G;
		gmix[2] = MG_G;
	}

	if (mask & FOURTH) {
		hmix[0] = (PB_G) + t*(PB_S);
		hmix[1] = (FE_G) + t*(FE_S);
		hmix[2] = (MG_G) + t*(MG_S);
	}

	if (mask & FIFTH) {
		smix[0] = PB_S;
		smix[1] = FE_S;
		smix[2] = MG_S;
	}

	if (mask & SIXTH) {
		double d2gdsdt[NA][NS], d2gds2[NA][NS], dsdt[NS];

		fillD2GDS2
		fillD2GDSDT

		[self pureOrder:SECOND t:t p:p s:NULL dt:dsdt dp:NULL dt2:NULL dtp:NULL dp2:NULL];

		cpmix[0] = D2PB_GDT2;
		cpmix[1] = D2FE_GDT2;
		cpmix[2] = D2MG_GDT2;

		for (i=0; i<NA; i++) {
			for (j=0; j<NS; j++)
				cpmix[i] += 2.0*d2gdsdt[i][j]*dsdt[j] + d2gds2[i][j]*SQUARE(dsdt[j]);
			cpmix[i] *= -t;
		}
	}

	if(mask & SEVENTH) {
		double d2gdsdt[NA][NS], d2gds2[NA][NS], d2gdt2[NA], d3gds3[NA][NS],
		d3gds2dt[NA][NS], d3gdsdt2[NA][NS], d3gdt3[NA], dsdt[NS], d2sdt2[NS],
		temp;

		fillD2GDT2
		fillD2GDSDT
		fillD2GDS2
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDSDT2
		fillD3GDT3

		[self pureOrder:SECOND | FOURTH t:t p:p s:NULL dt:dsdt dp:NULL dt2:d2sdt2 dtp:NULL dp2:NULL];

		for (i=0; i<NA; i++) {
			temp = d2gdt2[i];
			for (j=0; j<NS; j++)
				temp += 2.0*d2gdsdt[i][j]*dsdt[j] + d2gds2[i][j]*SQUARE(dsdt[j]);

			cpmixdt[i] = d3gdt3[i];
			for (j=0; j<NS; j++)
				cpmixdt[i] += 3.0*d3gdsdt2[i][j]*dsdt[j]
				+ 3.0*d2gdsdt[i][j]*d2sdt2[j]
				+ 3.0*d2gds2[i][j]*dsdt[j]*d2sdt2[j]
				+ 3.0*d3gds2dt[i][j]*dsdt[j]*dsdt[j]
				+ d3gds3[i][j]*dsdt[j]*dsdt[j]*dsdt[j];

			cpmixdt[i] = -t*cpmixdt[i] - temp;
		}
	}

	if (mask & EIGHTH) {
		vmix[0] = DPB_GDP;
		vmix[1] = DFE_GDP;
		vmix[2] = DMG_GDP;
	}

	if(mask & NINTH) {
		double d2gdsdt[NA][NS], d2gdsdp[NA][NS], d2gds2[NA][NS], d2gdtdp[NA],
		dsdt[NS], dsdp[NS];

		fillD2GDSDT
		fillD2GDS2
		fillD2GDSDP
		fillD2GDTDP

		[self pureOrder:SECOND | THIRD t:t p:p s:NULL dt:dsdt dp:dsdp dt2:NULL dtp:NULL dp2:NULL];

		for (i=0; i<NA; i++) {
			vmixdt[i] = d2gdtdp[i];
			for (j=0; j<NS; j++)
				vmixdt[i] += d2gdsdt[i][j]*dsdp[j] + d2gdsdp[i][j]*dsdt[j]
                + d2gds2[i][j]*dsdt[j]*dsdp[j];
		}
	}

	if(mask & TENTH) {
		double d2gdsdp[NA][NS], d2gds2[NA][NS], d2gdp2[NA], dsdp[NS];

		fillD2GDS2
		fillD2GDSDP
		fillD2GDP2

		[self pureOrder:THIRD t:t p:p s:NULL dt:NULL dp:dsdp dt2:NULL dtp:NULL dp2:NULL];

		for (i=0; i<NA; i++) {
			vmixdp[i] = d2gdp2[i];
			for (j=0; j<NS; j++)
				vmixdp[i] += 2.0*d2gdsdp[i][j]*dsdp[j] + d2gds2[i][j]*dsdp[j]*dsdp[j];
		}
	}

	if(mask & ELEVENTH) {
		double d2gdsdt[NA][NS], d2gdsdp[NA][NS], d2gds2[NA][NS], d3gds3[NA][NS],
		d3gds2dt[NA][NS], d3gdsdt2[NA][NS], d3gds2dp[NA][NS], d3gdsdtdp[NA][NS],
		d3gdt2dp[NA], dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS];

		fillD2GDSDT
		fillD2GDS2
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDSDT2
		fillD2GDSDP
		fillD3GDS2DP
		fillD3GDSDTDP
		fillD3GDT2DP

		[self pureOrder:SECOND | THIRD | FOURTH | FIFTH t:t p:p s:NULL dt:dsdt dp:dsdp dt2:d2sdt2 dtp:d2sdtdp dp2:NULL];

		for (i=0; i<NA; i++) {
			vmixdt2[i] = d3gdt2dp[i];
			for (j=0; j<NS; j++)
				vmixdt2[i] += d3gdsdt2[i][j]*dsdp[j]
				+ 2.0*d2gdsdt[i][j]*d2sdtdp[j]
				+ d2gdsdp[i][j]*d2sdt2[j] + 2.0*d3gdsdtdp[i][j]*dsdt[j]
				+ 2.0*d3gds2dt[i][j]*dsdt[j]*dsdp[j]
				+ d2gds2[i][j]*d2sdt2[j]*dsdp[j]
				+ 2.0*d2gds2[i][j]*dsdt[j]*d2sdtdp[j]
				+ d3gds2dp[i][j]*dsdt[j]*dsdt[j]
				+ d3gds3[i][j]*dsdt[j]*dsdt[j]*dsdp[j];
		}
	}

	if(mask & TWELFTH) {
		double d2gdsdt[NA][NS], d2gdsdp[NA][NS], d2gds2[NA][NS], d3gds3[NA][NS],
		d3gds2dt[NA][NS], d3gds2dp[NA][NS], d3gdsdtdp[NA][NS], d3gdsdp2[NA][NS],
		d3gdtdp2[NA], dsdt[NS], dsdp[NS], d2sdtdp[NS], d2sdp2[NS];

		fillD2GDSDT
		fillD2GDS2
		fillD3GDS3
		fillD3GDS2DT
		fillD2GDSDP
		fillD3GDS2DP
		fillD3GDSDTDP
		fillD3GDSDP2
		fillD3GDTDP2

		[self pureOrder:SECOND | THIRD | FIFTH | SIXTH t:t p:p s:NULL dt:dsdt dp:dsdp dt2:NULL dtp:d2sdtdp dp2:d2sdp2];

		for (i=0; i<NA; i++) {
			vmixdtdp[i] = d3gdtdp2[i];
			for (j=0; j<NS; j++)
				vmixdtdp[i] += 2.0*d3gdsdtdp[i][j]*dsdp[j] + d2gdsdt[i][j]*d2sdp2[j]
				+ 2.0*d2gdsdp[i][j]*d2sdtdp[j] + d3gdsdp2[i][j]*dsdt[j]
				+ 2.0*d3gds2dp[i][j]*dsdt[j]*dsdp[j]
				+ d2gds2[i][j]*dsdt[j]*d2sdp2[j]
				+ 2.0*d2gds2[i][j]*d2sdtdp[j]*dsdp[j]
				+ d3gds2dt[i][j]*dsdp[j]*dsdp[j]
				+ d3gds3[i][j]*dsdt[j]*dsdp[j]*dsdp[j];
		}
	}

	if(mask & THIRTEENTH) {
		double d2gdsdp[NA][NS], d2gds2[NA][NS], d3gds3[NA][NS], d3gds2dp[NA][NS],
		d3gdsdp2[NA][NS], d3gdp3[NA], dsdp[NS], d2sdp2[NS];

		fillD2GDS2
		fillD3GDS3
		fillD2GDSDP
		fillD3GDS2DP
		fillD3GDSDP2
		fillD3GDP3

		[self pureOrder:THIRD | SIXTH t:t p:p s:NULL dt:NULL dp:dsdp dt2:NULL dtp:NULL dp2:d2sdp2];

		for (i=0; i<NA; i++) {
			vmixdp2[i] = d3gdp3[i];
			for (j=0; j<NS; j++)
				vmixdp2[i] += 3.0*d3gdsdp2[i][j]*dsdp[j] + 3.0*d2gdsdp[i][j]*d2sdp2[j]
				+ 3.0*d2gds2[i][j]*dsdp[j]*d2sdp2[j]
				+ 3.0*d3gds2dp[i][j]*dsdp[j]*dsdp[j]
				+ d3gds3[i][j]*dsdp[j]*dsdp[j]*dsdp[j];
		}
	}

}

#undef fillD2GDS2
#undef fillD2GDSDT
#undef fillD2GDSDP
#undef fillD2GDT2
#undef fillD2GDTDP
#undef fillD2GDP2
#undef fillD3GDS3
#undef fillD3GDS2DT
#undef fillD3GDS2DP
#undef fillD3GDSDT2
#undef fillD3GDSDTDP
#undef fillD3GDSDP2
#undef fillD3GDT3
#undef fillD3GDT2DP
#undef fillD3GDTDP2
#undef fillD3GDP3

#undef PB_S
#undef PB_H
#undef PB_G
#undef DPB_GDS0
#undef DPB_GDT
#undef DPB_GDP
#undef D2PB_GDS0S0
#undef D2PB_GDS0DT
#undef D2PB_GDS0DP
#undef D2PB_GDT2
#undef D2PB_GDTP
#undef D2PB_GDP2
#undef D3PB_GDS0S0S0
#undef D3PB_GDS0S0DT
#undef D3PB_GDS0S0DP
#undef D3PB_GDS0DT2
#undef D3PB_GDS0DTDP
#undef D3PB_GDS0DP2
#undef D3PB_GDT3
#undef D3PB_GDT2DP
#undef D3PB_GDTDP2
#undef D3PB_GDP3

#undef FE_S
#undef FE_H
#undef FE_G
#undef DFE_GDS1
#undef DFE_GDT
#undef DFE_GDP
#undef D2FE_GDS1S1
#undef D2FE_GDS1DT
#undef D2FE_GDS1DP
#undef D2FE_GDT2
#undef D2FE_GDTP
#undef D2FE_GDP2
#undef D3FE_GDS1S1S1
#undef D3FE_GDS1S1DT
#undef D3FE_GDS1S1DP
#undef D3FE_GDS1DT2
#undef D3FE_GDS1DTDP
#undef D3FE_GDS1DP2
#undef D3FE_GDT3
#undef D3FE_GDT2DP
#undef D3FE_GDTDP2
#undef D3FE_GDP3

#undef MG_S
#undef MG_H
#undef MG_G
#undef DMG_GDS2
#undef DMG_GDT
#undef DMG_GDP
#undef D2MG_GDS2S2
#undef D2MG_GDS2DT
#undef D2MG_GDS2DP
#undef D2MG_GDT2
#undef D2MG_GDTP
#undef D2MG_GDP2
#undef D3MG_GDS2S2S2
#undef D3MG_GDS2S2DT
#undef D3MG_GDS2S2DP
#undef D3MG_GDS2DT2
#undef D3MG_GDS2DTDP
#undef D3MG_GDS2DP2
#undef D3MG_GDT3
#undef D3MG_GDT2DP
#undef D3MG_GDTDP2
#undef D3MG_GDP3


/*
 * Global (to this file): activity definitions and component transforms
 *    The function conRhm defines the conversion from m[i], to r[j]
 */
/* Order: X2, X3 */
#define FR0(i)      (i == 1) ? 1.0 - r[0] : - r[0]
#define FR1(i)      (i == 2) ? 1.0 - r[1] : - r[1]

/* Order: s1, s2, s3 */
#define GSS0(i)     (i == 0) ? 1.0 - s[0] : - s[0]
#define GSS1(i)     (i == 1) ? 1.0 - s[1] : - s[1]
#define GSS2(i)     (i == 2) ? 1.0 - s[2] : - s[2]

#define DFR0DR0(i)  - 1.0
#define DFR1DR1(i)  - 1.0

#define DGSS0DS0(i) - 1.0
#define DGSS1DS1(i) - 1.0
#define DGSS2DS2(i) - 1.0

#define ENDMEMBERS ((1.0-r[0]-r[1])*ends[0] + r[0]*ends[1] + r[1]*ends[2])

#define DENDDR0    (ends[1] - ends[0])
#define DENDDR1    (ends[2] - ends[0])

/*
 * Global (to this file): derivative definitions
 */

#define S  - R*(      xfe2M1*log(xfe2M1) +     xmg2M1*log(xmg2M1) \
+     xti4M1*log(xti4M1) +     xfe3M1*log(xfe3M1) \
+ 2.0*xfe2M2*log(xfe2M2) + 2.0*xmg2M2*log(xmg2M2) \
+ 2.0*xti4M2*log(xti4M2) + 2.0*xfe3M2*log(xfe3M2) )
#define H  (G0)+(GX2)*r[0]+(GX3)*r[1]+(GS1)*s[0]+(GS2)*s[1]+(GS3)*s[2] \
+(GX2X2)*r[0]*r[0]+(GX2X3)*r[0]*r[1]+(GX2S1)*r[0]*s[0] \
+(GX2S2)*r[0]*s[1]+(GX2S3)*r[0]*s[2]+(GX3X3)*r[1]*r[1] \
+(GX3S1)*r[1]*s[0]+(GX3S2)*r[1]*s[1]+(GX3S3)*r[1]*s[2] \
+(GS1S1)*s[0]*s[0]+(GS1S2)*s[0]*s[1]+(GS1S3)*s[0]*s[2] \
+(GS2S2)*s[1]*s[1]+(GS2S3)*s[1]*s[2]+(GS3S3)*s[2]*s[2]
#define G  H - t*(S)

/*----------------------------------------------------------------------------*/

#define DGDR0  R*t*(  0.5*log(xfe2M1) - 0.5*log(xfe3M1) \
+ 0.5*log(xfe2M2) + log(xti4M2) - 1.5*log(xfe3M2) ) \
+(GX2)+2.0*(GX2X2)*r[0]+(GX2X3)*r[1]+(GX2S1)*s[0]+(GX2S2)*s[1]+(GX2S3)*s[2]
#define DGDR1  R*t*(  0.5*log(xmg2M1) - 0.5*log(xfe3M1) \
+ 0.5*log(xmg2M2) + log(xti4M2) - 1.5*log(xfe3M2) ) \
+(GX3)+(GX2X3)*r[0]+2.0*(GX3X3)*r[1]+(GX3S1)*s[0]+(GX3S2)*s[1]+(GX3S3)*s[2]
#define DGDS0  R*t*(- 0.5*log(xti4M1) + 0.5*log(xfe3M1) \
+ 0.5*log(xti4M2) - 0.5*log(xfe3M2) ) \
+(GS1)+(GX2S1)*r[0]+(GX3S1)*r[1]+2.0*(GS1S1)*s[0]+(GS1S2)*s[1]+(GS1S3)*s[2]
#define DGDS1  R*t*(  0.5*log(xfe2M1) - 0.5*log(xti4M1) \
- 0.5*log(xfe2M2) + 0.5*log(xti4M2) ) \
+(GS2)+(GX2S2)*r[0]+(GX3S2)*r[1]+(GS1S2)*s[0]+2.0*(GS2S2)*s[1]+(GS2S3)*s[2]
#define DGDS2  R*t*(  0.5*log(xmg2M1) - 0.5*log(xti4M1) \
- 0.5*log(xmg2M2) + 0.5*log(xti4M2) ) \
+(GS3)+(GX2S3)*r[0]+(GX3S3)*r[1]+(GS1S3)*s[0]+(GS2S3)*s[1]+2.0*(GS3S3)*s[2]
#define DGDT   - (S)
#define DGDP   0.0

/*----------------------------------------------------------------------------*/

#define D2GDR0R0  R*t*(  0.25/xfe2M1  + 0.25/xfe3M1 \
+ 0.125/xfe2M2 + 0.5/xti4M2 + 1.125/xfe3M2 ) + 2.0*(GX2X2)
#define D2GDR0R1  R*t*(  0.25/xfe3M1 + 0.5/xti4M2 + 1.125/xfe3M2 ) + (GX2X3)
#define D2GDR0S0  R*t*(  - 0.25/xfe3M1 + 0.25/xti4M2 + 0.375/xfe3M2 ) + (GX2S1)
#define D2GDR0S1  R*t*(  0.25/xfe2M1 - 0.125/xfe2M2 + 0.25/xti4M2 ) +(GX2S2)
#define D2GDR0S2  R*t*(  0.25/xti4M2 ) + (GX2S3)
#define D2GDR0DT  R*  (  0.5*log(xfe2M1) - 0.5*log(xfe3M1) \
+ 0.5*log(xfe2M2) + log(xti4M2) - 1.5*log(xfe3M2) )
#define D2GDR0DP  0.0

#define D2GDR1R1  R*t*(  0.25/xmg2M1 + 0.25/xfe3M1 \
+ 0.125/xmg2M2 + 0.5/xti4M2 + 1.125/xfe3M2 ) + 2.0*(GX3X3)
#define D2GDR1S0  R*t*( - 0.25/xfe3M1 + 0.25/xti4M2 + 0.375/xfe3M2 ) + (GX3S1)
#define D2GDR1S1  R*t*( 0.25/xti4M2 ) + (GX3S2)
#define D2GDR1S2  R*t*(  0.25/xmg2M1 - 0.125/xmg2M2 + 0.25/xti4M2 ) + (GX3S3)
#define D2GDR1DT  R*  (  0.5*log(xmg2M1) - 0.5*log(xfe3M1) \
+ 0.5*log(xmg2M2) + log(xti4M2) - 1.5*log(xfe3M2) )
#define D2GDR1DP  0.0

#define D2GDS0S0  R*t*(  0.25/xti4M1  + 0.25/xfe3M1 \
+ 0.125/xti4M2 + 0.125/xfe3M2 ) +2.0*(GS1S1)
#define D2GDS0S1  R*t*(  0.25/xti4M1  + 0.125/xti4M2 ) + (GS1S2)
#define D2GDS0S2  R*t*(  0.25/xti4M1  + 0.125/xti4M2 ) + (GS1S3)
#define D2GDS0DT  R*  (- 0.5*log(xti4M1) + 0.5*log(xfe3M1) \
+ 0.5*log(xti4M2) - 0.5*log(xfe3M2) )
#define D2GDS0DP  0.0

#define D2GDS1S1  R*t*(  0.25/xfe2M1  + 0.25/xti4M1 \
+ 0.125/xfe2M2 + 0.125/xti4M2 ) + 2.0*(GS2S2)
#define D2GDS1S2  R*t*(  0.25/xti4M1  + 0.125/xti4M2 ) +(GS2S3)
#define D2GDS1DT  R*  (  0.5*log(xfe2M1) - 0.5*log(xti4M1) \
- 0.5*log(xfe2M2) + 0.5*log(xti4M2) )
#define D2GDS1DP  0.0

#define D2GDS2S2  R*t*(  0.25/xmg2M1  + 0.25/xti4M1 \
+ 0.125/xmg2M2 + 0.125/xti4M2 ) + 2.0*(GS3S3)
#define D2GDS2DT  R*  (  0.5*log(xmg2M1) - 0.5*log(xti4M1) \
- 0.5*log(xmg2M2) + 0.5*log(xti4M2) )
#define D2GDS2DP  0.0

#define D2GDT2    0.0
#define D2GDTDP   0.0
#define D2GDP2    0.0

/*----------------------------------------------------------------------------*/

#define D3GDR0R0R0  R*t*(- 0.125/SQUARE(xfe2M1)   + 0.125/SQUARE(xfe3M1) \
- 0.03125/SQUARE(xfe2M2) - 0.25/SQUARE(xti4M2)  + 0.84375/SQUARE(xfe3M2) )
#define D3GDR0R0R1  R*t*(  0.125/SQUARE(xfe3M1) - 0.25/SQUARE(xti4M2) + 0.84375/SQUARE(xfe3M2) )
#define D3GDR0R0S0  R*t*(- 0.125/SQUARE(xfe3M1) - 0.125/SQUARE(xti4M2) + 0.28125/SQUARE(xfe3M2) )
#define D3GDR0R0S1  R*t*(- 0.125/SQUARE(xfe2M1) + 0.03125/SQUARE(xfe2M2) - 0.125/SQUARE(xti4M2) )
#define D3GDR0R0S2  R*t*(- 0.125/SQUARE(xti4M2) )
#define D3GDR0R0DT  R*  (  0.25/xfe2M1  + 0.25/xfe3M1 + 0.125/xfe2M2 + 0.5/xti4M2 + 1.125/xfe3M2 )
#define D3GDR0R0DP  0.0

#define D3GDR0R1R1  R*t*(  0.125/SQUARE(xfe3M1) - 0.25/SQUARE(xti4M2) + 0.84375/SQUARE(xfe3M2) )
#define D3GDR0R1S0  R*t*(- 0.125/SQUARE(xfe3M1) - 0.125/SQUARE(xti4M2) + 0.28125/SQUARE(xfe3M2) )
#define D3GDR0R1S1  R*t*(- 0.125/SQUARE(xti4M2) )
#define D3GDR0R1S2  R*t*(- 0.125/SQUARE(xti4M2) )
#define D3GDR0R1DT  R*  (  0.25/xfe3M1 + 0.5/xti4M2 + 1.125/(xfe3M2) )
#define D3GDR0R1DP  0.0

#define D3GDR0S0S0  R*t*(  0.125/SQUARE(xfe3M1) - 0.0625/SQUARE(xti4M2) + 0.09375/SQUARE(xfe3M2) )
#define D3GDR0S0S1  R*t*(- 0.0625/SQUARE(xti4M2) )
#define D3GDR0S0S2  R*t*(- 0.0625/SQUARE(xti4M2) )
#define D3GDR0S0DT  R*  (  - 0.25/xfe3M1 + 0.25/xti4M2 + 0.375/xfe3M2 )
#define D3GDR0S0DP  0.0

#define D3GDR0S1S1  R*t*(- 0.125/SQUARE(xfe2M1) - 0.03125/SQUARE(xfe2M2) - 0.0625/SQUARE(xti4M2) )
#define D3GDR0S1S2  R*t*(- 0.0625/SQUARE(xti4M2) )
#define D3GDR0S1DT  R*  (  0.25/xfe2M1 - 0.125/xfe2M2 + 0.25/xti4M2 )
#define D3GDR0S1DP  0.0

#define D3GDR0S2S2  R*t*(- 0.0625/SQUARE(xti4M2) )
#define D3GDR0S2DT  R*  (  0.25/xti4M2 )
#define D3GDR0S2DP  0.0

#define D3GDR0DT2   0.0
#define D3GDR0DTDP  0.0
#define D3GDR0DP2   0.0

#define D3GDR1R1R1  R*t*(- 0.125/SQUARE(xmg2M1) + 0.125/SQUARE(xfe3M1) \
- 0.03125/SQUARE(xmg2M2) - 0.25/SQUARE(xti4M2) + 0.84375/SQUARE(xfe3M2) )
#define D3GDR1R1S0  R*t*(- 0.125/SQUARE(xfe3M1) - 0.125/SQUARE(xti4M2) + 0.28125/SQUARE(xfe3M2) )
#define D3GDR1R1S1  R*t*(- 0.125/SQUARE(xti4M2) )
#define D3GDR1R1S2  R*t*(- 0.125/SQUARE(xmg2M1) + 0.03125/SQUARE(xmg2M2) - 0.125/SQUARE(xti4M2) )
#define D3GDR1R1DT  R*  (  0.25/xmg2M1 + 0.25/xfe3M1 + 0.125/xmg2M2 + 0.5/xti4M2 + 1.125/xfe3M2 )
#define D3GDR1R1DP  0.0

#define D3GDR1S0S0  R*t*( 0.125/SQUARE(xfe3M1) - 0.0625/SQUARE(xti4M2) + 0.09375/SQUARE(xfe3M2) )
#define D3GDR1S0S1  R*t*(- 0.0625/SQUARE(xti4M2) )
#define D3GDR1S0S2  R*t*(- 0.0625/SQUARE(xti4M2) )
#define D3GDR1S0DT  R*  ( - 0.25/xfe3M1 + 0.25/xti4M2 + 0.375/xfe3M2 )
#define D3GDR1S0DP  0.0

#define D3GDR1S1S1  R*t*(- 0.0625/SQUARE(xti4M2) )
#define D3GDR1S1S2  R*t*(- 0.0625/SQUARE(xti4M2) )
#define D3GDR1S1DT  R*  ( 0.25/xti4M2 )
#define D3GDR1S1DP  0.0

#define D3GDR1S2S2  R*t*(- 0.125/SQUARE(xmg2M1) - 0.03125/SQUARE(xmg2M2) - 0.0625/SQUARE(xti4M2) )
#define D3GDR1S2DT  R*  (  0.25/xmg2M1 - 0.125/xmg2M2 + 0.25/xti4M2 )
#define D3GDR1S2DP  0.0

#define D3GDR1DT2   0.0
#define D3GDR1DTDP  0.0
#define D3GDR1DP2   0.0

#define D3GDS0S0S0  R*t*(  0.125/SQUARE(xti4M1)   - 0.125/SQUARE(xfe3M1) \
- 0.03125/SQUARE(xti4M2) + 0.03125/SQUARE(xfe3M2) )
#define D3GDS0S0S1  R*t*(  0.125/SQUARE(xti4M1) - 0.03125/SQUARE(xti4M2) )
#define D3GDS0S0S2  R*t*(  0.125/SQUARE(xti4M1) - 0.03125/SQUARE(xti4M2) )
#define D3GDS0S0DT  R*  (  0.25/xti4M1  + 0.25/xfe3M1 + 0.125/xti4M2 + 0.125/xfe3M2 )
#define D3GDS0S0DP  0.0

#define D3GDS0S1S1  R*t*(  0.125/SQUARE(xti4M1)  - 0.03125/SQUARE(xti4M2) )
#define D3GDS0S1S2  R*t*(  0.125/SQUARE(xti4M1)  - 0.03125/SQUARE(xti4M2) )
#define D3GDS0S1DT  R*  (  0.25/xti4M1  + 0.125/xti4M2 )
#define D3GDS0S1DP  0.0

#define D3GDS0S2S2  R*t*(  0.125/SQUARE(xti4M1)  - 0.03125/SQUARE(xti4M2) )
#define D3GDS0S2DT  R*  (  0.25/xti4M1  + 0.125/xti4M2 )
#define D3GDS0S2DP  0.0
#define D3GDS0DT2   0.0
#define D3GDS0DTDP  0.0
#define D3GDS0DP2   0.0

#define D3GDS1S1S1  R*t*(- 0.125/SQUARE(xfe2M1)   + 0.125/SQUARE(xti4M1) \
+ 0.03125/SQUARE(xfe2M2) - 0.03125/SQUARE(xti4M2) )
#define D3GDS1S1S2  R*t*(  0.125/SQUARE(xti4M1) - 0.03125/SQUARE(xti4M2) )
#define D3GDS1S1DT  R*  (  0.25/xfe2M1  + 0.25/xti4M1 + 0.125/xfe2M2 + 0.125/xti4M2 )
#define D3GDS1S1DP  0.0

#define D3GDS1S2S2  R*t*(  0.125/SQUARE(xti4M1)  - 0.03125/SQUARE(xti4M2) )
#define D3GDS1S2DT  R*  (  0.25/xti4M1  + 0.125/xti4M2 )
#define D3GDS1S2DP  0.0
#define D3GDS1DT2   0.0
#define D3GDS1DTDP  0.0
#define D3GDS1DP2   0.0

#define D3GDS2S2S2  R*t*(- 0.125/SQUARE(xmg2M1)   + 0.125/SQUARE(xti4M1) \
+ 0.03125/SQUARE(xmg2M2) - 0.03125/SQUARE(xti4M2) )
#define D3GDS2S2DT  R*  (  0.25/xmg2M1  + 0.25/xti4M1 + 0.125/xmg2M2 + 0.125/xti4M2 )
#define D3GDS2S2DP  0.0
#define D3GDS2DT2   0.0
#define D3GDS2DTDP  0.0
#define D3GDS2DP2   0.0

#define D3GDT3      0.0
#define D3GDT2DP    0.0
#define D3GDTDP2    0.0
#define D3GDP3      0.0

/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1; \
d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;

#define fillD2GDRDS \
d2gdrds[0][0] = D2GDR0S0; d2gdrds[0][1] = D2GDR0S1; d2gdrds[0][2] = D2GDR0S2; \
d2gdrds[1][0] = D2GDR1S0; d2gdrds[1][1] = D2GDR1S1; d2gdrds[1][2] = D2GDR1S2;

#define fillD2GDRDT \
d2gdrdt[0] = D2GDR0DT; d2gdrdt[1] = D2GDR1DT;

#define fillD2GDRDP \
d2gdrdp[0] = D2GDR0DP; d2gdrdp[1] = D2GDR1DP;

#define fillD2GDS2 \
d2gds2[0][0] = D2GDS0S0;     d2gds2[0][1] = D2GDS0S1;     d2gds2[0][2] = D2GDS0S2; \
d2gds2[1][0] = d2gds2[0][1]; d2gds2[1][1] = D2GDS1S1;     d2gds2[1][2] = D2GDS1S2; \
d2gds2[2][0] = d2gds2[0][2]; d2gds2[2][1] = d2gds2[1][2]; d2gds2[2][2] = D2GDS2S2;

#define fillD2GDSDT \
d2gdsdt[0] = D2GDS0DT;  d2gdsdt[1] = D2GDS1DT; d2gdsdt[2] = D2GDS2DT;

#define fillD2GDSDP \
d2gdsdp[0] = D2GDS0DP;  d2gdsdp[1] = D2GDS1DP; d2gdsdp[2] = D2GDS2DP;

#define fillD3GDR3 \
d3gdr3[0][0][0] = D3GDR0R0R0;      d3gdr3[0][0][1] = D3GDR0R0R1; \
d3gdr3[0][1][0] = d3gdr3[0][0][1]; d3gdr3[0][1][1] = D3GDR0R1R1; \
d3gdr3[1][0][0] = d3gdr3[0][0][1]; d3gdr3[1][0][1] = d3gdr3[0][1][1]; \
d3gdr3[1][1][0] = d3gdr3[0][1][1]; d3gdr3[1][1][1] = D3GDR1R1R1;

#define fillD3GDR2DS \
d3gdr2ds[0][0][0] = D3GDR0R0S0;        d3gdr2ds[0][0][1] = D3GDR0R0S1;        \
d3gdr2ds[0][0][2] = D3GDR0R0S2; \
d3gdr2ds[0][1][0] = D3GDR0R1S0;        d3gdr2ds[0][1][1] = D3GDR0R1S1;        \
d3gdr2ds[0][1][2] = D3GDR0R1S2; \
d3gdr2ds[1][0][0] = d3gdr2ds[0][1][0]; d3gdr2ds[1][0][1] = d3gdr2ds[0][1][1]; \
d3gdr2ds[1][0][2] = d3gdr2ds[0][1][2]; \
d3gdr2ds[1][1][0] = D3GDR1R1S0;        d3gdr2ds[1][1][1] = D3GDR1R1S1;        \
d3gdr2ds[1][1][2] = D3GDR1R1S2;

#define fillD3GDR2DT \
d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT; \
d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT;

#define fillD3GDR2DP \
d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP; \
d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP;

#define fillD3GDRDS2 \
d3gdrds2[0][0][0] = D3GDR0S0S0;        d3gdrds2[0][0][1] = D3GDR0S0S1; \
d3gdrds2[0][0][2] = D3GDR0S0S2; \
d3gdrds2[0][1][0] = d3gdrds2[0][0][1]; d3gdrds2[0][1][1] = D3GDR0S1S1; \
d3gdrds2[0][1][2] = D3GDR0S1S2; \
d3gdrds2[0][2][0] = d3gdrds2[0][0][2]; d3gdrds2[0][2][1] = d3gdrds2[0][1][2]; \
d3gdrds2[0][2][2] = D3GDR0S2S2; \
d3gdrds2[1][0][0] = D3GDR1S0S0;        d3gdrds2[1][0][1] = D3GDR1S0S1; \
d3gdrds2[1][0][2] = D3GDR1S0S2; \
d3gdrds2[1][1][0] = d3gdrds2[1][0][1]; d3gdrds2[1][1][1] = D3GDR1S1S1; \
d3gdrds2[1][1][2] = D3GDR1S1S2; \
d3gdrds2[1][2][0] = d3gdrds2[1][0][2]; d3gdrds2[1][2][1] = d3gdrds2[1][1][2]; \
d3gdrds2[1][2][2] = D3GDR1S2S2; \

#define fillD3GDRDSDT \
d3gdrdsdt[0][0] = D3GDR0S0DT;  d3gdrdsdt[0][1] = D3GDR0S1DT; \
d3gdrdsdt[0][2] = D3GDR0S2DT; \
d3gdrdsdt[1][0] = D3GDR1S0DT;  d3gdrdsdt[1][1] = D3GDR1S1DT; \
d3gdrdsdt[1][2] = D3GDR1S2DT; \

#define fillD3GDRDSDP \
d3gdrdsdp[0][0] = D3GDR0S0DP; d3gdrdsdp[0][1] = D3GDR0S1DP; \
d3gdrdsdp[0][2] = D3GDR0S2DP; \
d3gdrdsdp[1][0] = D3GDR1S0DP; d3gdrdsdp[1][1] = D3GDR1S1DP; \
d3gdrdsdp[1][2] = D3GDR1S2DP; \

#define fillD3GDS3 \
d3gds3[0][0][0] = D3GDS0S0S0;      d3gds3[0][0][1] = D3GDS0S0S1; \
d3gds3[0][0][2] = D3GDS0S0S2; \
d3gds3[0][1][0] = d3gds3[0][0][1]; d3gds3[0][1][1] = D3GDS0S1S1; \
d3gds3[0][1][2] = D3GDS0S1S2; \
d3gds3[0][2][0] = d3gds3[0][0][2]; d3gds3[0][2][1] = d3gds3[0][1][2]; \
d3gds3[0][2][2] = D3GDS0S2S2; \
d3gds3[1][0][0] = d3gds3[0][0][1]; d3gds3[1][0][1] = d3gds3[0][1][1]; \
d3gds3[1][0][2] = d3gds3[0][1][2]; \
d3gds3[1][1][0] = d3gds3[0][1][1]; d3gds3[1][1][1] = D3GDS1S1S1; \
d3gds3[1][1][2] = D3GDS1S1S2; \
d3gds3[1][2][0] = d3gds3[0][1][2]; d3gds3[1][2][1] = d3gds3[1][1][2]; \
d3gds3[1][2][2] = D3GDS1S2S2; \
d3gds3[2][0][0] = d3gds3[0][0][2]; d3gds3[2][0][1] = d3gds3[0][1][2]; \
d3gds3[2][0][2] = d3gds3[0][2][2]; \
d3gds3[2][1][0] = d3gds3[0][1][2]; d3gds3[2][1][1] = d3gds3[1][1][2]; \
d3gds3[2][1][2] = d3gds3[1][2][2]; \
d3gds3[2][2][0] = d3gds3[0][2][2]; d3gds3[2][2][1] = d3gds3[1][2][2]; \
d3gds3[2][2][2] = D3GDS2S2S2;

#define fillD3GDS2DT \
d3gds2dt[0][0] = D3GDS0S0DT;     d3gds2dt[0][1] = D3GDS0S1DT; \
d3gds2dt[0][2] = D3GDS0S2DT; \
d3gds2dt[1][0] = d3gds2dt[0][1]; d3gds2dt[1][1] = D3GDS1S1DT; \
d3gds2dt[1][2] = D3GDS1S2DT; \
d3gds2dt[2][0] = d3gds2dt[0][2]; d3gds2dt[2][1] = d3gds2dt[1][2]; \
d3gds2dt[2][2] = D3GDS2S2DT;

#define fillD3GDS2DP \
d3gds2dp[0][0] = D3GDS0S0DP;     d3gds2dp[0][1] = D3GDS0S1DP; \
d3gds2dp[0][2] = D3GDS0S2DP; \
d3gds2dp[1][0] = d3gds2dp[0][1]; d3gds2dp[1][1] = D3GDS1S1DP; \
d3gds2dp[1][2] = D3GDS1S2DP; \
d3gds2dp[2][0] = d3gds2dp[0][2]; d3gds2dp[2][1] = d3gds2dp[1][2]; \
d3gds2dp[2][2] = D3GDS2S2DP;

#define fillD3GDSDT2 \
d3gdsdt2[0] = D3GDS0DT2; d3gdsdt2[1] = D3GDS1DT2; d3gdsdt2[2] = D3GDS2DT2;

#define fillD3GDSDTDP \
d3gdsdtdp[0] = D3GDS0DTDP; d3gdsdtdp[1] = D3GDS1DTDP; \
d3gdsdtdp[2] = D3GDS2DTDP;

#define fillD3GDSDP2 \
d3gdsdp2[0] = D3GDS0DP2; d3gdsdp2[1] = D3GDS1DP2; d3gdsdp2[2] = D3GDS2DP2;

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
	int i, j, iter = 0;

	/* look-up or compute the current ordering state */
	if ( (t != tOld)       || (p != pOld) || (r[0] != rOld[0]) || (r[1] != rOld[1]) ) {
		double dgds[NS], sNew[NS];

		for (i=0; i<NS; i++) { sOld[i] = 2.0; dgds[i] = 0.0; }
		sNew[0] = 0.9*(1.0-r[0]-r[1]);
		sNew[1] = 0.9*r[0];
		sNew[2] = 0.9*r[1];

		while ( ((fabs(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) ||
				 (fabs(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON) ||
				 (fabs(sNew[2]-sOld[2]) > 10.0*DBL_EPSILON) ) && (iter < MAX_ITER)) {
			double s[NS];

			for (i=0; i<NS; i++) s[i] = sNew[i];

			xfe2M1 = (r[0] + s[1])/2.0;
			xmg2M1 = (r[1] + s[2])/2.0;
			xti4M1 = (1.0 - s[0] - s[1] - s[2])/2.0;
			xfe3M1 = (1.0 - r[0] - r[1] + s[0])/2.0;

			xfe2M2 = (r[0] - s[1])/4.0;
			xmg2M2 = (r[1] - s[2])/4.0;
			xti4M2 = (1.0 + 2.0*r[0] + 2.0*r[1] + s[0] + s[1] + s[2])/4.0;
			xfe3M2 = (3.0-3.0*r[0]-3.0*r[1]-s[0])/4.0;

			if (xfe2M1 <= 0.0) xfe2M1 = DBL_EPSILON; /* added in V1.0-7 */
			if (xmg2M1 <= 0.0) xmg2M1 = DBL_EPSILON; /* added in V1.0-7 */
			if (xti4M1 <= 0.0) xti4M1 = DBL_EPSILON; /* added in V1.0-7 */
			if (xfe3M1 <= 0.0) xfe3M1 = DBL_EPSILON; /* added in V1.0-7 */

			if (xfe2M2 <= 0.0) xfe2M2 = DBL_EPSILON; /* added in V1.0-7 */
			if (xmg2M2 <= 0.0) xmg2M2 = DBL_EPSILON; /* added in V1.0-7 */
			if (xti4M2 <= 0.0) xti4M2 = DBL_EPSILON; /* added in V1.0-7 */
			if (xfe3M2 <= 0.0) xfe3M2 = DBL_EPSILON; /* added in V1.0-7 */

			dgds[0] = DGDS0;
			dgds[1] = DGDS1;
			dgds[2] = DGDS2;

			invd2gds2[0][0] = D2GDS0S0;
			invd2gds2[0][1] = D2GDS0S1;
			invd2gds2[0][2] = D2GDS0S2;
			invd2gds2[1][0] = invd2gds2[0][1];
			invd2gds2[1][1] = D2GDS1S1;
			invd2gds2[1][2] = D2GDS1S2;
			invd2gds2[2][0] = invd2gds2[0][2];
			invd2gds2[2][1] = invd2gds2[1][2];
			invd2gds2[2][2] = D2GDS2S2;

			for (i=0; i<NS; i++) sOld[i] = s[i];

			[self gaussj:invd2gds2];

			for (i=0; i<NS; i++) {
				for(j=0; j<NS; j++) s[i] += - invd2gds2[i][j]*dgds[j];
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
					NSLog(@"ERROR in ORTHO-OXIDE.C (function ORDER). Failed to converge!\n");
					NSLog(@"  X Fe  = %13.6g, X2 Mg = %13.6g\n", r[0], r[1]);
					NSLog(@"  s1    = %13.6g, s2    = %13.6g, s3    = %13.6g\n", sOld[0], sOld[1], sOld[2]);
					NSLog(@"  dgds1 = %13.6g, dgds2 = %13.6g, dgds4 = %13.6g\n", dgds[0], dgds[1], dgds[2]);
					NSLog(@"  X Fe2+ M1: %13.6g  X Fe2+ M2: %13.6g\n", xfe2M1, xfe2M2);
					NSLog(@"  X Mg   M1: %13.6g  X MG   M2: %13.6g\n", xmg2M1, xmg2M2);
					NSLog(@"  X Ti   M1: %13.6g  X Ti   M2: %13.6g\n", xti4M1, xti4M2);
					NSLog(@"  X Fe3+ M1: %13.6g  X Fe3+ M2: %13.6g\n", xfe3M1, xfe3M2);
					break;
				}
			}
		}

	}

	if (mask & FIRST  ) {   /* return s        */
		for (i=0; i<NS; i++) s[i] = sOld[i];
	}
	if (mask & SECOND ) {   /* compute ds/dr:  */
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
		double d2gdsdt[NS];

		fillD2GDSDT

		for (i=0; i<NS; i++) {
			dt[i] = 0.0;
			for (j=0; j<NS; j++) dt[i] += - invd2gds2[i][j]*d2gdsdt[j];
		}
	}
	if (mask & FOURTH ) {   /* compute ds/dp:  */
		double d2gdsdp[NS];

		fillD2GDSDP

		for (i=0; i<NS; i++) {
			dp[i] = 0.0;
			for (j=0; j<NS; j++) dp[i] += - invd2gds2[i][j]*d2gdsdp[j];
		}
	}
	if (mask & FIFTH  ) {   /* compute d2s/dr2 */
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
	const char *NAMES[NA]    = { "pseudobrookite", "ferropseudobrookite", "karrooite" };
	const char *FORMULAS[NA] = { "Fe2TiO5", "FeTi2O5", "MgTi2O5" };
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
	 (2)  SECOND           THIRD | FOURTH | FIFTH | SIXTH | EIGHTH
	 (3)  THIRD            FOURTH | SEVENTH

	 (1) converts a vector of moles of elements into a vector of moles of
	 endmember rhm oxides components.
	 (2) calculates from a vector of moles of endmember components, one or
	 all of: r[], x[], dr[]/dm[] d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
	 (3) calculates from a vector of independent compositional variables
	 mole fractions of endmember components and/or the Jacobian matrix
	 dx[]/dr[]

	 In this routine it is assumed that the elements are in the order of atomic
	 numbers and that the order of rhm oxides components has been verified as:
	 m[0] = pseudobrookite      (Fe2TiO5) ,
	 m[1] = ferropseudobrookite (FeTi2O5),
	 m[2] = karrooite           (MgTi2O5),

	 ----------------------------------------------------------------------------*/

	int i, j, k;

	if (inpMask == FIRST && outMask == SECOND) {
		/* Converts a vector of moles of elements into a vector of moles of
		 end-member components.                                                 */
		double sumchg, fe2, fe3;
		static const int O  =  8;
		static const int Mg = 12;
		static const int Al = 13;
		static const int Ti = 22;
		static const int Cr = 24;
		static const int Mn = 25;
		static const int Fe = 26;
		static const int Co = 27;
		static const int Ni = 28;

		sumchg = 2.0*e[Mg] + 3.0*e[Al] + 4.0*e[Ti] + 3.0*e[Cr] + 2.0*e[Mn] + 2.0*e[Co] + 2.0*e[Ni];
		fe3 = 2.0*e[O] - sumchg - 2.0*e[Fe];
		if (fe3 < 0.0) fe3 = 0.0;
		fe2 = e[Fe] - fe3;

		m[0] = fe3/2.0; /* (fe3 + e[Al] + e[Cr])/2.0;   */  /* Moles of Fe2TiO5 */
		m[1] = fe2;     /* fe2 + e[Mn] + e[Co] + e[Ni]; */  /* Moles of FeTi2O5 */
		m[2] = e[Mg];                                       /* Moles of MgTi2O5 */

	} else if (inpMask == SECOND) {
		double sum;

		if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
			NSLog(@"Illegal call to conRhm with inpMask = %o and outMask = %o", inpMask, outMask);

		for (i=0, sum=0.0; i<NA; i++) sum += m[i];

		if (outMask & THIRD) {
			/* Converts a vector of moles of end-member components (m) into a vector
			 of independent compositional variables (r) required as input for the
			 remaining public functions.                                          */
			r[0] = (sum != 0.0) ? m[1]/sum : 0.0;  /* XFe = X FeTi2O5 */
			r[1] = (sum != 0.0) ? m[2]/sum : 0.0;  /* XMg = X MgTi2O5 */
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
					dm[0][j] = (j == 1) ? (1.0-m[1]/sum)/sum : -m[1]/SQUARE(sum);
					dm[1][j] = (j == 2) ? (1.0-m[2]/sum)/sum : -m[2]/SQUARE(sum);
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
						d2m[0][j][k]  = 2.0*m[1]/CUBE(sum);
						d2m[0][j][k] -= (j == 1) ? 1.0/SQUARE(sum) : 0.0;
						d2m[0][j][k] -= (k == 1) ? 1.0/SQUARE(sum) : 0.0;
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
							d3m[0][j][k][l]  = -6.0*m[1]/QUARTIC(sum);
							d3m[0][j][k][l] += (j == 1) ? 2.0/CUBE(sum) : 0.0;
							d3m[0][j][k][l] += (k == 1) ? 2.0/CUBE(sum) : 0.0;
							d3m[0][j][k][l] += (l == 1) ? 2.0/CUBE(sum) : 0.0;
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
			NSLog(@"Illegal call to conRhm with inpMask = %o and outMask = %o", inpMask, outMask);

		if (outMask & FOURTH) {
			/* Converts a vector of independent compositional variables (r) into a
			 vector of mole fractions of endmember components (x).                */
			x[0] = 1.0 - r[0] - r[1];
			x[1] = r[0];
			x[2] = r[1];
		}

		if (outMask & SEVENTH) {
			/* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
			for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
			dr[0][0] = -1.0; dr[0][1] = -1.0;
			dr[1][0] =  1.0;
			dr[2][1] =  1.0;
		}

	} else  {
		NSLog(@"Illegal call to conRhm with inpMask = %o and outMask = %o", inpMask, outMask);
	}

}

-(NSString *)displayFormula:(double)t
						  p:(double)p
						  r:(double [NA])r
{
	double totFe2, totMg, totFe3, totTi;

	totFe2 = r[0];
	totMg  = r[1];
	totFe3 = 2.0*(1.0-r[0]-r[1]);
	totTi  = 1.0+r[0]+r[1];

	return [NSString stringWithFormat:@"Fe''%4.2fMg%4.2fFe'''%4.2fTi%4.2fO5", totFe2, totMg, totFe3, totTi];
}

-(void)activity:(int)mask
			  t:(double)t
			  p:(double)p
			  r:(double [NA])r
			  a:(double [NA])a      // (pointer to a[]) activities              BINARY MASK: 0001
			 mu:(double [NA])mu     // (pointer to mu[]) chemical potentials    BINARY MASK: 0010
			 dx:(double [NA][NR])dx // (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100
{
	double s[NS], g, dgdr[NR];
	double fr[NA][NR];
	int i, j;

	for(i=0; i<NA; i++) {
		fr[i][0] = FR0(i); /* XFe */
		fr[i][1] = FR1(i); /* XMg */
	}

	[self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

	g       = G;
	dgdr[0] = DGDR0;
	dgdr[1] = DGDR1;

	if (mask & FIRST) {
		double a0[NA];

		[self pureEndMembers:FIRST t:(double)t	p:(double)p a:a0 mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		for(i=0; i<NA; i++) {
			for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
			a[i] = exp(a[i]/(R*t));
			if (a0[i] != 0.0) a[i] = a[i]/a0[i];
		}

	}

	if (mask & SECOND) {
		double mu0[NA];

		[self pureEndMembers:SECOND t:(double)t p:(double)p a:NULL mu:mu0 gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		for(i=0; i<NA; i++) {
			for (mu[i]=g-mu0[i], j=0; j<NR; j++) mu[i] += fr[i][j]*dgdr[j];
		}

	}

	if (mask & THIRD) {
		double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR],
		dfrdr[NA][NR], gs[NA][NS], dgsds[NA][NS], sum, a0[NA];
		int k, l;

		fillD2GDR2
		fillD2GDRDS
		fillD2GDS2

		[self pureEndMembers:FIRST t:(double)t	p:(double)p a:a0 mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		for(i=0; i<NA; i++) {
			gs[i][0] = GSS0(i);        /* s1  */
			gs[i][1] = GSS1(i);        /* s2  */
			gs[i][2] = GSS2(i);        /* s3  */
			dfrdr[i][0] = DFR0DR0(i);  /* XFe */
			dfrdr[i][1] = DFR1DR1(i);  /* XMg */
			dgsds[i][0] = DGSS0DS0(i); /* s1  */
			dgsds[i][1] = DGSS1DS1(i); /* s2  */
			dgsds[i][2] = DGSS2DS2(i); /* s3  */
		}

		[self order:SECOND t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

		for (i=0; i<NA; i++) {
			for (k=0; k<NR; k++) {
				/* compute activity of the i-th component */
				for (dx[i][k]=g, j=0; j<NR; j++) dx[i][k] += fr[i][j]*dgdr[j];
				dx[i][k] = exp(dx[i][k]/(R*t));
				if (a0[i] != 0.0) dx[i][k]  = dx[i][k]/a0[i];

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
		double ends[NA];

		*gmix = G;

		[self pureEndMembers:THIRD t:(double)t	p:(double)p a:NULL mu:NULL gmix:ends hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];
		*gmix -= ENDMEMBERS;
	}

	if(mask & SECOND) {
		double ends[NA];

		dx[0] = DGDR0;
		dx[1] = DGDR1;

		[self pureEndMembers:THIRD t:(double)t	p:(double)p a:NULL mu:NULL gmix:ends hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		dx[0] -= DENDDR0;
		dx[1] -= DENDDR1;
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
	double s[NS], ends[NA];

	[self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

	*hmix = (G) + t*(S);

	[self pureEndMembers:FOURTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:ends smix:NULL cpmix:NULL
				 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

	*hmix -= ENDMEMBERS;
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
		double ends[NA];

		*smix = S;

		[self pureEndMembers:FIFTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:ends cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		*smix -= ENDMEMBERS;
	}

	if(mask & SECOND) {
		double d2gdrds[NR][NS], d2gdrdt[NR], d2gds2[NS][NS], d2gdsdt[NS],
		dsdr[NS][NR], dsdt[NS], ends[NA];
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

		[self pureEndMembers:FIFTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:ends cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		dx[0] -= DENDDR0;
		dx[1] -= DENDDR1;
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
						+ d2gds2[k][l]*dsdr[k][j]*d2sdrdt[l][i];
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
	d2gdt2  = D2GDT2;

	if (mask & FIRST) {
		double ends[NA];

		*cpmix = d2gdt2;
		for (i=0; i<NS; i++) {
			*cpmix += 2.0*d2gdsdt[i]*dsdt[i];
			for (j=0; j<NS; j++) *cpmix += d2gds2[i][j]*dsdt[i]*dsdt[j];
		}
		*cpmix *= -t;

		[self pureEndMembers:SIXTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:ends
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		*cpmix -= ENDMEMBERS;
	}

	if(mask & SECOND) {
		double d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gdsdt2[NS], d2sdt2[NS],
		temp, ends[NA];
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

		[self pureEndMembers:SEVENTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:ends vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		*dt -= ENDMEMBERS;
	}

	if(mask & THIRD) {
		double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
		d3gds2dt[NS][NS], d2gdrds[NR][NS], d3gdrdt2[NR], d3gdsdt2[NS],
		dsdr[NS][NR], d2sdrdt[NS][NR], d2sdt2[NS], ends[NA];
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

		[self pureEndMembers:SIXTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:ends
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		dx[0] -= DENDDR0;
		dx[1] -= DENDDR1;
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
		double ends[NA];

		*vmix = DGDP;

		[self pureEndMembers:EIGHTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:ends vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		*vmix -= ENDMEMBERS;
	}

	if(mask & SECOND) {
		double d2gdrds[NR][NS], d2gdrdp[NR], d2gds2[NS][NS], d2gdsdp[NS],
		dsdr[NS][NR], dsdp[NS], ends[NA];
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

		[self pureEndMembers:EIGHTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:ends vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		dx[0] -= DENDDR0;
		dx[1] -= DENDDR1;
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
		double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], dsdt[NS], dsdp[NS],
		ends[NA];
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

		[self pureEndMembers:NINTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:ends vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		*dt -= ENDMEMBERS;
	}

	if(mask & FIFTH) {
		double d2gdp2 = D2GDP2;
		double d2gds2[NS][NS], d2gdsdp[NS], dsdp[NS], ends[NA];
		int i,j;

		fillD2GDS2
		fillD2GDSDP

		[self order:FOURTH t:t p:p r:r s:NULL dr:NULL dt:NULL dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

		*dp = d2gdp2;
		for (i=0; i<NS; i++) {
			*dp += 2.0*d2gdsdp[i]*dsdp[i];
			for (j=0; j<NS; j++) *dp += d2gds2[i][j]*dsdp[i]*dsdp[j];
		}

		[self pureEndMembers:TENTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:ends vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		*dp -= ENDMEMBERS;
	}

	if(mask & SIXTH) {
		double d3gdt2dp = D3GDT2DP;
		double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS],
		d3gds2dp[NS][NS], d3gds2dt[NS][NS], d3gdsdtdp[NS], d3gdsdt2[NS],
		dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS], ends[NA];
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

		[self pureEndMembers:ELEVENTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:ends vmixdtdp:NULL vmixdp2:NULL];

		*dt2 -= ENDMEMBERS;
	}

	if(mask & SEVENTH) {
		double d3gdtdp2 = D3GDTDP2;
		double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS],
		d3gds2dt[NS][NS], d3gds2dp[NS][NS], d3gdsdtdp[NS], d3gdsdp2[NS],
		dsdt[NS], dsdp[NS], d2sdtdp[NS], d2sdp2[NS], ends[NA];
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

		[self pureEndMembers:TWELFTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:ends vmixdp2:NULL];

		*dtdp -= ENDMEMBERS;
	}

	if(mask & EIGHTH) {
		double d3gdp3 = D3GDP3;
		double d2gds2[NS][NS], d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS],
		d3gdsdp2[NS], dsdp[NS], d2sdp2[NS], ends[NA];
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

		[self pureEndMembers:THIRTEENTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:ends];

		*dp2 -= ENDMEMBERS;
	}

	if(mask & NINTH) {
		double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
		d3gds2dp[NS][NS], d2gdrds[NR][NS], d3gdrdtdp[NR], d3gdsdtdp[NS],
		dsdt[NS], dsdp[NS], dsdr[NS][NR], d2sdrdt[NS][NR], d2sdrdp[NS][NR],
		d2gds2[NS][NS], d2gdsdt[NS], d3gdrdsdp[NR][NS], d2gdsdp[NS],
		d2sdtdp[NS], d3gds2dt[NS][NS], ends[NA];
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

		[self pureEndMembers:NINTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:ends vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		dxdt[0] -= DENDDR0;
		dxdt[1] -= DENDDR1;
	}

	if(mask & TENTH) {
		double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gds2dp[NS][NS],
		d2gdrds[NR][NS], dsdp[NS], dsdr[NS][NR], d2sdrdp[NS][NR], d2gds2[NS][NS],
		d3gdrdsdp[NR][NS], d3gdrdp2[NR], d3gdsdp2[NS], d2gdsdp[NS], d2sdp2[NS],
		ends[NA];
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

		[self pureEndMembers:TENTH t:(double)t p:(double)p a:NULL mu:NULL gmix:NULL hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:ends vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		dxdp[0] -= DENDDR0;
		dxdp[1] -= DENDDR1;
	}

}

-(NSUInteger)numberOfSolutionSpecies {
	return NA;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
	switch (index) {
		case 0:
			return @"pseudobrookite";
			break;
		case 1:
			return @"ferropseudobrookite";
			break;
		case 2:
			return @"karrooite";
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
	double deltaMu[NA], xNz[NA], x[NA], gamma[NA], xLast[NA], gammaLast[NA], affinity = 0.0;
	NSUInteger i, j, nz = 0, index[NA];

	BOOL debugS = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
	BOOL debugV = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
	if (debugV) NSLog(@"Entering [... affinityAndCompositionFromLiquidChemicalPotentialSum] ...");
    for (i=0; i<NA; i++) gamma[i] = 0.0;

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

		if (debugV) {
			NSLog(@"Iteration %lu", count);
			NSLog(@"%10.3g %10.3g %10.3g %10.3g", a[0], a[1], a[2], affinity);
			NSLog(@"%10.3g %10.3g %10.3g", x[0], x[1], x[2]);
            double g[NA];
			for (i=0, j=0; i<NA; i++) g[i] = (x[i] != 0.0) ? gamma[j++] : 0.0;
            NSLog(@"%13.6g %13.6g %13.6g", g[0], g[1], g[2]);
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
		NSLog(@"... Terminated (converged %@) for phase %@ in %lu iterations with delta affinity %f J for %f atoms.",
			  converged ? @"YES" : @"NO", [self phaseName], count, fabs(affinity-affinityLast), NATOMS);
		for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) NSLog(@"... ... Activity coefficient of component %@ is %f with mole fraction %f",
														 [[endmembers objectAtIndex:i] phaseName], gamma[j++], x[i]);
	}
	if (debugV) NSLog(@"Exiting [... affinityAndCompositionFromLiquidChemicalPotentialSum].");


	[results addObject:[NSNumber numberWithDouble:affinity]];                    // affinity in J
	for (i=0; i<NA; i++) [results addObject:[NSNumber numberWithDouble:x[i]]];   // composition in mole fraction of endmembers
	[results addObject:[NSNumber numberWithBool:converged]];                     // convergence flag
	[results addObject:[NSNumber numberWithUnsignedInteger:count]];                  // iteration count
	[results addObject:[NSNumber numberWithDouble:NATOMS]];                      // number of atoms used to scale affinity
	[results addObject:[NSNumber numberWithDouble:fabs(affinity-affinityLast)]]; // likely error in affinity

	return [NSArray arrayWithArray:results];
}

#import "SolutionPhase.h"

@end
