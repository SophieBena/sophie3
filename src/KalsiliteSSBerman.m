//
//  KalsiliteSSBerman.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/10/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import "KalsiliteSSBerman.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"

@implementation VacancyKalsilite

-(id)init {
	if ((self = [super init])) {
		highAlbite = [[BermanProperties alloc] initWithH:-3921618.0
													   S:224.412
													  k0:393.64
													  k1:-24.155E2
													  k2:-78.928E5
													  k3:107.064E7
													  v0:0.0
													  v1:0.0
													  v2:0.0
													  v3:0.0
													  v4:0.0];
		betaNepheline = [[BermanProperties alloc] initWithH:-2093004.0*4.0
														  S:124.641*4.0
														 k0:205.24*4.0
														 k1:-7.599E2*4.0
														 k2:-108.383E5*4.0
														 k3:208.182E7*4.0
														 l1:-50.249E-2*2.0
														 l2:165.95E-5*2.0
														 Tt:467.0
													 deltaH:241.0*4.0
														 v0:5.433181*8.0
														 v1:-2.0500e-6
														 v2:5.2000e-12
														 v3:31.802e-6
														 v4:213.0e-10];
	}
	return self;
}

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	return [highAlbite getGibbsFreeEnergyFromT:t andP:p] + [betaNepheline getGibbsFreeEnergyFromT:t andP:p]/2.0;
}
-(double)getEnthalpyFromT:(double)t andP:(double)p {
	return [highAlbite getEnthalpyFromT:t andP:p] + [betaNepheline getEnthalpyFromT:t andP:p]/2.0;
}
-(double)getEntropyFromT:(double)t andP:(double)p {
	return [highAlbite getEntropyFromT:t andP:p] + [betaNepheline getEntropyFromT:t andP:p]/2.0;
}
-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	return [highAlbite getHeatCapacityFromT:t andP:p] + [betaNepheline getHeatCapacityFromT:t andP:p]/2.0;
}
-(double)getDcpDtFromT:(double)t andP:(double)p {
	return [highAlbite getDcpDtFromT:t andP:p] + [betaNepheline getDcpDtFromT:t andP:p]/2.0;
}
-(double)getVolumeFromT:(double)t andP:(double)p {
	return [highAlbite getVolumeFromT:t andP:p] + [betaNepheline getVolumeFromT:t andP:p]/2.0;
}
-(double)getDvDtFromT:(double)t andP:(double)p {
	return [highAlbite getDvDtFromT:t andP:p] + [betaNepheline getDvDtFromT:t andP:p]/2.0;
}
-(double)getDvDpFromT:(double)t andP:(double)p {
	return [highAlbite getDvDpFromT:t andP:p] + [betaNepheline getDvDpFromT:t andP:p]/2.0;
}
-(double)getD2vDt2FromT:(double)t andP:(double)p {
	return [highAlbite getD2vDt2FromT:t andP:p] + [betaNepheline getD2vDt2FromT:t andP:p]/2.0;
}
-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	return [highAlbite getD2vDtDpFromT:t andP:p] + [betaNepheline getD2vDtDpFromT:t andP:p]/2.0;
}
-(double)getD2vDp2FromT:(double)t andP:(double)p {
	return [highAlbite getD2vDp2FromT:t andP:p] + [betaNepheline getD2vDp2FromT:t andP:p]/2.0;
}

@end

@implementation CalciumKalsilite

-(id)init {
	if ((self = [super init])) {
		anorthite = [[BermanProperties alloc] initWithH:-4228730.0+3.7*4184.0
													  S:200.186+3.7*4184.0/2200.0
													 k0:439.37
													 k1:-37.341E2
													 k2:0.0
													 k3:-31.702E7
													 v0:0.0
													 v1:0.0
													 v2:0.0
													 v3:0.0
													 v4:0.0];
		betaNepheline = [[BermanProperties alloc] initWithH:-2093004.0*4.0
														  S:124.641*4.0
														 k0:205.24*4.0
														 k1:-7.599E2*4.0
														 k2:-108.383E5*4.0
														 k3:208.182E7*4.0
														 l1:-50.249E-2*2.0
														 l2:165.95E-5*2.0
														 Tt:467.0
													 deltaH:241.0*4.0
														 v0:5.433181*8.0
														 v1:-2.0500e-6
														 v2:5.2000e-12
														 v3:31.802e-6
														 v4:213.0e-10];
	}
	return self;
}

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	return [anorthite getGibbsFreeEnergyFromT:t andP:p] + [betaNepheline getGibbsFreeEnergyFromT:t andP:p]/2.0 + 23096.0 -t*15.8765;
}
-(double)getEnthalpyFromT:(double)t andP:(double)p {
	return [anorthite getEnthalpyFromT:t andP:p] + [betaNepheline getEnthalpyFromT:t andP:p]/2.0 + 23096.0;
}
-(double)getEntropyFromT:(double)t andP:(double)p {
	return [anorthite getEntropyFromT:t andP:p] + [betaNepheline getEntropyFromT:t andP:p]/2.0 + 15.8765;
}
-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	return [anorthite getHeatCapacityFromT:t andP:p] + [betaNepheline getHeatCapacityFromT:t andP:p]/2.0;
}
-(double)getDcpDtFromT:(double)t andP:(double)p {
	return [anorthite getDcpDtFromT:t andP:p] + [betaNepheline getDcpDtFromT:t andP:p]/2.0;
}
-(double)getVolumeFromT:(double)t andP:(double)p {
	return [anorthite getVolumeFromT:t andP:p] + [betaNepheline getVolumeFromT:t andP:p]/2.0;
}
-(double)getDvDtFromT:(double)t andP:(double)p {
	return [anorthite getDvDtFromT:t andP:p] + [betaNepheline getDvDtFromT:t andP:p]/2.0;
}
-(double)getDvDpFromT:(double)t andP:(double)p {
	return [anorthite getDvDpFromT:t andP:p] + [betaNepheline getDvDpFromT:t andP:p]/2.0;
}
-(double)getD2vDt2FromT:(double)t andP:(double)p {
	return [anorthite getD2vDt2FromT:t andP:p] + [betaNepheline getD2vDt2FromT:t andP:p]/2.0;
}
-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	return [anorthite getD2vDtDpFromT:t andP:p] + [betaNepheline getD2vDtDpFromT:t andP:p]/2.0;
}
-(double)getD2vDp2FromT:(double)t andP:(double)p {
	return [anorthite getD2vDp2FromT:t andP:p] + [betaNepheline getD2vDp2FromT:t andP:p]/2.0;
}

@end

@implementation KalsiliteSSBerman

static NSArray *endmembers;

#define NR         3    /* Independent composition variables           */
#define NS         0    /* Ordering parameters s3 = X3                 */
#define NA         4    /* Endmember compositions                      */
#define NATOMS  28.0    /* Average number of atoms in the formula unit */

#pragma mark -
#pragma mark class methods

+(void)initialize {
	if (self == [KalsiliteSSBerman class]) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"Initialize(KalsiliteSSBerman) - entry ...");
		NSMutableArray *mutableEndmembers = [NSMutableArray arrayWithCapacity:NA];

		BermanProperties *nanepheline = [[BermanProperties alloc] initWithH:-2093004.0*4.0
																		  S:124.641*4.0
																		 k0:205.24*4.0
																		 k1:-7.599E2*4.0
																		 k2:-108.383E5*4.0
																		 k3:208.182E7*4.0
																		 l1:-50.249E-2*2.0
																		 l2:165.95E-5*2.0
																		 Tt:467.0
																	 deltaH:241.0*4.0
																		 v0:5.4131*4.0
																		 v1:-2.0500e-6
																		 v2:5.2000e-12
																		 v3:31.802e-6
																		 v4:213.0e-10];
		[nanepheline setPhaseFormula:@"Na4Al4Si4O16"];
		[nanepheline setPhaseName:@"na-nepheline"];
		[mutableEndmembers addObject:nanepheline];
		if (debug) NSLog(@"... allocated na-nepheline ...");

#define DH2     -1.35      * 1000.0 * 4.184 /* K4Al4Si4O16 kals-neph  joules */
#define DS2      0.0       * 1000.0 * 4.184 /* joules/T */
#define DV2     -0.00001   * 1000.0 * 4.184 /* joules/bar */

		BermanProperties *knepheline = [[BermanProperties alloc] initWithH:-2109563.55*4.0-(DH2)
																		 S:133.9653*4.0-(DS2)
																		k0:186.0*4.0
																		k1:0.0
																		k2:-131.067E5*4.0
																		k3:213.893E7*4.0
																		l1:-7.096454E-2*2.0
																		l2:21.682E-5*2.0
																		Tt:800.15
																	deltaH:1154.0*4.0
																		v0:6.043478*4.0-(DV2)
																		v1:-2.0500e-6
																		v2:5.2000e-12
																		v3:31.802e-6
																		v4:213.0e-10];
		[knepheline setPhaseFormula:@"K4Al4Si4O16"];
		[knepheline setPhaseName:@"k-nepheline"];
		[mutableEndmembers addObject:knepheline];
		if (debug) NSLog(@"... allocated k-nepheline ...");

#undef DH2
#undef DS2
#undef DV2

		VacancyKalsilite *vcnepheline = [[VacancyKalsilite alloc] init];
		[vcnepheline setPhaseFormula:@"Na3Al3Si5O16"];
		[vcnepheline setPhaseName:@"vc-nepheline"];
		[mutableEndmembers addObject:vcnepheline];
		if (debug) NSLog(@"... allocated vc-nepheline ...");

		CalciumKalsilite *canepheline = [[CalciumKalsilite alloc] init];
		[canepheline setPhaseFormula:@"CaNa2Al4Si4O16"];
		[canepheline setPhaseName:@"ca-nepheline"];
		[mutableEndmembers addObject:canepheline];
		if (debug) NSLog(@"... allocated ca-nepheline ...");

		endmembers = [NSArray arrayWithArray:mutableEndmembers];
	}
}

#pragma mark -
#pragma mark instance methods

-(id)init {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(KalsiliteSSBerman) ... entry ...");
		[self setPhaseName:@"Kalsilite ss"];
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

#define MAX_ITER 200    /* Maximum number of iterations allowed in order */

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Nepheline solution parameters:
 * Sack, R.O., Ghiorso, M.S. (1995)
 * Thermodynamics of Nepheline Solutions
 */

/* kals - neph */
#define DH1                  32145.67 /* Na4Al4Si4O16 joules */
#define DS1                     10.46 /* joules/T */
#define DV1                       0.0 /* joules/bar */
#define DH2                   -5648.4 /* K4Al4Si4O16 joules */
#define DS2                       0.0 /* joules/T */
#define DV2                     -0.04 /* joules/bar */
#define DH3                   14644.0 /* []Na3Al3Si5O16 joules */
#define DS3                     -18.7 /* joules/T */
#define DV3                       0.0 /* joules/bar */
#define DH4                   35184.0 /* []CaNa2Al2Si2O16 joules */
#define DS4                     -7.17 /* joules/T */
#define DV4                       0.0 /* joules/bar */

/* kalsilite structure */
#define WKALS                 29288.0 /* WNa4-K4 joules */
#define WKVKALS               30334.0 /* WK4-[]Na3 joules */
#define WNVKALS               14646.0 /* WNa4-[]Na3 joules */
#define WNCKALS                   0.0 /* WNa4-[]CaNa2 joules */
#define WKCKALS               14644.0 /* WK4-[]CaNa2 joules */
#define WVVKALS                   0.0 /* W[]Na3-[]CaNa2 joules */

/*
 * Global (to this file): variables
 */

#define R  8.3143

/* correction to zero point entropy of []CaNa2Al4Si4O16 */
#define S4      -15.8765
#define G4      -t*(S4)

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conSpn defines the conversion from m[i], to r[j]
 */
/* Order: X2, X3, X4 */
#define FR2(i)     (i == 1) ? 1.0 - r[0] : - r[0]
#define FR3(i)     (i == 2) ? 1.0 - r[1] : - r[1]
#define FR4(i)     (i == 3) ? 1.0 - r[2] : - r[2]

#define DFR2DR2(i) - 1.0
#define DFR3DR3(i) - 1.0
#define DFR4DR4(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */

#define S (DS1)*(1.0-r[0]-r[1]-r[2]) + (DS2)*r[0] + (DS3)*r[1] + ((DS4)+(S4))*r[2] \
- 4.0*R*(xk*log(xk) + (xvc-xca)*log(xvc-xca) + xna*log(xna) + xca*log(xca))
#define H (DH1)*(1.0-r[0]-r[1]-r[2]) + (DH2)*r[0] + (DH3)*r[1] + (DH4)*r[2] \
+ (WKALS)*r[0]*(1.0-r[0]-r[1]-r[2]) + (WNVKALS)*r[1]*(1.0-r[0]-r[1]-r[2]) \
+ (WNCKALS)*r[2]*(1.0-r[0]-r[1]-r[2]) + (WKVKALS)*r[0]*r[1] \
+ (WKCKALS)*r[0]*r[2] + (WVVKALS)*r[1]*r[2]
#define V (DV1)*(1.0-r[0]-r[1]-r[2]) + (DV2)*r[0] + (DV3)*r[1] + (DV4)*r[2]
#define G (H) - t*(S) + (p-1.0)*(V)

/*----------------------------------------------------------------------------*/

#define DGDR0 - ((DH1)-t*(DS1)+(p-1.0)*(DV1)) \
+ ((DH2)-t*(DS2)+(p-1.0)*(DV2)) + R*t*4.0*(log(xk) - log(xna)) \
+ (WKALS)*(1.0-2.0*r[0]-r[1]-r[2]) - (WNVKALS)*r[1] \
- (WNCKALS)*r[2] + (WKVKALS)*r[1] + (WKCKALS)*r[2]
#define DGDR1 - ((DH1)-t*(DS1)+(p-1.0)*(DV1)) \
+ ((DH3)-t*(DS3)+(p-1.0)*(DV3)) + R*t*(log(xvc-xca) - log(xna)) \
- (WKALS)*r[0] + (WNVKALS)*(1.0-r[0]-2.0*r[1]-r[2]) \
- (WNCKALS)*r[2] + (WKVKALS)*r[0] + (WVVKALS)*r[2]
#define DGDR2 - ((DH1)-t*(DS1)+(p-1.0)*(DV1)) \
+ ((DH4)-t*(DS4)+(p-1.0)*(DV4)) + (G4) \
+ R*t*(log(xca) - 2.0*log(xna) - 1.0) \
- (WKALS)*r[0] - (WNVKALS)*r[1] \
+ (WNCKALS)*(1.0-r[0]-r[1]-2.0*r[2]) \
+ (WKCKALS)*r[0] + (WVVKALS)*r[1]
#define DGDT  - (S)
#define DGDP  (V)

/*----------------------------------------------------------------------------*/

#define D2GDR0R0 R*t*4.0*(1.0/xk + 1.0/xna) - 2.0*(WKALS)
#define D2GDR0R1 R*t/xna - (WKALS) + (WKVKALS) - (WNVKALS)
#define D2GDR0R2 R*t*2.0/xna - (WKALS) - (WNCKALS) + (WKCKALS)
#define D2GDR0DT (DS1) - (DS2) + R*4.0*(log(xk) - log(xna))
#define D2GDR0DP - (DV1) + (DV2)

#define D2GDR1R1 R*t*0.25*(1.0/(xvc-xca) + 1.0/xna) - 2.0*(WNVKALS)
#define D2GDR1R2 R*t*(0.5/xna) - (WNVKALS) - (WNCKALS) + (WVVKALS)
#define D2GDR1DT (DS1) - (DS3) + R*(log(xvc-xca) - log(xna))
#define D2GDR1DP - (DV1) + (DV3)

#define D2GDR2R2 R*t*(0.25/xca + 1.0/xna) - 2.0*(WNCKALS)
#define D2GDR2DT (DS1) - (DS4) - (S4) + R*(log(xca) - 2.0*log(xna) - 1.0)
#define D2GDR2DP - (DV1) + (DV4)

#define D2GDT2   0.0
#define D2GDTDP  0.0
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/

#define D3GDR0R0R0 - R*t*4.0*(1.0/(xk*xk) - 1.0/(xna*xna))
#define D3GDR0R0R1 R*t/(xna*xna)

#define D3GDR0R0R2 R*t*2.0/(xna*xna)
#define D3GDR0R0DT R*4.0*(1.0/xk + 1.0/xna)
#define D3GDR0R0DP 0.0

#define D3GDR0R1R1 R*t*0.25/(xna*xna)
#define D3GDR0R1R2 R*t*0.5/(xna*xna)
#define D3GDR0R1DT R/xna
#define D3GDR0R1DP 0.0

#define D3GDR0R2R2 R*t/(xna*xna)
#define D3GDR0R2DT R*2.0/xna
#define D3GDR0R2DP 0.0

#define D3GDR1R1R1 - R*t*0.25*0.25*(1.0/((xvc-xca)*(xvc-xca)) - 1.0/(xna*xna))
#define D3GDR1R1R2 R*t*0.25*(0.5/(xna*xna))
#define D3GDR1R1DT R*0.25*(1.0/(xvc-xca) + 1.0/xna)
#define D3GDR1R1DP 0.0

#define D3GDR1R2R2 R*t*(0.5*0.5/(xna*xna))
#define D3GDR1R2DT R*(0.5/xna)
#define D3GDR1R2DP 0.0

#define D3GDR2R2R2 R*t*(- 0.25*0.25/(xca*xca) + 0.5/(xna*xna))
#define D3GDR2R2DT R*(0.25/xca + 1.0/xna)
#define D3GDR2R2DP 0.0

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
#define D3GDR2DT2  0.0
#define D3GDR2DTDP 0.0
#define D3GDR2DP2  0.0

/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1;     d2gdr2[0][2] = D2GDR0R2; \
d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;     d2gdr2[1][2] = D2GDR1R2; \
d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2]; d2gdr2[2][2] = D2GDR2R2;

#define fillD2GDRDT d2gdrdt[0] = D2GDR0DT;    d2gdrdt[1] = D2GDR1DT;    d2gdrdt[2] = D2GDR2DT;
#define fillD2GDRDP d2gdrdp[0] = D2GDR0DP;    d2gdrdp[1] = D2GDR1DP;    d2gdrdp[2] = D2GDR2DP;

#define fillD3GDR3 \
d3gdr3[0][0][0] = D3GDR0R0R0;		d3gdr3[0][0][1] = D3GDR0R0R1;      d3gdr3[0][0][2] = D3GDR0R0R2; \
d3gdr3[0][1][0] = d3gdr3[0][0][1];	d3gdr3[0][1][1] = D3GDR0R1R1;      d3gdr3[0][1][2] = D3GDR0R1R2; \
d3gdr3[0][2][0] = d3gdr3[0][0][2];     d3gdr3[0][2][1] = d3gdr3[0][1][2]; d3gdr3[0][2][2] = D3GDR0R2R2; \
d3gdr3[1][0][0] = d3gdr3[0][0][1];	d3gdr3[1][0][1] = d3gdr3[0][1][1]; d3gdr3[1][0][2] = d3gdr3[0][1][2]; \
d3gdr3[1][1][0] = d3gdr3[0][1][1];	d3gdr3[1][1][1] = D3GDR1R1R1;      d3gdr3[1][1][2] = D3GDR1R1R2; \
d3gdr3[1][2][0] = d3gdr3[0][1][2];	d3gdr3[1][2][1] = d3gdr3[1][1][2]; d3gdr3[1][2][2] = D3GDR1R2R2; \
d3gdr3[2][0][0] = d3gdr3[0][0][2];	d3gdr3[2][0][1] = d3gdr3[0][1][2]; d3gdr3[2][0][2] = d3gdr3[0][2][2]; \
d3gdr3[2][1][0] = d3gdr3[0][1][2];	d3gdr3[2][1][1] = d3gdr3[1][1][2]; d3gdr3[2][1][2] = d3gdr3[1][2][2]; \
d3gdr3[2][2][0] = d3gdr3[0][2][2];	d3gdr3[2][2][1] = d3gdr3[1][2][2]; d3gdr3[2][2][2] = D3GDR2R2R2; \

#define fillD3GDR2DT \
d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT;     d3gdr2dt[0][2] = D3GDR0R2DT; \
d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT;     d3gdr2dt[1][2] = D3GDR1R2DT; \
d3gdr2dt[2][0] = d3gdr2dt[0][2]; d3gdr2dt[2][1] = d3gdr2dt[1][2]; d3gdr2dt[2][2] = D3GDR2R2DT;

#define fillD3GDR2DP \
d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP;     d3gdr2dp[0][2] = D3GDR0R2DP; \
d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP;     d3gdr2dp[1][2] = D3GDR1R2DP; \
d3gdr2dp[2][0] = d3gdr2dp[0][2]; d3gdr2dp[2][1] = d3gdr2dp[1][2]; d3gdr2dp[2][2] = D3GDR2R2DP;

#define fillD3GDRDT2  d3gdrdt2[0] = D3GDR0DT2;   d3gdrdt2[1] = D3GDR1DT2;   d3gdrdt2[2] = D3GDR2DT2;
#define fillD3GDRDTDP d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP; d3gdrdtdp[2] = D3GDR2DTDP;
#define fillD3GDRDP2  d3gdrdp2[0] = D3GDR0DP2;   d3gdrdp2[1] = D3GDR1DP2;   d3gdrdp2[2] = D3GDR2DP2;

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
	const char *NAMES[NA]    = { "na-nepheline", "k-nepheline", "vc-nepheline", "ca-nepheline" };
	const char *FORMULAS[NA] = { "Na4Al4Si4O16", "K4Al4Si4O16", "Na3Al3Si5O16", "CaNa2Al4Si4O16" };
	BOOL result = YES;
	int i;
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

		result = result && (r[0] >= 0.0) && (r[0] <= 1.0-r[1]/4.0-r[2]/2.0);
		result = result && (r[1] >= 0.0) && (r[1] <= 1.0);
		result = result && (r[2] >= 0.0) && (r[2] <= 1.0);

		result = result && (4.0*r[0] >= 0.0) && (4.0*r[0] <= 4.0+sqrt(DBL_EPSILON));               /* tot K  */
		result = result && (4.0-4.0*r[0]-r[1]-2.0*r[2] >= 0.0)                                     /* tot Na */
		&& (4.0-4.0*r[0]-r[1]-2.0*r[2] <= 4.0+sqrt(DBL_EPSILON));
		result = result && (r[2] >= 0.0) && (r[2] <= 1.0+sqrt(DBL_EPSILON));                       /* tot Ca */
		result = result && (r[1]+r[2] >= 0.0) && (r[1]+r[2] <= 1.0+sqrt(DBL_EPSILON));             /* tot Vc */
		result = result && (4.0-r[1] >= 3.0-sqrt(DBL_EPSILON)) && (4.0-r[1] <= 5.0+sqrt(DBL_EPSILON));   /* tot Al */
		result = result && (4.0+r[1] >= 3.0-sqrt(DBL_EPSILON)) && (4.0+r[1] <= 5.0+sqrt(DBL_EPSILON));   /* tot Si */
	}
	/* Check bounds on moles of endmember components */
	if (mask & SIXTH) {
		for (i=0, sum=0.0; i<NA; i++) sum += m[i];
		result = result && (sum >= 0.0);
		if (sum > 0.0) {
			double totalNa = 4.0*m[0] + 3.0*m[2] + 2.0*m[3];
			double totalK  = 4.0*m[1];
			double totalCa = m[3];
			double totalVc = m[2] + m[3];
			double totalAl = 4.0*m[0] + 4.0*m[1] + 3.0*m[2] + 4.0*m[3];
			double totalSi = 4.0*m[0] + 4.0*m[1] + 5.0*m[2] + 4.0*m[3];
			// result = result && (m[0] >= -3.0*m[2]/4.0-m[3]/2.0-DBL_EPSILON) && (m[0] <= sum+DBL_EPSILON);
			// result = result && (m[1] >= 0.0) && (m[1] <= sum-m[2]/4.0-m[3]/2.0+DBL_EPSILON);
			// result = result && (m[2] >= 0.0) && (m[2] <= sum+DBL_EPSILON);
			// result = result && (m[3] >= 0.0) && (m[3] <= sum+DBL_EPSILON);
			result = result && (totalK  >= 0.0) && (totalK  <= 4.0*sum);   /* tot K  */
			result = result && (totalNa >= 0.0) && (totalNa <= 4.0*sum);   /* tot Na */
			result = result && (totalCa >= 0.0) && (totalCa <= sum);       /* tot Ca */
			result = result && (totalVc >= 0.0) && (totalVc <= sum);       /* tot Vc */
			result = result && (totalAl >= 3.0*sum) && (totalAl <= 4.0*sum); /* tot Al */
			result = result && (totalSi >= 4.0*sum) && (totalSi <= 5.0*sum); /* tot Si */
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
	 endmember spinel components.
	 (2) calculates from a vector of moles of endmember components, one or
	 all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
	 (3) calculates from a vector of independent compositional variables
	 mole fractions of endmember components and/or the Jacobian matrix
	 dx[]/dr[]

	 In this routine it is assumed that the elements are in the order of atomic
	 numbers and that the order of spinel components has been verified as:
	 m[0] = na-nepheline   (Na4Al4Si4O16) ,
	 m[1] = k-nepheline    (K4Al4Si4O16),
	 m[2] = vc-nepheline   (Na3Al3Si5O16),
	 m[3] = ca-nepheline   (CaNa2Al4Si4O16),

	 ----------------------------------------------------------------------------*/

	int i, j, k;

	if (inpMask == FIRST && outMask == SECOND) {
		/* Converts a vector of moles of elements into a vector of moles of
		 end-member components.                                                 */
		static const int O  =  8;
		static const int Na = 11;
		static const int Al = 13;
		static const int Si = 14;
		static const int K  = 19;
		static const int Ca = 20;
		double vacancy  = 3.0*e[O]/4.0 - e[Si] - e[Al] - e[Na] - e[K] - e[Ca];

		if (fabs(vacancy-e[Ca]) < 10.0*DBL_EPSILON) vacancy = e[Ca];

		m[0] = (e[Na] - 3.0*vacancy + e[Ca])/4.0;  /* Moles of Na4Al4Si4O16 */
		m[1] = e[K]/4.0;                           /* Moles of K4Al4Si4O16  */
		m[2] = vacancy-e[Ca];                      /* Moles of Na3Al3Si5O16 */
		m[3] = e[Ca];                              /* Moles of CaNa2Al4Si4O16 */

	} else if (inpMask == SECOND) {
		double sum;

		if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
			NSLog(@"Illegal call to conNph with inpMask = %o and outMask = %o\n",
				  inpMask, outMask);

		for (i=0, sum=0.0; i<NA; i++) sum += m[i];

		if (outMask & THIRD) {
			/* Converts a vector of moles of end-member components (m) into a vector
			 of independent compositional variables (r) required as input for the
			 remaining public functions.                                          */
			r[0] = (sum != 0.0) ? m[1]/sum : 0.0;  /* X2 = X K4Al4Si4O16    */
			r[1] = (sum != 0.0) ? m[2]/sum : 0.0;  /* X3 = X Na3Al3Si5O16   */
			r[2] = (sum != 0.0) ? m[3]/sum : 0.0;  /* X4 = X CaNa2Al4Si4O16 */
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
					dm[0][j] = (j == 1) ? (1.0 - m[1]/sum)/sum : -m[1]/SQUARE(sum);
					dm[1][j] = (j == 2) ? (1.0 - m[2]/sum)/sum : -m[2]/SQUARE(sum);
					dm[2][j] = (j == 3) ? (1.0 - m[3]/sum)/sum : -m[3]/SQUARE(sum);
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
						d2m[2][j][k]  = 2.0*m[3]/CUBE(sum);
						d2m[2][j][k] -= (j == 3) ? 1.0/SQUARE(sum) : 0.0;
						d2m[2][j][k] -= (k == 3) ? 1.0/SQUARE(sum) : 0.0;
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
							d3m[2][j][k][l]  = -6.0*m[3]/QUARTIC(sum);
							d3m[2][j][k][l] += (j == 3) ? 2.0/CUBE(sum) : 0.0;
							d3m[2][j][k][l] += (k == 3) ? 2.0/CUBE(sum) : 0.0;
							d3m[2][j][k][l] += (l == 3) ? 2.0/CUBE(sum) : 0.0;
						}
					}
				}
			}
		}

	} else if (inpMask == THIRD) {

		if (outMask & ~(FOURTH | SEVENTH))
			NSLog(@"Illegal call to conNph with inpMask = %o and outMask = %o\n",
				  inpMask, outMask);

		if (outMask & FOURTH) {
			/* Converts a vector of independent compositional variables (r) into a
			 vector of mole fractions of endmember components (x).                */
			x[0] = 1.0 - r[0] - r[1] - r[2];
			x[1] = r[0];
			x[2] = r[1];
			x[3] = r[2];
		}

		if (outMask & SEVENTH) {
			/* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
			for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
			dr[0][0] = -1.0; dr[0][1] = -1.0; dr[0][2] = -1.0;
			dr[1][0] =  1.0;
			dr[2][1] =  1.0;
			dr[3][2] =  1.0;
		}

	} else  {
		NSLog(@"Illegal call to conNph with inpMask = %o and outMask = %o\n",
			  inpMask, outMask);
	}

}

-(NSString *)displayFormula:(double)t
						  p:(double)p
						  r:(double [NA])r
{
	double totNa, totK, totCa, totVc, totAl, totSi;

	totNa = 4.0*(1.0-r[0]-r[1]-r[2]) + 3.0*r[1] + 2.0*r[2];
	totK  = 4.0*r[0];
	totCa = r[2];
	totVc = r[1] + r[2];
	totAl = 4.0*(1.0-r[0]-r[1]-r[2]) + 4.0*r[0] + 3.0*r[1] + 4.0*r[2];
	totSi = 4.0*(1.0-r[0]-r[1]-r[2]) + 4.0*r[0] + 5.0*r[1] + 4.0*r[2];

	return [NSString stringWithFormat:@"Na%4.2fK%4.2fCa%4.2f[]%4.2fAl%4.2fSi%4.2fO16",
			totNa, totK, totCa, totVc, totAl, totSi];
}

-(void)activity:(int)mask
			  t:(double)t
			  p:(double)p
			  r:(double [NA])r
			  a:(double [NA])a      // (pointer to a[]) activities              BINARY MASK: 0001
			 mu:(double [NA])mu     // (pointer to mu[]) chemical potentials    BINARY MASK: 0010
			 dx:(double [NA][NR])dx // (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100
{
	double g, dgdr[NR];
	double fr[NA][NR];
	int i, j;
	double xk  = r[0];
	double xvc = (r[1]+r[2])/4.0;
	double xna = 1.0 - r[0] - r[1]/4.0 - r[2]/2.0;
	double xca = r[2]/4.0;

	if (xk  <= DBL_EPSILON) xk  = DBL_EPSILON;
	if (xvc <= DBL_EPSILON) xvc = DBL_EPSILON;
	if (xna <= DBL_EPSILON) xna = DBL_EPSILON;
	if (xca <= DBL_EPSILON) xca = DBL_EPSILON;
	if (xk  >= 1.0 - DBL_EPSILON) xk  = 1.0 - DBL_EPSILON;
	if (xvc >= 1.0 - DBL_EPSILON) xvc = 1.0 - DBL_EPSILON;
	if (xna >= 1.0 - DBL_EPSILON) xna = 1.0 - DBL_EPSILON;
	if (xca >= 1.0 - DBL_EPSILON) xca = 1.0 - DBL_EPSILON;
	if ((xvc-xca) <= 0.0) xvc = xca + DBL_EPSILON;

	for(i=0; i<NA; i++) {
		fr[i][0] = FR2(i); /* X2 */
		fr[i][1] = FR3(i); /* X3 */
		fr[i][2] = FR4(i); /* X4 */
	}

	g       = G;
	dgdr[0] = DGDR0;
	dgdr[1] = DGDR1;
	dgdr[2] = DGDR2;

	if (mask & FIRST) {
		for(i=0; i<NA; i++) {
			for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
			a[i] = exp(a[i]/(R*t));
		}
	}

	if (mask & SECOND) {
		for(i=0; i<NA; i++)
			for (mu[i]=g, j=0; j<NR; j++) mu[i] += fr[i][j]*dgdr[j];
	}

	if (mask & THIRD) {
		double d2gdr2[NR][NR], dfrdr[NA][NR], sum;
		int k;

		fillD2GDR2

		for(i=0; i<NA; i++) {
			dfrdr[i][0] = DFR2DR2(i); /* X2 */
			dfrdr[i][1] = DFR3DR3(i); /* X3 */
			dfrdr[i][2] = DFR4DR4(i); /* X4 */
		}

		for (i=0; i<NA; i++) {
			for (k=0; k<NR; k++) {
				/* compute activity of the i-th component */
				for (dx[i][k]=g, j=0; j<NR; j++) dx[i][k] += fr[i][j]*dgdr[j];
				dx[i][k] = exp(dx[i][k]/(R*t));

				/* compute derivative of i-th activity with respect to r(k) */
				sum = (1.0+dfrdr[i][k])*dgdr[k];
				for (j=0; j<NR; j++) sum += fr[i][j]*d2gdr2[j][k];
				dx[i][k] *= sum/(R*t);
			}
		}
	}

	if (mask & FOURTH) {
		/* implement exclusion criteria on quantities for preclb routines         */
		static const double exclusion[NA] = {
			0.05,  /* exclusion criteria on the mole fraction of Na+  */
			0.05,  /* exclusion criteria on the mole fraction of K+   */
			0.05,  /* exclusion criteria on the mole fraction of Ca++ */
			0.05,  /* exclusion criteria on the mole fraction of Vc   */
		};
		double x[NA], totNa, totK, totCa, totVc;

		totNa = 4.0*(1.0-r[0]-r[1]-r[2]) + 3.0*r[1] + 2.0*r[2];
		totK  = 4.0*r[0];
		totCa = r[2];
		totVc = r[1] + r[2];

		x[0] = totNa/4.0;   /* Na+ averaged on both sites */
		x[1] = totK/4.0;    /* K+ averaged on both sites  */
		x[2] = totVc-totCa; /* Vc on large site           */
		x[3] = totCa;       /* Ca on small site           */

		for (i=0; i<NA; i++) {
			if (x[i] < exclusion[i]) {
				if (mask & FIRST)  a[i]  = 0.0;
				if (mask & SECOND) mu[i] = 0.0;
				if (mask & THIRD)  for (j=0; j<NR; j++) dx[i][j] = 0.0;
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
	double xk  = r[0];
	double xvc = (r[1]+r[2])/4.0;
	double xna = 1.0 - r[0] - r[1]/4.0 - r[2]/2.0;
	double xca = r[2]/4.0;

	if (xk  <= DBL_EPSILON) xk  = DBL_EPSILON;
	if (xvc <= DBL_EPSILON) xvc = DBL_EPSILON;
	if (xna <= DBL_EPSILON) xna = DBL_EPSILON;
	if (xca <= DBL_EPSILON) xca = DBL_EPSILON;
	if (xk  >= 1.0 - DBL_EPSILON) xk  = 1.0 - DBL_EPSILON;
	if (xvc >= 1.0 - DBL_EPSILON) xvc = 1.0 - DBL_EPSILON;
	if (xna >= 1.0 - DBL_EPSILON) xna = 1.0 - DBL_EPSILON;
	if (xca >= 1.0 - DBL_EPSILON) xca = 1.0 - DBL_EPSILON;
	if ((xvc-xca) <= 0.0) xvc = xca + DBL_EPSILON;

	if (mask & FIRST) {
		*gmix = G;
	}

	if(mask & SECOND) {
		dx[0] = DGDR0;
		dx[1] = DGDR1;
		dx[2] = DGDR2;
	}

	if(mask & THIRD) {
		double d2gdr2[NR][NR];
		int i, j;

		fillD2GDR2
		for (i=0; i<NR; i++) for (j=0; j<NR; j++) dx2[i][j] = d2gdr2[i][j];
	}

	if(mask & FOURTH) {
		double d3gdr3[NR][NR][NR];
		int i, j, k;

		fillD3GDR3
		for (i=0; i<NR; i++) for (j=0; j<NR; j++) for (k=0; k<NR; k++) dx3[i][j][k] = d3gdr3[i][j][k];
	}

}

-(void)hmix:(int)mask
		  t:(double)t
		  p:(double)p
		  r:(double [NR])r
	   hmix:(double *)hmix // Enthalpy of mixing BINARY MASK: 1
/*  This function calculates enthalpy of mixing corrected to 1 bar. */
{
	double xk  = r[0];
	double xvc = (r[1]+r[2])/4.0;
	double xna = 1.0 - r[0] - r[1]/4.0 - r[2]/2.0;
	double xca = r[2]/4.0;

	if (xk  <= DBL_EPSILON) xk  = DBL_EPSILON;
	if (xvc <= DBL_EPSILON) xvc = DBL_EPSILON;
	if (xna <= DBL_EPSILON) xna = DBL_EPSILON;
	if (xca <= DBL_EPSILON) xca = DBL_EPSILON;
	if (xk  >= 1.0 - DBL_EPSILON) xk  = 1.0 - DBL_EPSILON;
	if (xvc >= 1.0 - DBL_EPSILON) xvc = 1.0 - DBL_EPSILON;
	if (xna >= 1.0 - DBL_EPSILON) xna = 1.0 - DBL_EPSILON;
	if (xca >= 1.0 - DBL_EPSILON) xca = 1.0 - DBL_EPSILON;
	if ((xvc-xca) <= 0.0) xvc = xca + DBL_EPSILON;

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
	double xk  = r[0];
	double xvc = (r[1]+r[2])/4.0;
	double xna = 1.0 - r[0] - r[1]/4.0 - r[2]/2.0;
	double xca = r[2]/4.0;

	if (xk  <= DBL_EPSILON) xk  = DBL_EPSILON;
	if (xvc <= DBL_EPSILON) xvc = DBL_EPSILON;
	if (xna <= DBL_EPSILON) xna = DBL_EPSILON;
	if (xca <= DBL_EPSILON) xca = DBL_EPSILON;
	if (xk  >= 1.0 - DBL_EPSILON) xk  = 1.0 - DBL_EPSILON;
	if (xvc >= 1.0 - DBL_EPSILON) xvc = 1.0 - DBL_EPSILON;
	if (xna >= 1.0 - DBL_EPSILON) xna = 1.0 - DBL_EPSILON;
	if (xca >= 1.0 - DBL_EPSILON) xca = 1.0 - DBL_EPSILON;
	if ((xvc-xca) <= 0.0) xvc = xca + DBL_EPSILON;

	if (mask & FIRST) {
		*smix = S;
	}

	if(mask & SECOND) {
		double d2gdrdt[NR];
		int i;

		fillD2GDRDT
		for (i=0; i<NR; i++) dx[i] = -d2gdrdt[i];
	}

	if(mask & THIRD) {
		double d3gdr2dt[NR][NR];
		int i, j;

		fillD3GDR2DT

		for (i=0; i<NR; i++) for (j=0; j<NR; j++) dx2[i][j] = -d3gdr2dt[i][j];
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
	double d2gdt2;
	int i;

	d2gdt2  = D2GDT2;

	if (mask & FIRST) {
		*cpmix = -t*d2gdt2;
	}

	if(mask & SECOND) {
		double d3gdt3 = D3GDT3;

		*dt = -t*d3gdt3 - d2gdt2;
	}

	if(mask & THIRD) {
		double d3gdrdt2[NR];

		fillD3GDRDT2
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

		fillD2GDRDP
		for (i=0; i<NR; i++) dx[i] = d2gdrdp[i];
	}

	if(mask & THIRD) {
		double d3gdr2dp[NR][NR];
		int i, j;

		fillD3GDR2DP
		for (i=0; i<NR; i++) for (j=0; j<NR; j++) dx2[i][j] = d3gdr2dp[i][j];
	}

	if(mask & FOURTH) {
		double d2gdtdp = D2GDTDP;

		*dt = d2gdtdp;
	}

	if(mask & FIFTH) {
		double d2gdp2 = D2GDP2;

		*dp = d2gdp2;
	}

	if(mask & SIXTH) {
		double d3gdt2dp = D3GDT2DP;

		*dt2 = d3gdt2dp;
	}

	if(mask & SEVENTH) {
		double d3gdtdp2 = D3GDTDP2;

		*dtdp = d3gdtdp2;
	}

	if(mask & EIGHTH) {
		double d3gdp3 = D3GDP3;

		*dp2 = d3gdp3;
	}

	if(mask & NINTH) {
		double d3gdrdtdp[NR];
		int i;

		fillD3GDRDTDP
		for (i=0; i<NR; i++) dxdt[i]=d3gdrdtdp[i];
	}

	if(mask & TENTH) {
		double d3gdrdp2[NR];
		int i;

		fillD3GDRDP2
		for (i=0; i<NR; i++) dxdp[i]=d3gdrdp2[i];
	}

}

#define NAS 6

-(NSUInteger)numberOfSolutionSpecies {
	return NAS;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
	switch (index) {
		case 0:
			return @"na-nepheline";
			break;
		case 1:
			return @"k-nepheline";
			break;
		case 2:
			return @"vc-nepheline";
			break;
		case 3:
			return @"ca-nepheline";
			break;
		case 4:
			return @"vc-k-nepheline";
			break;
		case 5:
			return @"ca-k-nepheline";
			break;
		default:
			return @"";
			break;
	}
}

-(DoubleVector *)convertMolesOfSpeciesToMolesOfComponents:(double *)mSpecies {
	DoubleVector *mComponentWrapper =[[DoubleVector alloc] initWithSize:NA];
	double *mComponents = [mComponentWrapper pointerToDouble];
	mComponents[0] = mSpecies[0] - 3.0*mSpecies[4]/4.0 - mSpecies[5]/2.0;
	mComponents[1] = mSpecies[1] + 3.0*mSpecies[4]/4.0 + mSpecies[5]/2.0;
	mComponents[2] = mSpecies[2] + mSpecies[4];
	mComponents[3] = mSpecies[3] + mSpecies[5];
	return mComponentWrapper;
}

-(DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
	DoubleVector *muSpeciesWrapper = [[DoubleVector alloc] initWithSize:NAS];
	double *muSpecies = [muSpeciesWrapper pointerToDouble];
	double *muComponents = [[self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p] pointerToDouble];
	for (NSUInteger i=0; i<NA; i++) muSpecies[i] = muComponents[i];
	// []K3Al3Si5O16
	muSpecies[4] = muSpecies[2] + 3.0*(muSpecies[1] - muSpecies[0])/4.0;
	// []CaK2Al4Si4O16
	muSpecies[5] = muSpecies[3] +     (muSpecies[1] - muSpecies[0])/2.0;
	return muSpeciesWrapper;
}

-(DoubleVector *)elementalCompositionOfSpeciesAtIndex:(NSUInteger)index {
	DoubleVector *speciesElementArrayWrapper = nil;

	if (index < NA) speciesElementArrayWrapper = [[endmembers objectAtIndex:index] formulaAsElementArray];
	else if (index == 4) {
		speciesElementArrayWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
		double *speciesElementArray = [speciesElementArrayWrapper pointerToDouble];
        double *component = [[[endmembers objectAtIndex:2] formulaAsElementArray] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] += component[i];
		component = [[[endmembers objectAtIndex:1] formulaAsElementArray] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] += 3.0*component[i]/4.0;
		component = [[[endmembers objectAtIndex:0] formulaAsElementArray] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] -= 3.0*component[i]/4.0;
	} else if (index == 5) {
		speciesElementArrayWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
		double *speciesElementArray = [speciesElementArrayWrapper pointerToDouble];
        double *component = [[[endmembers objectAtIndex:3] formulaAsElementArray] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] += component[i];
		component = [[[endmembers objectAtIndex:1] formulaAsElementArray] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] += component[i]/2.0;
		component = [[[endmembers objectAtIndex:0] formulaAsElementArray] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] -= component[i]/2.0;
	}
	return speciesElementArrayWrapper;
}

-(void)correctActivityCoefficients:(double [NAS])gamma forComposition:(double [NAS])x { }

// --> SolutionPhaseProtocol public function
-(NSArray *)affinityAndCompositionFromLiquidChemicalPotentialSum:(double *)chemicalPotentials andT:(double)t andP:(double)p {
	NSMutableArray *results = [NSMutableArray arrayWithCapacity:NA+1];
	double mu0[NAS], deltaMu[NAS], xNz[NAS], x[NAS], gamma[NAS], xLast[NAS], gammaLast[NAS], affinity = 0.0;
	NSUInteger i, j, nz = 0, index[NAS];
	double adjustmentCoefficient[NAS] = { 4.0, 4.0, 4.0, 4.0, 4.0, 4.0 }, reducedAdjustmentCoefficient[NAS];

	BOOL debugS = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
	BOOL debugV = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
	if (debugV) NSLog(@"Entering [... affinityAndCompositionFromLiquidChemicalPotentialSum] ...");

	// Compute solid -> liquid delta mus and deflate composition space
	for (i=0; i<NA; i++) {
		if (chemicalPotentials[i] != 0.0) {
			mu0[i] = [[endmembers objectAtIndex:i] getGibbsFreeEnergyFromT:t andP:p];
			deltaMu[nz] = chemicalPotentials[i] - mu0[i];
			index[nz] = i;
			gamma[nz] = 1.0;
			reducedAdjustmentCoefficient[nz] = adjustmentCoefficient[i];
			nz++;
		} else mu0[i] = 0.0;
		x[i] = 0.0;
		xLast[i] = 0.0;
	}

	// []K3Al3Si5O16
	if ((chemicalPotentials[0] != 0.0) && (chemicalPotentials[1] != 0.0) && (chemicalPotentials[2] != 0.0)) {
		double r[NR] = { 0.75, 1.0, 0.0};
		mu0[4] = mu0[2] + 3.0*(mu0[1] - mu0[0])/4.0 + (H);
		deltaMu[nz] = chemicalPotentials[2] + 3.0*(chemicalPotentials[1] - chemicalPotentials[0])/4.0 - mu0[4];
		index[nz] = 4;
		gamma[nz] = 1.0;
		reducedAdjustmentCoefficient[nz] = adjustmentCoefficient[4];
		nz++;
	} else mu0[4] = 0.0;
	x[4] = 0.0;
	xLast[4] = 0.0;

	// []CaK2Al4Si4O16
	if ((chemicalPotentials[0] != 0.0) && (chemicalPotentials[1] != 0.0) && (chemicalPotentials[3] != 0.0)) {
		double r[NR] = { 0.5, 0.0, 1.0};
		mu0[5] = mu0[3] + (mu0[1] - mu0[0])/2.0 + (H);
		deltaMu[nz] = chemicalPotentials[3] + (chemicalPotentials[1] - chemicalPotentials[0])/2.0 - mu0[5];
		index[nz] = 5;
		gamma[nz] = 1.0;
		reducedAdjustmentCoefficient[nz] = adjustmentCoefficient[5];
		nz++;
	} else mu0[5] = 0.0;
	x[5] = 0.0;
	xLast[5] = 0.0;

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
				xNz[i] = exp(((deltaMu[i]-R*t*reducedAdjustmentCoefficient[i]*log(gamma[i]))
							  -(deltaMu[nz-2]-R*t*reducedAdjustmentCoefficient[nz-2]*log(gamma[nz-2])))/(reducedAdjustmentCoefficient[i]*R*t));
				sum += xNz[i];
			}
			xNz[nz-2] = exp(((deltaMu[nz-2]-R*t*reducedAdjustmentCoefficient[nz-2]*log(gamma[nz-2]))
							 -(deltaMu[nz-1]-R*t*reducedAdjustmentCoefficient[nz-1]*log(gamma[nz-1])))/(reducedAdjustmentCoefficient[nz-2]*R*t));

			xNz[nz-2] /= 1.0 + xNz[nz-2]*sum;
			xNz[nz-1] = 1.0 - xNz[nz-2];
			if (nz > 2) for (i=0; i<(nz-2); i++) {
				xNz[i] *= xNz[nz-2];
				xNz[nz-1] -= xNz[i];
			}
			for (i=0; i<nz; i++) if (xNz[i] <= DBL_EPSILON) xNz[i] = DBL_EPSILON;

			/* compute the chemical affinity (choice of mu[] is arbitrary) */
			affinity = -(deltaMu[0]-R*t*reducedAdjustmentCoefficient[0]*log(gamma[0])) + R*t*reducedAdjustmentCoefficient[0]*log(xNz[0]);
		}

		// Reinflate the solution
		for (i=0; i<nz; i++) x[index[i]] = xNz[i];

		// Determine activity coefficients
		double a[NAS], mu[NA], r[NR];
		xReduced[0] = x[0] - 3.0*x[4]/4.0 - x[5]/2.0;
		xReduced[1] = x[1] + 3.0*x[4]/4.0 + x[5]/2.0;
		xReduced[2] = x[2] + x[4];
		xReduced[3] = x[3] + x[5];
		[self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:xReduced r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
		[self activity:FIRST | SECOND t:t p:p r:r a:a mu:mu dx:NULL];
		if (![self testPermissibleValuesOfComponents:xReduced]) {
			NSLog(@"Composition estimate is infeasible.");
			NSLog(@"species X0 %@ = %g", [[self nameOfSolutionSpeciesAtIndex:0] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[0]);
			NSLog(@"species X1 %@ = %g", [[self nameOfSolutionSpeciesAtIndex:1] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[1]);
			NSLog(@"species X2 %@ = %g", [[self nameOfSolutionSpeciesAtIndex:2] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[2]);
			NSLog(@"species X3 %@ = %g", [[self nameOfSolutionSpeciesAtIndex:3] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[3]);
			NSLog(@"species X4 %@ = %g", [[self nameOfSolutionSpeciesAtIndex:4] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[4]);
			NSLog(@"species X5 %@ = %g", [[self nameOfSolutionSpeciesAtIndex:5] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[5]);
			NSLog(@"%@ X0 = %13.6g, a = %13.6g, mu = %13.6g", [[self nameOfSolutionSpeciesAtIndex:0] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
				  xReduced[0], a[0], mu[0]);
			NSLog(@"%@ X1 = %13.6g, a = %13.6g, mu = %13.6g", [[self nameOfSolutionSpeciesAtIndex:1] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
				  xReduced[1], a[1], mu[1]);
			NSLog(@"%@ X2 = %13.6g, a = %13.6g, mu = %13.6g", [[self nameOfSolutionSpeciesAtIndex:2] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
				  xReduced[2], a[2], mu[2]);
			NSLog(@"%@ X3 = %13.6g, a = %13.6g, mu = %13.6g", [[self nameOfSolutionSpeciesAtIndex:3] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
				  xReduced[3], a[3], mu[3]);
		}

		for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) gamma[j++] = pow(a[i], 1.0/adjustmentCoefficient[i])/x[i];

		if ((chemicalPotentials[0] != 0.0) && (chemicalPotentials[1] != 0.0) && (chemicalPotentials[2] != 0.0)) {
			a[4] = exp((mu[2] + 3.0*(mu[1] - mu[0])/4.0 + mu0[2] + 3.0*(mu0[1] - mu0[0])/4.0 - mu0[4])/(R*t));
			gamma[j++] = pow(a[4], 1.0/adjustmentCoefficient[4])/x[4];
		}
		if ((chemicalPotentials[0] != 0.0) && (chemicalPotentials[1] != 0.0) && (chemicalPotentials[3] != 0.0)) {
			a[5] = exp((mu[3] + (mu[1] - mu[0])/2.0 + mu0[3] + (mu0[1] - mu0[0])/2.0 - mu0[5])/(R*t));
			gamma[j]   = pow(a[5], 1.0/adjustmentCoefficient[5])/x[5];
		}
		[self correctActivityCoefficients:gamma forComposition:x];

		if (debugV) {
			NSLog(@"Iteration %lu", count);
			NSLog(@"%13.6g %13.6g %13.6g %13.6g %13.6g %13.6g %13.6g", a[0], a[1], a[2], a[3], a[4], a[5], affinity);
			NSLog(@"%13.6g %13.6g %13.6g %13.6g %13.6g %13.6g", x[0], x[1], x[2], x[3], x[4], x[5]);
			double g[NAS];
			for (i=0, j=0; i<NAS; i++) g[i] = (x[i] != 0.0) ? gamma[j++] : 0.0;
			NSLog(@"%13.6g %13.6g %13.6g %13.6g %13.6g %13.6g", g[0], g[1], g[2], g[3], g[4], g[5]);
		}

		if (count > 0) { // pure empiricism
			for (i=0; i<nz; i++) gamma[i] = (gamma[i]+gammaLast[i])/2.0;
		}
		for (i=0; i<nz; i++) gammaLast[i] = gamma[i];

		converged = (fabs(affinity-affinityLast) < 0.1);
		count++;

	} while (count < 100 && !converged);

	if (debugS) {
		NSLog(@"... Terminated (converged %@) for phase %@ in %lu iterations with affinity %f J (delta %f) for %f atoms.",
			  converged ? @"YES" : @"NO", [self phaseName], count, affinity, fabs(affinity-affinityLast), NATOMS);
		for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) NSLog(@"... ... Activity coefficient of %@ is %13.6g with mole fraction %13.6g (reduced %13.6g)",
														 [[[endmembers objectAtIndex:i] phaseName] stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
														 gamma[j++], x[i], xReduced[i]);
		NSLog(@"... ... Activity coefficient of %@ is %13.6g with mole fraction %13.6g", @"vc-k-nepheline ", gamma[j++], x[4]);
		NSLog(@"... ... Activity coefficient of %@ is %13.6g with mole fraction %13.6g", @"ca-k-nepheline ", gamma[j],   x[5]);
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
