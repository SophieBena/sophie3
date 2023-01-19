//
//  NephelineSSBerman.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/7/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import "NephelineSSBerman.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"
#import "IntegerVector.h"
#import "MathSupport.h"

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

@implementation VacancyNepheline

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

@implementation CalciumNepheline

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

@implementation NephelineSSBerman

static NSArray *endmembers;
static NSCountedSet *instanceSet;

#define NR         3    /* Four independent composition variables      */
#define NS         1    /* Three ordering parameters s3 = X3           */
#define NA         4    /* Five endmember compositions                 */
#define NATOMS  28.0    /* Average number of atoms in the formula unit */

#pragma mark -
#pragma mark class methods

+(void)initialize {
	if (self == [NephelineSSBerman class]) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"Initialize(NephelineSSBerman) - entry ...");
		instanceSet = [[NSCountedSet alloc] initWithCapacity:1];
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

		BermanProperties *knepheline = [[BermanProperties alloc] initWithH:-2109563.55*4.0-(DH2)-3572.76
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

		VacancyNepheline *vcnepheline = [[VacancyNepheline alloc] init];
		[vcnepheline setPhaseFormula:@"Na3Al3Si5O16"];
		[vcnepheline setPhaseName:@"vc-nepheline"];
		[mutableEndmembers addObject:vcnepheline];
		if (debug) NSLog(@"... allocated vc-nepheline ...");

		CalciumNepheline *canepheline = [[CalciumNepheline alloc] init];
		[canepheline setPhaseFormula:@"CaNa2Al4Si4O16"];
		[canepheline setPhaseName:@"ca-nepheline"];
		[mutableEndmembers addObject:canepheline];
		if (debug) NSLog(@"... allocated ca-nepheline ...");

		endmembers = [NSArray arrayWithArray:mutableEndmembers];
	}
}

#pragma mark -
#pragma mark instance methods

@synthesize operationParent;

-(id)initWithCompositionConstraint:(NSString *)name {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(NephelineSSBerman) with %@ composition ... entry ...", name);

		[self setPhaseName:[NSString stringWithString:name]];
		computeMixingQuantities = NO;
		operationParent = @"";

		tOld = -9999.0;
		pOld = -9999.0;
		for (NSUInteger i=0; i<NR; i++) rOld[i] = -9999.0;
		for (NSUInteger i=0; i<NS; i++) sOld[i] = 2.0;
		xkls  = 0.0;
		xvcls = 0.0;
		xnals = 0.0;
		xkss  = 0.0;
		xcass = 0.0;
		xnass = 0.0;
		if (debug) NSLog(@"... exiting.");
	}
	return self;
}

-(id)init {
	return [self initWithCompositionConstraint:@"Feldspathoid"];
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

#define HEX                 -31744.63 /* joules */
#define SEX                    -20.92 /* joules/T */
#define VEX                     -1.05 /* joules/bar */
#define HX                  -13893.18 /* joules */
#define SX                      12.55 /* joules/T */
#define VX                        0.0 /* joules/bar */

#define H23                   73520.0 /* joules */
#define H24                   42560.0 /* joules */
#define S23                       0.0 /* joules/K */
#define S24                       0.0 /* joules/K */
/* nepheline structure */
#define WHNAKLS               6861.76 /* W(LS)Na-K joules */
#define WVNAKLS                  0.33 /* joules/bar */
#define WHNAKSS              51002.96 /* W(SS)Na-K joules */
#define WVNAKSS                  0.54 /* joules/bar */
#define WVN                   14644.0 /* W[]Si-NaAl joules */
#define WVK                    8368.0 /* W[]Si-KAl joules */
#define WCANA                     0.0 /* W[]Ca-Na2 joules */
#define WCAK                      0.0 /* W[]Ca-K2 joules */
#define WPLAG                     0.0 /* WCaAl-NaSi joules */
/* kalsilite structure */
#define WKALS                 29288.0 /* WNa4-K4 joules */
#define WKVKALS               30334.0 /* WK4-[]Na3 joules */
#define WNVKALS               14646.0 /* WNa4-[]Na3 joules */
#define WNCKALS                   0.0 /* WNa4-[]CaNa2 joules */
#define WKCKALS               14644.0 /* WK4-[]CaNa2 joules */
#define WVVKALS                   0.0 /* W[]Na3-[]CaNa2 joules */

#define DWKNALS                   0.0 /* joules */
#define DWKNASS                   0.0 /* joules */
#define DWVNASS                   0.0 /* joules */
#define DWVNALS                   0.0 /* joules */
#define DWVKLS                    0.0 /* joules */
#define AWLS                      0.0 /* joules */
#define AWSS                      0.0 /* joules */
#define AWVNA                     0.0 /* joules */
#define AWVK                      0.0 /* joules */

#define R  8.3143

/* correction to zero point entropy of []CaNa2Al4Si4O16 */
#define S4      -15.8765
#define G4      -t*(S4)

#define GEX     (HEX)-t*(SEX)+(p-1.0)*(VEX)
#define GX      (HX)-t*(SX)+(p-1.0)*(VX)
#define WNAKLS  (WHNAKLS)+(p-1.0)*(WVNAKLS)
#define WNAKSS  (WHNAKSS)+(p-1.0)*(WVNAKSS)
#define G23     (H23)-t*(S23)
#define G24     (H24)-t*(S24)

#define MAX_ITER 200    /* Maximum number of iterations allowed in order */

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conSpn defines the conversion from m[i], to r[j]
 */
/* Order: X2, X3, X4 */
#define FR2(i)     (i == 1) ? 1.0 - r[0] : - r[0]
#define FR3(i)     (i == 2) ? 1.0 - r[1] : - r[1]
#define FR4(i)     (i == 3) ? 1.0 - r[2] : - r[2]

/* Order: S1 */
#define GS1(i)     - s[0]

#define DFR2DR2(i) - 1.0
#define DFR3DR3(i) - 1.0
#define DFR4DR4(i) - 1.0

#define DGS1DS1(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */

#define S (S4)*r[2] - R*(xkls*log(xkls) + xvcls*log(xvcls) + \
xnals*log(xnals) - (1.0-3.0*xcass)*log(1.0-3.0*xcass) + \
3.0*xkss*log(xkss) + 3.0*xcass*log(xcass) + 3.0*xnass*log(xnass)) \
+ 0.25*(2.0*(SEX)+(SX))*s[0] + (SX)*r[0]*(1.0-r[0]) \
+ (3.0/16.0)*(SX)*s[0]*s[0] - 0.5*(SX)*r[0]*s[0] \
+ 0.5*(2.0*(S23)+(SEX)-(SX))*r[0]*r[1] + 0.125*((SX)-(SEX)-2.0*(S23))*r[1]*s[0] \
+ (6.0*(S23)+3.0*(SEX)-15.0*(SX))*r[0]*r[2]/18.0 \
+ (-9.0*(SEX)-3.0*(SX)+6.0*(S23)-36.0*(S24))*r[2]*s[0]/24.0
#define H 0.25*(2.0*(HEX)+(HX)+3.0*(WHNAKLS)-(WHNAKSS))*s[0] \
+ ((HX)+(WHNAKLS)+(WHNAKSS))*r[0]*(1.0-r[0]) \
+ (1.0/16.0)*(3.0*(HX)-9.0*(WHNAKLS)-(WHNAKSS))*s[0]*s[0] \
- 0.5*((HX)+3.0*(WHNAKLS)-(WHNAKSS))*r[0]*s[0] + (WVN)*r[1]*(1.0-r[1]) \
+ 0.5*(2.0*(H23)+(HEX)-(HX)-2.0*(WHNAKLS)-2.0*(WVN)+2.0*(WVK))*r[0]*r[1] \
+ 0.125*((HX)-(HEX)-2.0*(H23)-6.0*(WHNAKLS)-6.0*(WVN)+6.0*(WVK))*r[1]*s[0] \
+ (DWKNALS)*(r[0]*r[0]+0.5*r[0]*s[0]-3.0*s[0]*s[0]/16.0-r[0]*r[0]*r[0] \
-5.0*r[0]*r[0]*s[0]/4.0-3.0*r[0]*s[0]*s[0]/16.0 \
+9.0*s[0]*s[0]*s[0]/64.0-r[0]*r[1]/3.0-r[1]*s[0]/4.0 \
-2.0*r[0]*r[0]*r[1]/3.0+r[0]*r[1]*r[1]/3.0+r[1]*r[1]*s[0]/4.0 \
+3.0*r[1]*s[0]*s[0]/8.0) \
+ (DWKNASS)*(r[0]*r[0]+0.5*r[0]*s[0]-3.0*s[0]*s[0]/16.0-r[0]*r[0]*r[0] \
-r[0]*r[0]*s[0]/4.0+5.0*r[0]*s[0]*s[0]/16.0 \
-3.0*s[0]*s[0]*s[0]/64.0) \
+ (DWVNASS)*(4.0*r[0]*r[1]/3.0-4.0*r[0]*r[0]*r[1]/3.0-r[0]*r[1]*r[1]/3.0 \
-r[1]*r[1]*s[0]/4.0-r[1]*s[0]*s[0]/4.0) \
+ (DWVNALS)*(2.0*r[0]*r[1]/3.0-r[1]*s[0]/2.0-2.0*r[0]*r[0]*r[1]/3.0 \
-2.0*r[0]*r[1]*r[1]/3.0+0.5*r[1]*r[1]*s[0] \
+3.0*r[1]*s[0]*s[0]/8.0) \
+ (DWVKLS)*(r[0]*r[1]/3.0+r[1]*s[0]/4.0+2.0*r[0]*r[0]*r[1]/3.0 \
-r[0]*r[1]*r[1]/3.0-r[1]*r[1]*s[0]/4.0 \
-3.0*r[1]*s[0]*s[0]/8.0) \
+ (AWLS)*(r[0]+3.0*s[0]/4.0-3.0*r[0]*r[0]-9.0*r[0]*s[0]/2.0 \
-27.0*s[0]*s[0]/16.0+2.0*r[0]*r[0]*r[0]+9.0*r[0]*r[0]*s[0]/2.0 \
+27.0*r[0]*s[0]*s[0]/8.0+27.0*s[0]*s[0]*s[0]/32.0 + r[0]*r[1] \
+3.0*r[1]*s[0]/4.0-2.0*r[0]*r[1]*r[1]-3.0*r[1]*r[1]*s[0]/2.0) \
+ (AWSS)*(r[0]-s[0]/4.0-3.0*r[0]*r[0]+3.0*r[0]*s[0]/2.0-3.0*s[0]*s[0]/16.0 \
+2.0*r[0]*r[0]*r[0]-3.0*r[0]*r[0]*s[0]/2.0+3.0*r[0]*s[0]*s[0]/8.0 \
-s[0]*s[0]*s[0]/32.0) \
+ (AWVNA)*(r[1]-r[0]*r[1]-3.0*r[1]*s[0]/4.0-3.0*r[1]*r[1]+2.0*r[0]*r[1]*r[1] \
+ 3.0*r[1]*r[1]*s[0]/2.0+2.0*r[1]*r[1]*r[1]) \
+ (AWVK)*(r[0]*r[1]+3.0*r[1]*s[0]/4.0-2.0*r[0]*r[1]*r[1] \
-3.0*r[1]*r[1]*s[0]/2.0) \
+ (WCANA)*r[2]*(1.0-r[2]) \
+ (6.0*(H23)+3.0*(HEX)-15.0*(HX)-18.0*(WHNAKLS)-4.0*(WHNAKSS)+36.0*(WCAK) \
-36.0*(WCANA)+18.0*(WVN)-18.0*(WVK))*r[0]*r[2]/18.0 \
+ ((WPLAG)-(WCANA)-(WVN))*r[1]*r[2] \
+ (18.0*(WVN)-18.0*(WVK)+36.0*(WCAK)-36.0*(WCANA)-18.0*(WHNAKLS) \
+4.0*(WHNAKSS)-9.0*(HEX)-3.0*(HX)+6.0*(H23)-36.0*(H24))*r[2]*s[0]/24.0
#define V 0.25*(2.0*(VEX)+(VX)+3.0*(WVNAKLS)-(WVNAKSS))*s[0] \
+ ((VX)+(WVNAKLS)+(WVNAKSS))*r[0]*(1.0-r[0]) \
+ (1.0/16.0)*(3.0*(VX)-9.0*(WVNAKLS)-(WVNAKSS))*s[0]*s[0] \
- 0.5*((VX)+3.0*(WVNAKLS)-(WVNAKSS))*r[0]*s[0] \
+ 0.5*((VEX)-(VX)-2.0*(WVNAKLS))*r[0]*r[1] \
+ 0.125*((VX)-(VEX)-6.0*(WVNAKLS))*r[1]*s[0] \
+ (3.0*(VEX)-15.0*(VX)-18.0*(WVNAKLS)-4.0*(WVNAKSS))*r[0]*r[2]/18.0 \
+ (-18.0*(WVNAKLS)+4.0*(WVNAKSS)-9.0*(VEX)-3.0*(VX))*r[2]*s[0]/24.0
#define G (H) - t*(S) + (p-1.0)*(V)

/*----------------------------------------------------------------------------*/

#define DGDR0 R*t*(log(xkls) - log(xnals) + 3.0*log(xkss) - 3.0*log(xnass)) \
+ ((GX)+(WNAKLS)+(WNAKSS))*(1.0-2.0*r[0]) \
- 0.5*((GX)+3.0*(WNAKLS)-(WNAKSS))*s[0] \
+ 0.5*(2.0*(G23)+(GEX)-(GX)-2.0*(WNAKLS)-2.0*(WVN)+2.0*(WVK))*r[1] \
+ (DWKNALS)*(2.0*r[0]+0.5*s[0]-3.0*r[0]*r[0]-5.0*r[0]*s[0]/2.0 \
-3.0*s[0]*s[0]/16.0-r[1]/3.0-4.0*r[0]*r[1]/3.0 \
+r[1]*r[1]/3.0+3.0*r[1]*s[0]*s[0]/8.0) \
+ (DWKNASS)*(2.0*r[0]+0.5*s[0]-3.0*r[0]*r[0]-r[0]*s[0]/2.0 \
+5.0*s[0]*s[0]/16.0) \
+ (DWVNASS)*(4.0*r[1]/3.0-8.0*r[0]*r[1]/3.0-r[1]*r[1]/3.0) \
+ (DWVNALS)*(2.0*r[1]/3.0-4.0*r[0]*r[1]/3.0-2.0*r[1]*r[1]/3.0) \
+ (DWVKLS)*(r[1]/3.0+4.0*r[0]*r[1]/3.0-r[1]*r[1]/3.0) \
+ (AWLS)*(1.0-6.0*r[0]-9.0*s[0]/2.0+6.0*r[0]*r[0]+9.0*r[0]*s[0] \
+27.0*s[0]*s[0]/8.0+r[1]-2.0*r[1]*r[1]) \
+ (AWSS)*(1.0-6.0*r[0]+3.0*s[0]/2.0+6.0*r[0]*r[0] \
-6.0*r[0]*s[0]/2.0+3.0*s[0]*s[0]/8.0) \
+ (AWVNA)*(-r[1]+2.0*r[1]*r[1]) + (AWVK)*(r[1]-2.0*r[1]*r[1]) \
+ (6.0*(G23)+3.0*(GEX)-15.0*(GX)-18.0*(WNAKLS)-4.0*(WNAKSS)+36.0*(WCAK) \
-36.0*(WCANA)+18.0*(WVN)-18.0*(WVK))*r[2]/18.0
#define DGDR1 R*t*(log(xvcls) - log(xnals)) + (WVN)*(1.0-2.0*r[1]) \
+ 0.5*(2.0*(G23)+(GEX)-(GX)-2.0*(WNAKLS)-2.0*(WVN)+2.0*(WVK))*r[0] \
+ 0.125*((GX)-(GEX)-2.0*(G23)-6.0*(WNAKLS)-6.0*(WVN)+6.0*(WVK))*s[0] \
+ (DWKNALS)*(-r[0]/3.0-s[0]/4.0-2.0*r[0]*r[0]/3.0+2.0*r[0]*r[1]/3.0 \
+r[1]*s[0]/2.0+3.0*s[0]*s[0]/8.0) \
+ (DWVNASS)*(4.0*r[0]/3.0-4.0*r[0]*r[0]/3.0-2.0*r[0]*r[1]/3.0 \
-r[1]*s[0]/2.0-s[0]*s[0]/4.0) \
+ (DWVNALS)*(2.0*r[0]/3.0-s[0]/2.0-2.0*r[0]*r[0]/3.0 \
-4.0*r[0]*r[1]/3.0+r[1]*s[0]+3.0*s[0]*s[0]/8.0) \
+ (DWVKLS)*(r[0]/3.0+s[0]/4.0+2.0*r[0]*r[0]/3.0 \
-2.0*r[0]*r[1]/3.0-r[1]*s[0]/2.0-3.0*s[0]*s[0]/8.0) \
+ (AWLS)*(r[0]+3.0*s[0]/4.0-4.0*r[0]*r[1]-6.0*r[1]*s[0]/2.0) \
+ (AWVNA)*(1.0-r[0]-3.0*s[0]/4.0-6.0*r[1]+4.0*r[0]*r[1] \
+ 6.0*r[1]*s[0]/2.0+6.0*r[1]*r[1]) \
+ (AWVK)*(r[0]+3.0*s[0]/4.0-4.0*r[0]*r[1]-6.0*r[1]*s[0]/2.0) \
+ ((WPLAG)-(WCANA)-(WVN))*r[2]
#define DGDR2 (G4) + R*t*(log(xvcls) - log(xnals) + log(1.0-3.0*xcass) \
+ log(xcass) - log(xnass) + 1.0) \
+ (WCANA)*(1.0-2.0*r[2]) \
+ (6.0*(G23)+3.0*(GEX)-15.0*(GX)-18.0*(WNAKLS)-4.0*(WNAKSS)+36.0*(WCAK) \
-36.0*(WCANA)+18.0*(WVN)-18.0*(WVK))*r[0]/18.0 \
+ ((WPLAG)-(WCANA)-(WVN))*r[1] \
+ (18.0*(WVN)-18.0*(WVK)+36.0*(WCAK)-36.0*(WCANA)-18.0*(WNAKLS) \
+4.0*(WNAKSS)-9.0*(GEX)-3.0*(GX)+6.0*(G23)-36.0*(G24))*s[0]/24.0
#define DGDS0 R*t*(0.75*log(xkls) - 0.75*log(xnals) \
- 3.0*0.25*log(xkss) + 3.0*0.25*log(xnass)) \
+ 0.25*(2.0*(GEX)+(GX)+3.0*(WNAKLS)-(WNAKSS)) \
+ (1.0/8.0)*(3.0*(GX)-9.0*(WNAKLS)-(WNAKSS))*s[0] \
- 0.5*((GX)+3.0*(WNAKLS)-(WNAKSS))*r[0] \
+ 0.125*((GX)-(GEX)-2.0*(G23)-6.0*(WNAKLS)-6.0*(WVN)+6.0*(WVK))*r[1] \
+ (DWKNALS)*(0.5*r[0]-3.0*s[0]/8.0-5.0*r[0]*r[0]/4.0 \
-3.0*r[0]*s[0]/8.0+27.0*s[0]*s[0]/64.0-r[1]/4.0 \
+r[1]*r[1]/4.0+3.0*r[1]*s[0]/4.0) \
+ (DWKNASS)*(0.5*r[0]-3.0*s[0]/8.0-r[0]*r[0]/4.0+5.0*r[0]*s[0]/8.0 \
-9.0*s[0]*s[0]/64.0) \
+ (DWVNASS)*(-r[1]*r[1]/4.0-r[1]*s[0]/2.0) \
+ (DWVNALS)*(-r[1]/2.0+0.5*r[1]*r[1]+3.0*r[1]*s[0]/4.0) \
+ (DWVKLS)*(r[1]/4.0-r[1]*r[1]/4.0-3.0*r[1]*s[0]/4.0) \
+ (AWLS)*(3.0/4.0-9.0*r[0]/2.0-27.0*s[0]/8.0+9.0*r[0]*r[0]/2.0 \
+27.0*r[0]*s[0]/4.0+81.0*s[0]*s[0]/32.0+3.0*r[1]/4.0 \
-3.0*r[1]*r[1]/2.0) \
+ (AWSS)*(-1.0/4.0+3.0*r[0]/2.0-3.0*s[0]/8.0-3.0*r[0]*r[0]/2.0 \
+3.0*r[0]*s[0]/4.0-3.0*s[0]*s[0]/32.0) \
+ (AWVNA)*(-3.0*r[1]/4.0+3.0*r[1]*r[1]/2.0) \
+ (AWVK)*(3.0*r[1]/4.0-3.0*r[1]*r[1]/2.0) \
+ (18.0*(WVN)-18.0*(WVK)+36.0*(WCAK)-36.0*(WCANA)-18.0*(WNAKLS) \
+4.0*(WNAKSS)-9.0*(GEX)-3.0*(GX)+6.0*(G23)-36.0*(G24))*r[2]/24.0
#define DGDT  - (S)
#define DGDP  (V)

/*----------------------------------------------------------------------------*/

#define D2GDR0R0 R*t*(1.0/xkls + 1.0/xnals + 3.0/xkss + 3.0/xnass) \
- 2.0*((GX)+(WNAKLS)+(WNAKSS)) \
+ (DWKNALS)*(2.0-6.0*r[0]-5.0*s[0]/2.0-4.0*r[1]/3.0) \
+ (DWKNASS)*(2.0-6.0*r[0]-s[0]/2.0) + (DWVNASS)*(-8.0*r[1]/3.0) \
- (DWVNALS)*(4.0*r[1]/3.0) + (DWVKLS)*(4.0*r[1]/3.0) \
+ (AWLS)*(-6.0+12.0*r[0]+9.0*s[0]) \
+ (AWSS)*(-6.0+12.0*r[0]-6.0*s[0]/2.0)
#define D2GDR0R1 R*t/xnals \
+ 0.5*(2.0*(G23)+(GEX)-(GX)-2.0*(WNAKLS)-2.0*(WVN)+2.0*(WVK)) \
+ (DWKNALS)*(-1.0/3.0-4.0*r[0]/3.0+2.0*r[1]/3.0+3.0*s[0]*s[0]/8.0) \
+ (DWVNASS)*(4.0/3.0-8.0*r[0]/3.0-2.0*r[1]/3.0) \
+ (DWVNALS)*(2.0/3.0-4.0*r[0]/3.0-4.0*r[1]/3.0) \
+ (DWVKLS)*(1.0/3.0+4.0*r[0]/3.0-2.0*r[1]/3.0) \
+ (AWLS)*(1.0-4.0*r[1]) + (AWVNA)*(-1.0+4.0*r[1]) \
+ (AWVK)*(1.0-4.0*r[1])
#define D2GDR0R2 R*t*(1.0/xnals + 1.0/xnass) \
+ (6.0*(G23)+3.0*(GEX)-15.0*(GX)-18.0*(WNAKLS)-4.0*(WNAKSS)+36.0*(WCAK) \
-36.0*(WCANA)+18.0*(WVN)-18.0*(WVK))/18.0
#define D2GDR0S0 R*t*(0.75/xkls + 0.75/xnals - 3.0*0.25/xkss \
- 3.0*0.25/xnass) \
- 0.5*((GX)+3.0*(WNAKLS)-(WNAKSS)) \
+ (DWKNALS)*(0.5-5.0*r[0]/2.0-6.0*s[0]/16.0+6.0*r[1]*s[0]/8.0) \
+ (DWKNASS)*(0.5-r[0]/2.0+10.0*s[0]/16.0) \
+ (AWLS)*(-9.0/2.0+9.0*r[0]+27.0*s[0]/4.0+r[1]) \
+ (AWSS)*(3.0/2.0-6.0*r[0]/2.0+3.0*s[0]/4.0)
#define D2GDR0DT R*(log(xkls) - log(xnals) + 3.0*log(xkss) \
- 3.0*log(xnass)) - (SX)*(1.0-2.0*r[0]) \
+ 0.5*(SX)*s[0] - 0.5*(2.0*(S23)+(SEX)-(SX))*r[1] \
+ (-6.0*(S23)-3.0*(SEX)+15.0*(SX))*r[2]/18.0
#define D2GDR0DP ((VX)+(WVNAKLS)+(WVNAKSS))*(1.0-2.0*r[0]) \
- 0.5*((VX)+3.0*(WVNAKLS)-(WVNAKSS))*s[0] \
+ 0.5*((VEX)-(VX)-2.0*(WVNAKLS))*r[1] \
+ (3.0*(VEX)-15.0*(VX)-18.0*(WVNAKLS)-4.0*(WVNAKSS))*r[2]/18.0

#define D2GDR1R1 R*t*(1.0/xvcls + 1.0/xnals) - 2.0*(WVN) \
+ (DWKNALS)*(2.0*r[0]/3.0+s[0]/2.0) \
+ (DWVNASS)*(-2.0*r[0]/3.0-s[0]/2.0) \
+ (DWVNALS)*(-4.0*r[0]/3.0+s[0]) \
+ (DWVKLS)*(-2.0*r[0]/3.0-s[0]/2.0) \
+ (AWLS)*(-4.0*r[0]-6.0*s[0]/2.0) \
+ (AWVNA)*(-6.0+4.0*r[0]+6.0*s[0]/2.0+12.0*r[1]) \
+ (AWVK)*(-4.0*r[0]-6.0*s[0]/2.0)
#define D2GDR1R2 R*t*(1.0/xvcls + 1.0/xnals) + ((WPLAG)-(WCANA)-(WVN))
#define D2GDR1S0 R*t*0.75/xnals \
+ 0.125*((GX)-(GEX)-2.0*(G23)-6.0*(WNAKLS)-6.0*(WVN)+6.0*(WVK)) \
+ (DWKNALS)*(-1.0/4.0+r[1]/2.0+3.0*s[0]/4.0) \
+ (DWVNASS)*(-r[1]/2.0-s[0]/2.0) \
+ (DWVNALS)*(-1.0/2.0+r[1]+3.0*s[0]/4.0) \
+ (DWVKLS)*(1.0/4.0-r[1]/2.0-3.0*s[0]/4.0) \
+ (AWLS)*(3.0/4.0-6.0*r[1]/2.0) \
+ (AWVNA)*(-3.0/4.0+6.0*r[1]/2.0) \
+ (AWVK)*(3.0/4.0-6.0*r[1]/2.0)
#define D2GDR1DT R*(log(xvcls) - log(xnals)) \
- 0.5*(2.0*(S23)+(SEX)-(SX))*r[0] - 0.125*((SX)-(SEX)-2.0*(S23))*s[0]
#define D2GDR1DP 0.5*((VEX)-(VX)-2.0*(WVNAKLS))*r[0] \
+ 0.125*((VX)-(VEX)-6.0*(WVNAKLS))*s[0]

#define D2GDR2R2 R*t*(1.0/xvcls + 1.0/xnals - 1.0/(1.0-3.0*xcass) \
+ 1.0/(3.0*xcass) + 1.0/(3.0*xnass)) - 2.0*(WCANA)
#define D2GDR2S0 R*t*(0.75/xnals - 0.25/xnass) \
+ (18.0*(WVN)-18.0*(WVK)+36.0*(WCAK)-36.0*(WCANA)-18.0*(WNAKLS) \
+4.0*(WNAKSS)-9.0*(GEX)-3.0*(GX)+6.0*(G23)-36.0*(G24))/24.0
#define D2GDR2DT - (S4) + R*(log(xvcls) - log(xnals) + log(1.0-3.0*xcass) \
+ log(xcass) - log(xnass) + 1.0) \
+ (-6.0*(S23)-3.0*(SEX)+15.0*(SX))*r[0]/18.0 \
+ (9.0*(SEX)+3.0*(SX)-6.0*(S23)+36.0*(S24))*s[0]/24.0
#define D2GDR2DP (3.0*(VEX)-15.0*(VX)-18.0*(WVNAKLS)-4.0*(WVNAKSS))*r[0]/18.0 \
+ (-18.0*(WVNAKLS)+4.0*(WVNAKSS)-9.0*(VEX)-3.0*(VX))*s[0]/24.0

#define D2GDS0S0 R*t*(0.75*0.75/xkls + 0.75*0.75/xnals \
+ 3.0*0.25*0.25/xkss + 3.0*0.25*0.25/xnass) \
+ (1.0/8.0)*(3.0*(GX)-9.0*(WNAKLS)-(WNAKSS)) \
+ (DWKNALS)*(-3.0/8.0-3.0*r[0]/8.0+27.0*s[0]/32.0 \
+3.0*r[1]/4.0) \
+ (DWKNASS)*(-3.0/8.0+5.0*r[0]/8.0-9.0*s[0]/32.0) \
+ (DWVNASS)*(-r[1]/2.0) + (DWVNALS)*(3.0*r[1]/4.0) \
+ (DWVKLS)*(-3.0*r[1]/4.0) \
+ (AWLS)*(-27.0/8.0+27.0*r[0]/4.0+81.0*s[0]/16.0) \
+ (AWSS)*(-3.0/8.0+3.0*r[0]/4.0-3.0*s[0]/16.0)
#define D2GDS0DT R*(0.75*log(xkls) - 0.75*log(xnals) \
- 3.0*0.25*log(xkss) + 3.0*0.25*log(xnass)) \
- 0.25*(2.0*(SEX)+(SX)) - (1.0/8.0)*(3.0*(SX))*s[0] \
+ 0.5*((SX))*r[0] - 0.125*((SX)-(SEX)-2.0*(S23))*r[1] \
+ (9.0*(SEX)+3.0*(SX)-6.0*(S23)+36.0*(S24))*r[2]/24.0
#define D2GDS0DP 0.25*(2.0*(VEX)+(VX)+3.0*(WVNAKLS)-(WVNAKSS)) \
+ (1.0/8.0)*(3.0*(VX)-9.0*(WVNAKLS)-(WVNAKSS))*s[0] \
- 0.5*((VX)+3.0*(WVNAKLS)-(WVNAKSS))*r[0] \
+ 0.125*((VX)-(VEX)-6.0*(WVNAKLS))*r[1] \
+ (-18.0*(WVNAKLS)+4.0*(WVNAKSS)-9.0*(VEX)-3.0*(VX))*r[2]/24.0

#define D2GDT2   0.0
#define D2GDTDP  0.0
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/

#define D3GDR0R0R0 - R*t*(1.0/(xkls*xkls) - 1.0/(xnals*xnals) \
+ 3.0/(xkss*xkss) - 3.0/(xnass*xnass)) - 6.0*(DWKNALS) \
- 6.0*(DWKNASS) + (AWLS)*(12.0) + (AWSS)*(12.0)
#define D3GDR0R0R1 R*t/(xnals*xnals) \
+ (DWKNALS)*(-4.0/3.0) + (DWVNASS)*(-8.0/3.0) \
- (DWVNALS)*(4.0/3.0) + (DWVKLS)*(4.0/3.0)

#define D3GDR0R0R2 R*t*(1.0/(xnals*xnals) + 1.0/(xnass*xnass))
#define D3GDR0R0S0 - R*t*(0.75/(xkls*xkls) - 0.75/(xnals*xnals) \
- 3.0*0.25/(xkss*xkss) + 3.0*0.25/(xnass*xnass)) \
+ (DWKNALS)*(-5.0/2.0) + (DWKNASS)*(-1.0/2.0) \
+ (AWLS)*(9.0) + (AWSS)*(-6.0/2.0)
#define D3GDR0R0DT R*(1.0/xkls + 1.0/xnals + 3.0/xkss + 3.0/xnass) \
+ 2.0*(SX)
#define D3GDR0R0DP - 2.0*((VX)+(WVNAKLS)+(WVNAKSS))

#define D3GDR0R1R1 R*t/(xnals*xnals) \
+ (DWKNALS)*(2.0/3.0) + (DWVNASS)*(-2.0/3.0) \
+ (DWVNALS)*(-4.0/3.0) + (DWVKLS)*(-2.0/3.0) \
+ (AWLS)*(-4.0) + (AWVNA)*(4.0) + (AWVK)*(-4.0)
#define D3GDR0R1R2 R*t/(xnals*xnals)
#define D3GDR0R1S0 R*t*0.75/(xnals*xnals) \
+ (DWKNALS)*(3.0*s[0]/4.0)
#define D3GDR0R1DT R/xnals - 0.5*(2.0*(S23)+(SEX)-(SX))
#define D3GDR0R1DP 0.5*((VEX)-(VX)-2.0*(WVNAKLS))

#define D3GDR0R2R2 R*t*(1.0/(xnals*xnals) + 1.0/(3.0*xnass*xnass))
#define D3GDR0R2S0 R*t*(0.75/(xnals*xnals) - 0.25/(xnass*xnass))
#define D3GDR0R2DT R*(1.0/xnals + 1.0/xnass) \
+ (-6.0*(S23)-3.0*(SEX)+15.0*(SX))/18.0
#define D3GDR0R2DP (3.0*(VEX)-15.0*(VX)-18.0*(WVNAKLS)-4.0*(WVNAKSS))/18.0

#define D3GDR0S0S0 - R*t*(0.75*0.75/(xkls*xkls) - 0.75*0.75/(xnals*xnals) \
+ 3.0*0.25*0.25/(xkss*xkss) - 3.0*0.25*0.25/(xnass*xnass)) \
+ (DWKNALS)*(-6.0/16.0+6.0*r[1]/8.0) \
+ (DWKNASS)*(10.0/16.0) + (AWLS)*(27.0/4.0) + (AWSS)*(3.0/4.0)
#define D3GDR0S0DT R*(0.75/xkls + 0.75/xnals - 3.0*0.25/xkss \
- 3.0*0.25/xnass) + 0.5*(SX)
#define D3GDR0S0DP - 0.5*((VX)+3.0*(WVNAKLS)-(WVNAKSS))

#define D3GDR1R1R1 - R*t*(1.0/(xvcls*xvcls) - 1.0/(xnals*xnals)) \
+ (AWVNA)*(12.0)
#define D3GDR1R1R2 R*t*(-1.0/(xvcls*xvcls) + 1.0/(xnals*xnals))
#define D3GDR1R1S0 R*t*0.75/(xnals*xnals) \
+ (DWKNALS)*(1.0/2.0) + (DWVNASS)*(-1.0/2.0) \
+ (DWVNALS) + (DWVKLS)*(-1.0/2.0) + (AWLS)*(-6.0/2.0) \
+ (AWVNA)*(6.0/2.0)+ (AWVK)*(-6.0/2.0)
#define D3GDR1R1DT R*(1.0/xvcls + 1.0/xnals)
#define D3GDR1R1DP 0.0

#define D3GDR1R2R2 R*t*(-1.0/(xvcls*xvcls) + 1.0/(xnals*xnals))
#define D3GDR1R2S0 R*t*0.75/(xnals*xnals)
#define D3GDR1R2DT  R*(1.0/xvcls + 1.0/xnals)
#define D3GDR1R2DP 0.0

#define D3GDR1S0S0 R*t*0.75*0.75/(xnals*xnals) \
+ (DWKNALS)*(3.0/4.0) + (DWVNASS)*(-1.0/2.0) \
+ (DWVNALS)*(3.0/4.0) + (DWVKLS)*(-3.0/4.0)
#define D3GDR1S0DT R*0.75/xnals - 0.125*(-2.0*(S23)+(SX)-(SEX))
#define D3GDR1S0DP 0.125*((VX)-(VEX)-6.0*(WVNAKLS))

#define D3GDR2R2R2 R*t*(-1.0/(xvcls*xvcls) + 1.0/(xnals*xnals) \
- 1.0/SQUARE(1.0-3.0*xcass) \
- 1.0/SQUARE(3.0*xcass) + 1.0/SQUARE(3.0*xnass))
#define D3GDR2R2S0 R*t*(0.75/(xnals*xnals) - 0.25/(3.0*xnass*xnass))
#define D3GDR2R2DT R*(1.0/xvcls + 1.0/xnals - 1.0/(1.0-3.0*xcass) \
+ 1.0/(3.0*xcass) + 1.0/(3.0*xnass))
#define D3GDR2R2DP 0.0

#define D3GDR2S0S0 R*t*(0.75*0.75/(xnals*xnals) + 0.25*0.25/(xnass*xnass))
#define D3GDR2S0DT R*(0.75/xnals - 0.25/xnass) \
+ (9.0*(SEX)+3.0*(SX)-6.0*(S23)+36.0*(S24))/24.0
#define D3GDR2S0DP (-18.0*(WVNAKLS)+4.0*(WVNAKSS)-9.0*(VEX)-3.0*(VX))/24.0

#define D3GDS0S0S0 - R*t*(0.75*0.75*0.75/(xkls*xkls) \
- 0.75*0.75*0.75/(xnals*xnals) - 3.0*0.25*0.25*0.25/(xkss*xkss) \
+ 3.0*0.25*0.25*0.25/(xnass*xnass)) + (DWKNALS)*(27.0/32.0) \
+ (DWKNASS)*(-9.0/32.0) + (AWLS)*(81.0/16.0) \
+ (AWSS)*(-3.0/16.0)
#define D3GDS0S0DT R*(0.75*0.75/xkls + 0.75*0.75/xnals \
+ 3.0*0.25*0.25/xkss + 3.0*0.25*0.25/xnass) \
- (1.0/8.0)*(3.0*(SX))
#define D3GDS0S0DP  (1.0/8.0)*(3.0*(VX)-9.0*(WVNAKLS)-(WVNAKSS))
#define D3GDS0DT2   0.0
#define D3GDS0DTDP  0.0
#define D3GDS0DP2   0.0

#define D3GDT3      0.0
#define D3GDT2DP    0.0
#define D3GDTDP2    0.0
#define D3GDP3      0.0

#define D3GDR0DT2   0.0
#define D3GDR0DTDP  0.0
#define D3GDR0DP2   0.0
#define D3GDR1DT2   0.0
#define D3GDR1DTDP  0.0
#define D3GDR1DP2   0.0
#define D3GDR2DT2   0.0
#define D3GDR2DTDP  0.0
#define D3GDR2DP2   0.0

/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1;     d2gdr2[0][2] = D2GDR0R2; \
d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;     d2gdr2[1][2] = D2GDR1R2; \
d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2]; d2gdr2[2][2] = D2GDR2R2;

#define fillD2GDRDS d2gdrds[0][0] = D2GDR0S0; d2gdrds[1][0] = D2GDR1S0; d2gdrds[2][0] = D2GDR2S0;
#define fillD2GDRDT d2gdrdt[0] = D2GDR0DT;    d2gdrdt[1] = D2GDR1DT;    d2gdrdt[2] = D2GDR2DT;
#define fillD2GDRDP d2gdrdp[0] = D2GDR0DP;    d2gdrdp[1] = D2GDR1DP;    d2gdrdp[2] = D2GDR2DP;
#define fillD2GDS2  d2gds2[0][0] = D2GDS0S0;
#define fillD2GDSDT d2gdsdt[0] = D2GDS0DT;
#define fillD2GDSDP d2gdsdp[0] = D2GDS0DP;

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

#define fillD3GDR2DS \
d3gdr2ds[0][0][0] = D3GDR0R0S0;        d3gdr2ds[0][1][0] = D3GDR0R1S0;        d3gdr2ds[0][2][0] = D3GDR0R2S0; \
d3gdr2ds[1][0][0] = d3gdr2ds[0][1][0]; d3gdr2ds[1][1][0] = D3GDR1R1S0;        d3gdr2ds[1][2][0] = D3GDR1R2S0; \
d3gdr2ds[2][0][0] = d3gdr2ds[0][2][0]; d3gdr2ds[2][1][0] = d3gdr2ds[1][2][0]; d3gdr2ds[2][2][0] = D3GDR2R2S0;

#define fillD3GDR2DT \
d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT;     d3gdr2dt[0][2] = D3GDR0R2DT; \
d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT;     d3gdr2dt[1][2] = D3GDR1R2DT; \
d3gdr2dt[2][0] = d3gdr2dt[0][2]; d3gdr2dt[2][1] = d3gdr2dt[1][2]; d3gdr2dt[2][2] = D3GDR2R2DT;

#define fillD3GDR2DP \
d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP;     d3gdr2dp[0][2] = D3GDR0R2DP; \
d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP;     d3gdr2dp[1][2] = D3GDR1R2DP; \
d3gdr2dp[2][0] = d3gdr2dp[0][2]; d3gdr2dp[2][1] = d3gdr2dp[1][2]; d3gdr2dp[2][2] = D3GDR2R2DP;

#define fillD3GDRDS2 \
d3gdrds2[0][0][0] = D3GDR0S0S0;  d3gdrds2[1][0][0] = D3GDR1S0S0;  d3gdrds2[2][0][0] = D3GDR2S0S0;

#define fillD3GDRDSDT \
d3gdrdsdt[0][0] = D3GDR0S0DT;    d3gdrdsdt[1][0] = D3GDR1S0DT;    d3gdrdsdt[2][0] = D3GDR2S0DT;

#define fillD3GDRDSDP \
d3gdrdsdp[0][0] = D3GDR0S0DP;    d3gdrdsdp[1][0] = D3GDR1S0DP;    d3gdrdsdp[2][0] = D3GDR2S0DP;

#define fillD3GDS3    d3gds3[0][0][0] = D3GDS0S0S0;
#define fillD3GDS2DT  d3gds2dt[0][0] = D3GDS0S0DT;
#define fillD3GDS2DP  d3gds2dp[0][0] = D3GDS0S0DP;
#define fillD3GDSDT2  d3gdsdt2[0] = D3GDS0DT2;
#define fillD3GDSDTDP d3gdsdtdp[0] = D3GDS0DTDP;
#define fillD3GDSDP2  d3gdsdp2[0] = D3GDS0DP2;
#define fillD3GDRDT2  d3gdrdt2[0] = D3GDR0DT2;   d3gdrdt2[1] = D3GDR1DT2;   d3gdrdt2[2] = D3GDR2DT2;
#define fillD3GDRDTDP d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP; d3gdrdtdp[2] = D3GDR2DTDP;
#define fillD3GDRDP2  d3gdrdp2[0] = D3GDR0DP2;   d3gdrdp2[1] = D3GDR1DP2;   d3gdrdp2[2] = D3GDR2DP2;

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
	if ( (t != tOld)       || (p != pOld) ||
		(r[0] != rOld[0]) || (r[1] != rOld[1]) || (r[2] != rOld[2]) ) {
		double dgds[NS], sNew[NS];
		int skipCheck = FALSE;
		double xk  = r[0];
		double xvc = (r[1]+r[2])/4.0;
		double xna = 1.0 - r[0] - r[1]/4.0 - r[2]/2.0;

		for (i=0; i<NS; i++) { sOld[i] = 2.0; sNew[i] = 0.0; dgds[i] = 0.0; }

		/* Initial guess assumes random distribution */
		sNew[0] = -4.0*r[0]*(r[1]+2.0*r[2]/3.0)/(4.0-r[1]-2.0*r[2]);

		xkls  = r[0] + 3.0*sNew[0]/4.0;
		xvcls = r[1] + r[2];
		xnals = 1.0 - r[0] - r[1] - r[2] - 3.0*sNew[0]/4.0;
		xkss  = r[0] - sNew[0]/4.0;
		xcass = r[2]/3.0;
		xnass = 1.0 - r[0] - r[2]/3.0 + sNew[0]/4.0;

		if (xkls  <= DBL_EPSILON) xkls  = DBL_EPSILON;
		if (xvcls <= DBL_EPSILON) xvcls = DBL_EPSILON;
		if (xnals <= DBL_EPSILON) xnals = DBL_EPSILON;
		if (xkss  <= DBL_EPSILON) xkss  = DBL_EPSILON;
		if (xcass <= DBL_EPSILON) xcass = DBL_EPSILON;
		if (xnass <= DBL_EPSILON) xnass = DBL_EPSILON;

		if (xkls  >= 1.0     - DBL_EPSILON) xkls  = 1.0     - DBL_EPSILON;
		if (xvcls >= 1.0     - DBL_EPSILON) xvcls = 1.0     - DBL_EPSILON;
		if (xnals >= 1.0     - DBL_EPSILON) xnals = 1.0     - DBL_EPSILON;
		if (xkss  >= 1.0     - DBL_EPSILON) xkss  = 1.0     - DBL_EPSILON;
		if (xcass >= 1.0/3.0 - DBL_EPSILON) xcass = 1.0/3.0 - DBL_EPSILON;
		if (xnass >= 1.0     - DBL_EPSILON) xnass = 1.0     - DBL_EPSILON;

		sNew[0] = xkls - xkss;

		/* Check for irrelevant compositions */
		//if (xk  < sqrt(DBL_EPSILON) || xna < sqrt(DBL_EPSILON) ||
		//	xvc > ((double) (1.0/4.0) - sqrt(DBL_EPSILON)) ) {
		//	for (i=0; i<NS; i++) sOld[i] = sNew[i];
		//	skipCheck = TRUE;
		//}

		while ( (ABS(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) && (iter < MAX_ITER)) {
			double s[NS], sCorr[NS], lambda = 1.0;

			for (i=0; i<NS; i++) s[i] = sNew[i];

			dgds[0] = DGDS0;
			invd2gds2[0][0] = D2GDS0S0;

			for (i=0; i<NS; i++) sOld[i] = s[i];

            invd2gds2[0][0] = (invd2gds2[0][0] != 0.0) ? 1.0/invd2gds2[0][0] : 0.0;

			for (i=0; i<NS; i++) {
				for(j=0, sCorr[i]=0.0; j<NS; j++) sCorr[i] += - invd2gds2[i][j]*dgds[j];
				s[i] = sOld[i] + lambda*sCorr[i];
			}

			xkls  = r[0] + 3.0*s[0]/4.0;
			xvcls = r[1] + r[2];
			xnals = 1.0 - r[0] - r[1] - r[2] - 3.0*s[0]/4.0;
			xkss  = r[0] - s[0]/4.0;
			xcass = r[2]/3.0;
			xnass = 1.0 - r[0] - r[2]/3.0 + s[0]/4.0;

			while ((   xkls  < 0.0 || xkls  > 1.0     || xvcls < 0.0 || xvcls > 1.0
					|| xnals < 0.0 || xnals > 1.0     || xkss  < 0.0 || xkss  > 1.0
					|| xcass < 0.0 || xcass > 1.0/3.0 || xnass < 0.0 || xnass > 1.0)
				   && lambda > DBL_EPSILON ) {
				lambda /= 2.0;
				for (i=0; i<NS; i++) s[i] = sOld[i] + lambda*sCorr[i];
				xkls  = r[0] + 3.0*s[0]/4.0;
				xvcls = r[1] + r[2];
				xnals = 1.0 - r[0] - r[1] - r[2] - 3.0*s[0]/4.0;
				xkss  = r[0] - s[0]/4.0;
				xcass = r[2]/3.0;
				xnass = 1.0 - r[0] - r[2]/3.0 + s[0]/4.0;
			}

			if (xkls  <= DBL_EPSILON) xkls  = DBL_EPSILON;
			if (xvcls <= DBL_EPSILON) xvcls = DBL_EPSILON;
			if (xnals <= DBL_EPSILON) xnals = DBL_EPSILON;
			if (xkss  <= DBL_EPSILON) xkss  = DBL_EPSILON;
			if (xcass <= DBL_EPSILON) xcass = DBL_EPSILON;
			if (xnass <= DBL_EPSILON) xnass = DBL_EPSILON;

			if (xkls  >= 1.0     - DBL_EPSILON) xkls  = 1.0     - DBL_EPSILON;
			if (xvcls >= 1.0     - DBL_EPSILON) xvcls = 1.0     - DBL_EPSILON;
			if (xnals >= 1.0     - DBL_EPSILON) xnals = 1.0     - DBL_EPSILON;
			if (xkss  >= 1.0     - DBL_EPSILON) xkss  = 1.0     - DBL_EPSILON;
			if (xcass >= 1.0/3.0 - DBL_EPSILON) xcass = 1.0/3.0 - DBL_EPSILON;
			if (xnass >= 1.0     - DBL_EPSILON) xnass = 1.0     - DBL_EPSILON;

			s[0] = xkls - xkss;

			for (i=0; i<NS; i++) sNew[i] = s[i];
			iter++;
		}
		tOld = t;
		pOld = p;
		for (i=0; i<NR; i++) rOld[i] = r[i];

        BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
		if (debug) {
			for (i=0; i<NS; i++) {
				if (!skipCheck && ABS(dgds[i]) > sqrt(DBL_EPSILON)
					&& (ABS(          xk *dgds[i]) > sqrt(DBL_EPSILON) &&
						ABS(          xna*dgds[i]) > sqrt(DBL_EPSILON) &&
						ABS((1.0/4.0-xvc)*dgds[i]) > sqrt(DBL_EPSILON) )
					&& ABS(sOld[i]) > DBL_EPSILON) {
					NSLog(@"ERROR in NEPHELINE.C (function ORDER). Failed to converge!\n");
					if (iter >= MAX_ITER)
						NSLog(@"  Iteration limit (%4d) exceeded.\n", iter);
					NSLog(@"  X2    = %13.6g, X3    = %13.6g\n", r[0], r[1]);
					NSLog(@"  X4    = %13.6g, s1    = %13.6g\n", r[2], sOld[0]);
					NSLog(@"  dgds1 = %13.6g\n", dgds[0]);
					NSLog(@"  X K  ls: %13.6g  X K  ss: %13.6g\n", xkls,  xkss);
					NSLog(@"  X Na ls: %13.6g  X Na ss: %13.6g\n", xnals, xnass);
					NSLog(@"  X Vc ls: %13.6g  X Ca ss: %13.6g\n", xvcls, xcass);
					break;
				}
			}
		}

	}

	if (mask & FIRST  ) {   /* return s        */
		for (i=0; i<NS; i++) s[i] = sOld[i];
	}

	if (mask & SECOND) {   /* compute ds/dr:  */
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

	if (mask & THIRD) {   /* compute ds/dt:  */
		double *s = sOld;
		double d2gdsdt[NS];

		fillD2GDSDT

		for (i=0; i<NS; i++) {
			dt[i] = 0.0;
			for (j=0; j<NS; j++) dt[i] += - invd2gds2[i][j]*d2gdsdt[j];
		}
	}

	if (mask & FOURTH) {   /* compute ds/dp:  */
		double *s = sOld;
		double d2gdsdp[NS];

		fillD2GDSDP

		for (i=0; i<NS; i++) {
			dp[i] = 0.0;
			for (j=0; j<NS; j++) dp[i] += - invd2gds2[i][j]*d2gdsdp[j];
		}
	}

	if (mask & FIFTH) {   /* compute d2s/dr2 */
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

	if (mask & SIXTH) {   /* compute d2s/drt */
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

	if (mask & EIGHTH) {   /* compute d2s/dt2 */
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

	if (mask & NINTH) {   /* compute d2s/dtp */
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

	if (mask & TENTH) {   /* compute d2s/dp2 */
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

		result = result && (4.0*r[0] >= 0.0) && (4.0*r[0] <= 4.0+sqrt(DBL_EPSILON));                   /* tot K  */
		result = result && (4.0-4.0*r[0]-r[1]-2.0*r[2] >= 0.0)                                         /* tot Na */
		&& (4.0-4.0*r[0]-r[1]-2.0*r[2] <= 4.0+sqrt(DBL_EPSILON));
		result = result && (r[2] >= 0.0) && (r[2] <= 1.0+sqrt(DBL_EPSILON));                           /* tot Ca */
		result = result && (r[1]+r[2] >= 0.0) && (r[1]+r[2] <= 1.0+sqrt(DBL_EPSILON));                 /* tot Vc */
		result = result && (4.0-r[1] >= 3.0-sqrt(DBL_EPSILON)) && (4.0-r[1] <= 5.0+sqrt(DBL_EPSILON)); /* tot Al */
		result = result && (4.0+r[1] >= 3.0-sqrt(DBL_EPSILON)) && (4.0+r[1] <= 5.0+sqrt(DBL_EPSILON)); /* tot Si */
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
	double s[NS], g, dgdr[NR];
	double fr[NA][NR];
	int i, j;

	for(i=0; i<NA; i++) {
		fr[i][0] = FR2(i); /* X2 */
		fr[i][1] = FR3(i); /* X3 */
		fr[i][2] = FR4(i); /* X4 */
	}

	[self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

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
		double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR],
		dfrdr[NA][NR], gs[NA][NS], dgsds[NA][NS], sum;
		int k, l;

		fillD2GDR2
		fillD2GDRDS
		fillD2GDS2

		for(i=0; i<NA; i++) {
			gs[i][0] = GS1(i); /* s1 */
			dfrdr[i][0] = DFR2DR2(i); /* X2 */
			dfrdr[i][1] = DFR3DR3(i); /* X3 */
			dfrdr[i][2] = DFR4DR4(i); /* X4 */
			dgsds[i][0] = DGS1DS1(i); /* s1 */
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
		dx[2] = DGDR2;
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
		double r[NR] = { 0.75, 1.0, 0.0}, s[NS] = { -1.0 };
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
		double r[NR] = { 0.5, 0.0, 1.0}, s[NS] = { -2.0/3.0 };
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

	double projectedNaOverNaPlusK = ((xReduced[0]+xReduced[1]) != 0.0) ? xReduced[0]/(xReduced[0]+xReduced[1]) : 1.0;
	if      ((projectedNaOverNaPlusK < 0.35) && [[self phaseName] isEqualToString:@"Nepheline ss"]) affinity = 0.0;
	else if ((projectedNaOverNaPlusK > 0.35) && [[self phaseName] isEqualToString:@"Panunzite"])    affinity = 0.0;

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

	if ((refMoles[1] == 0.0) || (refMoles[2] == 0.0) ||
		([[self phaseName] isEqualToString:@"Feldspathoid"]) ||
		(([[self phaseName] isEqualToString:@"Nepheline ss"]) &&
		 ([instanceSet countForObject:[@"Panunzite" stringByAppendingString:[self operationParent]]] > 0)) ||
		(([[self phaseName] isEqualToString:@"Panunzite"]) &&
		 ([instanceSet countForObject:[@"Nepheline ss" stringByAppendingString:[self operationParent]]] > 0))
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

	// Specific to feldspathoid ...
	moles[0] = refMoles[1]/refTotalMoles; // switch in projected K for Na
	moles[1] = refMoles[0]/refTotalMoles; // switch in projected Na for K
	moles[2] = refMoles[2]/refTotalMoles; // Na-Vc
	moles[3] = refMoles[3]/refTotalMoles; // Ca-Vc

	if (debugS) {
		NSLog(@"... d = %d, r = %d, rNorm = %10.3e Ne = %10.3e Ks = %10.3e Vc = %10.3e Ca = %10.3e Affinity/RT = %10.3e Ne = %10.3e Ks = %10.3e Vc = %10.3e Ca = %10.3e",
			  NA, NA, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, refMoles[0]/refTotalMoles, refMoles[1]/refTotalMoles, refMoles[2]/refTotalMoles, refMoles[3]/refTotalMoles);
		NSLog(@"... d = %d, r = %d, rNorm = %10.3e Ne = %10.3e Ks = %10.3e Vc = %10.3e Ca = %10.3e Affinity/RT = %10.3e Ne = %10.3e Ks = %10.3e Vc = %10.3e Ca = %10.3e",
			  NA, NA, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, moles[0], moles[1], moles[2], moles[3]);
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

		// Specific to feldspathoid ...
		if (debugS) {
			NSLog(@"... d = %d, r = %ld, rNorm = %10.3e Ne = %10.3e Ks = %10.3e Vc = %10.3e Ca = %10.3e Affinity/RT = %10.3e Ne = %10.3e Ks = %10.3e Vc = %10.3e Ca = %10.3e sL = %10.3e",
				  NA, pseudorank, rNorm, dMatrix[0][0], dMatrix[1][0], dMatrix[2][0], dMatrix[3][0], affinityScaledByRT, moles[0], moles[1], moles[2], moles[3], stepLength);
		}
		// ... end

		iters++;
	}

	if (debugS) {
		NSLog(@"... Soln: Reference composition is a %@", [self phaseName]);
		NSLog(@"... There are %lu instances of Nepheline ss instantiated for this class",
			  [instanceSet countForObject:[@"Nepheline ss" stringByAppendingString:[self operationParent]]]);
		NSLog(@"... There are %lu instances of Panunzite    instantiated for this class",
			  [instanceSet countForObject:[@"Panunzite" stringByAppendingString:[self operationParent]]]);
		NSLog(@"... There are %lu instances of Feldspathoid instantiated for this class",
			  [instanceSet countForObject:[@"Feldspathoid" stringByAppendingString:[self operationParent]]]);
		NSLog(@"... Soln: affinity = %g", affinityScaledByRT*8.314472*t);
		for (NSUInteger i=0; i<NA; i++) NSLog(@"... Soln: %6.3f Ref: %6.3f delta: %13.6e", moles[i], refMoles[i]/refTotalMoles, moles[i]-refMoles[i]/refTotalMoles);
	}

	NSMutableArray *composition = [NSMutableArray arrayWithCapacity:NA];
	for (NSUInteger i=0; i<NA; i++) [composition addObject:[NSNumber numberWithDouble:moles[i]]];
	BOOL additionalPhaseDetected = (affinityScaledByRT > tolerance) ? YES : NO;
	if (!converged) additionalPhaseDetected = NO;

	NSString *nameOfCoexistingPhase = @"";
	if (additionalPhaseDetected) {
		if      ([[self phaseName] isEqualToString:@"Nepheline ss"]) nameOfCoexistingPhase = @"Panunzite";
		else if ([[self phaseName] isEqualToString:@"Panunzite"])    nameOfCoexistingPhase = @"Nepheline ss";
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
	// order Na4Al4Si4O16, K4Al4Si4O16, Na3[]Al3Si5O16, CaNa2[]Al4Si4O16
	double totProjectedNa = refMoles[0];
	double totProjectedK  = refMoles[1];

	if ((totProjectedNa+totProjectedK) == 0.0) return [NSString stringWithString:[self phaseName]];
	if ((totProjectedNa/(totProjectedNa+totProjectedK)) < 0.35) return @"Panunzite";
	else                                                        return @"Nepheline ss";
}

#import "SolutionPhase.h"

@end
