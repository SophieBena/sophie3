//
//  RhombohedralBerman.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 12/4/10.
//  Copyright 2010 OFM Research Inc. All rights reserved.
//

#import "RhombohedralBerman.h"
#import "BermanProperties.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"

@implementation RhombohedralBerman

static NSArray *endmembers;

#define NR       4  /* Four independent composition variables      */
#define NS       3  /* Three ordering parameters                   */
#define NA       5  /* Five endmember compositions                 */
#define NATOMS 5.0  /* Average number of atoms in the formula unit */

#pragma mark -
#pragma mark class methods

+(void)initialize {
	if (self == [RhombohedralBerman class]) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"Initialize(RhombohedralBerman) - entry ...");
		NSMutableArray *mutableEndmembers = [NSMutableArray arrayWithCapacity:NA];

		BermanProperties *geikielite = [[BermanProperties alloc] initWithH:-1572560.0
																		 S:74.56
																		k0:146.20
																		k1:-4.160E2
																		k2:-39.998E5
																		k3:40.233E7
																		v0:3.086
																		v1:-0.584e-6
																		v2:1.230e-12
																		v3:27.248e-6
																		v4:29.968e-10];
		[geikielite setPhaseFormula:@"MgTiO3"];
		[geikielite setPhaseName:@"geikielite"];
		[mutableEndmembers addObject:geikielite];
		if (debug) NSLog(@"... allocated geikielite ...");

		BermanProperties *hematite = [[BermanProperties alloc] initWithH:-825627.0
																	   S:87.437
																	  k0:146.86
																	  k1:0.0
																	  k2:-55.768e5
																	  k3:52.563e7
																	  l1:-7.403e-2
																	  l2:27.921e-5
																	  Tt:955.0
																  deltaH:1287.0
																	  v0:3.027
																	  v1:-0.479e-6
																	  v2:0.304e-12
																	  v3:38.310e-6
																	  v4:1.650e-10];
		[hematite setPhaseFormula:@"Fe2O3"];
		[hematite setPhaseName:@"hematite"];
		[mutableEndmembers addObject:hematite];
		if (debug) NSLog(@"... allocated hematite ...");

		BermanProperties *ilmenite = [[BermanProperties alloc] initWithH:-1231947.0
																	   S:108.628
																	  k0:150.0
																	  k1:-4.416e2
																	  k2:-33.237e5
																	  k3:34.815e7
																	  v0:3.170
																	  v1:-0.584e-6
																	  v2:1.230e-12
																	  v3:27.248e-6
																	  v4:29.968e-10];
		[ilmenite setPhaseFormula:@"FeTiO3"];
		[ilmenite setPhaseName:@"ilmenite"];
		[mutableEndmembers addObject:ilmenite];
		if (debug) NSLog(@"... allocated ilmenite ...");

		BermanProperties *pyrophanite = [[BermanProperties alloc] initWithH:-1350707.0
																		  S:104.935
																		 k0:150.00
																		 k1:-4.416E2
																		 k2:-33.237E5
																		 k3:34.815E7
																		 v0:2.8859
																		 v1:-0.584e-6
																		 v2:1.230e-12
																		 v3:27.248e-6
																		 v4:29.968e-10];
		[pyrophanite setPhaseFormula:@"MnTiO3"];
		[pyrophanite setPhaseName:@"pyrophanite"];
		[mutableEndmembers addObject:pyrophanite];
		if (debug) NSLog(@"... allocated pyrophanite ...");

		// the adjustment of 20 kJ is arbitrary just to keep Al in the 1 - 3% level
		BermanProperties *corundum = [[BermanProperties alloc] initWithH:-1675700.0+20000.0
																	   S:50.820
																	  k0:155.02
																	  k1:-8.284E2
																	  k2:-38.614E5
																	  k3:40.908E7
																	  v0:2.558
																	  v1:-0.385E-6
																	  v2:0.375E-12
																	  v3:21.342E-6
																	  v4:47.180E-10];
		[corundum setPhaseFormula:@"Al2O3"];
		[corundum setPhaseName:@"corundum"];
		[mutableEndmembers addObject:corundum];
		if (debug) NSLog(@"... allocated corundum ...");

		endmembers = [NSArray arrayWithArray:mutableEndmembers];
	}
}

#pragma mark -
#pragma mark instance methods

-(id)init {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(RhombohedralBerman) ... entry ...");
		[self setPhaseName:@"Ilmenite ss"];
		computeMixingQuantities = NO;
		tOld = -9999.0;
		pOld = -9999.0;
		tOldPure = -9999.0;
		pOldPure = -9999.0;
		for (NSUInteger i=0; i<NR; i++) rOld[i] = -9999.0;
		for (NSUInteger i=0; i<NS; i++) sOld[i] = 2.0;
		for (NSUInteger i=0; i<NS; i++) sOldPure[i] = 2.0;
		xmg2a = 0.0;
		xfe2a = 0.0;
		xmn2a = 0.0;
		xti4a = 0.0;
		xal3a = 0.0;
		xfe3a = 0.0;
		xmg2b = 0.0;
		xfe2b = 0.0;
		xmn2b = 0.0;
		xti4b = 0.0;
		xal3b = 0.0;
		xfe3b = 0.0;
		xmg2ID = 0.0;
		xfe2ID = 0.0;
		xmn2ID = 0.0;
		xti4ID = 0.0;
		xal3ID = 0.0;
		xfe3ID = 0.0;
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

#define MAX_ITER 100    /* Maximum number of iterations allowed in order */

#define SQUARE(x)  ((x)*(x))
#define CUBE(x)    ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Rhombohedral solution parameters:
 * Ghiorso, M.S., Evans, B.W. (2008)
 */

static const double dvilm     =  0.010758; /* joules/bar */
static const double dvgei     =  0.010758; /* joules/bar */
static const double dvpyr     =  0.010758; /* joules/bar */

static const double wvilm     =  0.035089; /* joules/bar */
static const double wvgei     =  0.035089; /* joules/bar */
static const double wvpyr     =  0.035089; /* joules/bar */

static const double dwvhmilm  =  0.013701; /* joules/bar */
static const double dwvhmgei  =  0.013701; /* joules/bar */
static const double dwvhmpyr  =  0.013701; /* joules/bar */

static const double wvhmilm2  =  0.000000; /* joules/bar */
static const double wvhmgei2  =  0.000000; /* joules/bar */
static const double wvhmpyr2  =  0.000000; /* joules/bar */

static const double dwvcrnilm =  0.013701; /* joules/bar - guess */
static const double dwvcrngei =  0.013701; /* joules/bar - guess */
static const double dwvcrnpyr =  0.013701; /* joules/bar - guess */

static const double wvhmilm   = -0.11764;  /* joules/bar */
static const double wvhmgei   = -0.11764;  /* joules/bar - guess */
static const double wvhmpyr   = -0.11764;  /* joules/bar - guess */
static const double wvilmgei  =  0.000000; /* joules/bar - guess */
static const double wvilmpyr  =  0.000000; /* joules/bar - guess */
static const double wvgeipyr  =  0.000000; /* joules/bar - guess */

static const double dhilm     =  17477.0; /* joules ordering analysis     	- fixed    */
static const double dhgei     =  17477.0; /* joules ordering analysis     	- fixed    */
static const double dhpyr     =  17477.0; /* joules ordering analysis     	- fixed    */

static const double whilm     =   3189.0; /* joules ordering analysis     	- fixed    */
static const double whgei     =   3189.0; /* joules ordering analysis     	- fixed    */
static const double whpyr     =   3189.0; /* joules ordering analysis     	- fixed    */

static const double dwhhmilm  =  -5626.63; /* -1510.8 joules ordering analysis - dwhhmilm/2+whhmilm2/4 is fixed at -3021.6 */
static const double dwhhmgei  =  -5626.63; /* -1510.8 joules ordering analysis - dwhhmgei/2+whhmgei2/4 is fixed at -3021.6 */
static const double dwhhmpyr  =  -5626.63; /* -1510.8 joules ordering analysis - dwhhmpyr/2+whhmpyr2/4 is fixed at -3021.6 */

static const double whhmilm2  =   -833.14; /* -9064.8 joules ordering analysis - dwhhmilm/2+whhmilm2/4 is fixed at -3021.6 */
static const double whhmgei2  =   -833.14; /* -9064.8 joules ordering analysis - dwhhmgei/2+whhmgei2/4 is fixed at -3021.6 */
static const double whhmpyr2  =   -833.14; /* -9064.8 joules ordering analysis - dwhhmpyr/2+whhmpyr2/4 is fixed at -3021.6 */

static const double dwhcrnilm =      0.0; /* joules - guess                   - fixed    */
static const double dwhcrngei =      0.0; /* joules - guess                   - fixed    */
static const double dwhcrnpyr =      0.0; /* joules - guess                   - fixed    */

static const double whmcrn    =  69000.0; /* 69000.0 joules Majzlan et al., 2002      - fixed   */
static const double whhmilm   =  22535.6; /* 21200.0 joules solvus analysis	       - variable */
static const double whhmgei   =  22535.6; /* 21200.0 joules guess		       - variable */
static const double whhmpyr   =  22535.6; /* 21200.0 joules guess		       - variable */
static const double wcrnilm   =  22535.6; /* 21200.0 joules guess		       - variable */
static const double wcrngei   =  22535.6; /* 21200.0 joules guess		       - variable */
static const double wcrnpyr   =  22535.6; /* 21200.0 joules guess		       - variable */
static const double whilmgei  =   2600.0; /*  2600.0 joules Pownceby and Fisher-White - fixed   */
static const double whilmpyr  =   2200.0; /*  2200.0 joules O'Neill et al.	       - fixed    */
static const double whgeipyr  =   2600.0; /*  2600.0 joules analogy with ilm-gei      - fixed   */

static const double whilmgeiT =  88099.9; /* 63694.0 guess, equal to whilmgei zeroes (G*)st	  */
static const double whilmpyrT =  30244.0; /* 45214.0 guess, equal to whilmpyr zeroes (G*)su	  */
static const double whgeipyrT =   2600.0; /* 54454.0 guess, equal to whilmpyr zeroes (G*)tu	  */

static const double whilmilmgei =    0.0; /* joules - guess */
static const double whilmgeigei =    0.0; /* joules - guess */
static const double whilmilmpyr =    0.0; /* joules - guess */
static const double whilmpyrpyr =    0.0; /* joules - guess */
static const double whgeigeipyr =    0.0; /* joules - guess */
static const double whgeipyrpyr =    0.0; /* joules - guess */

static const double SROconst = 0.0730205; /* Scaling factor for SRO */

#define R  8.3143

/*
 *=============================================================================
 * Local function to compute ordering state and associated derivatives of
 * pure component endmembers
 */

#define IL_S          -R*((1.0+s[0])*log(1.0+s[0]) + (1.0-s[0])*log(1.0-s[0]) - 2.0*log(2.0))
#define IL_H          (dhilm+(p-1.0)*dvilm)*(1.0 - s[0]*s[0]) + (whilm+(p-1.0)*wvilm)*s[0]*s[0]*(1.0 - s[0]*s[0])
#define IL_G          IL_H - t*(IL_S)

#define DIL_GDS0      R*t*(log(1.0+s[0]) - log(1.0-s[0])) \
- 2.0*(dhilm+(p-1.0)*dvilm)*s[0] + (whilm+(p-1.0)*wvilm)*(2.0*s[0] - 4.0*s[0]*s[0]*s[0])
#define DIL_GDT       -(IL_S)
#define DIL_GDP       dvilm*(1.0 - s[0]*s[0]) + wvilm*s[0]*s[0]*(1.0 - s[0]*s[0])

#define D2IL_GDS0S0   R*t*(1.0/(1.0+s[0]) + 1.0/(1.0-s[0])) - 2.0*(dhilm+(p-1.0)*dvilm) + (whilm+(p-1.0)*wvilm)*(2.0 - 12.0*s[0]*s[0])
#define D2IL_GDS0DT   R*(log(1.0+s[0]) - log(1.0-s[0]))
#define D2IL_GDS0DP   - 2.0*dvilm*s[0] + wvilm*(2.0*s[0] - 4.0*s[0]*s[0]*s[0])
#define D2IL_GDT2     0.0
#define D2IL_GDTP     0.0
#define D2IL_GDP2     0.0

#define D3IL_GDS0S0S0 R*t*(-1.0/SQUARE(1.0+s[0]) + 1.0/SQUARE(1.0-s[0])) - 24.0*(whilm+(p-1.0)*wvilm)*s[0]
#define D3IL_GDS0S0DT R*(1.0/(1.0+s[0]) + 1.0/(1.0-s[0]))
#define D3IL_GDS0S0DP - 2.0*dvilm + wvilm*(2.0 - 12.0*s[0]*s[0])
#define D3IL_GDS0DT2  0.0
#define D3IL_GDS0DTDP 0.0
#define D3IL_GDS0DP2  0.0
#define D3IL_GDT3     0.0
#define D3IL_GDT2DP   0.0
#define D3IL_GDTDP2   0.0
#define D3IL_GDP3     0.0

#define GK_S          -R*((1.0+s[1])*log(1.0+s[1]) + (1.0-s[1])*log(1.0-s[1]) - 2.0*log(2.0))
#define GK_H          (dhgei+(p-1.0)*dvgei)*(1.0 - s[1]*s[1]) + (whgei+(p-1.0)*wvgei)*s[1]*s[1]*(1.0 - s[1]*s[1])
#define GK_G          GK_H - t*(GK_S)

#define DGK_GDS1      R*t*(log(1.0+s[1]) - log(1.0-s[1])) \
- 2.0*(dhgei+(p-1.0)*dvgei)*s[1] + (whgei+(p-1.0)*wvgei)*(2.0*s[1] - 4.0*s[1]*s[1]*s[1])
#define DGK_GDT       -(GK_S)
#define DGK_GDP       dvgei*(1.0 - s[1]*s[1]) + wvgei*s[1]*s[1]*(1.0 - s[1]*s[1])

#define D2GK_GDS1S1   R*t*(1.0/(1.0+s[1]) + 1.0/(1.0-s[1])) - 2.0*(dhgei+(p-1.0)*dvgei) + (whgei+(p-1.0)*wvgei)*(2.0 - 12.0*s[1]*s[1])
#define D2GK_GDS1DT   R*(log(1.0+s[1]) - log(1.0-s[1]))
#define D2GK_GDS1DP   - 2.0*dvgei*s[1] + wvgei*(2.0*s[1] - 4.0*s[1]*s[1]*s[1])
#define D2GK_GDT2     0.0
#define D2GK_GDTP     0.0
#define D2GK_GDP2     0.0

#define D3GK_GDS1S1S1 R*t*(-1.0/SQUARE(1.0+s[1]) + 1.0/SQUARE(1.0-s[1])) - 24.0*(whgei+(p-1.0)*wvgei)*s[1]
#define D3GK_GDS1S1DT R*(1.0/(1.0+s[1]) + 1.0/(1.0-s[1]))
#define D3GK_GDS1S1DP - 2.0*dvgei + wvgei*(2.0 - 12.0*s[1]*s[1])
#define D3GK_GDS1DT2  0.0
#define D3GK_GDS1DTDP 0.0
#define D3GK_GDS1DP2  0.0
#define D3GK_GDT3     0.0
#define D3GK_GDT2DP   0.0
#define D3GK_GDTDP2   0.0
#define D3GK_GDP3     0.0

#define PY_S          -R*((1.0+s[2])*log(1.0+s[2]) + (1.0-s[2])*log(1.0-s[2]) - 2.0*log(2.0))
#define PY_H          (dhpyr+(p-1.0)*dvpyr)*(1.0 - s[2]*s[2]) + (whpyr+(p-1.0)*wvpyr)*s[2]*s[2]*(1.0 - s[2]*s[2])
#define PY_G          PY_H - t*(PY_S)

#define DPY_GDS2      R*t*(log(1.0+s[2]) - log(1.0-s[2])) \
- 2.0*(dhpyr+(p-1.0)*dvpyr)*s[2] + (whpyr+(p-1.0)*wvpyr)*(2.0*s[2] - 4.0*s[2]*s[2]*s[2])
#define DPY_GDT       -(PY_S)
#define DPY_GDP       dvpyr*(1.0 - s[2]*s[2]) + wvpyr*s[2]*s[2]*(1.0 - s[2]*s[2])

#define D2PY_GDS2S2   R*t*(1.0/(1.0+s[2]) + 1.0/(1.0-s[2])) - 2.0*(dhpyr+(p-1.0)*dvpyr) + (whpyr+(p-1.0)*wvpyr)*(2.0 - 12.0*s[2]*s[2])
#define D2PY_GDS2DT   R*(log(1.0+s[2]) - log(1.0-s[2]))
#define D2PY_GDS2DP   - 2.0*dvpyr*s[2] + wvpyr*(2.0*s[2] - 4.0*s[2]*s[2]*s[2])
#define D2PY_GDT2     0.0
#define D2PY_GDTP     0.0
#define D2PY_GDP2     0.0

#define D3PY_GDS2S2S2 R*t*(-1.0/SQUARE(1.0+s[2]) + 1.0/SQUARE(1.0-s[2])) - 24.0*(whpyr+(p-1.0)*wvpyr)*s[2]
#define D3PY_GDS2S2DT R*(1.0/(1.0+s[2]) + 1.0/(1.0-s[2]))
#define D3PY_GDS2S2DP - 2.0*dvpyr + wvpyr*(2.0 - 12.0*s[2]*s[2])
#define D3PY_GDS2DT2  0.0
#define D3PY_GDS2DTDP 0.0
#define D3PY_GDS2DP2  0.0
#define D3PY_GDT3     0.0
#define D3PY_GDT2DP   0.0
#define D3PY_GDTDP2   0.0
#define D3PY_GDP3     0.0

#define HM_S          SROconst*2.0*R*log(2.0)
#define HM_H          0.0
#define HM_G          HM_H - t*(HM_S)
#define DHM_GDT       - SROconst*2.0*R*log(2.0)
#define DHM_GDP       0.0
#define D2HM_GDT2     0.0
#define D2HM_GDTP     0.0
#define D2HM_GDP2     0.0
#define D3HM_GDT3     0.0
#define D3HM_GDT2DP   0.0
#define D3HM_GDTDP2   0.0
#define D3HM_GDP3     0.0

#define CR_S          SROconst*2.0*R*log(2.0)
#define CR_H          0.0
#define CR_G          CR_H - t*(CR_S)
#define DCR_GDT       - SROconst*2.0*R*log(2.0)
#define DCR_GDP       0.0
#define D2CR_GDT2     0.0
#define D2CR_GDTP     0.0
#define D2CR_GDP2     0.0
#define D3CR_GDT3     0.0
#define D3CR_GDT2DP   0.0
#define D3CR_GDTDP2   0.0
#define D3CR_GDP3     0.0

#define fillD2GDSDT   d2gdsdt[0]   = D2IL_GDS0DT;   d2gdsdt[1]   = D2GK_GDS1DT;   d2gdsdt[2]   = D2PY_GDS2DT;
#define fillD2GDSDP   d2gdsdp[0]   = D2IL_GDS0DP;   d2gdsdp[1]   = D2GK_GDS1DP;   d2gdsdp[2]   = D2PY_GDS2DP;
#define fillD3GDS3    d3gds3[0]    = D3IL_GDS0S0S0; d3gds3[1]    = D3GK_GDS1S1S1; d3gds3[2]    = D3PY_GDS2S2S2;
#define fillD3GDS2DT  d3gds2dt[0]  = D3IL_GDS0S0DT; d3gds2dt[1]  = D3GK_GDS1S1DT; d3gds2dt[2]  = D3PY_GDS2S2DT;
#define fillD3GDS2DP  d3gds2dp[0]  = D3IL_GDS0S0DP; d3gds2dp[1]  = D3GK_GDS1S1DP; d3gds2dp[2]  = D3PY_GDS2S2DP;
#define fillD3GDSDT2  d3gdsdt2[0]  = D3IL_GDS0DT2;  d3gdsdt2[1]  = D3GK_GDS1DT2;  d3gdsdt2[2]  = D3PY_GDS2DT2;
#define fillD3GDSDTDP d3gdsdtdp[0] = D3IL_GDS0DTDP; d3gdsdtdp[1] = D3GK_GDS1DTDP; d3gdsdtdp[2] = D3PY_GDS2DTDP;
#define fillD3GDSDP2  d3gdsdp2[0]  = D3IL_GDS0DP2;  d3gdsdp2[1]  = D3GK_GDS1DP2;  d3gdsdp2[2]  = D3PY_GDS2DP2;

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
	const int iterMax = 1000;
	int i;

	if ((t != tOldPure) || (p != pOldPure)) {
		double dgds[NS], sNew[NS];
		int iter = 0;
		for (i=0; i<NS; i++) { sOldPure[i] = 2.0; sNew[i] = 0.98; }
		while ((iter < iterMax) && (
									(fabs(sNew[0]-sOldPure[0]) > 10.0*DBL_EPSILON) ||
									(fabs(sNew[1]-sOldPure[1]) > 10.0*DBL_EPSILON) ||
									(fabs(sNew[2]-sOldPure[2]) > 10.0*DBL_EPSILON) )) {
			double s[NS];

			for (i=0; i<NS; i++) s[i] = sNew[i];

			dgds[0] = DIL_GDS0;
			dgds[1] = DGK_GDS1;
			dgds[2] = DPY_GDS2;

			d2gds2Pure[0] = D2IL_GDS0S0;
			d2gds2Pure[1] = D2GK_GDS1S1;
			d2gds2Pure[2] = D2PY_GDS2S2;

			for (i=0; i<NS; i++) sOldPure[i] = s[i];

			for (i=0; i<NS; i++) {
				s[i] += - dgds[i]/d2gds2Pure[i];
				s[i] = MIN(s[i],  1.0 - 10.0*DBL_EPSILON);
				s[i] = MAX(s[i],  0.0);
			}

			for (i=0; i<NS; i++) sNew[i] = s[i];
			iter++;
		}
		tOldPure = t;
		pOldPure = p;
		if (iter == iterMax) {
			// Condition arises at a T above the 1st order transition when Newton's method has
			// arrived at a local minimum.  In this case, the global minimum is at zero.
			for (i=0; i<NS; i++) sOldPure[i] = 0.0;
		}
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
		double *s = sOldPure;
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
		for (i=0; i<NS; i++) dt2[i] = - (d3gdsdt2[i] + 2.0*d3gds2dt[i]*dsdt[i] + d3gds3[i]*dsdt[i]*dsdt[i])/d2gds2Pure[i];
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
			dtp[i] = - (d3gdsdtdp[i] + d3gds2dt[i]*dsdp[i] + d3gds2dp[i]*dsdt[i] + d3gds3[i]*dsdt[i]*dsdp[i])/d2gds2Pure[i];

	}

	if (mask & SIXTH  ) {   /* compute d2s/dp2 */
		double *s = sOldPure;
		double d2gdsdp[NS], d3gds3[NS], d3gds2dp[NS], d3gdsdp2[NS], dsdp[NS];

		fillD2GDSDP
		fillD3GDS3
		fillD3GDS2DP
		fillD3GDSDP2

		for (i=0; i<NS; i++) dsdp[i] = - d2gdsdp[i]/d2gds2Pure[i];

		for (i=0; i<NS; i++) dp2[i] = - (d3gdsdp2[i] + 2.0*d3gds2dp[i]*dsdp[i] + d3gds3[i]*dsdp[i]*dsdp[i])/d2gds2Pure[i];
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
d2gds2[0][0] = 0.0;         d2gds2[0][1] = D2GK_GDS1S1;  d2gds2[0][2] = 0.0;         \
d2gds2[1][0] = 0.0;         d2gds2[1][1] = 0.0;          d2gds2[1][2] = 0.0;         \
d2gds2[2][0] = D2IL_GDS0S0; d2gds2[2][1] = 0.0;          d2gds2[2][2] = 0.0;         \
d2gds2[3][0] = 0.0;         d2gds2[3][1] = 0.0;          d2gds2[3][2] = D2PY_GDS2S2; \
d2gds2[4][0] = 0.0;         d2gds2[4][1] = 0.0;          d2gds2[4][2] = 0.0;

#define fillD2GDSDT \
d2gdsdt[0][0] = 0.0;         d2gdsdt[0][1] = D2GK_GDS1DT;  d2gdsdt[0][2] =  0.0;        \
d2gdsdt[1][0] = 0.0;         d2gdsdt[1][1] = 0.0;          d2gdsdt[1][2] =  0.0;        \
d2gdsdt[2][0] = D2IL_GDS0DT; d2gdsdt[2][1] = 0.0;          d2gdsdt[2][2] =  0.0;        \
d2gdsdt[3][0] = 0.0;         d2gdsdt[3][1] = 0.0;          d2gdsdt[3][2] = D2PY_GDS2DT; \
d2gdsdt[4][0] = 0.0;         d2gdsdt[4][1] = 0.0;          d2gdsdt[4][2] =  0.0;

#define fillD2GDSDP \
d2gdsdp[0][0] = 0.0;         d2gdsdp[0][1] = D2GK_GDS1DP;  d2gdsdp[0][2] = 0.0;         \
d2gdsdp[1][0] = 0.0;         d2gdsdp[1][1] = 0.0;          d2gdsdp[1][2] = 0.0;         \
d2gdsdp[2][0] = D2IL_GDS0DP; d2gdsdp[2][1] = 0.0;          d2gdsdp[2][2] = 0.0;         \
d2gdsdp[3][0] = 0.0;         d2gdsdp[3][1] = 0.0;          d2gdsdp[3][2] = D2PY_GDS2DP; \
d2gdsdp[4][0] = 0.0;         d2gdsdp[4][1] = 0.0;          d2gdsdp[4][2] = 0.0;

#define fillD2GDT2 \
d2gdt2[0]  = D2GK_GDT2; d2gdt2[1]  = D2HM_GDT2; d2gdt2[2]  = D2IL_GDT2; d2gdt2[3]  = D2PY_GDT2; d2gdt2[4]  = D2CR_GDT2;

#define fillD2GDTDP \
d2gdtdp[0] = D2GK_GDTP; d2gdtdp[1] = D2HM_GDTP; d2gdtdp[2] = D2IL_GDTP; d2gdtdp[3] = D2PY_GDTP; d2gdtdp[4] = D2CR_GDTP;

#define fillD2GDP2 \
d2gdp2[0]  = D2GK_GDP2; d2gdp2[1]  = D2HM_GDP2; d2gdp2[2]  = D2IL_GDP2; d2gdp2[3]  = D2PY_GDP2; d2gdp2[4]  = D2CR_GDP2;

#define fillD3GDS3 \
d3gds3[0][0] = 0.0;           d3gds3[0][1] = D3GK_GDS1S1S1; d3gds3[0][2] = 0.0;           \
d3gds3[1][0] = 0.0;           d3gds3[1][1] = 0.0;           d3gds3[1][2] = 0.0;           \
d3gds3[2][0] = D3IL_GDS0S0S0; d3gds3[2][1] = 0.0;           d3gds3[2][2] = 0.0;           \
d3gds3[3][0] = 0.0;           d3gds3[3][1] = 0.0;           d3gds3[3][2] = D3PY_GDS2S2S2; \
d3gds3[4][0] = 0.0;           d3gds3[4][1] = 0.0;           d3gds3[4][2] = 0.0;

#define fillD3GDS2DT \
d3gds2dt[0][0] = 0.0;           d3gds2dt[0][1] = D3GK_GDS1S1DT; d3gds2dt[0][2] = 0.0;           \
d3gds2dt[1][0] = 0.0;           d3gds2dt[1][1] = 0.0;           d3gds2dt[1][2] = 0.0;           \
d3gds2dt[2][0] = D3IL_GDS0S0DT; d3gds2dt[2][1] = 0.0;           d3gds2dt[2][2] = 0.0;           \
d3gds2dt[3][0] = 0.0;           d3gds2dt[3][1] = 0.0;           d3gds2dt[3][2] = D3PY_GDS2S2DT; \
d3gds2dt[4][0] = 0.0;           d3gds2dt[4][1] = 0.0;           d3gds2dt[4][2] = 0.0;

#define fillD3GDS2DP \
d3gds2dp[0][0] = 0.0;           d3gds2dp[0][1] = D3GK_GDS1S1DP; d3gds2dp[0][2] = 0.0;           \
d3gds2dp[1][0] = 0.0;           d3gds2dp[1][1] = 0.0;           d3gds2dp[1][2] = 0.0;           \
d3gds2dp[2][0] = D3IL_GDS0S0DP; d3gds2dp[2][1] = 0.0;           d3gds2dp[2][2] = 0.0;           \
d3gds2dp[3][0] = 0.0;           d3gds2dp[3][1] = 0.0;           d3gds2dp[3][2] = D3PY_GDS2S2DP; \
d3gds2dp[4][0] = 0.0;           d3gds2dp[4][1] = 0.0;           d3gds2dp[4][2] = 0.0;

#define fillD3GDSDT2 \
d3gdsdt2[0][0] = 0.0;          d3gdsdt2[0][1] = D3GK_GDS1DT2; d3gdsdt2[0][2] = 0.0;           \
d3gdsdt2[1][0] = 0.0;          d3gdsdt2[1][1] = 0.0;          d3gdsdt2[1][2] = 0.0;           \
d3gdsdt2[2][0] = D3IL_GDS0DT2; d3gdsdt2[2][1] = 0.0;          d3gdsdt2[2][2] = 0.0;           \
d3gdsdt2[3][0] = 0.0;          d3gdsdt2[3][1] = 0.0;          d3gdsdt2[3][2] =  D3PY_GDS2DT2; \
d3gdsdt2[4][0] = 0.0;          d3gdsdt2[4][1] = 0.0;          d3gdsdt2[4][2] = 0.0;

#define fillD3GDSDTDP \
d3gdsdtdp[0][0] = 0.0;           d3gdsdtdp[0][1] = D3GK_GDS1DTDP; d3gdsdtdp[0][2] = 0.0;           \
d3gdsdtdp[1][0] = 0.0;           d3gdsdtdp[1][1] = 0.0;           d3gdsdtdp[1][2] = 0.0;           \
d3gdsdtdp[2][0] = D3IL_GDS0DTDP; d3gdsdtdp[2][1] = 0.0;           d3gdsdtdp[2][2] = 0.0;           \
d3gdsdtdp[3][0] = 0.0;           d3gdsdtdp[3][1] = 0.0;           d3gdsdtdp[3][2] = D3PY_GDS2DTDP; \
d3gdsdtdp[4][0] = 0.0;           d3gdsdtdp[4][1] = 0.0;           d3gdsdtdp[4][2] = 0.0;

#define fillD3GDSDP2 \
d3gdsdp2[0][0] = 0.0;          d3gdsdp2[0][1] = D3GK_GDS1DP2; d3gdsdp2[0][2] = 0.0;          \
d3gdsdp2[1][0] = 0.0;          d3gdsdp2[1][1] = 0.0;          d3gdsdp2[1][2] = 0.0;          \
d3gdsdp2[2][0] = D3IL_GDS0DP2; d3gdsdp2[2][1] = 0.0;          d3gdsdp2[2][2] = 0.0;          \
d3gdsdp2[3][0] = 0.0;          d3gdsdp2[3][1] = 0.0;          d3gdsdp2[3][2] = D3PY_GDS2DP2; \
d3gdsdp2[4][0] = 0.0;          d3gdsdp2[4][1] = 0.0;          d3gdsdp2[4][2] = 0.0;

#define fillD3GDT3 \
d3gdt3[0]   = D3GK_GDT3;   d3gdt3[1]   = D3HM_GDT3;   d3gdt3[2]   = D3IL_GDT3;   d3gdt3[3]   = D3PY_GDT3;   d3gdt3[4]   = D3CR_GDT3;

#define fillD3GDT2DP \
d3gdt2dp[0] = D3GK_GDT2DP; d3gdt2dp[1] = D3HM_GDT2DP; d3gdt2dp[2] = D3IL_GDT2DP; d3gdt2dp[3] = D3PY_GDT2DP; d3gdt2dp[4] = D3CR_GDT2DP;

#define fillD3GDTDP2 \
d3gdtdp2[0] = D3GK_GDTDP2; d3gdtdp2[1] = D3HM_GDTDP2; d3gdtdp2[2] = D3IL_GDTDP2; d3gdtdp2[3] = D3PY_GDTDP2; d3gdtdp2[4] = D3CR_GDTDP2;

#define fillD3GDP3 \
d3gdp3[0]   = D3GK_GDP3;   d3gdp3[1]   = D3HM_GDP3;   d3gdp3[2]   = D3IL_GDP3;   d3gdp3[3]   = D3PY_GDP3;   d3gdp3[4]   = D3CR_GDP3;

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
		a[0] = GK_G;
		a[0] = exp(a[0]/(R*t));
		a[1] = HM_G;
		a[1] = exp(a[1]/(R*t));
		a[2] = IL_G;
		a[2] = exp(a[2]/(R*t));
		a[3] = PY_G;
		a[3] = exp(a[3]/(R*t));
		a[4] = CR_G;
		a[4] = exp(a[4]/(R*t));
	}

	if (mask & SECOND) {
		mu[0] = GK_G;
		mu[1] = HM_G;
		mu[2] = IL_G;
		mu[3] = PY_G;
		mu[4] = CR_G;
	}

	if (mask & THIRD) {
		gmix[0] = GK_G;
		gmix[1] = HM_G;
		gmix[2] = IL_G;
		gmix[3] = PY_G;
		gmix[4] = CR_G;
	}

	if (mask & FOURTH) {
		hmix[0] = (GK_G) + t*(GK_S);
		hmix[1] = (HM_G) - t*(DHM_GDT);
		hmix[2] = (IL_G) + t*(IL_S);
		hmix[3] = (PY_G) + t*(PY_S);
		hmix[4] = (CR_G) - t*(DCR_GDT);
	}

	if (mask & FIFTH) {
		smix[0] = GK_S;
		smix[1] = -(DHM_GDT);
		smix[2] = IL_S;
		smix[3] = PY_S;
		smix[4] = -(DCR_GDT);
	}

	if (mask & SIXTH) {
		double d2gdsdt[NA][NS], d2gds2[NA][NS], dsdt[NS];

		fillD2GDS2
		fillD2GDSDT

		[self pureOrder:SECOND t:t p:p s:NULL dt:dsdt dp:NULL dt2:NULL dtp:NULL dp2:NULL];

		cpmix[0] = D2GK_GDT2;
		cpmix[1] = D2HM_GDT2;
		cpmix[2] = D2IL_GDT2;
		cpmix[3] = D2PY_GDT2;
		cpmix[4] = D2CR_GDT2;

		for (i=0; i<NA; i++) {
			for (j=0; j<NS; j++) cpmix[i] += 2.0*d2gdsdt[i][j]*dsdt[j] + d2gds2[i][j]*SQUARE(dsdt[j]);
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
		vmix[0] = DGK_GDP;
		vmix[1] = DHM_GDP;
		vmix[2] = DIL_GDP;
		vmix[3] = DPY_GDP;
		vmix[4] = DCR_GDP;
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
			for (j=0; j<NS; j++) vmixdt[i] += d2gdsdt[i][j]*dsdp[j] + d2gdsdp[i][j]*dsdt[j]
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
			for (j=0; j<NS; j++) vmixdp[i] += 2.0*d2gdsdp[i][j]*dsdp[j] + d2gds2[i][j]*dsdp[j]*dsdp[j];
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
			for (j=0; j<NS; j++) vmixdtdp[i] += 2.0*d3gdsdtdp[i][j]*dsdp[j] + d2gdsdt[i][j]*d2sdp2[j]
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
			for (j=0; j<NS; j++) vmixdp2[i] += 3.0*d3gdsdp2[i][j]*dsdp[j] + 3.0*d2gdsdp[i][j]*d2sdp2[j]
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

#undef IL_S
#undef IL_H
#undef IL_G
#undef DIL_GDS0
#undef DIL_GDT
#undef DIL_GDP
#undef D2IL_GDS0S0
#undef D2IL_GDS0DT
#undef D2IL_GDS0DP
#undef D2IL_GDT2
#undef D2IL_GDTP
#undef D2IL_GDP2
#undef D3IL_GDS0S0S0
#undef D3IL_GDS0S0DT
#undef D3IL_GDS0S0DP
#undef D3IL_GDS0DT2
#undef D3IL_GDS0DTDP
#undef D3IL_GDS0DP2
#undef D3IL_GDT3
#undef D3IL_GDT2DP
#undef D3IL_GDTDP2
#undef D3IL_GDP3

#undef GK_S
#undef GK_H
#undef GK_G
#undef DGK_GDS1
#undef DGK_GDT
#undef DGK_GDP
#undef D2GK_GDS1S1
#undef D2GK_GDS1DT
#undef D2GK_GDS1DP
#undef D2GK_GDT2
#undef D2GK_GDTP
#undef D2GK_GDP2
#undef D3GK_GDS1S1S1
#undef D3GK_GDS1S1DT
#undef D3GK_GDS1S1DP
#undef D3GK_GDS1DT2
#undef D3GK_GDS1DTDP
#undef D3GK_GDS1DP2
#undef D3GK_GDT3
#undef D3GK_GDT2DP
#undef D3GK_GDTDP2
#undef D3GK_GDP3

#undef PY_S
#undef PY_H
#undef PY_G
#undef DPY_GDS2
#undef DPY_GDT
#undef DPY_GDP
#undef D2PY_GDS2S2
#undef D2PY_GDS2DT
#undef D2PY_GDS2DP
#undef D2PY_GDT2
#undef D2PY_GDTP
#undef D2PY_GDP2
#undef D3PY_GDS2S2S2
#undef D3PY_GDS2S2DT
#undef D3PY_GDS2S2DP
#undef D3PY_GDS2DT2
#undef D3PY_GDS2DTDP
#undef D3PY_GDS2DP2
#undef D3PY_GDT3
#undef D3PY_GDT2DP
#undef D3PY_GDTDP2
#undef D3PY_GDP3

#undef HM_S
#undef HM_H
#undef HM_G
#undef DHM_GDT
#undef DHM_GDP
#undef D2HM_GDT2
#undef D2HM_GDTP
#undef D2HM_GDP2
#undef D3HM_GDT3
#undef D3HM_GDT2DP
#undef D3HM_GDTDP2
#undef D3HM_GDP3

#undef CR_S
#undef CR_H
#undef CR_G
#undef DCR_GDT
#undef DCR_GDP
#undef D2CR_GDT2
#undef D2CR_GDTP
#undef D2CR_GDP2
#undef D3CR_GDT3
#undef D3CR_GDT2DP
#undef D3CR_GDTDP2
#undef D3CR_GDP3

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conRhm defines the conversion from m[i], to r[j]
 */
/* Order: Xil, Xgk, Xpy, Xcn */
#define FR0(i)     (i == 2) ? 1.0 - r[0] : - r[0]
#define FR1(i)     (i == 0) ? 1.0 - r[1] : - r[1]
#define FR2(i)     (i == 3) ? 1.0 - r[2] : - r[2]
#define FR3(i)     (i == 4) ? 1.0 - r[3] : - r[3]

/* Order: s, t, u */
#define GS0(i)     (i == 2) ? 1.0 - s[0] : - s[0]
#define GS1(i)     (i == 0) ? 1.0 - s[1] : - s[1]
#define GS2(i)     (i == 3) ? 1.0 - s[2] : - s[2]

#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0
#define DFR2DR2(i) - 1.0
#define DFR3DR3(i) - 1.0

#define DGS0DS0(i) - 1.0
#define DGS1DS1(i) - 1.0
#define DGS2DS2(i) - 1.0

#define ENDMEMBERS (r[1]*ends[0] + (1.0-r[0]-r[1]-r[2]-r[3])*ends[1] + r[0]*ends[2] + r[2]*ends[3] + r[3]*ends[4])

#define DENDDR0    (ends[2] - ends[1])
#define DENDDR1    (ends[0] - ends[1])
#define DENDDR2    (ends[3] - ends[1])
#define DENDDR3    (ends[4] - ends[1])

/*
 * Global (to this file): derivative definitions
 */

#define S  -R*(  xfe2a*log(xfe2a) + xmg2a*log(xmg2a) + xmn2a*log(xmn2a) + xti4a*log(xti4a) \
+ xfe2b*log(xfe2b) + xmg2b*log(xmg2b) + xmn2b*log(xmn2b) + xti4b*log(xti4b) \
- 2.0*xfe2ID*log(xfe2ID) - 2.0*xmg2ID*log(xmg2ID) - 2.0*xmn2ID*log(xmn2ID) - 2.0*xti4ID*log(xti4ID) \
) \
-(1.0-SROconst)*2.0*R*(  xfe3ID*log(xfe3ID) + xal3ID*log(xal3ID) + xfe2ID*log(xfe2ID) \
+ xmg2ID*log(xmg2ID) + xmn2ID*log(xmn2ID) + xti4ID*log(xti4ID) ) + SROconst*2.0*R*log(2.0)

#define H    ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
+ ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
+ ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
+  whmcrn*r[3]*(1.0-r[0]-r[1]-r[2]-r[3]) \
+ (whilmgei+(p-1.0)*wvilmgei+whilmgeiT)*r[0]*r[1]/2.0 + (whilmgei-whilmgeiT)*s[0]*s[1]/2.0 \
+ (whilmpyr+(p-1.0)*wvilmpyr+whilmpyrT)*r[0]*r[2]/2.0 + (whilmpyr-whilmpyrT)*s[0]*s[2]/2.0 \
+ (whgeipyr+(p-1.0)*wvgeipyr+whgeipyrT)*r[1]*r[2]/2.0 + (whgeipyr-whgeipyrT)*s[1]*s[2]/2.0 \
+ (wcrnilm+(dwhcrnilm+(p-1.0)*dwvcrnilm)*(r[0]-r[3]))*r[0]*r[3] \
+ (wcrngei+(dwhcrngei+(p-1.0)*dwvcrngei)*(r[1]-r[3]))*r[1]*r[3] \
+ (wcrnpyr+(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*(r[2]-r[3]))*r[2]*r[3] \
+ (dhilm+(p-1.0)*dvilm+((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
-(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[3]/2.0)*(r[0]*r[0]-s[0]*s[0]) \
+ (dhgei+(p-1.0)*dvgei+((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
-(dwhcrngei+(p-1.0)*dwvcrngei)*r[3]/2.0)*(r[1]*r[1]-s[1]*s[1]) \
+ (dhpyr+(p-1.0)*dvpyr+((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
-(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[3]/2.0)*(r[2]*r[2]-s[2]*s[2]) \
+ (whilm+(p-1.0)*wvilm)*s[0]*s[0]*(r[0]*r[0]-s[0]*s[0]) \
+ (whgei+(p-1.0)*wvgei)*s[1]*s[1]*(r[1]*r[1]-s[1]*s[1]) \
+ (whpyr+(p-1.0)*wvpyr)*s[2]*s[2]*(r[2]*r[2]-s[2]*s[2]) \
+ (whilmilmgei)*(r[0]*r[0]-s[0]*s[0])*r[1]/4.0 + (whilmgeigei)*(r[1]*r[1]-s[1]*s[1])*r[0]/4.0 \
+ (whilmilmpyr)*(r[0]*r[0]-s[0]*s[0])*r[2]/4.0 + (whilmpyrpyr)*(r[2]*r[2]-s[2]*s[2])*r[0]/4.0 \
+ (whgeigeipyr)*(r[1]*r[1]-s[1]*s[1])*r[2]/4.0 + (whgeipyrpyr)*(r[2]*r[2]-s[2]*s[2])*r[1]/4.0

#define G  H - t*(S)

#define DGDR0  R*t*(  0.5*log(xfe2a) + 0.5*log(xti4a) + 0.5*log(xfe2b) + 0.5*log(xti4b) - log(xfe2ID) - log(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*(- log(xfe3ID) + 0.5*log(xfe2ID) + 0.5*log(xti4ID) ) \
+ ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
- ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1] \
- ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2] \
- 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
-     (dwhhmgei+(p-1.0)*dwvhmgei)*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
-     (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
-  whmcrn*r[3] \
+ (whilmgei+(p-1.0)*wvilmgei+whilmgeiT)*r[1]/2.0 \
+ (whilmpyr+(p-1.0)*wvilmpyr+whilmpyrT)*r[2]/2.0 \
+ (wcrnilm+(dwhcrnilm+(p-1.0)*dwvcrnilm)*(r[0]-r[3]))*r[3] + (dwhcrnilm+(p-1.0)*dwvcrnilm)*r[0]*r[3] \
+ 2.0*(dhilm+(p-1.0)*dvilm+((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
-(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[3]/2.0)*r[0] \
- ((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(r[0]*r[0]-s[0]*s[0]) \
- ((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(r[1]*r[1]-s[1]*s[1]) \
- ((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(r[2]*r[2]-s[2]*s[2]) \
+ 2.0*(whilm+(p-1.0)*wvilm)*s[0]*s[0]*r[0] \
+ (whilmilmgei)*r[0]*r[1]/2.0 + (whilmgeigei)*(r[1]*r[1]-s[1]*s[1])/4.0 \
+ (whilmilmpyr)*r[0]*r[2]/2.0 + (whilmpyrpyr)*(r[2]*r[2]-s[2]*s[2])/4.0

#define DGDR1  R*t*(  0.5*log(xmg2a) + 0.5*log(xti4a) + 0.5*log(xmg2b) + 0.5*log(xti4b) - log(xmg2ID) - log(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*(- log(xfe3ID) + 0.5*log(xmg2ID) + 0.5*log(xti4ID) ) \
- ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0] \
+ ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
- ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2] \
-     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
- 2.0*(dwhhmgei+(p-1.0)*dwvhmgei)*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
-     (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
-  whmcrn*r[3] \
+ (whilmgei+(p-1.0)*wvilmgei+whilmgeiT)*r[0]/2.0 \
+ (whgeipyr+(p-1.0)*wvgeipyr+whgeipyrT)*r[2]/2.0 \
+ (wcrngei+(dwhcrngei+(p-1.0)*dwvcrngei)*(r[1]-r[3]))*r[3] +(dwhcrngei+(p-1.0)*dwvcrngei)*r[1]*r[3] \
- ((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(r[0]*r[0]-s[0]*s[0]) \
+ 2.0*(dhgei+(p-1.0)*dvgei+((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
-(dwhcrngei+(p-1.0)*dwvcrngei)*r[3]/2.0)*r[1] \
- ((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(r[1]*r[1]-s[1]*s[1]) \
- ((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(r[2]*r[2]-s[2]*s[2]) \
+ 2.0*(whgei+(p-1.0)*wvgei)*s[1]*s[1]*r[1] \
+ (whilmilmgei)*(r[0]*r[0]-s[0]*s[0])/4.0 + (whilmgeigei)*r[1]*r[0]/2.0 \
+ (whgeigeipyr)*r[1]*r[2]/2.0 + (whgeipyrpyr)*(r[2]*r[2]-s[2]*s[2])/4.0

#define DGDR2  R*t*(  0.5*log(xmn2a) + 0.5*log(xti4a) + 0.5*log(xmn2b) + 0.5*log(xti4b) - log(xmn2ID) - log(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*(- log(xfe3ID) + 0.5*log(xmn2ID) + 0.5*log(xti4ID) ) \
- ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0] \
- ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1] \
+ ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
-     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
-     (dwhhmgei+(p-1.0)*dwvhmgei)*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
- 2.0*(dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
-  whmcrn*r[3] \
+ (whilmpyr+(p-1.0)*wvilmpyr+whilmpyrT)*r[0]/2.0 \
+ (whgeipyr+(p-1.0)*wvgeipyr+whgeipyrT)*r[1]/2.0 \
+ (wcrnpyr+(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*(r[2]-r[3]))*r[3] + (dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[2]*r[3] \
- ((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(r[0]*r[0]-s[0]*s[0]) \
- ((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(r[1]*r[1]-s[1]*s[1]) \
+ 2.0*(dhpyr+(p-1.0)*dvpyr+((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
-(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[3]/2.0)*r[2] \
- ((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(r[2]*r[2]-s[2]*s[2]) \
+ 2.0*(whpyr+(p-1.0)*wvpyr)*s[2]*s[2]*r[2] \
+ (whilmilmpyr)*(r[0]*r[0]-s[0]*s[0])/4.0 + (whilmpyrpyr)*r[2]*r[0]/2.0 \
+ (whgeigeipyr)*(r[1]*r[1]-s[1]*s[1])/4.0 + (whgeipyrpyr)*r[2]*r[1]/2.0

#define DGDR3  (1.0-SROconst)*2.0*R*t*(- log(xfe3ID) + log(xal3ID) ) \
- ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0] \
- ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1] \
- ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2] \
- (dwhhmilm+(p-1.0)*dwvhmilm)*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
- (dwhhmgei+(p-1.0)*dwvhmgei)*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
- (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
+  whmcrn*(1.0-r[0]-r[1]-r[2]-2.0*r[3]) \
+ (wcrnilm+(dwhcrnilm+(p-1.0)*dwvcrnilm)*(r[0]-r[3]))*r[0] - (dwhcrnilm+(p-1.0)*dwvcrnilm)*r[0]*r[3] \
+ (wcrngei+(dwhcrngei+(p-1.0)*dwvcrngei)*(r[1]-r[3]))*r[1] - (dwhcrngei+(p-1.0)*dwvcrngei)*r[1]*r[3] \
+ (wcrnpyr+(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*(r[2]-r[3]))*r[2] - (dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[2]*r[3] \
- ((dwhhmilm+(p-1.0)*dwvhmilm)/2.0 + (whhmilm2+(p-1.0)*wvhmilm2)/4.0 + (dwhcrnilm+(p-1.0)*dwvcrnilm)/2.0)*(r[0]*r[0]-s[0]*s[0]) \
- ((dwhhmgei+(p-1.0)*dwvhmgei)/2.0 + (whhmgei2+(p-1.0)*wvhmgei2)/4.0 + (dwhcrngei+(p-1.0)*dwvcrngei)/2.0)*(r[1]*r[1]-s[1]*s[1]) \
- ((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0 + (whhmpyr2+(p-1.0)*wvhmpyr2)/4.0 + (dwhcrnpyr+(p-1.0)*dwvcrnpyr)/2.0)*(r[2]*r[2]-s[2]*s[2])

#define DGDS0  0.5*R*t*( log(xfe2a) - log(xti4a) - log(xfe2b) + log(xti4b) ) \
+ (whilmgei-whilmgeiT)*s[1]/2.0 + (whilmpyr-whilmpyrT)*s[2]/2.0 \
- 2.0*(dhilm+(p-1.0)*dvilm+((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
-(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[3]/2.0)*s[0] \
+ 2.0*(whilm+(p-1.0)*wvilm)*s[0]*(r[0]*r[0]-2.0*s[0]*s[0]) \
- (whilmilmgei)*s[0]*r[1]/2.0 - (whilmilmpyr)*s[0]*r[2]/2.0

#define DGDS1  0.5*R*t*( log(xmg2a) - log(xti4a) - log(xmg2b) + log(xti4b) ) \
+ (whilmgei-whilmgeiT)*s[0]/2.0 + (whgeipyr-whgeipyrT)*s[2]/2.0 \
- 2.0*(dhgei+(p-1.0)*dvgei+((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
-(dwhcrngei+(p-1.0)*dwvcrngei)*r[3]/2.0)*s[1] \
+ 2.0*(whgei+(p-1.0)*wvgei)*s[1]*(r[1]*r[1]-2.0*s[1]*s[1]) \
- (whilmgeigei)*s[1]*r[0]/2.0 - (whgeigeipyr)*s[1]*r[2]/2.0

#define DGDS2  0.5*R*t*( log(xmn2a) - log(xti4a) - log(xmn2b) + log(xti4b) ) \
+ (whilmpyr-whilmpyrT)*s[0]/2.0 + (whgeipyr-whgeipyrT)*s[1]/2.0 \
- 2.0*(dhpyr+(p-1.0)*dvpyr+((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
-(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[3]/2.0)*s[2] \
+ 2.0*(whpyr+(p-1.0)*wvpyr)*s[2]*(r[2]*r[2]-2.0*s[2]*s[2]) \
- (whilmpyrpyr)*s[2]*r[0]/2.0 - (whgeipyrpyr)*s[2]*r[1]/2.0

#define DGDT   - (S)

#define DGDP     (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
+ (wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
+ (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
+  wvilmgei*r[0]*r[1] \
+  wvilmpyr*r[0]*r[2] \
+  wvgeipyr*r[1]*r[2] \
+  dwvcrnilm*(r[0]-r[3])*r[0]*r[3] \
+  dwvcrngei*(r[1]-r[3])*r[1]*r[3] \
+  dwvcrnpyr*(r[2]-r[3])*r[2]*r[3] \
+ (dvilm + (dwvhmilm/2.0+wvhmilm2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) - dwvcrnilm*r[3]/2.0)*(r[0]*r[0]-s[0]*s[0]) \
+ (dvgei + (dwvhmgei/2.0+wvhmgei2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) - dwvcrngei*r[3]/2.0)*(r[1]*r[1]-s[1]*s[1]) \
+ (dvpyr + (dwvhmpyr/2.0+wvhmpyr2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) - dwvcrnpyr*r[3]/2.0)*(r[2]*r[2]-s[2]*s[2]) \
+  wvilm*s[0]*s[0]*(r[0]*r[0]-s[0]*s[0]) \
+  wvgei*s[1]*s[1]*(r[1]*r[1]-s[1]*s[1]) \
+  wvpyr*s[2]*s[2]*(r[2]*r[2]-s[2]*s[2])

#define D2GDR0R0  R*t*(  0.25/xfe2a + 0.25/xti4a + 0.25/xfe2b + 0.25/xti4b - 0.5/xfe2ID - 0.5/xti4ID ) \
+ (1.0-SROconst)*2.0*R*t*( 1.0/xfe3ID + 0.25/xfe2ID + 0.25/xti4ID ) \
- 2.0*((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
- 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
+ (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
+ (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
- 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
+	(dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
+	(dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
+ 2.0*(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[3] \
+ 2.0*(dhilm+(p-1.0)*dvilm+((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
-(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[3]/2.0) \
- 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*r[0] \
+ 2.0*(whilm+(p-1.0)*wvilm)*s[0]*s[0] \
+ (whilmilmgei)*r[1]/2.0 + (whilmilmpyr)*r[2]/2.0

#define D2GDR0R1  R*t*( 0.25/xti4a + 0.25/xti4b - 0.5/xti4ID ) \
+ (1.0-SROconst)*2.0*R*t*( 1.0/xfe3ID + 0.25/xti4ID ) \
- ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
- (dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
- ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-4.0*r[1]-r[2]-r[3])) \
+ (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
+ 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
-	(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
+	(dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
+ (whilmgei+(p-1.0)*wvilmgei+whilmgeiT)/2.0 \
- 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*r[0] \
- 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*r[1] \
+ (whilmilmgei)*r[0]/2.0 + (whilmgeigei)*r[1]/2.0

#define D2GDR0R2  R*t*( 0.25/xti4a + 0.25/xti4b - 0.5/xti4ID ) \
+ (1.0-SROconst)*2.0*R*t*( 1.0/xfe3ID + 0.25/xti4ID ) \
- ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
- (dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
+ (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
- ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-4.0*r[2]-r[3])) \
+ 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
+	(dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
-	(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
+ (whilmpyr+(p-1.0)*wvilmpyr+whilmpyrT)/2.0 \
- 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*r[0] \
- 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*r[2] \
+ (whilmilmpyr)*r[0]/2.0 + (whilmpyrpyr)*r[2]/2.0

#define D2GDR0R3  (1.0-SROconst)*2.0*R*t*( 1.0/xfe3ID ) \
- ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
- (dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
+ (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
+ (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
+ 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
+	(dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
+	(dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
-  whmcrn \
+ (wcrnilm+(dwhcrnilm+(p-1.0)*dwvcrnilm)*(r[0]-2.0*r[3])) + (dwhcrnilm+(p-1.0)*dwvcrnilm)*r[0] \
- 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0+(dwhcrnilm+(p-1.0)*dwvcrnilm)/2.0)*r[0]

#define D2GDR0S0  0.25*R*t*( 1.0/xfe2a - 1.0/xti4a - 1.0/xfe2b + 1.0/xti4b ) \
+ 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*s[0] \
+ 4.0*(whilm+(p-1.0)*wvilm)*s[0]*r[0]

#define D2GDR0S1  0.25*R*t*( - 1.0/xti4a + 1.0/xti4b ) \
+ 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*s[1] \
- (whilmgeigei)*s[1]/2.0

#define D2GDR0S2  0.25*R*t*( - 1.0/xti4a + 1.0/xti4b ) \
+ 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*s[2] \
- (whilmpyrpyr)*s[2]/2.0

#define D2GDR0DT  R*(  0.5*log(xfe2a) + 0.5*log(xti4a) + 0.5*log(xfe2b) + 0.5*log(xti4b) - log(xfe2ID) - log(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*(- log(xfe3ID) + 0.5*log(xfe2ID) + 0.5*log(xti4ID) )

#define D2GDR0DP (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
- (wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1] \
- (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2] \
- 2.0*dwvhmilm*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
-     dwvhmgei*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
-     dwvhmpyr*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
+ wvilmgei*r[1] \
+ wvilmpyr*r[2] \
+ dwvcrnilm*(r[0]-r[3])*r[3] + dwvcrnilm*r[0]*r[3] \
+ 2.0*(dvilm+(dwvhmilm/2.0+wvhmilm2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrnilm*r[3]/2.0)*r[0] \
- (dwvhmilm/2.0+wvhmilm2/4.0)*(r[0]*r[0]-s[0]*s[0]) \
- (dwvhmgei/2.0+wvhmgei2/4.0)*(r[1]*r[1]-s[1]*s[1]) \
- (dwvhmpyr/2.0+wvhmpyr2/4.0)*(r[2]*r[2]-s[2]*s[2]) \
+ 2.0*wvilm*s[0]*s[0]*r[0]

#define D2GDR1R1  R*t*(  0.25/xmg2a + 0.25/xti4a + 0.25/xmg2b + 0.25/xti4b - 0.5/xmg2ID - 0.5/xti4ID ) \
+ (1.0-SROconst)*2.0*R*t*( 1.0/xfe3ID + 0.25/xmg2ID + 0.25/xti4ID ) \
+ (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
- 2.0*((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3])) \
- 2.0*(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
+ (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
+     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
- 2.0*(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
+     (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
+ 2.0*(dwhcrngei+(p-1.0)*dwvcrngei)*r[3] \
+ 2.0*(dhgei+(p-1.0)*dvgei+((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
-(dwhcrngei+(p-1.0)*dwvcrngei)*r[3]/2.0) \
- 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*r[1] \
+ 2.0*(whgei+(p-1.0)*wvgei)*s[1]*s[1] \
+ (whilmgeigei)*r[0]/2.0 + (whgeigeipyr)*r[2]/2.0

#define D2GDR1R2  R*t*(  0.25/xti4a + 0.25/xti4b - 0.5/xti4ID ) \
+ (1.0-SROconst)*2.0*R*t*( 1.0/xfe3ID + 0.25/xti4ID ) \
+ (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
- ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3])) \
- (dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
- ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-4.0*r[2]-r[3])) \
+     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
+ 2.0*(dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
-     (dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
+ (whgeipyr+(p-1.0)*wvgeipyr+whgeipyrT)/2.0 \
- 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*r[1] \
- 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*r[2] \
+ (whgeigeipyr)*r[1]/2.0 + (whgeipyrpyr)*r[2]/2.0

#define D2GDR1R3  (1.0-SROconst)*2.0*R*t*( 1.0/xfe3ID ) \
+ (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
- ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3])) \
- (dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
+ (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
+     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
+ 2.0*(dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
+     (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
-  whmcrn \
+ (wcrngei+(dwhcrngei+(p-1.0)*dwvcrngei)*(r[1]-2.0*r[3])) +(dwhcrngei+(p-1.0)*dwvcrngei)*r[1] \
- 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0+(dwhcrngei+(p-1.0)*dwvcrngei)/2.0)*r[1]

#define D2GDR1S0  0.25*R*t*( - 1.0/xti4a + 1.0/xti4b ) \
+ 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*s[0] \
- (whilmilmgei)*s[0]/2.0

#define D2GDR1S1  0.25*R*t*( 1.0/xmg2a - 1.0/xti4a - 1.0/xmg2b + 1.0/xti4b ) \
+ 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*s[1] \
+ 4.0*(whgei+(p-1.0)*wvgei)*s[1]*r[1]

#define D2GDR1S2  0.25*R*t*( - 1.0/xti4a + 1.0/xti4b ) \
+ 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*s[2] \
- (whgeipyrpyr)*s[2]/2.0

#define D2GDR1DT  R*(  0.5*log(xmg2a) + 0.5*log(xti4a) + 0.5*log(xmg2b) + 0.5*log(xti4b) - log(xmg2ID) - log(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*(- log(xfe3ID) + 0.5*log(xmg2ID) + 0.5*log(xti4ID) )

#define D2GDR1DP - (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0] \
+ (wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
- (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2] \
-     dwvhmilm*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
- 2.0*dwvhmgei*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
-     dwvhmpyr*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
+ wvilmgei*r[0] \
+ wvgeipyr*r[2] \
+ dwvcrngei*(r[1]-r[3])*r[3] +dwvcrngei*r[1]*r[3] \
- (dwvhmilm/2.0+wvhmilm2/4.0)*(r[0]*r[0]-s[0]*s[0]) \
+ 2.0*(dvgei+(dwvhmgei/2.0+wvhmgei2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrngei*r[3]/2.0)*r[1] \
- (dwvhmgei/2.0+wvhmgei2/4.0)*(r[1]*r[1]-s[1]*s[1]) \
- (dwvhmpyr/2.0+wvhmpyr2/4.0)*(r[2]*r[2]-s[2]*s[2]) \
+ 2.0*wvgei*s[1]*s[1]*r[1]

#define D2GDR2R2  R*t*(  0.25/xmn2a + 0.25/xti4a + 0.25/xmn2b + 0.25/xti4b - 0.5/xmn2ID - 0.5/xti4ID ) \
+ (1.0-SROconst)*2.0*R*t*( 1.0/xfe3ID + 0.25/xmn2ID + 0.25/xti4ID ) \
+ (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
+ (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
- 2.0*((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3])) \
- 2.0*(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
+     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
+     (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
- 2.0*(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
+ 2.0*(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[3] \
+ 2.0*(dhpyr+(p-1.0)*dvpyr+((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
-(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[3]/2.0) \
- 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*r[2] \
+ 2.0*(whpyr+(p-1.0)*wvpyr)*s[2]*s[2] \
+ (whilmpyrpyr)*r[0]/2.0 + (whgeipyrpyr)*r[1]/2.0

#define D2GDR2R3  (1.0-SROconst)*2.0*R*t*( 1.0/xfe3ID ) \
+ (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
+ (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
- ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3])) \
- (dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
+     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
+     (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
+ 2.0*(dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
-  whmcrn \
+ (wcrnpyr+(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*(r[2]-2.0*r[3])) + (dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[2] \
- 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0+(dwhcrnpyr+(p-1.0)*dwvcrnpyr)/2.0)*r[2]

#define D2GDR2S0  0.25*R*t*( - 1.0/xti4a + 1.0/xti4b ) \
+ 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*s[0] \
- (whilmilmpyr)*s[0]/2.0

#define D2GDR2S1  0.25*R*t*( - 1.0/xti4a + 1.0/xti4b ) \
+ 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*s[1] \
- (whgeigeipyr)*s[1]/2.0

#define D2GDR2S2  0.25*R*t*( 1.0/xmn2a - 1.0/xti4a - 1.0/xmn2b + 1.0/xti4b ) \
+ 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*s[2] \
+ 4.0*(whpyr+(p-1.0)*wvpyr)*s[2]*r[2]

#define D2GDR2DT  R*(  0.5*log(xmn2a) + 0.5*log(xti4a) + 0.5*log(xmn2b) + 0.5*log(xti4b) - log(xmn2ID) - log(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*(- log(xfe3ID) + 0.5*log(xmn2ID) + 0.5*log(xti4ID) )

#define D2GDR2DP  - (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0] \
- (wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1] \
+ (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
-     dwvhmilm*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
-     dwvhmgei*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
- 2.0*dwvhmpyr*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
+ wvilmpyr*r[0] \
+ wvgeipyr*r[1] \
+ dwvcrnpyr*(r[2]-r[3])*r[3] + dwvcrnpyr*r[2]*r[3] \
- (dwvhmilm/2.0+wvhmilm2/4.0)*(r[0]*r[0]-s[0]*s[0]) \
- (dwvhmgei/2.0+wvhmgei2/4.0)*(r[1]*r[1]-s[1]*s[1]) \
+ 2.0*(dvpyr+(dwvhmpyr/2.0+wvhmgei2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrnpyr*r[3]/2.0)*r[2] \
- (dwvhmpyr/2.0+wvhmpyr2/4.0)*(r[2]*r[2]-s[2]*s[2]) \
+ 2.0*wvpyr*s[2]*s[2]*r[2]

#define D2GDR3R3  (1.0-SROconst)*2.0*R*t*( 1.0/xfe3ID + 1.0/xal3ID ) \
+ 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
+ 2.0*(dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
+ 2.0*(dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
- 2.0*whmcrn \
- 2.0*(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[0] \
- 2.0*(dwhcrngei+(p-1.0)*dwvcrngei)*r[1] \
- 2.0*(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[2]

#define D2GDR3S0 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0 + (whhmilm2+(p-1.0)*wvhmilm2)/4.0 + (dwhcrnilm+(p-1.0)*dwvcrnilm)/2.0)*s[0]

#define D2GDR3S1 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0 + (whhmgei2+(p-1.0)*wvhmgei2)/4.0 + (dwhcrngei+(p-1.0)*dwvcrngei)/2.0)*s[1]

#define D2GDR3S2 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0 + (whhmpyr2+(p-1.0)*wvhmpyr2)/4.0 + (dwhcrnpyr+(p-1.0)*dwvcrnpyr)/2.0)*s[2]

#define D2GDR3DT  (1.0-SROconst)*2.0*R*(- log(xfe3ID) + log(xal3ID) )

#define D2GDR3DP - (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0] \
- (wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1] \
- (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2] \
- dwvhmilm*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
- dwvhmgei*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
- dwvhmpyr*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
+ (dwvcrnilm*(r[0]-r[3]))*r[0] - dwvcrnilm*r[0]*r[3] \
+ (dwvcrngei*(r[1]-r[3]))*r[1] - dwvcrngei*r[1]*r[3] \
+ (dwvcrnpyr*(r[2]-r[3]))*r[2] - dwvcrnpyr*r[2]*r[3] \
- (dwvhmilm/2.0+wvhmilm2/4.0+dwvcrnilm/2.0)*(r[0]*r[0]-s[0]*s[0]) \
- (dwvhmgei/2.0+wvhmgei2/4.0+dwvcrngei/2.0)*(r[1]*r[1]-s[1]*s[1]) \
- (dwvhmpyr/2.0+wvhmpyr2/4.0+dwvcrnpyr/2.0)*(r[2]*r[2]-s[2]*s[2])

#define D2GDS0S0  0.25*R*t*( 1.0/xfe2a + 1.0/xti4a + 1.0/xfe2b + 1.0/xti4b ) \
- 2.0*(dhilm+(p-1.0)*dvilm+((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
-(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[3]/2.0) \
+ 2.0*(whilm+(p-1.0)*wvilm)*(r[0]*r[0]-6.0*s[0]*s[0]) \
- (whilmilmgei)*r[1]/2.0 - (whilmilmpyr)*r[2]/2.0

#define D2GDS0S1  0.25*R*t*( 1.0/xti4a + 1.0/xti4b ) + (whilmgei-whilmgeiT)/2.0

#define D2GDS0S2  0.25*R*t*( 1.0/xti4a + 1.0/xti4b ) + (whilmpyr-whilmpyrT)/2.0

#define D2GDS0DT  0.5*R*( log(xfe2a) - log(xti4a) - log(xfe2b) + log(xti4b) )

#define D2GDS0DP - 2.0*(dvilm+(dwvhmilm/2.0+wvhmilm2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrnilm*r[3]/2.0)*s[0] \
+ 2.0*wvilm*s[0]*(r[0]*r[0]-2.0*s[0]*s[0])

#define D2GDS1S1  0.25*R*t*( 1.0/xmg2a + 1.0/xti4a + 1.0/xmg2b + 1.0/xti4b ) \
- 2.0*(dhgei+(p-1.0)*dvgei+((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
-(dwhcrngei+(p-1.0)*dwvcrngei)*r[3]/2.0) \
+ 2.0*(whgei+(p-1.0)*wvgei)*(r[1]*r[1]-6.0*s[1]*s[1]) \
- (whilmgeigei)*r[0]/2.0 - (whgeigeipyr)*r[2]/2.0

#define D2GDS1S2  0.25*R*t*( 1.0/xti4a + 1.0/xti4b ) + (whgeipyr-whgeipyrT)/2.0

#define D2GDS1DT  0.5*R*( log(xmg2a) - log(xti4a) - log(xmg2b) + log(xti4b) )

#define D2GDS1DP - 2.0*(dvgei+(dwvhmgei/2.0+wvhmgei2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrngei*r[3]/2.0)*s[1] \
+ 2.0*wvgei*s[1]*(r[1]*r[1]-2.0*s[1]*s[1])

#define D2GDS2S2  0.25*R*t*( 1.0/xmn2a + 1.0/xti4a + 1.0/xmn2b + 1.0/xti4b ) \
- 2.0*(dhpyr+(p-1.0)*dvpyr+((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
-(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[3]/2.0) \
+ 2.0*(whpyr+(p-1.0)*wvpyr)*(r[2]*r[2]-6.0*s[2]*s[2]) \
- (whilmpyrpyr)*r[0]/2.0 - (whgeipyrpyr)*r[1]/2.0

#define D2GDS2DT  0.5*R*( log(xmn2a) - log(xti4a) - log(xmn2b) + log(xti4b) )

#define D2GDS2DP - 2.0*(dvpyr+(dwvhmpyr/2.0+wvhmpyr2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrnpyr*r[3]/2.0)*s[2] \
+ 2.0*wvpyr*s[2]*(r[2]*r[2]-2.0*s[2]*s[2])

#define D2GDT2    0.0
#define D2GDTDP   0.0
#define D2GDP2    0.0

/*----------------------------------------------------------------------------*/
/* Excess free energy modifications incomplete below this point               */
/*----------------------------------------------------------------------------*/


#define D3GDR0R0R0 R*t*(- 0.125/SQUARE(xfe2a) - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xfe2b) - 0.125/SQUARE(xti4b) \
+ 0.25/SQUARE(xfe2ID) + 0.25/SQUARE(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*(1.0/SQUARE(xfe3ID) - 0.125/SQUARE(xfe2ID) - 0.125/SQUARE(xti4ID) ) \
+ 12.0*(dwhhmilm+(p-1.0)*dwvhmilm) \
- 6.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)

#define D3GDR0R0R1 R*t*(  0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR0R0R2 R*t*(  0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR0R0R3 (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR0R0S0 R*t*( - 0.125/SQUARE(xfe2a) + 0.125/SQUARE(xti4a) + 0.125/SQUARE(xfe2b) - 0.125/SQUARE(xti4b) )  \
+ 4.0*(whilm+(p-1.0)*wvilm)*s[0]
#define D3GDR0R0S1 R*t*( 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR0R0S2 R*t*( 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR0R0DT R*(  0.25/xfe2a + 0.25/xti4a + 0.25/xfe2b + 0.25/xti4b - 0.5/xfe2ID - 0.5/xti4ID ) \
+ (1.0-SROconst)*2.0*R*( 1.0/xfe3ID + 0.25/xfe2ID + 0.25/xti4ID )
#define D3GDR0R0DP - 2.0*(wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
- 2.0*dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
+ dwvhmgei*r[1] \
+ dwvhmpyr*r[2] \
- 2.0*dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
+	dwvhmgei*r[1] \
+	dwvhmpyr*r[2] \
+ 2.0*dwvcrnilm*r[3] \
+ 2.0*(dvilm+(dwvhmilm/2.0+wvhmilm2/4.0)*(1.0-2.0*r[0]-r[1]-r[2]-r[3])-dwvcrnilm*r[3]/2.0) \
- 2.0*(dwvhmilm/2.0+wvhmilm2/4.0)*r[0] \
+ 2.0*wvilm*s[0]*s[0]
/*end-OK*/
#define D3GDR0R1R1 R*t*( 0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR0R1R2 R*t*( 0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR0R1R3 (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR0R1S0 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR0R1S1 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR0R1S2 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR0R1DT R*( 0.25/xti4a + 0.25/xti4b - 0.5/xti4ID ) \
+ (1.0-SROconst)*2.0*R*( 1.0/xfe3ID + 0.25/xti4ID )
#define D3GDR0R1DP - (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
- (dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
- (wvhmgei+dwvhmgei*(1.0-r[0]-4.0*r[1]-r[2]-r[3])) \
+ dwvhmpyr*r[2] \
+ 2.0*dwvhmilm*r[0] \
-	dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
+	dwvhmpyr*r[2] \
+ wvilmgei \
- 2.0*(dwvhmilm/2.0+wvhmilm2/4.0)*r[0] \
- 2.0*(dwvhmgei/2.0+wvhmgei2/4.0)*r[1]
/*end-OK*/
#define D3GDR0R2R2 R*t*( 0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR0R2R3 (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR0R2S0 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR0R2S1 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR0R2S2 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR0R2DT R*( 0.25/xti4a + 0.25/xti4b - 0.5/xti4ID ) \
+ (1.0-SROconst)*2.0*R*( 1.0/xfe3ID + 0.25/xti4ID )
#define D3GDR0R2DP - (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
- dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
+ dwvhmgei*r[1] \
- (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-4.0*r[2]-r[3])) \
+ 2.0*dwvhmilm*r[0] \
+	dwvhmgei*r[1] \
-	dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
+ wvilmpyr \
- 2.0*(dwvhmilm/2.0+wvhmilm2/4.0)*r[0] \
- 2.0*(dwvhmpyr/2.0+wvhmpyr2/4.0)*r[2]
/*end-OK*/
#define D3GDR0R3R3 (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR0R3S0 0.0
#define D3GDR0R3S1 0.0
#define D3GDR0R3S2 0.0
#define D3GDR0R3DT (1.0-SROconst)*2.0*R*( 1.0/xfe3ID )
#define D3GDR0R3DP - (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
- dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
+ dwvhmgei*r[1] \
+ dwvhmpyr*r[2] \
+ 2.0*dwvhmilm*r[0] \
+	dwvhmgei*r[1] \
+	dwvhmpyr*r[2] \
+ (dwvcrnilm*(r[0]-2.0*r[3])) + dwvcrnilm*r[0] \
- 2.0*(dwvhmilm/2.0+wvhmilm2/4.0+dwvcrnilm/2.0)*r[0]

#define D3GDR0S0S0 R*t*( - 0.125/SQUARE(xfe2a) - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xfe2b) - 0.125/SQUARE(xti4b) ) \
+ 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0) + 4.0*(whilm+(p-1.0)*wvilm)*r[0]
#define D3GDR0S0S1 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR0S0S2 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR0S0DT R*( 0.25/xfe2a - 0.25/xti4a - 0.25/xfe2b + 0.25/xti4b )
#define D3GDR0S0DP 2.0*(dwvhmilm/2.0+wvhmilm2/4.0)*s[0] + 4.0*wvilm*s[0]*r[0]

#define D3GDR0S1S1 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) ) \
+ 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0) - (whilmgeigei)/2.0
#define D3GDR0S1S2 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR0S1DT R*( - 0.25/xti4a + 0.25/xti4b )
#define D3GDR0S1DP 2.0*(dwvhmgei/2.0+wvhmgei2/4.0)*s[1]

#define D3GDR0S2S2 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) ) \
+ 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0) - (whilmpyrpyr)/2.0
#define D3GDR0S2DT R*( - 0.25/xti4a + 0.25/xti4b )
#define D3GDR0S2DP 2.0*(dwvhmpyr/2.0+wvhmpyr2/4.0)*s[2]

#define D3GDR0DT2  0.0
#define D3GDR0DTDP 0.0
#define D3GDR0DP2  0.0
/*end-OK*/
#define D3GDR1R1R1 R*t*(  0.125/SQUARE(xmg2a) + 0.125/SQUARE(xti4a) + 0.125/SQUARE(xmg2b) + 0.125/SQUARE(xti4b) \
- 0.25/SQUARE(xmg2ID) - 0.25/SQUARE(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xmg2ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR1R1R2 R*t*(  0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID)  + 0.125/SQUARE(xti4ID) )
#define D3GDR1R1R3 (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR1R1S0 R*t*(   0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR1R1S1 R*t*( - 0.125/SQUARE(xmg2a) + 0.125/SQUARE(xti4a) + 0.125/SQUARE(xmg2b) - 0.125/SQUARE(xti4b) ) \
+ 4.0*(whgei+(p-1.0)*wvgei)*s[1]
#define D3GDR1R1S2 R*t*(   0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR1R1DT R*(  0.25/xmg2a + 0.25/xti4a + 0.25/xmg2b + 0.25/xti4b - 0.5/xmg2ID - 0.5/xti4ID ) \
+ (1.0-SROconst)*2.0*R*( 1.0/xfe3ID + 0.25/xmg2ID + 0.25/xti4ID )
#define D3GDR1R1DP dwvhmilm*r[0] \
- 2.0*(wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3])) \
- 2.0*dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
+ dwvhmpyr*r[2] \
+ dwvhmilm*r[0] \
- 2.0*dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
+ dwvhmpyr*r[2] \
+ 2.0*dwvcrngei*r[3] \
+ 2.0*(dvgei+(dwvhmgei/2.0+wvhmgei2/4.0)*(1.0-r[0]-2.0*r[1]-r[2]-r[3])-dwvcrngei*r[3]/2.0) \
- 2.0*(dwvhmgei/2.0+wvhmgei2/4.0)*r[1] \
+ 2.0*wvgei*s[1]*s[1]
/*end-OK*/
#define D3GDR1R2R2 R*t*(  0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR1R2R3 (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR1R2S0 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR1R2S1 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR1R2S2 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR1R2DT R*(  0.25/xti4a + 0.25/xti4b - 0.5/xti4ID ) \
+ (1.0-SROconst)*2.0*R*( 1.0/xfe3ID + 0.25/xti4ID )
#define D3GDR1R2DP (dwvhmilm)*r[0] \
- ((wvhmgei)+(dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3])) \
- (dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
- ((wvhmpyr)+(dwvhmpyr)*(1.0-r[0]-r[1]-4.0*r[2]-r[3])) \
+     (dwvhmilm)*r[0] \
+ 2.0*(dwvhmgei)*r[1] \
-     (dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
+ (wvgeipyr) \
- 2.0*((dwvhmgei)/2.0+(wvhmgei2)/4.0)*r[1] \
- 2.0*((dwvhmpyr)/2.0+(wvhmpyr2)/4.0)*r[2]
/*end-OK*/
#define D3GDR1R3R3 (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR1R3S0 0.0
#define D3GDR1R3S1 0.0
#define D3GDR1R3S2 0.0
#define D3GDR1R3DT (1.0-SROconst)*2.0*R*( 1.0/xfe3ID )
#define D3GDR1R3DP dwvhmilm*r[0] \
- (wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3])) \
- dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
+ dwvhmpyr*r[2] \
+     dwvhmilm*r[0] \
+ 2.0*dwvhmgei*r[1] \
+     dwvhmpyr*r[2] \
+ dwvcrngei*(r[1]-2.0*r[3])+dwvcrngei*r[1] \
- 2.0*(dwvhmgei/2.0 + wvhmgei2/4.0 + dwvcrngei/2.0)*r[1]

#define D3GDR1S0S0 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) ) \
+ 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0) - (whilmilmgei)/2.0
#define D3GDR1S0S1 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR1S0S2 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR1S0DT 0.25*R*( - 1.0/xti4a + 1.0/xti4b )
#define D3GDR1S0DP 2.0*(dwvhmilm/2.0+wvhmilm2/4.0)*s[0]

#define D3GDR1S1S1 R*t*( - 0.125/SQUARE(xmg2a) - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xmg2b) - 0.125/SQUARE(xti4b) ) \
+ 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0) + 4.0*(whgei+(p-1.0)*wvgei)*r[1]
#define D3GDR1S1S2 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR1S1DT R*( 0.25/xmg2a - 0.25/xti4a - 0.25/xmg2b + 0.25/xti4b )
#define D3GDR1S1DP 2.0*(dwvhmgei/2.0+wvhmgei2/4.0)*s[1] + 4.0*wvgei*s[1]*r[1]

#define D3GDR1S2S2 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) ) \
+ 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0) - (whgeipyrpyr)/2.0
#define D3GDR1S2DT 0.25*R*( - 1.0/xti4a + 1.0/xti4b )
#define D3GDR1S2DP 2.0*(dwvhmpyr/2.0+wvhmpyr2/4.0)*s[2]

#define D3GDR1DT2  0.0
#define D3GDR1DTDP 0.0
#define D3GDR1DP2  0.0
/*end-OK*/
#define D3GDR2R2R2 R*t*(  0.125/SQUARE(xmn2a) + 0.125/SQUARE(xti4a) + 0.125/SQUARE(xmn2b) + 0.125/SQUARE(xti4b) \
- 0.25/SQUARE(xmn2ID) - 0.25/SQUARE(xti4ID) ) \
+ (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xmn2ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR2R2R3 (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR2R2S0 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR2R2S1 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR2R2S2 - 0.125*R*t*( 1.0/SQUARE(xmn2a) - 1.0/SQUARE(xti4a) - 1.0/SQUARE(xmn2b) + 1.0/SQUARE(xti4b) ) \
+ 4.0*(whpyr+(p-1.0)*wvpyr)*s[2]
#define D3GDR2R2DT R*(  0.25/xmn2a + 0.25/xti4a + 0.25/xmn2b + 0.25/xti4b - 0.5/xmn2ID - 0.5/xti4ID ) \
+ (1.0-SROconst)*2.0*R*( 1.0/xfe3ID + 0.25/xmn2ID + 0.25/xti4ID )
#define D3GDR2R2DP dwvhmilm*r[0] \
+ dwvhmgei*r[1] \
- 2.0*(wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3])) \
- 2.0*dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
+     dwvhmilm*r[0] \
+     dwvhmgei*r[1] \
- 2.0*dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
+ 2.0*dwvcrnpyr*r[3] \
+ 2.0*(dvpyr+(dwvhmpyr/2.0+wvhmpyr2/4.0)*(1.0-r[0]-r[1]-2.0*r[2]-r[3])-dwvcrnpyr*r[3]/2.0) \
- 2.0*(dwvhmpyr/2.0+wvhmpyr2/4.0)*r[2] \
+ 2.0*wvpyr*s[2]*s[2]
/*end-OK*/
#define D3GDR2R3R3 (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR2R3S0 0.0
#define D3GDR2R3S1 0.0
#define D3GDR2R3S2 0.0
#define D3GDR2R3DT (1.0-SROconst)*2.0*R*( 1.0/xfe3ID )
#define D3GDR2R3DP dwvhmilm*r[0] \
+ dwvhmgei*r[1] \
- (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3])) \
- dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
+     dwvhmilm*r[0] \
+     dwvhmgei*r[1] \
+ 2.0*dwvhmpyr*r[2] \
+ dwvcrnpyr*(r[2]-2.0*r[3]) + dwvcrnpyr*r[2] \
- 2.0*(dwvhmpyr/2.0+wvhmpyr2/4.0+dwvcrnpyr/2.0)*r[2]

#define D3GDR2S0S0 - 0.125*R*t*( 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) ) \
+ 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0) - (whilmilmpyr)/2.0
#define D3GDR2S0S1 - 0.125*R*t*( 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR2S0S2 - 0.125*R*t*( 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR2S0DT R*( - 0.25/xti4a + 0.25/xti4b )
#define D3GDR2S0DP 2.0*(dwvhmilm/2.0+wvhmilm2/4.0)*s[0]

#define D3GDR2S1S1 - 0.125*R*t*( 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) ) \
+ 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0) - (whgeigeipyr)/2.0
#define D3GDR2S1S2 - 0.125*R*t*( 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR2S1DT 0.25*R*( - 1.0/xti4a + 1.0/xti4b )
#define D3GDR2S1DP 2.0*(dwvhmgei/2.0+wvhmgei2/4.0)*s[1]

#define D3GDR2S2S2 R*t*( - 0.125/SQUARE(xmn2a) - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xmn2b) - 0.125/SQUARE(xti4b) ) \
+ 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0) + 4.0*(whpyr+(p-1.0)*wvpyr)*r[2]
#define D3GDR2S2DT 0.25*R*( 1.0/xmn2a - 1.0/xti4a - 1.0/xmn2b + 1.0/xti4b )
#define D3GDR2S2DP 2.0*(dwvhmpyr/2.0+wvhmpyr2/4.0)*s[2] + 4.0*wvpyr*s[2]*r[2]

#define D3GDR2DT2  0.0
#define D3GDR2DTDP 0.0
#define D3GDR2DP2  0.0

#define D3GDR3R3R3 (1.0-SROconst)*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 1.0/SQUARE(xal3ID) )
#define D3GDR3R3S0 0.0
#define D3GDR3R3S1 0.0
#define D3GDR3R3S2 0.0
#define D3GDR3R3DT (1.0-SROconst)*2.0*R*( 1.0/xfe3ID + 1.0/xal3ID )
#define D3GDR3R3DP   2.0*dwvhmilm*r[0] \
+ 2.0*dwvhmgei*r[1] \
+ 2.0*dwvhmpyr*r[2] \
- 2.0*dwvcrnilm*r[0] \
- 2.0*dwvcrngei*r[1] \
- 2.0*dwvcrnpyr*r[2]

#define D3GDR3S0S0 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0 + (whhmilm2+(p-1.0)*wvhmilm2)/4.0 + (dwhcrnilm+(p-1.0)*dwvcrnilm)/2.0)
#define D3GDR3S0S1 0.0
#define D3GDR3S0S2 0.0
#define D3GDR3S0DT 0.0
#define D3GDR3S0DP 2.0*(dwvhmilm/2.0 + wvhmilm2/4.0 + dwvcrnilm/2.0)*s[0]

#define D3GDR3S1S1 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0 + (whhmgei2+(p-1.0)*wvhmgei2)/4.0 + (dwhcrngei+(p-1.0)*dwvcrngei)/2.0)
#define D3GDR3S1S2 0.0
#define D3GDR3S1DT 0.0
#define D3GDR3S1DP 2.0*(dwvhmgei/2.0 + wvhmgei2/4.0 + dwvcrngei/2.0)*s[1]

#define D3GDR3S2S2 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0 + (whhmpyr2+(p-1.0)*wvhmpyr2)/4.0 + (dwhcrnpyr+(p-1.0)*dwvcrnpyr)/2.0)
#define D3GDR3S2DT 0.0
#define D3GDR3S2DP 2.0*(dwvhmpyr/2.0 + wvhmpyr2/4.0 + dwvcrnpyr/2.0)*s[2]

#define D3GDR3DT2  0.0
#define D3GDR3DTDP 0.0
#define D3GDR3DP2  0.0

/*----------------------------------------------------------------------------*/
/* Excess free energy modifications are complete below this point             */
/*----------------------------------------------------------------------------*/

#define D3GDS0S0S0  - 0.125*R*t*(   1.0/SQUARE(xfe2a) - 1.0/SQUARE(xti4a) - 1.0/SQUARE(xfe2b) + 1.0/SQUARE(xti4b) ) \
- 24.0*(whilm+(p-1.0)*wvilm)*s[0]
#define D3GDS0S0S1  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS0S0S2  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS0S0DT  0.25*R*( 1.0/xfe2a + 1.0/xti4a + 1.0/xfe2b + 1.0/xti4b )
#define D3GDS0S0DP  - 2.0*(dvilm+(dwvhmilm/2.0+wvhmilm2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrnilm*r[3]/2.0) + 2.0*wvilm*(r[0]*r[0]-6.0*s[0]*s[0])
#define D3GDS0S1S1  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS0S1S2  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS0S1DT  0.25*R*( 1.0/xti4a + 1.0/xti4b )
#define D3GDS0S1DP  0.0
#define D3GDS0S2S2  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS0S2DT  0.25*R*( 1.0/xti4a + 1.0/xti4b )
#define D3GDS0S2DP  0.0
#define D3GDS0DT2   0.0
#define D3GDS0DTDP  0.0
#define D3GDS0DP2   0.0

#define D3GDS1S1S1  - 0.125*R*t*( 1.0/SQUARE(xmg2a) - 1.0/SQUARE(xti4a) - 1.0/SQUARE(xmg2b) + 1.0/SQUARE(xti4b) ) \
- 24.0*(whgei+(p-1.0)*wvgei)*s[1]
#define D3GDS1S1S2  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS1S1DT  0.25*R*( 1.0/xmg2a + 1.0/xti4a + 1.0/xmg2b + 1.0/xti4b )
#define D3GDS1S1DP  - 2.0*(dvgei+(dwvhmgei/2.0+wvhmgei2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrngei*r[3]/2.0) + 2.0*wvgei*(r[1]*r[1]-6.0*s[1]*s[1])
#define D3GDS1S2S2  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS1S2DT  0.25*R*( 1.0/xti4a + 1.0/xti4b )
#define D3GDS1S2DP  0.0
#define D3GDS1DT2   0.0
#define D3GDS1DTDP  0.0
#define D3GDS1DP2   0.0

#define D3GDS2S2S2  - 0.125*R*t*( 1.0/SQUARE(xmn2a) - 1.0/SQUARE(xti4a) - 1.0/SQUARE(xmn2b) + 1.0/SQUARE(xti4b) ) \
- 24.0*(whpyr+(p-1.0)*wvpyr)*s[2]
#define D3GDS2S2DT  0.25*R*( 1.0/xmn2a + 1.0/xti4a + 1.0/xmn2b + 1.0/xti4b )
#define D3GDS2S2DP  - 2.0*(dvpyr+(dwvhmpyr/2.0+wvhmpyr2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrnpyr*r[3]/2.0) + 2.0*wvpyr*(r[2]*r[2]-6.0*s[2]*s[2])
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
d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1;      d2gdr2[0][2] = D2GDR0R2;     d2gdr2[0][3] = D2GDR0R3; \
d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;      d2gdr2[1][2] = D2GDR1R2;     d2gdr2[1][3] = D2GDR1R3; \
d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2];  d2gdr2[2][2] = D2GDR2R2;     d2gdr2[2][3] = D2GDR2R3; \
d2gdr2[3][0] = d2gdr2[0][3]; d2gdr2[3][1] = d2gdr2[1][3];  d2gdr2[3][2] = d2gdr2[2][3]; d2gdr2[3][3] = D2GDR3R3;

#define fillD2GDRDS \
d2gdrds[0][0] = D2GDR0S0; d2gdrds[0][1] = D2GDR0S1; d2gdrds[0][2] = D2GDR0S2; \
d2gdrds[1][0] = D2GDR1S0; d2gdrds[1][1] = D2GDR1S1; d2gdrds[1][2] = D2GDR1S2; \
d2gdrds[2][0] = D2GDR2S0; d2gdrds[2][1] = D2GDR2S1; d2gdrds[2][2] = D2GDR2S2; \
d2gdrds[3][0] = D2GDR3S0; d2gdrds[3][1] = D2GDR3S1; d2gdrds[3][2] = D2GDR3S2;

#define fillD2GDRDT \
d2gdrdt[0] = D2GDR0DT; d2gdrdt[1] = D2GDR1DT; d2gdrdt[2] = D2GDR2DT; d2gdrdt[3] = D2GDR3DT;

#define fillD2GDRDP \
d2gdrdp[0] = D2GDR0DP; d2gdrdp[1] = D2GDR1DP; d2gdrdp[2] = D2GDR2DP; d2gdrdp[3] = D2GDR3DP;

#define fillD2GDS2 \
d2gds2[0][0] = D2GDS0S0;     d2gds2[0][1] = D2GDS0S1;     d2gds2[0][2] = D2GDS0S2; \
d2gds2[1][0] = d2gds2[0][1]; d2gds2[1][1] = D2GDS1S1;     d2gds2[1][2] = D2GDS1S2; \
d2gds2[2][0] = d2gds2[0][2]; d2gds2[2][1] = d2gds2[1][2]; d2gds2[2][2] = D2GDS2S2;

#define fillD2GDSDT \
d2gdsdt[0] = D2GDS0DT;  d2gdsdt[1] = D2GDS1DT; d2gdsdt[2] = D2GDS2DT;

#define fillD2GDSDP \
d2gdsdp[0] = D2GDS0DP;  d2gdsdp[1] = D2GDS1DP; d2gdsdp[2] = D2GDS2DP;

#define fillD3GDR3 \
d3gdr3[0][0][0] = D3GDR0R0R0;	     d3gdr3[0][0][1] = D3GDR0R0R1;       d3gdr3[0][0][2] = D3GDR0R0R2;      d3gdr3[0][0][3] = D3GDR0R0R3;      \
d3gdr3[0][1][0] = d3gdr3[0][0][1];  d3gdr3[0][1][1] = D3GDR0R1R1;       d3gdr3[0][1][2] = D3GDR0R1R2;      d3gdr3[0][1][3] = D3GDR0R1R3;      \
d3gdr3[0][2][0] = d3gdr3[0][0][2];  d3gdr3[0][2][1] = d3gdr3[0][1][2];  d3gdr3[0][2][2] = D3GDR0R2R2;      d3gdr3[0][2][3] = D3GDR0R2R3;      \
d3gdr3[0][3][0] = d3gdr3[0][0][3];  d3gdr3[0][3][1] = d3gdr3[0][1][3];  d3gdr3[0][3][2] = d3gdr3[0][2][3]; d3gdr3[0][3][3] = D3GDR0R3R3;      \
d3gdr3[1][0][0] = d3gdr3[0][0][1];  d3gdr3[1][0][1] = d3gdr3[0][1][1];	 d3gdr3[1][0][2] = d3gdr3[0][1][2]; d3gdr3[1][0][3] = d3gdr3[0][1][3]; \
d3gdr3[1][1][0] = d3gdr3[0][1][1];  d3gdr3[1][1][1] = D3GDR1R1R1;       d3gdr3[1][1][2] = D3GDR1R1R2;      d3gdr3[1][1][3] = D3GDR1R1R3;      \
d3gdr3[1][2][0] = d3gdr3[0][1][2];  d3gdr3[1][2][1] = d3gdr3[1][1][2];  d3gdr3[1][2][2] = D3GDR1R2R2;      d3gdr3[1][2][3] = D3GDR1R2R3;      \
d3gdr3[1][3][0] = d3gdr3[0][1][3];  d3gdr3[1][3][1] = d3gdr3[1][1][3];  d3gdr3[1][3][2] = d3gdr3[1][2][3]; d3gdr3[1][3][3] = D3GDR1R3R3;      \
d3gdr3[2][0][0] = d3gdr3[0][0][2];  d3gdr3[2][0][1] = d3gdr3[0][1][2];  d3gdr3[2][0][2] = d3gdr3[0][2][2]; d3gdr3[2][0][3] = d3gdr3[0][2][3]; \
d3gdr3[2][1][0] = d3gdr3[0][1][2];  d3gdr3[2][1][1] = d3gdr3[1][1][2];	 d3gdr3[2][1][2] = d3gdr3[1][2][2]; d3gdr3[2][1][3] = d3gdr3[1][2][3]; \
d3gdr3[2][2][0] = d3gdr3[0][2][2];  d3gdr3[2][2][1] = d3gdr3[1][2][2];  d3gdr3[2][2][2] = D3GDR2R2R2;      d3gdr3[2][2][3] = D3GDR2R2R3;      \
d3gdr3[2][3][0] = d3gdr3[0][2][3];  d3gdr3[2][3][1] = d3gdr3[1][2][3];  d3gdr3[2][3][2] = d3gdr3[2][2][3]; d3gdr3[2][3][3] = D3GDR2R3R3;      \
d3gdr3[3][0][0] = d3gdr3[0][0][3];  d3gdr3[3][0][1] = d3gdr3[0][1][3];  d3gdr3[3][0][2] = d3gdr3[0][2][3]; d3gdr3[3][0][3] = d3gdr3[0][3][3]; \
d3gdr3[3][1][0] = d3gdr3[0][1][3];  d3gdr3[3][1][1] = d3gdr3[1][1][3];  d3gdr3[3][1][2] = d3gdr3[1][2][3]; d3gdr3[3][1][3] = d3gdr3[1][3][3]; \
d3gdr3[3][2][0] = d3gdr3[0][2][3];  d3gdr3[3][2][1] = d3gdr3[1][2][3];  d3gdr3[3][2][2] = d3gdr3[2][2][3]; d3gdr3[3][2][3] = d3gdr3[2][3][3]; \
d3gdr3[3][3][0] = d3gdr3[0][3][3];  d3gdr3[3][3][1] = d3gdr3[1][3][3];  d3gdr3[3][3][2] = d3gdr3[2][3][3]; d3gdr3[3][3][3] = D3GDR3R3R3;

#define fillD3GDR2DS \
d3gdr2ds[0][0][0] = D3GDR0R0S0;        d3gdr2ds[0][0][1] = D3GDR0R0S1;        d3gdr2ds[0][0][2] = D3GDR0R0S2;	      \
d3gdr2ds[0][1][0] = D3GDR0R1S0;        d3gdr2ds[0][1][1] = D3GDR0R1S1;        d3gdr2ds[0][1][2] = D3GDR0R1S2;	      \
d3gdr2ds[0][2][0] = D3GDR0R2S0;        d3gdr2ds[0][2][1] = D3GDR0R2S1;        d3gdr2ds[0][2][2] = D3GDR0R2S2;	      \
d3gdr2ds[0][3][0] = D3GDR0R3S0;        d3gdr2ds[0][3][1] = D3GDR0R3S1;        d3gdr2ds[0][3][2] = D3GDR0R3S2;	      \
d3gdr2ds[1][0][0] = d3gdr2ds[0][1][0]; d3gdr2ds[1][0][1] = d3gdr2ds[0][1][1]; d3gdr2ds[1][0][2] = d3gdr2ds[0][1][2]; \
d3gdr2ds[1][1][0] = D3GDR1R1S0;        d3gdr2ds[1][1][1] = D3GDR1R1S1;        d3gdr2ds[1][1][2] = D3GDR1R1S2;	      \
d3gdr2ds[1][2][0] = D3GDR1R2S0;        d3gdr2ds[1][2][1] = D3GDR1R2S1;        d3gdr2ds[1][2][2] = D3GDR1R2S2;	      \
d3gdr2ds[1][3][0] = D3GDR1R3S0;        d3gdr2ds[1][3][1] = D3GDR1R3S1;        d3gdr2ds[1][3][2] = D3GDR1R3S2;	      \
d3gdr2ds[2][0][0] = d3gdr2ds[0][2][0]; d3gdr2ds[2][0][1] = d3gdr2ds[0][2][1]; d3gdr2ds[2][0][2] = d3gdr2ds[0][2][2]; \
d3gdr2ds[2][1][0] = d3gdr2ds[1][2][0]; d3gdr2ds[2][1][1] = d3gdr2ds[1][2][1]; d3gdr2ds[2][1][2] = d3gdr2ds[1][2][2]; \
d3gdr2ds[2][2][0] = D3GDR2R2S0;        d3gdr2ds[2][2][1] = D3GDR2R2S1;        d3gdr2ds[2][2][2] = D3GDR2R2S2;	      \
d3gdr2ds[2][3][0] = D3GDR2R3S0;        d3gdr2ds[2][3][1] = D3GDR2R3S1;        d3gdr2ds[2][3][2] = D3GDR2R3S2;	      \
d3gdr2ds[3][0][0] = d3gdr2ds[0][3][0]; d3gdr2ds[3][0][1] = d3gdr2ds[0][3][1]; d3gdr2ds[3][0][2] = d3gdr2ds[0][3][2]; \
d3gdr2ds[3][1][0] = d3gdr2ds[1][3][0]; d3gdr2ds[3][1][1] = d3gdr2ds[1][3][1]; d3gdr2ds[3][1][2] = d3gdr2ds[1][3][2]; \
d3gdr2ds[3][2][0] = d3gdr2ds[2][3][0]; d3gdr2ds[3][2][1] = d3gdr2ds[2][3][1]; d3gdr2ds[3][2][2] = d3gdr2ds[2][3][2]; \
d3gdr2ds[3][3][0] = D3GDR3R3S0;	d3gdr2ds[3][3][1] = D3GDR3R3S1;        d3gdr2ds[3][3][2] = D3GDR3R3S2;

#define fillD3GDR2DT \
d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT;      d3gdr2dt[0][2] = D3GDR0R2DT;     d3gdr2dt[0][3] = D3GDR0R3DT; \
d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT;      d3gdr2dt[1][2] = D3GDR1R2DT;     d3gdr2dt[1][3] = D3GDR1R3DT; \
d3gdr2dt[2][0] = d3gdr2dt[0][2]; d3gdr2dt[2][1] = d3gdr2dt[1][2];  d3gdr2dt[2][2] = D3GDR2R2DT;     d3gdr2dt[2][3] = D3GDR2R3DT; \
d3gdr2dt[3][0] = d3gdr2dt[0][3]; d3gdr2dt[3][1] = d3gdr2dt[1][3];  d3gdr2dt[3][2] = d3gdr2dt[2][3]; d3gdr2dt[3][3] = D3GDR3R3DT;

#define fillD3GDR2DP \
d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP;      d3gdr2dp[0][2] = D3GDR0R2DP;     d3gdr2dp[0][3] = D3GDR0R3DP; \
d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP;      d3gdr2dp[1][2] = D3GDR1R2DP;     d3gdr2dp[1][3] = D3GDR1R3DP; \
d3gdr2dp[2][0] = d3gdr2dp[0][2]; d3gdr2dp[2][1] = d3gdr2dp[1][2];  d3gdr2dp[2][2] = D3GDR2R2DP;     d3gdr2dp[2][3] = D3GDR2R3DP; \
d3gdr2dp[3][0] = d3gdr2dp[0][3]; d3gdr2dp[3][1] = d3gdr2dp[1][3];  d3gdr2dp[3][2] = d3gdr2dp[2][3]; d3gdr2dp[3][3] = D3GDR3R3DP;

#define fillD3GDRDS2 \
d3gdrds2[0][0][0] = D3GDR0S0S0;        d3gdrds2[0][0][1] = D3GDR0S0S1;        d3gdrds2[0][0][2] = D3GDR0S0S2; \
d3gdrds2[0][1][0] = d3gdrds2[0][0][1]; d3gdrds2[0][1][1] = D3GDR0S1S1;        d3gdrds2[0][1][2] = D3GDR0S1S2; \
d3gdrds2[0][2][0] = d3gdrds2[0][0][2]; d3gdrds2[0][2][1] = d3gdrds2[0][1][2]; d3gdrds2[0][2][2] = D3GDR0S2S2; \
d3gdrds2[1][0][0] = D3GDR1S0S0;        d3gdrds2[1][0][1] = D3GDR1S0S1;        d3gdrds2[1][0][2] = D3GDR1S0S2; \
d3gdrds2[1][1][0] = d3gdrds2[1][0][1]; d3gdrds2[1][1][1] = D3GDR1S1S1;        d3gdrds2[1][1][2] = D3GDR1S1S2; \
d3gdrds2[1][2][0] = d3gdrds2[1][0][2]; d3gdrds2[1][2][1] = d3gdrds2[1][1][2]; d3gdrds2[1][2][2] = D3GDR1S2S2; \
d3gdrds2[2][0][0] = D3GDR2S0S0;        d3gdrds2[2][0][1] = D3GDR2S0S1;        d3gdrds2[2][0][2] = D3GDR2S0S2; \
d3gdrds2[2][1][0] = d3gdrds2[2][0][1]; d3gdrds2[2][1][1] = D3GDR2S1S1;        d3gdrds2[2][1][2] = D3GDR2S1S2; \
d3gdrds2[2][2][0] = d3gdrds2[2][0][2]; d3gdrds2[2][2][1] = d3gdrds2[2][1][2]; d3gdrds2[2][2][2] = D3GDR2S2S2; \
d3gdrds2[3][0][0] = D3GDR3S0S0;        d3gdrds2[3][0][1] = D3GDR3S0S1;        d3gdrds2[3][0][2] = D3GDR3S0S2; \
d3gdrds2[3][1][0] = d3gdrds2[3][0][1]; d3gdrds2[3][1][1] = D3GDR3S1S1;        d3gdrds2[3][1][2] = D3GDR3S1S2; \
d3gdrds2[3][2][0] = d3gdrds2[3][0][2]; d3gdrds2[3][2][1] = d3gdrds2[3][1][2]; d3gdrds2[3][2][2] = D3GDR3S2S2;

#define fillD3GDRDSDT \
d3gdrdsdt[0][0] = D3GDR0S0DT; d3gdrdsdt[0][1] = D3GDR0S1DT; d3gdrdsdt[0][2] = D3GDR0S2DT; \
d3gdrdsdt[1][0] = D3GDR1S0DT; d3gdrdsdt[1][1] = D3GDR1S1DT; d3gdrdsdt[1][2] = D3GDR1S2DT; \
d3gdrdsdt[2][0] = D3GDR2S0DT; d3gdrdsdt[2][1] = D3GDR2S1DT; d3gdrdsdt[2][2] = D3GDR2S2DT; \
d3gdrdsdt[3][0] = D3GDR3S0DT; d3gdrdsdt[3][1] = D3GDR3S1DT; d3gdrdsdt[3][2] = D3GDR3S2DT;

#define fillD3GDRDSDP \
d3gdrdsdp[0][0] = D3GDR0S0DP; d3gdrdsdp[0][1] = D3GDR0S1DP; d3gdrdsdp[0][2] = D3GDR0S2DP; \
d3gdrdsdp[1][0] = D3GDR1S0DP; d3gdrdsdp[1][1] = D3GDR1S1DP; d3gdrdsdp[1][2] = D3GDR1S2DP; \
d3gdrdsdp[2][0] = D3GDR2S0DP; d3gdrdsdp[2][1] = D3GDR2S1DP; d3gdrdsdp[2][2] = D3GDR2S2DP; \
d3gdrdsdp[3][0] = D3GDR3S0DP; d3gdrdsdp[3][1] = D3GDR3S1DP; d3gdrdsdp[3][2] = D3GDR3S2DP;

#define fillD3GDS3 \
d3gds3[0][0][0] = D3GDS0S0S0;      d3gds3[0][0][1] = D3GDS0S0S1;      d3gds3[0][0][2] = D3GDS0S0S2; \
d3gds3[0][1][0] = d3gds3[0][0][1]; d3gds3[0][1][1] = D3GDS0S1S1;      d3gds3[0][1][2] = D3GDS0S1S2; \
d3gds3[0][2][0] = d3gds3[0][0][2]; d3gds3[0][2][1] = d3gds3[0][1][2]; d3gds3[0][2][2] = D3GDS0S2S2; \
d3gds3[1][0][0] = d3gds3[0][0][1]; d3gds3[1][0][1] = d3gds3[0][1][1]; d3gds3[1][0][2] = d3gds3[0][1][2]; \
d3gds3[1][1][0] = d3gds3[0][1][1]; d3gds3[1][1][1] = D3GDS1S1S1;      d3gds3[1][1][2] = D3GDS1S1S2; \
d3gds3[1][2][0] = d3gds3[0][1][2]; d3gds3[1][2][1] = d3gds3[1][1][2]; d3gds3[1][2][2] = D3GDS1S2S2; \
d3gds3[2][0][0] = d3gds3[0][0][2]; d3gds3[2][0][1] = d3gds3[0][1][2]; d3gds3[2][0][2] = d3gds3[0][2][2]; \
d3gds3[2][1][0] = d3gds3[0][1][2]; d3gds3[2][1][1] = d3gds3[1][1][2]; d3gds3[2][1][2] = d3gds3[1][2][2]; \
d3gds3[2][2][0] = d3gds3[0][2][2]; d3gds3[2][2][1] = d3gds3[1][2][2]; d3gds3[2][2][2] = D3GDS2S2S2;

#define fillD3GDS2DT \
d3gds2dt[0][0] = D3GDS0S0DT;     d3gds2dt[0][1] = D3GDS0S1DT;     d3gds2dt[0][2] = D3GDS0S2DT; \
d3gds2dt[1][0] = d3gds2dt[0][1]; d3gds2dt[1][1] = D3GDS1S1DT;     d3gds2dt[1][2] = D3GDS1S2DT; \
d3gds2dt[2][0] = d3gds2dt[0][2]; d3gds2dt[2][1] = d3gds2dt[1][2]; d3gds2dt[2][2] = D3GDS2S2DT;

#define fillD3GDS2DP \
d3gds2dp[0][0] = D3GDS0S0DP;     d3gds2dp[0][1] = D3GDS0S1DP;     d3gds2dp[0][2] = D3GDS0S2DP; \
d3gds2dp[1][0] = d3gds2dp[0][1]; d3gds2dp[1][1] = D3GDS1S1DP;     d3gds2dp[1][2] = D3GDS1S2DP; \
d3gds2dp[2][0] = d3gds2dp[0][2]; d3gds2dp[2][1] = d3gds2dp[1][2]; d3gds2dp[2][2] = D3GDS2S2DP;

#define fillD3GDSDT2 \
d3gdsdt2[0] = D3GDS0DT2; d3gdsdt2[1] = D3GDS1DT2; d3gdsdt2[2] = D3GDS2DT2;

#define fillD3GDSDTDP \
d3gdsdtdp[0] = D3GDS0DTDP; d3gdsdtdp[1] = D3GDS1DTDP; d3gdsdtdp[2] = D3GDS2DTDP;

#define fillD3GDSDP2 \
d3gdsdp2[0] = D3GDS0DP2; d3gdsdp2[1] = D3GDS1DP2; d3gdsdp2[2] = D3GDS2DP2;

#define fillD3GDRDT2 \
d3gdrdt2[0] = D3GDR0DT2; d3gdrdt2[1] = D3GDR1DT2; d3gdrdt2[2] = D3GDR2DT2; d3gdrdt2[3] = D3GDR3DT2;

#define fillD3GDRDTDP \
d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP; d3gdrdtdp[2] = D3GDR2DTDP; d3gdrdtdp[3] = D3GDR3DTDP;

#define fillD3GDRDP2 \
d3gdrdp2[0] = D3GDR0DP2; d3gdrdp2[1] = D3GDR1DP2; d3gdrdp2[2] = D3GDR2DP2; d3gdrdp2[3] = D3GDR3DP2;

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

	xfe2ID = (r[0]/2.0 > 0.0)		         ? r[0]/2.0		           : DBL_EPSILON;
	xmg2ID = (r[1]/2.0 > 0.0)		         ? r[1]/2.0		           : DBL_EPSILON;
	xmn2ID = (r[2]/2.0 > 0.0)		         ? r[2]/2.0		           : DBL_EPSILON;
	xti4ID = ((r[0]+r[1]+r[2])/2.0 > 0.0)    ? (r[0]+r[1]+r[2])/2.0    : DBL_EPSILON;
	xal3ID = (r[3] > 0.0) 		             ? r[3]		               : DBL_EPSILON;
	xfe3ID = (1.0-r[0]-r[1]-r[2]-r[3] > 0.0) ? 1.0-r[0]-r[1]-r[2]-r[3] : DBL_EPSILON;

	/* look-up or compute the current ordering state */
	if ( (t != tOld)       || (p != pOld) ||
		(r[0] != rOld[0]) || (r[1] != rOld[1]) || (r[2] != rOld[2]) || (r[3] != rOld[3]) ) {
		double dgds[NS], sNew[NS];
		for (i=0; i<NS; i++) { sOld[i] = 2.0; sNew[i] = 0.9*r[i]; }

		while ( ((fabs(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) ||
				 (fabs(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON) ||
				 (fabs(sNew[2]-sOld[2]) > 10.0*DBL_EPSILON) ) && (iter < MAX_ITER)) {
			double s[NS];

			for (i=0; i<NS; i++) s[i] = sNew[i];

			xfe2a = (r[0] + s[0])/2.0;
			xmg2a = (r[1] + s[1])/2.0;
			xmn2a = (r[2] + s[2])/2.0;
			xti4a = (r[0] - s[0] + r[1] - s[1] + r[2] - s[2])/2.0;
			xal3a = r[3];
			xfe3a = 1.0 - r[0] - r[1] - r[2] - r[3];

			xfe2b = (r[0] - s[0])/2.0;
			xmg2b = (r[1] - s[1])/2.0;
			xmn2b = (r[2] - s[2])/2.0;
			xti4b = (r[0] + s[0] + r[1] + s[1] + r[2] + s[2])/2.0;
			xal3b = r[3];
			xfe3b = 1.0 - r[0] - r[1] - r[2] - r[3];

			if (xfe2a <= 0.0) xfe2a = DBL_EPSILON; /* added in V1.0-7 */
			if (xmg2a <= 0.0) xmg2a = DBL_EPSILON; /* added in V1.0-7 */
			if (xmn2a <= 0.0) xmn2a = DBL_EPSILON; /* added in V1.0-7 */
			if (xti4a <= 0.0) xti4a = DBL_EPSILON; /* added in V1.0-7 */
			if (xal3a <= 0.0) xal3a = DBL_EPSILON; /* added in V1.0-7 */
			if (xfe3a <= 0.0) xfe3a = DBL_EPSILON; /* added in V1.0-7 */

			if (xfe2b <= 0.0) xfe2b = DBL_EPSILON; /* added in V1.0-7 */
			if (xmg2b <= 0.0) xmg2b = DBL_EPSILON; /* added in V1.0-7 */
			if (xmn2b <= 0.0) xmn2b = DBL_EPSILON; /* added in V1.0-7 */
			if (xti4b <= 0.0) xti4b = DBL_EPSILON; /* added in V1.0-7 */
			if (xal3b <= 0.0) xal3b = DBL_EPSILON; /* added in V1.0-7 */
			if (xfe3b <= 0.0) xfe3b = DBL_EPSILON; /* added in V1.0-7 */

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
				s[i] = MIN(s[i],  r[i] - DBL_EPSILON);
				s[i] = MAX(s[i], -r[i] + DBL_EPSILON);
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
				if (dgds[i]*r[i] > sqrt(DBL_EPSILON) && fabs(sOld[i]) > 10.0*DBL_EPSILON) {
					NSLog(@"ERROR in RHOMBOHEDRAL.C (function ORDER). Failed to converge!");
					NSLog(@"  t,p  = %13.6g, %13.6g", t, p);
					NSLog(@"  X Il = %13.6g, X Gk = %13.6g, X Py = %13.6g, X Crn = %13.6g", r[0], r[1], r[2], r[3]);
					NSLog(@"  s    = %13.6g, %13.6g, %13.6g", sOld[0], sOld[1], sOld[2]);
					NSLog(@"  dgds = %13.6g, %13.6g, %13.6g", dgds[0], dgds[1], dgds[2]);
					NSLog(@"  X Fe2+ A: %13.6g  X Fe2+ B: %13.6g", xfe2a, xfe2b);
					NSLog(@"  X Mg   A: %13.6g  X MG   B: %13.6g", xmg2a, xmg2b);
					NSLog(@"  X Mn   A: %13.6g  X Mn   B: %13.6g", xmn2a, xmn2b);
					NSLog(@"  X Ti   A: %13.6g  X Ti   B: %13.6g", xti4a, xti4b);
					NSLog(@"  X Al3+ A: %13.6g  X Al3+ B: %13.6g", xal3a, xal3b);
					NSLog(@"  X Fe3+ A: %13.6g  X Fe3+ B: %13.6g", xfe3a, xfe3b);
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
	const char *NAMES[NA]    = { "geikielite", "hematite", "ilmenite", "pyrophanite", "corundum" };
	const char *FORMULAS[NA] = { "MgTiO3", "Fe2O3", "FeTiO3", "MnTiO3", "Al2O3" };
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
	 (3)  THIRD            FOURTH | SEVENTH | EIGHTH

	 (1) converts a vector of moles of elements into a vector of moles of
	 endmember rhm oxides components.
	 (2) calculates from a vector of moles of endmember components, one or
	 all of: r[], x[], dr[]/dm[] d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
	 (3) calculates from a vector of independent compositional variables
	 mole fractions of endmember components and/or the Jacobian matrix
	 dx[]/dr[] or a vector of ordering parameters (returned in *x)

	 In this routine it is assumed that the elements are in the order of atomic
	 numbers and that the order of rhm oxides components has been verified as:
	 m[0] = geikielite  (MgTiO3) ,
	 m[1] = hematite    (Fe2O3),
	 m[2] = ilmenite    (FeTiO3),
	 m[3] = pyrophanite (MnTiO3),
	 m[4] = corundum    (Al2O3)

	 ----------------------------------------------------------------------------*/

	int i, j, k;

	if (inpMask == FIRST && outMask == SECOND) {
		/* Converts a vector of moles of elements into a vector of moles of
		 end-member components.                                                 */
		double sumcat, sumchg, fe2, fe3;
        static const int Mg = 12;
        static const int Al = 13;
        static const int Ti = 22;
        static const int Cr = 24;
        static const int Mn = 25;
        static const int Fe = 26;
        static const int Co = 27;
        static const int Ni = 28;

        sumcat = e[Mg] + e[Al] + e[Ti] + e[Cr] + e[Mn] + e[Fe] + e[Co] + e[Ni];

        sumchg = 2.0*e[Mg] + 3.0*e[Al] + 4.0*e[Ti] + 3.0*e[Cr] + 2.0*e[Mn]
               + 2.0*e[Co] + 2.0*e[Ni];
        fe3 = 6.0*sumcat/2.0 - sumchg - 2.0*e[Fe];
        fe2 = e[Fe] - fe3;

        m[0] = e[Mg];                     /* Moles of MgTiO3 */
        m[1] = fe3/2.0;                   /* Moles of Fe2O3  */
        m[2] = fe2;                       /* Moles of FeTiO3 */
        m[3] = e[Mn] + e[Co] + e[Ni];     /* Moles of MnTiO3 */
        m[4] = (e[Al] + e[Cr])/2.0;       /* Moles of Al2O3 */

        if (m[1] < 0.0) m[1] = 0.0;
        if (m[2] < 0.0) m[2] = 0.0;

	} else if (inpMask == SECOND) {
		double sum;

		if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
			NSLog(@"Illegal call to conRhm with inpMask = %o and outMask = %o", inpMask, outMask);

		for (i=0, sum=0.0; i<NA; i++) sum += m[i];

		if (outMask & THIRD) {
			/* Converts a vector of moles of end-member components (m) into a vector
			 of independent compositional variables (r) required as input for the
			 remaining public functions.                                          */
			r[0] = (sum != 0.0) ? m[2]/sum : 0.0;  /* Xil = X FeTiO3 */
			r[1] = (sum != 0.0) ? m[0]/sum : 0.0;  /* Xgk = X MgTiO3 */
			r[2] = (sum != 0.0) ? m[3]/sum : 0.0;  /* Xpy = X MnTiO3 */
			r[3] = (sum != 0.0) ? m[4]/sum : 0.0;  /* XCn = X Al2O3  */
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
					dm[0][j] = (j == 2) ? (1.0-m[2]/sum)/sum : -m[2]/SQUARE(sum);
					dm[1][j] = (j == 0) ? (1.0-m[0]/sum)/sum : -m[0]/SQUARE(sum);
					dm[2][j] = (j == 3) ? (1.0-m[3]/sum)/sum : -m[3]/SQUARE(sum);
					dm[3][j] = (j == 4) ? (1.0-m[4]/sum)/sum : -m[4]/SQUARE(sum);
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
						d2m[0][j][k]  = 2.0*m[2]/CUBE(sum);
						d2m[0][j][k] -= (j == 2) ? 1.0/SQUARE(sum) : 0.0;
						d2m[0][j][k] -= (k == 2) ? 1.0/SQUARE(sum) : 0.0;
						d2m[1][j][k]  = 2.0*m[0]/CUBE(sum);
						d2m[1][j][k] -= (j == 0) ? 1.0/SQUARE(sum) : 0.0;
						d2m[1][j][k] -= (k == 0) ? 1.0/SQUARE(sum) : 0.0;
						d2m[2][j][k]  = 2.0*m[3]/CUBE(sum);
						d2m[2][j][k] -= (j == 3) ? 1.0/SQUARE(sum) : 0.0;
						d2m[2][j][k] -= (k == 3) ? 1.0/SQUARE(sum) : 0.0;
						d2m[3][j][k]  = 2.0*m[4]/CUBE(sum);
						d2m[3][j][k] -= (j == 4) ? 1.0/SQUARE(sum) : 0.0;
						d2m[3][j][k] -= (k == 4) ? 1.0/SQUARE(sum) : 0.0;
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
							d3m[0][j][k][l]  = -6.0*m[2]/QUARTIC(sum);
							d3m[0][j][k][l] += (j == 2) ? 2.0/CUBE(sum) : 0.0;
							d3m[0][j][k][l] += (k == 2) ? 2.0/CUBE(sum) : 0.0;
							d3m[0][j][k][l] += (l == 2) ? 2.0/CUBE(sum) : 0.0;
							d3m[1][j][k][l]  = -6.0*m[0]/QUARTIC(sum);
							d3m[1][j][k][l] += (j == 0) ? 2.0/CUBE(sum) : 0.0;
							d3m[1][j][k][l] += (k == 0) ? 2.0/CUBE(sum) : 0.0;
							d3m[1][j][k][l] += (l == 0) ? 2.0/CUBE(sum) : 0.0;
							d3m[2][j][k][l]  = -6.0*m[3]/QUARTIC(sum);
							d3m[2][j][k][l] += (j == 3) ? 2.0/CUBE(sum) : 0.0;
							d3m[2][j][k][l] += (k == 3) ? 2.0/CUBE(sum) : 0.0;
							d3m[2][j][k][l] += (l == 3) ? 2.0/CUBE(sum) : 0.0;
							d3m[3][j][k][l]  = -6.0*m[4]/QUARTIC(sum);
							d3m[3][j][k][l] += (j == 4) ? 2.0/CUBE(sum) : 0.0;
							d3m[3][j][k][l] += (k == 4) ? 2.0/CUBE(sum) : 0.0;
							d3m[3][j][k][l] += (l == 4) ? 2.0/CUBE(sum) : 0.0;
						}
					}
				}
			}
		}

	} else if (inpMask == THIRD) {

		if (outMask & ~(FOURTH | SEVENTH | EIGHTH))
			NSLog(@"Illegal call to conRhm with inpMask = %o and outMask = %o", inpMask, outMask);

		if (outMask & FOURTH) {
			/* Converts a vector of independent compositional variables (r) into a
			 vector of mole fractions of endmember components (x).                */
			x[0] = r[1];
			x[1] = 1.0 - r[0] - r[1] - r[2] - r[3];
			x[2] = r[0];
			x[3] = r[2];
			x[4] = r[3];
		}

		if (outMask & SEVENTH) {
			/* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
			for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
			dr[0][1] =  1.0;
			dr[1][0] = -1.0; dr[1][1] = -1.0; dr[1][2] = -1.0; dr[1][3] = -1.0;
			dr[2][0] = 1.0;
			dr[3][2] = 1.0;
			dr[4][3] = 1.0;
		}

	} else  {
		NSLog(@"Illegal call to conRhm with inpMask = %o and outMask = %o", inpMask, outMask);
	}

}

-(NSString *)displayFormula:(double)t
						  p:(double)p
						  r:(double [NA])r
{

	double totMn, totFe2, totMg, totFe3, totAl, totTi;

	totMn  = r[2];
	totFe2 = r[0];
	totMg  = r[1];
	totFe3 = 2.0*(1.0-r[0]-r[1]-r[2]-r[3]);
	totAl  = 2.0*r[3];
	totTi  = r[0]+r[1]+r[2];

	return [NSString stringWithFormat:@"Mn%4.2fFe''%4.2fMg%4.2fFe'''%4.2fAl%4.2fTi%4.2fO3",
			totMn, totFe2, totMg, totFe3, totAl, totTi];
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
		fr[i][0] = FR0(i); /* Xil */
		fr[i][1] = FR1(i); /* Xgk */
		fr[i][2] = FR2(i); /* Xpy */
		fr[i][3] = FR3(i); /* Xcn */
	}

	[self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

	g       = G;
	dgdr[0] = DGDR0;
	dgdr[1] = DGDR1;
	dgdr[2] = DGDR2;
	dgdr[3] = DGDR3;

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
			gs[i][0] = GS0(i);        /* s   */
			gs[i][1] = GS1(i);        /* t   */
			gs[i][2] = GS2(i);        /* u   */
			dfrdr[i][0] = DFR0DR0(i); /* Xil */
			dfrdr[i][1] = DFR1DR1(i); /* Xgk */
			dfrdr[i][2] = DFR2DR2(i); /* Xpy */
			dfrdr[i][3] = DFR3DR3(i); /* Xcn */
			dgsds[i][0] = DGS0DS0(i); /* s   */
			dgsds[i][1] = DGS1DS1(i); /* t   */
			dgsds[i][2] = DGS2DS2(i); /* u   */
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
		dx[2] = DGDR2;
		dx[3] = DGDR3;

		[self pureEndMembers:THIRD t:(double)t	p:(double)p a:NULL mu:NULL gmix:ends hmix:NULL smix:NULL cpmix:NULL
					 cpmixdt:NULL vmix:NULL vmixdt:NULL vmixdp:NULL vmixdt2:NULL vmixdtdp:NULL vmixdp2:NULL];

		dx[0] -= DENDDR0;
		dx[1] -= DENDDR1;
		dx[2] -= DENDDR2;
		dx[3] -= DENDDR3;
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

	*hmix = (G) - t*(DGDT);

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

		*smix = -(DGDT);

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
		dx[2] -= DENDDR2;
		dx[3] -= DENDDR3;
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
		dx[2] -= DENDDR2;
		dx[3] -= DENDDR3;
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
		dx[2] -= DENDDR2;
		dx[3] -= DENDDR3;
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
		dxdt[2] -= DENDDR2;
		dxdt[3] -= DENDDR3;
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
		dxdp[2] -= DENDDR2;
		dxdp[3] -= DENDDR3;
	}

}

-(NSUInteger)numberOfSolutionSpecies {
	return NA;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
	switch (index) {
		case 0:
			return @"geikielite";
			break;
		case 1:
			return @"hematite";
			break;
		case 2:
			return @"ilmenite";
			break;
		case 3:
			return @"pyrophanite";
			break;
		case 4:
			return @"corundum";
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

		if (count > 25) { // pure empiricism
			for (i=0; i<nz; i++) gamma[i] = (gamma[i]+gammaLast[i])/2.0;
		}
		for (i=0; i<nz; i++) gammaLast[i] = gamma[i];

		[self correctActivityCoefficients:gamma forComposition:x];

		if (debugV) {
			NSLog(@"Iteration %lu", count);
			NSLog(@"%13.6g %13.6g %13.6g %13.6g %13.6g %13.6g", a[0], a[1], a[2], a[3], a[4], affinity);
			NSLog(@"%13.6g %13.6g %13.6g %13.6g %13.6g", x[0], x[1], x[2], x[3], x[4]);
            double g[NA];
			for (i=0, j=0; i<NA; i++) g[i] = (x[i] != 0.0) ? gamma[j++] : 0.0;
            NSLog(@"%13.6g %13.6g %13.6g %13.6g %13.6g", g[0], g[1], g[2], g[3], g[4]);
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
