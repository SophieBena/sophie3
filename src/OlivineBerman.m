//
//  OlivineBerman.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/15/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "OlivineBerman.h"
#import "BermanProperties.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"


@implementation OlivineBerman

static NSArray *endmembers;

#define NR         5    /* Five independent composition variables      */
#define NS         4    /* Four ordering parameters                    */
#define NA         6    /* Six endmember compositions                  */
#define NATOMS   7.0    /* Average number of atoms in the formula unit */

#pragma mark -
#pragma mark class methods

+(void)initialize {
	if (self == [OlivineBerman class]) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"Initialize(OlivineBerman) - entry ...");
		NSMutableArray *mutableEndmembers = [NSMutableArray arrayWithCapacity:NA];

		BermanProperties *tephroite = [[BermanProperties alloc] initWithH:-1732000.0
																		S:155.900
																	   k0:219.89
																	   k1:-12.710e2
																	   k2:-20.496e5
																	   k3:17.652e7
																	   v0:4.889
																	   v1:-0.784E-6
																	   v2:0.0
																	   v3:25.14E-6
																	   v4:64.7E-10];
		[tephroite setPhaseFormula:@"Mn2SiO4"];
		[tephroite setPhaseName:@"tephroite"];
		[mutableEndmembers addObject:tephroite];
		if (debug) NSLog(@"... allocated tephroite ...");

		BermanProperties *fayalite = [[BermanProperties alloc] initWithH:-1479360.0
																	   S:150.930
																	  k0:248.93
																	  k1:-19.239e2
																	  k2:0.0
																	  k3:-13.910e7
																	  v0:4.630
																	  v1:-0.730E-6
																	  v2:0.0
																	  v3:26.546E-6
																	  v4:79.482E-10];
		[fayalite setPhaseFormula:@"Fe2SiO4"];
		[fayalite setPhaseName:@"fayalite"];
		[mutableEndmembers addObject:fayalite];
		if (debug) NSLog(@"... allocated fayalite ...");

		BermanProperties *coOlivine = [[BermanProperties alloc] initWithH:-1414100.0
																		S:142.600
																	   k0:201.048
																	   k1:-0.369E2
																	   k2:-71.81E5
																	   k3:90.05E7
																	   v0:4.459
																	   v1:-0.64E-6
																	   v2:0.0E-12
																	   v3:28.422E-6
																	   v4:35.355E-10];
		[coOlivine setPhaseFormula:@"Co2SiO4"];
		[coOlivine setPhaseName:@"co-olivine"];
		[mutableEndmembers addObject:coOlivine];
		if (debug) NSLog(@"... allocated co-olivine ...");

		BermanProperties *niOlivine = [[BermanProperties alloc] initWithH:-1395300.0
																		S:128.100
																	   k0:214.997
																	   k1:-10.308E2
																	   k2:-49.445E5
																	   k3:62.370E7
																	   v0:4.259
																	   v1:-0.671E-6
																	   v2:0.0E-12
																	   v3:28.422E-6
																	   v4:35.355E-10];
		[niOlivine setPhaseFormula:@"Ni2SiO4"];
		[niOlivine setPhaseName:@"ni-olivine"];
		[mutableEndmembers addObject:niOlivine];
		if (debug) NSLog(@"... allocated ni-olivine ...");

		BermanProperties *monticellite = [[BermanProperties alloc] initWithH:-2250027.0
																		   S:108.300
																		  k0:226.34
																		  k1:-15.427E2
																		  k2:-11.797E5
																		  k3:-2.329E7
																		  v0:5.148
																		  v1:-0.904E-6
																		  v2:2.00E-12
																		  v3:27.863E-6
																		  v4:76.339E-10];
		[monticellite setPhaseFormula:@"CaMgSiO4"];
		[monticellite setPhaseName:@"monticellite"];
		[mutableEndmembers addObject:monticellite];
		if (debug) NSLog(@"... allocated monticellite ...");

		BermanProperties *forsterite = [[BermanProperties alloc] initWithH:-2174420.0
																		 S:94.010
																		k0:238.64
																		k1:-20.013E2
																		k2:0.0
																		k3:-11.624E7
																		v0:4.366
																		v1:-0.791E-6
																		v2:1.351E-12
																		v3:29.464E-6
																		v4:88.633E-10];
		[forsterite setPhaseFormula:@"Mg2SiO4"];
		[forsterite setPhaseName:@"forsterite"];
		[mutableEndmembers addObject:forsterite];
		if (debug) NSLog(@"... allocated forsterite and exiting.");

		endmembers = [NSArray arrayWithArray:mutableEndmembers];
	}
}

#pragma mark -
#pragma mark instance methods

-(id)init {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(OlivineBerman) ... entry ...");
		[self setPhaseName:@"Olivine"];
		computeMixingQuantities = NO;
		tOld = -9999.0;
		pOld = -9999.0;
		for (NSUInteger i=0; i<NR; i++) rOld[i] = -9999.0;
		for (NSUInteger i=0; i<NS; i++) sOld[i] = 2.0;
		xm1mg = 0.0;
		xm1mn = 0.0;
		xm1fe = 0.0;
		xm1co = 0.0;
		xm1ni = 0.0;
		xm2mg = 0.0;
		xm2mn = 0.0;
		xm2fe = 0.0;
		xm2co = 0.0;
		xm2ni = 0.0;
		xm2ca = 0.0;
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
#define CUBE(x) ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/**
 ============================================================================
 Olivine solution parameters:

 Hirschmann, M. (1991)  Thermodynamics of multicomponent olivines and the
   solution properties of (Ni,Mg,Fe)2SiO4 and (Ca,Mg,Fe)SiO4 olivines
   American Mineralogist 77:1232-1248.

 Sack, R.O., Ghiorso, M.S. (1989)
   Importance of considerations of mixing properties in establishing
   an internally consistent thermodyanmic database:  Thermochemistry of
   minerals in the system Mg2SIO4-Fe2SiO4-SiO2 Contributions to Mineralogy
   and Petrology 102: 41-68

 Properties of Mn- and Co- bearing olivines - Hirschmann and Ghiorso
   1994 GCA
 */

#define R  8.3143

#define HEXMGMN  15.80 *1000.0 /* joules */
#define HEXMGFE  00.00 *1000.0 /* joules */
#define HEXMGCO -15.00 *1000.0 /* joules */
#define HEXMGNI -19.75 *1000.0 /* joules */
#define HEXMNFE -11.80 *1000.0 /* joules */
#define HEXMNCO  00.00 *1000.0 /* joules */
#define HEXMNNI  00.00 *1000.0 /* joules */
#define HEXFECO  00.00 *1000.0 /* joules */
#define HEXFENI -20.00 *1000.0 /* joules */
#define HEXCONI   0.00 *1000.0 /* joules */

#define VEXMGMN  00.00         /* joules/bar */
#define VEXMGFE  00.00         /* joules/bar */
#define VEXMGCO  00.00         /* joules/bar */
#define VEXMGNI  00.000        /* joules/bar */
#define VEXMNFE  00.00         /* joules/bar */
#define VEXMNCO  00.00         /* joules/bar */
#define VEXMNNI  00.00         /* joules/bar */
#define VEXFECO  00.00         /* joules/bar */
#define VEXFENI  00.000        /* joules/bar */
#define VEXCONI  00.00         /* joules/bar */

#define HXMGMN    8.75 *1000.0 /* joules */
#define HXMGFE   10.15 *1000.0 /* joules */
#define HXMGCO   03.00 *1000.0 /* joules */
#define HXMGNI    2.20 *1000.0 /* joules */
#define HXMNFE    0.50 *1000.0 /* joules */
#define HXMNCO    0.00 *1000.0 /* joules */
#define HXMNNI   00.00 *1000.0 /* joules */
#define HXFECO   03.00 *1000.0 /* joules */
#define HXFENI   10.00 *1000.0 /* joules */
#define HXCONI    0.00 *1000.0 /* joules */

#define VXMGMN   00.00         /* joules/bar */
#define VXMGFE   00.015        /* joules/bar */
#define VXMGCO   00.00         /* joules/bar */
#define VXMGNI   00.000        /* joules/bar */
#define VXMNFE   00.00         /* joules/bar */
#define VXMNCO   00.00         /* joules/bar */
#define VXMNNI   00.00         /* joules/bar */
#define VXFECO   00.00         /* joules/bar */
#define VXFENI   00.000        /* joules/bar */
#define VXCONI   00.00         /* joules/bar */

#define WH1MGMN    6.625*1000.0 /* joules */
#define WH2MGMN    6.625*1000.0 /* joules */
#define WH1MGFE    5.075*1000.0 /* joules */
#define WH2MGFE    5.075*1000.0 /* joules */
#define WH1MGCO    1.50 *1000.0 /* joules */
#define WH2MGCO    1.50 *1000.0 /* joules */
#define WH1MGNI    -.600*1000.0 /* joules */
#define WH2MGNI    2.800*1000.0 /* joules */
#define WH1MNFE    1.75 *1000.0 /* joules */
#define WH2MNFE    1.75 *1000.0 /* joules */
#define WH1MNCO    0.00 *1000.0 /* joules */
#define WH2MNCO    0.00 *1000.0 /* joules */
#define WH1MNNI    0.00 *1000.0 /* joules */
#define WH2MNNI    0.00 *1000.0 /* joules */
#define WH1FECO    1.50 *1000.0 /* joules */
#define WH2FECO    1.50 *1000.0 /* joules */
#define WH1FENI    5.000*1000.0 /* joules */
#define WH2FENI    5.000*1000.0 /* joules */
#define WH1CONI    0.00 *1000.0 /* joules */
#define WH2CONI    0.00 *1000.0 /* joules */

#define WV1MGMN   00.00         /* joules/bar */
#define WV2MGMN   00.00         /* joules/bar */
#define WV1MGFE   00.0000       /* joules/bar */
#define WV2MGFE   00.0000       /* joules/bar */
#define WV1MGCO   00.00         /* joules/bar */
#define WV2MGCO   00.00         /* joules/bar */
#define WV1MGNI   00.0000       /* joules/bar */
#define WV2MGNI   00.0000       /* joules/bar */
#define WV1MNFE   00.00         /* joules/bar */
#define WV2MNFE   00.00         /* joules/bar */
#define WV1MNCO   00.00         /* joules/bar */
#define WV2MNCO   00.00         /* joules/bar */
#define WV1MNNI   00.00         /* joules/bar */
#define WV2MNNI   00.00         /* joules/bar */
#define WV1FECO   00.00         /* joules/bar */
#define WV2FECO   00.00         /* joules/bar */
#define WV1FENI   00.0000       /* joules/bar */
#define WV2FENI   00.0000       /* joules/bar */
#define WV1CONI   00.00         /* joules/bar */
#define WV2CONI   00.00         /* joules/bar */

#define WH2CAMG   34.50 *1000.0 /* joules */
#define WH2CAMN   16.00 *1000.0 /* joules */
#define WH2CAFE   21.90 *1000.0 /* joules */
#define WH2CACO   30.00 *1000.0 /* joules */
#define WH2CANI   40.00 *1000.0 /* joules */

#define WV2CAMG   00.35         /* joules/bar */
#define WV2CAMN   00.00         /* joules/bar */
#define WV2CAFE   00.00         /* joules/bar */
#define WV2CACO   00.00         /* joules/bar */
#define WV2CANI   00.00         /* joules/bar */

#define F_MN     09.50 *1000.0 /* joules */
#define F_FE     09.50 *1000.0 /* joules */
#define F_CO     00.00 *1000.0 /* joules */
#define F_NI     00.00 *1000.0 /* joules */

#define GEXMGMN  (HEXMGMN) + (p-1.0)*(VEXMGMN)
#define GEXMGFE  (HEXMGFE) + (p-1.0)*(VEXMGFE)
#define GEXMGCO  (HEXMGCO) + (p-1.0)*(VEXMGCO)
#define GEXMGNI  (HEXMGNI) + (p-1.0)*(VEXMGNI)
#define GEXMNFE  (HEXMNFE) + (p-1.0)*(VEXMNFE)
#define GEXMNCO  (HEXMNCO) + (p-1.0)*(VEXMNCO)
#define GEXMNNI  (HEXMNNI) + (p-1.0)*(VEXMNNI)
#define GEXFECO  (HEXFECO) + (p-1.0)*(VEXFECO)
#define GEXFENI  (HEXFENI) + (p-1.0)*(VEXFENI)
#define GEXCONI  (HEXCONI) + (p-1.0)*(VEXCONI)

#define GXMGMN  (HXMGMN) + (p-1.0)*(VXMGMN)
#define GXMGFE  (HXMGFE) + (p-1.0)*(VXMGFE)
#define GXMGCO  (HXMGCO) + (p-1.0)*(VXMGCO)
#define GXMGNI  (HXMGNI) + (p-1.0)*(VXMGNI)
#define GXMNFE  (HXMNFE) + (p-1.0)*(VXMNFE)
#define GXMNCO  (HXMNCO) + (p-1.0)*(VXMNCO)
#define GXMNNI  (HXMNNI) + (p-1.0)*(VXMNNI)
#define GXFECO  (HXFECO) + (p-1.0)*(VXFECO)
#define GXFENI  (HXFENI) + (p-1.0)*(VXFENI)
#define GXCONI  (HXCONI) + (p-1.0)*(VXCONI)

#define W1MGMN  (WH1MGMN) + (p-1.0)*(WV1MGMN)
#define W2MGMN  (WH2MGMN) + (p-1.0)*(WV2MGMN)
#define W1MGFE  (WH1MGFE) + (p-1.0)*(WV1MGFE)
#define W2MGFE  (WH2MGFE) + (p-1.0)*(WV2MGFE)
#define W1MGCO  (WH1MGCO) + (p-1.0)*(WV1MGCO)
#define W2MGCO  (WH2MGCO) + (p-1.0)*(WV2MGCO)
#define W1MGNI  (WH1MGNI) + (p-1.0)*(WV1MGNI)
#define W2MGNI  (WH2MGNI) + (p-1.0)*(WV2MGNI)
#define W1MNFE  (WH1MNFE) + (p-1.0)*(WV1MNFE)
#define W2MNFE  (WH2MNFE) + (p-1.0)*(WV2MNFE)
#define W1MNCO  (WH1MNCO) + (p-1.0)*(WV1MNCO)
#define W2MNCO  (WH2MNCO) + (p-1.0)*(WV2MNCO)
#define W1MNNI  (WH1MNNI) + (p-1.0)*(WV1MNNI)
#define W2MNNI  (WH2MNNI) + (p-1.0)*(WV2MNNI)
#define W1FECO  (WH1FECO) + (p-1.0)*(WV1FECO)
#define W2FECO  (WH2FECO) + (p-1.0)*(WV2FECO)
#define W1FENI  (WH1FENI) + (p-1.0)*(WV1FENI)
#define W2FENI  (WH2FENI) + (p-1.0)*(WV2FENI)
#define W1CONI  (WH1CONI) + (p-1.0)*(WV1CONI)
#define W2CONI  (WH2CONI) + (p-1.0)*(WV2CONI)

#define W2CAMG  (WH2CAMG) + (p-1.0)*(WV2CAMG)
#define W2CAMN  (WH2CAMN) + (p-1.0)*(WV2CAMN)
#define W2CAFE  (WH2CAFE) + (p-1.0)*(WV2CAFE)
#define W2CACO  (WH2CACO) + (p-1.0)*(WV2CACO)
#define W2CANI  (WH2CANI) + (p-1.0)*(WV2CANI)


/* Definitions of Taylor expansion coefficients in terms of solution
 * parameters. Independent variables are r1,r2,r3,r4,r5,s1,s2,s3,s4
 */

#define G0  0.25*(     ((GXMNFE)+(W1MNFE)+(W2MNFE)) +((GXMNCO)+(W1MNCO)+(W2MNCO)) +((GXMNNI)+(W1MNNI)+(W2MNNI)) +((GXFECO)+(W1FECO)+(W2FECO)) \
+((GXFENI)+(W1FENI)+(W2FENI)) +((GXCONI)+(W1CONI)+(W2CONI)) -2.0*((GXMGMN)+(W1MGMN)+(W2MGMN)) -2.0*((GXMGFE)+(W1MGFE)+(W2MGFE)) \
-2.0*((GXMGCO)+(W1MGCO)+(W2MGCO)) -2.0*((GXMGNI)+(W1MGNI)+(W2MGNI)))
#define GR1 0.25*(-3.0*((GXMGMN)+(W1MGMN)+(W2MGMN)) \
-((GXMGFE)+(W1MGFE)+(W2MGFE)) \
-((GXMGCO)+(W1MGCO)+(W2MGCO)) \
-((GXMGNI)+(W1MGNI)+(W2MGNI)) \
+((GXMNFE)+(W1MNFE)+(W2MNFE)) \
+((GXMNCO)+(W1MNCO)+(W2MNCO)) \
+((GXMNNI)+(W1MNNI)+(W2MNNI)))
#define GR2  0.25*(   -((GXMGMN)+(W1MGMN)+(W2MGMN)) \
-3.0*((GXMGFE)+(W1MGFE)+(W2MGFE)) \
-((GXMGCO)+(W1MGCO)+(W2MGCO)) \
-((GXMGNI)+(W1MGNI)+(W2MGNI)) \
+((GXMNFE)+(W1MNFE)+(W2MNFE)) \
+((GXFECO)+(W1FECO)+(W2FECO)) \
+((GXFENI)+(W1FENI)+(W2FENI)))
#define GR3  0.25*(   -((GXMGMN)+(W1MGMN)+(W2MGMN)) \
-((GXMGFE)+(W1MGFE)+(W2MGFE)) \
-3.0*((GXMGCO)+(W1MGCO)+(W2MGCO)) \
-((GXMGNI)+(W1MGNI)+(W2MGNI)) \
+((GXMNCO)+(W1MNCO)+(W2MNCO)) \
+((GXFECO)+(W1FECO)+(W2FECO)) \
+((GXCONI)+(W1CONI)+(W2CONI)))
#define GR4 0.25*(    -((GXMGMN)+(W1MGMN)+(W2MGMN)) \
-((GXMGFE)+(W1MGFE)+(W2MGFE)) \
-((GXMGCO)+(W1MGCO)+(W2MGCO)) \
-3.0*((GXMGNI)+(W1MGNI)+(W2MGNI)) \
+((GXMNNI)+(W1MNNI)+(W2MNNI)) \
+((GXFENI)+(W1FENI)+(W2FENI)) \
+((GXCONI)+(W1CONI)+(W2CONI)))
#define GR5        -(W2CAMG) \
-0.25*(((F_MN)+(GEXMGMN)+(GXMGMN) \
-2.0*(W2CAMN)+2.0*(W2MGMN)) \
+((F_FE)+(GEXMGFE)+(GXMGFE) \
-2.0*(W2CAFE)+2.0*(W2MGFE)) \
+((F_CO)+(GEXMGCO)+(GXMGCO) \
-2.0*(W2CACO)+2.0*(W2MGCO)) \
+((F_NI)+(GEXMGNI)+(GXMGNI) \
-2.0*(W2CANI)+2.0*(W2MGNI)))
#define GR1R5             -0.25*((F_MN)+(GEXMGMN)+(GXMGMN) +2.0*(W2CAMG)-2.0*(W2CAMN)+2.0*(W2MGMN))
#define GR2R5             -0.25*((F_FE)+(GEXMGFE)+(GXMGFE) +2.0*(W2CAMG)-2.0*(W2CAFE)+2.0*(W2MGFE))
#define GR3R5             -0.25*((F_CO)+(GEXMGCO)+(GXMGCO) +2.0*(W2CAMG)-2.0*(W2CACO)+2.0*(W2MGCO))
#define GS1   0.25*(( (GEXMGMN)+3.0*(W1MGMN)-3.0*(W2MGMN)) \
-((GEXMGFE)-(W1MGFE)+(W2MGFE)) \
-((GEXMGCO)-(W1MGCO)+(W2MGCO)) \
-((GEXMGNI)-(W1MGNI)+(W2MGNI)) \
+((GEXMNFE)-(W1MNFE)+(W2MNFE)) \
+((GEXMNCO)-(W1MNCO)+(W2MNCO)) \
+((GEXMNNI)-(W1MNNI)+(W2MNNI)))
#define GS2           0.25*((-(GEXMGMN)+(W1MGMN)-(W2MGMN)) \
+((GEXMGFE)+3.0*(W1MGFE)-3.0*(W2MGFE)) \
-((GEXMGCO)-(W1MGCO)+(W2MGCO)) \
-((GEXMGNI)-(W1MGNI)+(W2MGNI)) \
-((GEXMNFE)+(W1MNFE)-(W2MNFE)) \
+((GEXFECO)-(W1FECO)+(W2FECO)) \
+((GEXFENI)-(W1FENI)+(W2FENI)))
#define GS3           0.25*(( (GEXMGMN)-(W1MGMN)+(W2MGMN)) \
+((GEXMGFE)-(W1MGFE)+(W2MGFE)) \
-((GEXMGCO)+3.0*(W1MGCO)-3.0*(W2MGCO)) \
+((GEXMGNI)-(W1MGNI)+(W2MGNI)) \
+((GEXMNCO)+(W1MNCO)-(W2MNCO)) \
+((GEXFECO)+(W1FECO)-(W2FECO)) \
-((GEXCONI)-(W1CONI)+(W2CONI)))
#define GS4           0.25*(( (GEXMGMN)-(W1MGMN)+(W2MGMN)) \
+((GEXMGFE)-(W1MGFE)+(W2MGFE)) \
+((GEXMGCO)-(W1MGCO)+(W2MGCO)) \
-((GEXMGNI)+3.0*(W1MGNI)-3.0*(W2MGNI)) \
+((GEXMNNI)+(W1MNNI)-(W2MNNI)) \
+((GEXFENI)+(W1FENI)-(W2FENI)) \
+((GEXCONI)+(W1CONI)-(W2CONI)))
#define GR1R1         -0.25*(  (GXMGMN)+(W1MGMN)+(W2MGMN))
#define GR1R2          0.25*(( (GXMNFE)+(W1MNFE)+(W2MNFE)) -((GXMGMN)+(W1MGMN)+(W2MGMN)) -((GXMGFE)+(W1MGFE)+(W2MGFE)))
#define GR1R3          0.25*(( (GXMNCO)+(W1MNCO)+(W2MNCO)) -((GXMGMN)+(W1MGMN)+(W2MGMN)) -((GXMGCO)+(W1MGCO)+(W2MGCO)))
#define GR1R4          0.25*(( (GXMNNI)+(W1MNNI)+(W2MNNI)) -((GXMGMN)+(W1MGMN)+(W2MGMN)) -((GXMGNI)+(W1MGNI)+(W2MGNI)))
#define GR1R5             -0.25*((F_MN)+(GEXMGMN)+(GXMGMN) +2.0*(W2CAMG)-2.0*(W2CAMN)+2.0*(W2MGMN))
#define GR1S1                       0.5*((W1MGMN)-(W2MGMN))
#define GR1S2         0.25*((-(GEXMGMN)-(W1MGMN)+(W2MGMN)) +((GEXMGFE)+(W1MGFE)-(W2MGFE)) -((GEXMNFE)+(W1MNFE)-(W2MNFE)))
#define GR1S3         0.25*(( (GEXMGMN)-(W1MGMN)+(W2MGMN)) -((GEXMGCO)+(W1MGCO)-(W2MGCO)) +((GEXMNCO)+(W1MNCO)-(W2MNCO)))
#define GR1S4         0.25*(( (GEXMGMN)-(W1MGMN)+(W2MGMN)) -((GEXMGNI)+(W1MGNI)-(W2MGNI)) +((GEXMNNI)+(W1MNNI)-(W2MNNI)))
#define GR2R2         -0.25*(( (GXMGFE)+(W1MGFE)+(W2MGFE)))
#define GR2R3          0.25*(( (GXFECO)+(W1FECO)+(W2FECO)) -((GXMGFE)+(W1MGFE)+(W2MGFE)) -((GXMGCO)+(W1MGCO)+(W2MGCO)))
#define GR2R4          0.25*(( (GXFENI)+(W1FENI)+(W2FENI)) -((GXMGFE)+(W1MGFE)+(W2MGFE)) -((GXMGNI)+(W1MGNI)+(W2MGNI)))
#define GR2R5             -0.25*((F_FE)+(GEXMGFE)+(GXMGFE) +2.0*(W2CAMG)-2.0*(W2CAFE)+2.0*(W2MGFE))
#define GR2S1         0.25*(( (GEXMGMN)+(W1MGMN)-(W2MGMN)) -((GEXMGFE)-(W1MGFE)+(W2MGFE)) +((GEXMNFE)-(W1MNFE)+(W2MNFE)))
#define GR2S2                       0.5*((W1MGFE)-(W2MGFE))
#define GR2S3         0.25*(( (GEXMGFE)-(W1MGFE)+(W2MGFE)) -((GEXMGCO)+(W1MGCO)-(W2MGCO)) +((GEXFECO)+(W1FECO)-(W2FECO)))
#define GR2S4         0.25*(( (GEXMGFE)-(W1MGFE)+(W2MGFE)) -((GEXMGNI)+(W1MGNI)-(W2MGNI)) +((GEXFENI)+(W1FENI)-(W2FENI)))
#define GR3R3          -0.25*(  (GXMGCO)+(W1MGCO)+(W2MGCO))
#define GR3R4          0.25*(( (GXCONI)+(W1CONI)+(W2CONI)) -((GXMGCO)+(W1MGCO)+(W2MGCO)) -((GXMGNI)+(W1MGNI)+(W2MGNI)))
#define GR3R5             -0.25*((F_CO)+(GEXMGCO)+(GXMGCO) +2.0*(W2CAMG)-2.0*(W2CACO)+2.0*(W2MGCO))
#define GR3S1         0.25*(( (GEXMGMN)+(W1MGMN)-(W2MGMN)) -((GEXMGCO)-(W1MGCO)+(W2MGCO)) +((GEXMNCO)-(W1MNCO)+(W2MNCO)))
#define GR3S2         0.25*(( (GEXMGFE)+(W1MGFE)-(W2MGFE)) -((GEXMGCO)-(W1MGCO)+(W2MGCO)) +((GEXFECO)-(W1FECO)+(W2FECO)))
#define GR3S3                      0.5*(-(W1MGCO)+(W2MGCO))
#define GR3S4          0.25*(( (GEXMGCO)+(W1MGCO)-(W2MGCO)) -((GEXMGNI)-(W1MGNI)+(W2MGNI)) +((GEXCONI)-(W1CONI)+(W2CONI)))
#define GR4R4           -0.25*(  (GXMGNI)+(W1MGNI)+(W2MGNI))
#define GR4R5             -0.25*((F_NI)+(GEXMGNI)+(GXMGNI) +2.0*(W2CAMG)-2.0*(W2CANI)+2.0*(W2MGNI))
#define GR4S1          0.25*(( (GEXMGMN)+(W1MGMN)-(W2MGMN)) -((GEXMGNI)-(W1MGNI)+(W2MGNI)) +((GEXMNNI)-(W1MNNI)+(W2MNNI)))
#define GR4S2          0.25*(( (GEXMGFE)+(W1MGFE)-(W2MGFE)) -((GEXMGNI)-(W1MGNI)+(W2MGNI)) +((GEXFENI)-(W1FENI)+(W2FENI)))
#define GR4S3          0.25*((-(GEXMGCO)+(W1MGCO)-(W2MGCO)) +((GEXMGNI)-(W1MGNI)+(W2MGNI)) -((GEXCONI)-(W1CONI)+(W2CONI)))
#define GR4S4                       0.5*(-(W1MGNI)+(W2MGNI))
#define GR5R5        -1.0*(W2CAMG)
#define GR5S1                0.25*((F_MN)+(GEXMGMN)+(GXMGMN) -2.0*(W2CAMG)+2.0*(W2CAMN)-2.0*(W2MGMN))
#define GR5S2                0.25*((F_FE)+(GEXMGFE)+(GXMGFE) -2.0*(W2CAMG)+2.0*(W2CAFE)-2.0*(W2MGFE))
#define GR5S3               -0.25*((F_CO)+(GEXMGCO)+(GXMGCO) -2.0*(W2CAMG)+2.0*(W2CACO)-2.0*(W2MGCO))
#define GR5S4               -0.25*((F_NI)+(GEXMGNI)+(GXMGNI) -2.0*(W2CAMG)+2.0*(W2CANI)-2.0*(W2MGNI))
#define GS1S1               0.25*((GXMGMN)-(W1MGMN)-(W2MGMN))
#define GS1S2            0.25*(( (GXMGMN)-(W1MGMN)-(W2MGMN)) +((GXMGFE)-(W1MGFE)-(W2MGFE)) -((GXMNFE)-(W1MNFE)-(W2MNFE)))
#define GS1S3            0.25*(-((GXMGMN)-(W1MGMN)-(W2MGMN)) -((GXMGCO)-(W1MGCO)-(W2MGCO)) +((GXMNCO)-(W1MNCO)-(W2MNCO)))
#define GS1S4            0.25*(-((GXMGMN)-(W1MGMN)-(W2MGMN)) -((GXMGNI)-(W1MGNI)-(W2MGNI)) +((GXMNNI)-(W1MNNI)-(W2MNNI)))
#define GS2S2             0.25*( (GXMGFE)-(W1MGFE)-(W2MGFE))
#define GS2S3            0.25*(-((GXMGFE)-(W1MGFE)-(W2MGFE)) -((GXMGCO)-(W1MGCO)-(W2MGCO)) +((GXFECO)-(W1FECO)-(W2FECO)))
#define GS2S4            0.25*(-((GXMGFE)-(W1MGFE)-(W2MGFE)) -((GXMGNI)-(W1MGNI)-(W2MGNI)) +((GXFENI)-(W1FENI)-(W2FENI)))
#define GS3S3             0.25*( (GXMGCO)-(W1MGCO)-(W2MGCO))
#define GS3S4            0.25*( ((GXMGCO)-(W1MGCO)-(W2MGCO)) +((GXMGNI)-(W1MGNI)-(W2MGNI)) -((GXCONI)-(W1CONI)-(W2CONI)))
#define GS4S4             0.25*( (GXMGNI)-(W1MGNI)-(W2MGNI))


#define V0  0.25*(     ((VXMNFE)+(WV1MNFE)+(WV2MNFE)) \
+((VXMNCO)+(WV1MNCO)+(WV2MNCO)) \
+((VXMNNI)+(WV1MNNI)+(WV2MNNI)) \
+((VXFECO)+(WV1FECO)+(WV2FECO)) \
+((VXFENI)+(WV1FENI)+(WV2FENI)) \
+((VXCONI)+(WV1CONI)+(WV2CONI)) \
-2.0*((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
-2.0*((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
-2.0*((VXMGCO)+(WV1MGCO)+(WV2MGCO)) \
-2.0*((VXMGNI)+(WV1MGNI)+(WV2MGNI)))
#define VR1 0.25*(-3.0*((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
-((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
-((VXMGCO)+(WV1MGCO)+(WV2MGCO)) \
-((VXMGNI)+(WV1MGNI)+(WV2MGNI)) \
+((VXMNFE)+(WV1MNFE)+(WV2MNFE)) \
+((VXMNCO)+(WV1MNCO)+(WV2MNCO)) \
+((VXMNNI)+(WV1MNNI)+(WV2MNNI)))
#define VR2  0.25*(   -((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
-3.0*((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
-((VXMGCO)+(WV1MGCO)+(WV2MGCO)) \
-((VXMGNI)+(WV1MGNI)+(WV2MGNI)) \
+((VXMNFE)+(WV1MNFE)+(WV2MNFE)) \
+((VXFECO)+(WV1FECO)+(WV2FECO)) \
+((VXFENI)+(WV1FENI)+(WV2FENI)))
#define VR3  0.25*(   -((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
-((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
-3.0*((VXMGCO)+(WV1MGCO)+(WV2MGCO)) \
-((VXMGNI)+(WV1MGNI)+(WV2MGNI)) \
+((VXMNCO)+(WV1MNCO)+(WV2MNCO)) \
+((VXFECO)+(WV1FECO)+(WV2FECO)) \
+((VXCONI)+(WV1CONI)+(WV2CONI)))
#define VR4 0.25*(    -((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
-((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
-((VXMGCO)+(WV1MGCO)+(WV2MGCO)) \
-3.0*((VXMGNI)+(WV1MGNI)+(WV2MGNI)) \
+((VXMNNI)+(WV1MNNI)+(WV2MNNI)) \
+((VXFENI)+(WV1FENI)+(WV2FENI)) \
+((VXCONI)+(WV1CONI)+(WV2CONI)))
#define VR5        -(WV2CAMG) \
-0.25*(((VEXMGMN)+(VXMGMN) \
-2.0*(WV2CAMN)+2.0*(WV2MGMN)) \
+((VEXMGFE)+(VXMGFE) \
-2.0*(WV2CAFE)+2.0*(WV2MGFE)) \
+((VEXMGCO)+(VXMGCO) \
-2.0*(WV2CACO)+2.0*(WV2MGCO)) \
+((VEXMGNI)+(VXMGNI) \
-2.0*(WV2CANI)+2.0*(WV2MGNI)))
#define VS1   0.25*(( (VEXMGMN)+3.0*(WV1MGMN)-3.0*(WV2MGMN)) \
-((VEXMGFE)-(WV1MGFE)+(WV2MGFE)) \
-((VEXMGCO)-(WV1MGCO)+(WV2MGCO)) \
-((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
+((VEXMNFE)-(WV1MNFE)+(WV2MNFE)) \
+((VEXMNCO)-(WV1MNCO)+(WV2MNCO)) \
+((VEXMNNI)-(WV1MNNI)+(WV2MNNI)))
#define VS2           0.25*((-(VEXMGMN)+(WV1MGMN)-(WV2MGMN)) \
+((VEXMGFE)+3.0*(WV1MGFE)-3.0*(WV2MGFE)) \
-((VEXMGCO)-(WV1MGCO)+(WV2MGCO)) \
-((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
-((VEXMNFE)+(WV1MNFE)-(WV2MNFE)) \
+((VEXFECO)-(WV1FECO)+(WV2FECO)) \
+((VEXFENI)-(WV1FENI)+(WV2FENI)))
#define VS3           0.25*(( (VEXMGMN)-(WV1MGMN)+(WV2MGMN)) \
+((VEXMGFE)-(WV1MGFE)+(WV2MGFE)) \
-((VEXMGCO)+3.0*(WV1MGCO)-3.0*(WV2MGCO)) \
+((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
+((VEXMNCO)+(WV1MNCO)-(WV2MNCO)) \
+((VEXFECO)+(WV1FECO)-(WV2FECO)) \
-((VEXCONI)-(WV1CONI)+(WV2CONI)))
#define VS4           0.25*(( (VEXMGMN)-(WV1MGMN)+(WV2MGMN)) \
+((VEXMGFE)-(WV1MGFE)+(WV2MGFE)) \
+((VEXMGCO)-(WV1MGCO)+(WV2MGCO)) \
-((VEXMGNI)+3.0*(WV1MGNI)-3.0*(WV2MGNI)) \
+((VEXMNNI)+(WV1MNNI)-(WV2MNNI)) \
+((VEXFENI)+(WV1FENI)-(WV2FENI)) \
+((VEXCONI)+(WV1CONI)-(WV2CONI)))
#define VR1R1         -0.25*(  (VXMGMN)+(WV1MGMN)+(WV2MGMN))
#define VR1R2          0.25*(( (VXMNFE)+(WV1MNFE)+(WV2MNFE)) \
-((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
-(( VXMGFE)+(WV1MGFE)+(WV2MGFE)))
#define VR1R3          0.25*(( (VXMNCO)+(WV1MNCO)+(WV2MNCO)) \
-((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
-((VXMGCO)+(WV1MGCO)+(WV2MGCO)))
#define VR1R4          0.25*(( (VXMNNI)+(WV1MNNI)+(WV2MNNI)) \
-((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
-((VXMGNI)+(WV1MGNI)+(WV2MGNI)))
#define VR1R5             -0.25*((VEXMGMN)+(VXMGMN) \
+2.0*(WV2CAMG)-2.0*(WV2CAMN)+2.0*(WV2MGMN))
#define VR1S1                       0.5*((WV1MGMN)-(WV2MGMN))
#define VR1S2         0.25*((-(VEXMGMN)-(WV1MGMN)+(WV2MGMN)) \
+((VEXMGFE)+(WV1MGFE)-(WV2MGFE)) \
-((VEXMNFE)+(WV1MNFE)-(WV2MNFE)))
#define VR1S3         0.25*(( (VEXMGMN)-(WV1MGMN)+(WV2MGMN)) \
-((VEXMGCO)+(WV1MGCO)-(WV2MGCO)) \
+((VEXMNCO)+(WV1MNCO)-(WV2MNCO)))
#define VR1S4         0.25*(( (VEXMGMN)-(WV1MGMN)+(WV2MGMN)) \
-((VEXMGNI)+(WV1MGNI)-(WV2MGNI)) \
+((VEXMNNI)+(WV1MNNI)-(WV2MNNI)))
#define VR2R2         -0.25*(( (VXMGFE)+(WV1MGFE)+(WV2MGFE)))
#define VR2R3          0.25*(( (VXFECO)+(WV1FECO)+(WV2FECO)) \
-((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
-((VXMGCO)+(WV1MGCO)+(WV2MGCO)))
#define VR2R4          0.25*(( (VXFENI)+(WV1FENI)+(WV2FENI)) \
-((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
-((VXMGNI)+(WV1MGNI)+(WV2MGNI)))
#define VR2R5             -0.25*((VEXMGFE)+(VXMGFE) \
+2.0*(WV2CAMG)-2.0*(WV2CAFE)+2.0*(WV2MGFE))
#define VR2S1         0.25*(( (VEXMGMN)+(WV1MGMN)-(WV2MGMN)) \
-((VEXMGFE)-(WV1MGFE)+(WV2MGFE)) \
+((VEXMNFE)-(WV1MNFE)+(WV2MNFE)))
#define VR2S2                       0.5*((WV1MGFE)-(WV2MGFE))
#define VR2S3         0.25*(( (VEXMGFE)-(WV1MGFE)+(WV2MGFE)) \
-((VEXMGCO)+(WV1MGCO)-(WV2MGCO)) \
+((VEXFECO)+(WV1FECO)-(WV2FECO)))
#define VR2S4         0.25*(( (VEXMGFE)-(WV1MGFE)+(WV2MGFE)) \
-((VEXMGNI)+(WV1MGNI)-(WV2MGNI)) \
+((VEXFENI)+(WV1FENI)-(WV2FENI)))
#define VR3R3          -0.25*(  (VXMGCO)+(WV1MGCO)+(WV2MGCO))
#define VR3R4          0.25*(( (VXCONI)+(WV1CONI)+(WV2CONI)) \
-((VXMGCO)+(WV1MGCO)+(WV2MGCO)) \
-((VXMGNI)+(WV1MGNI)+(WV2MGNI)))
#define VR3R5             -0.25*((VEXMGCO)+(VXMGCO) \
+2.0*(WV2CAMG)-2.0*(WV2CACO)+2.0*(WV2MGCO))
#define VR3S1         0.25*(( (VEXMGMN)+(WV1MGMN)-(WV2MGMN)) \
-((VEXMGCO)-(WV1MGCO)+(WV2MGCO)) \
+((VEXMNCO)-(WV1MNCO)+(WV2MNCO)))
#define VR3S2         0.25*(( (VEXMGFE)+(WV1MGFE)-(WV2MGFE)) \
-((VEXMGCO)-(WV1MGCO)+(WV2MGCO)) \
+((VEXFECO)-(WV1FECO)+(WV2FECO)))
#define VR3S3                      0.5*(-(WV1MGCO)+(WV2MGCO))
#define VR3S4          0.25*(( (VEXMGCO)+(WV1MGCO)-(WV2MGCO)) \
-((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
+((VEXCONI)-(WV1CONI)+(WV2CONI)))
#define VR4R4           -0.25*(  (VXMGNI)+(WV1MGNI)+(WV2MGNI))
#define VR4R5              -0.25*((VEXMGNI)+(VXMGNI) \
+2.0*(WV2CAMG)-2.0*(WV2CANI)+2.0*(WV2MGNI))
#define VR4S1          0.25*(( (VEXMGMN)+(WV1MGMN)-(WV2MGMN)) \
-((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
+((VEXMNNI)-(WV1MNNI)+(WV2MNNI)))
#define VR4S2          0.25*(( (VEXMGFE)+(WV1MGFE)-(WV2MGFE)) \
-((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
+((VEXFENI)-(WV1FENI)+(WV2FENI)))
#define VR4S3          0.25*((-(VEXMGCO)+(WV1MGCO)-(WV2MGCO)) \
+((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
-((VEXCONI)-(WV1CONI)+(WV2CONI)))
#define VR4S4                       0.5*(-(WV1MGNI)+(WV2MGNI))
#define VR5R5                                      -(WV2CAMG)
#define VR5S1                0.25*((VEXMGMN)+(VXMGMN) \
-2.0*(WV2CAMG)+2.0*(WV2CAMN)-2.0*(WV2MGMN))
#define VR5S2                0.25*((VEXMGFE)+(VXMGFE) \
-2.0*(WV2CAMG)+2.0*(WV2CAFE)-2.0*(WV2MGFE))
#define VR5S3               -0.25*((VEXMGCO)+(VXMGCO) \
-2.0*(WV2CAMG)+2.0*(WV2CACO)-2.0*(WV2MGCO))
#define VR5S4               -0.25*((VEXMGNI)+(VXMGNI) \
-2.0*(WV2CAMG)+2.0*(WV2CANI)-2.0*(WV2MGNI))
#define VS1S1               0.25*((VXMGMN)-(WV1MGMN)-(WV2MGMN))
#define VS1S2            0.25*(( (VXMGMN)-(WV1MGMN)-(WV2MGMN)) \
+((VXMGFE)-(WV1MGFE)-(WV2MGFE)) \
-((VXMNFE)-(WV1MNFE)-(WV2MNFE)))
#define VS1S3            0.25*(-((VXMGMN)-(WV1MGMN)-(WV2MGMN)) \
-((VXMGCO)-(WV1MGCO)-(WV2MGCO)) \
+((VXMNCO)-(WV1MNCO)-(WV2MNCO)))
#define VS1S4            0.25*(-((VXMGMN)-(WV1MGMN)-(WV2MGMN)) \
-((VXMGNI)-(WV1MGNI)-(WV2MGNI)) \
+((VXMNNI)-(WV1MNNI)-(WV2MNNI)))
#define VS2S2             0.25*(  (VXMGFE)-(WV1MGFE)-(WV2MGFE))
#define VS2S3            0.25*(-((VXMGFE)-(WV1MGFE)-(WV2MGFE)) \
-((VXMGCO)-(WV1MGCO)-(WV2MGCO)) \
+((VXFECO)-(WV1FECO)-(WV2FECO)))
#define VS2S4            0.25*(-((VXMGFE)-(WV1MGFE)-(WV2MGFE)) \
-((VXMGNI)-(WV1MGNI)-(WV2MGNI)) \
+((VXFENI)-(WV1FENI)-(WV2FENI)))
#define VS3S3             0.25*(  (VXMGCO)-(WV1MGCO)-(WV2MGCO))
#define VS3S4            0.25*( ((VXMGCO)-(WV1MGCO)-(WV2MGCO)) \
+((VXMGNI)-(WV1MGNI)-(WV2MGNI)) \
-((VXCONI)-(WV1CONI)-(WV2CONI)))
#define VS4S4             0.25*(  (VXMGNI)-(WV1MGNI)-(WV2MGNI))

/*
 * "Darken Equation" coefficients -
 * Global (to this file): activity definitions and component transforms
 *    The function conOlv defines the conversion from m[i], to r[j]
 */
/* Order R1,R2,R3,R4,R5 */
#define FR0(i)     (i == 0) ? 1.0 - r[0] : -(1.0 + r[0])
#define FR1(i)     (i == 1) ? 1.0 - r[1] : -(1.0 + r[1])
#define FR2(i)     (i == 2) ? 1.0 - r[2] : -(1.0 + r[2])
#define FR3(i)     (i == 3) ? 1.0 - r[3] : -(1.0 + r[3])
#define FR4(i)     (i == 4) ? 1.0 - r[4] : -r[4]

/* Order: S1, S2, S3, S4 */
#define GGS0(i)     -s[0]
#define GGS1(i)     -s[1]
#define GGS2(i)     -s[2]
#define GGS3(i)     -s[3]

#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0
#define DFR2DR2(i) - 1.0
#define DFR3DR3(i) - 1.0
#define DFR4DR4(i) - 1.0

#define DGS0DS0(i) - 1.0
#define DGS1DS1(i) - 1.0
#define DGS2DS2(i) - 1.0
#define DGS3DS3(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */

#define S -R*( xm1mn*log(xm1mn) + xm2mn*log(xm2mn) + xm1fe*log(xm1fe)+ \
xm2fe*log(xm2fe) + xm1co*log(xm1co) + xm2co*log(xm2co)+ \
xm1ni*log(xm1ni) + xm2ni*log(xm2ni) + xm2ca*log(xm2ca)+ \
xm1mg*log(xm1mg) + xm2mg*log(xm2mg))

/*  enthalpy here is enthalpy at P of interest */
#define H     (G0) + \
(GR1)*r[0] + (GR2)*r[1] + (GR3)*r[2] + (GR4)*r[3] + \
(GR5)*r[4] + (GS1)*s[0] + (GS2)*s[1] + (GS3)*s[2] + \
(GS4)*s[3] + \
(GR1R1)*r[0]*r[0] + (GR1R2)*r[0]*r[1] + (GR1R3)*r[0]*r[2] + \
(GR1R4)*r[0]*r[3] + (GR1R5)*r[0]*r[4] + (GR1S1)*r[0]*s[0] + \
(GR1S2)*r[0]*s[1] + (GR1S3)*r[0]*s[2] + (GR1S4)*r[0]*s[3] + \
(GR2R2)*r[1]*r[1] + (GR2R3)*r[1]*r[2] + (GR2R4)*r[1]*r[3] + \
(GR2R5)*r[1]*r[4] + (GR2S1)*r[1]*s[0] + (GR2S2)*r[1]*s[1] + \
(GR2S3)*r[1]*s[2] + (GR2S4)*r[1]*s[3] + (GR3R3)*r[2]*r[2] + \
(GR3R4)*r[2]*r[3] + (GR3R5)*r[2]*r[4] + (GR3S1)*r[2]*s[0] + \
(GR3S2)*r[2]*s[1] + (GR3S3)*r[2]*s[2] + (GR3S4)*r[2]*s[3] + \
(GR4R4)*r[3]*r[3] + (GR4R5)*r[3]*r[4] + (GR4S1)*r[3]*s[0] + \
(GR4S2)*r[3]*s[1] + (GR4S3)*r[3]*s[2] + (GR4S4)*r[3]*s[3] + \
(GR5R5)*r[4]*r[4] + (GR5S1)*r[4]*s[0] + (GR5S2)*r[4]*s[1] + \
(GR5S3)*r[4]*s[2] + (GR5S4)*r[4]*s[3] + (GS1S1)*s[0]*s[0] + \
(GS1S2)*s[0]*s[1] + (GS1S3)*s[0]*s[2] + (GS1S4)*s[0]*s[3] + \
(GS2S2)*s[1]*s[1] + (GS2S3)*s[1]*s[2] + (GS2S4)*s[1]*s[3] + \
(GS3S3)*s[2]*s[2] + (GS3S4)*s[2]*s[3] + (GS4S4)*s[3]*s[3]
#define V    (V0) + \
(VR1)*r[0] + (VR2)*r[1] + (VR3)*r[2] + (VR4)*r[3] + \
(VR5)*r[4] + (VS1)*s[0] + (VS2)*s[1] + (VS3)*s[2] + \
(VS4)*s[3] + \
(VR1R1)*r[0]*r[0] + (VR1R2)*r[0]*r[1] + (VR1R3)*r[0]*r[2] + \
(VR1R4)*r[0]*r[3] + (VR1R5)*r[0]*r[4] + (VR1S1)*r[0]*s[0] + \
(VR1S2)*r[0]*s[1] + (VR1S3)*r[0]*s[2] + (VR1S4)*r[0]*s[3] + \
(VR2R2)*r[1]*r[1] + (VR2R3)*r[1]*r[2] + (VR2R4)*r[1]*r[3] + \
(VR2R5)*r[1]*r[4] + (VR2S1)*r[1]*s[0] + (VR2S2)*r[1]*s[1] + \
(VR2S3)*r[1]*s[2] + (VR2S4)*r[1]*s[3] + (VR3R3)*r[2]*r[2] + \
(VR3R4)*r[2]*r[3] + (VR3R5)*r[2]*r[4] + (VR3S1)*r[2]*s[0] + \
(VR3S2)*r[2]*s[1] + (VR3S3)*r[2]*s[2] + (VR3S4)*r[2]*s[3] + \
(VR4R4)*r[3]*r[3] + (VR4R5)*r[3]*r[4] + (VR4S1)*r[3]*s[0] + \
(VR4S2)*r[3]*s[1] + (VR4S3)*r[3]*s[2] + (VR4S4)*r[3]*s[3] + \
(VR5R5)*r[4]*r[4] + (VR5S1)*r[4]*s[0] + (VR5S2)*r[4]*s[1] + \
(VR5S3)*r[4]*s[2] + (VR5S4)*r[4]*s[3] + (VS1S1)*s[0]*s[0] + \
(VS1S2)*s[0]*s[1] + (VS1S3)*s[0]*s[2] + (VS1S4)*s[0]*s[3] + \
(VS2S2)*s[1]*s[1] + (VS2S3)*s[1]*s[2] + (VS2S4)*s[1]*s[3] + \
(VS3S3)*s[2]*s[2] + (VS3S4)*s[2]*s[3] + (VS4S4)*s[3]*s[3]

#define G    (H) - t*(S)

/*----------------------------------------------------------------------------*/

#define DGDR0  (GR1) + 2.0*(GR1R1)*r[0] + \
(GR1R2)*r[1] + (GR1R3)*r[2] + (GR1R4)*r[3] + (GR1R5)*r[4] + \
(GR1S1)*s[0] + (GR1S2)*s[1] + (GR1S3)*s[2] + (GR1S4)*s[3] + \
0.5*R*t*(log(xm1mn*xm2mn/xm1mg/xm2mg))
#define DGDR1  (GR2) + 2.0*(GR2R2)*r[1] + \
(GR1R2)*r[0] + (GR2R3)*r[2] + (GR2R4)*r[3] + (GR2R5)*r[4] + \
(GR2S1)*s[0] + (GR2S2)*s[1] + (GR2S3)*s[2] + (GR2S4)*s[3] + \
0.5*R*t*(log(xm1fe*xm2fe/xm1mg/xm2mg))
#define DGDR2  (GR3) + 2.0*(GR3R3)*r[2] + \
(GR1R3)*r[0] + (GR2R3)*r[1] + (GR3R4)*r[3] +(GR3R5)*r[4] + \
(GR3S1)*s[0] + (GR3S2)*s[1] + (GR3S3)*s[2] +(GR3S4)*s[3] + \
0.5*R*t*(log(xm1co*xm2co/xm1mg/xm2mg))
#define DGDR3  (GR4) + 2.0*(GR4R4)*r[3] + \
(GR1R4)*r[0] + (GR2R4)*r[1] + (GR3R4)*r[2] +(GR4R5)*r[4] + \
(GR4S1)*s[0] + (GR4S2)*s[1] + (GR4S3)*s[2] +(GR4S4)*s[3] + \
0.5*R*t*(log(xm1ni*xm2ni/xm1mg/xm2mg))
#define DGDR4  (GR5) + 2.0*(GR5R5)*r[4] + \
(GR1R5)*r[0] + (GR2R5)*r[1] + (GR3R5)*r[2] + (GR4R5)*r[3] + \
(GR5S1)*s[0] + (GR5S2)*s[1] + (GR5S3)*s[2] + (GR5S4)*s[3] + \
R*t*(log(xm2ca/xm2mg))
#define DGDS0  (GS1) + 2.0*(GS1S1)*s[0] + \
(GR1S1)*r[0] + (GR2S1)*r[1] + (GR3S1)*r[2] + (GR4S1)*r[3] + \
(GR5S1)*r[4] + (GS1S2)*s[1] + (GS1S3)*s[2] + (GS1S4)*s[3] + \
0.5*R*t*(log(xm2mn*xm1mg/xm1mn/xm2mg))
#define DGDS1  (GS2) + 2.0*(GS2S2)*s[1] + \
(GR1S2)*r[0] + (GR2S2)*r[1] + (GR3S2)*r[2] + (GR4S2)*r[3] + \
(GR5S2)*r[4] + (GS1S2)*s[0] + (GS2S3)*s[2] + (GS2S4)*s[3] + \
0.5*R*t*(log(xm2fe*xm1mg/xm1fe/xm2mg))
#define DGDS2  (GS3) + 2.0*(GS3S3)*s[2] + \
(GR1S3)*r[0] + (GR2S3)*r[1] + (GR3S3)*r[2] + (GR4S3)*r[3] + \
(GR5S3)*r[4] + (GS1S3)*s[0] + (GS2S3)*s[1] + (GS3S4)*s[3] + \
0.5*R*t*(log(xm1co*xm2mg/xm2co/xm1mg))
#define DGDS3  (GS4) + 2.0*(GS4S4)*s[3] + \
(GR1S4)*r[0] + (GR2S4)*r[1] + (GR3S4)*r[2] + (GR4S4)*r[3] + \
(GR5S4)*r[4] + (GS1S4)*s[0] + (GS2S4)*s[1] + (GS3S4)*s[2] + \
0.5*R*t*(log(xm1ni*xm2mg/xm2ni/xm1mg))
#define DGDT  (S)
#define DGDP  (V)

/*----------------------------------------------------------------------------*/

#define D2GDR0R0 2.0*(GR1R1) + 0.25*R*t*(\
1.0/xm1mn + 1.0/xm2mn + 1.0/xm1mg + 1.0/xm2mg)
#define D2GDR0R1 (GR1R2) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDR0R2 (GR1R3) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDR0R3 (GR1R4) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDR0R4 (GR1R5) + 0.50*R*t*(1.0/xm2mg)
#define D2GDR0S0 (GR1S1) +  0.25*R*t*( \
-1.0/xm1mn + 1.0/xm2mn - 1.0/xm1mg + 1.0/xm2mg)
#define D2GDR0S1 (GR1S2) + 0.25*R*t*(- 1.0/xm1mg + 1.0/xm2mg)
#define D2GDR0S2 (GR1S3) + 0.25*R*t*(  1.0/xm1mg - 1.0/xm2mg)
#define D2GDR0S3 (GR1S4) + 0.25*R*t*(  1.0/xm1mg - 1.0/xm2mg)
#define D2GDR0DT 0.5*R*(log(xm1mn*xm2mn/xm1mg/xm2mg))
#define D2GDR0DP (VR1) + 2.0*(VR1R1)*r[0] + \
(VR1R2)*r[1] + (VR1R3)*r[2] + (VR1R4)*r[3] + (VR1R5)*r[4] + \
(VR1S1)*s[0] + (VR1S2)*s[1] + (VR1S3)*s[2] + (VR1S4)*s[3]

#define D2GDR1R1 2.0*(GR2R2) + 0.25*R*t* \
(1.0/xm1fe + 1.0/xm2fe + 1.0/xm1mg + 1.0/xm2mg)
#define D2GDR1R2 (GR2R3) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDR1R3 (GR2R4) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDR1R4 (GR2R5) + 0.5*R*t*(1.0/xm2mg)
#define D2GDR1S0 (GR2S1) + 0.25*R*t*(1.0/xm2mg - 1.0/xm1mg)
#define D2GDR1S1 (GR2S2) + 0.25*R*t*( \
-1.0/xm1fe + 1.0/xm2fe - 1.0/xm1mg + 1.0/xm2mg)
#define D2GDR1S2 (GR2S3) + 0.25*R*t*(  1.0/xm1mg - 1.0/xm2mg)
#define D2GDR1S3 (GR2S4) + 0.25*R*t*(  1.0/xm1mg - 1.0/xm2mg)
#define D2GDR1DT 0.5*R*(log(xm1fe*xm2fe/xm1mg/xm2mg))
#define D2GDR1DP (VR2) + 2.0*(VR2R2)*r[1] + \
(VR1R2)*r[0] + (VR2R3)*r[2] + (VR2R4)*r[3] + (VR2R5)*r[4] + \
(VR2S1)*s[0] + (VR2S2)*s[1] + (VR2S3)*s[2] + (VR2S4)*s[3]

#define D2GDR2R2 2.0*(GR3R3) + 0.25*R*t*( \
1.0/xm1co + 1.0/xm2co + 1.0/xm1mg + 1.0/xm2mg)
#define D2GDR2R3 (GR3R4) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDR2R4 (GR3R5) + 0.5*R*t*(1.0/xm2mg)
#define D2GDR2S0 (GR3S1) + 0.25*R*t*(1.0/xm2mg-1.0/xm1mg)
#define D2GDR2S1 (GR3S2) + 0.25*R*t*(1.0/xm2mg-1.0/xm1mg)
#define D2GDR2S2 (GR3S3) + 0.25*R*t*( \
1.0/xm1co - 1.0/xm2co + 1.0/xm1mg - 1.0/xm2mg)
#define D2GDR2S3 (GR3S4) + 0.25*R*t*(1.0/xm1mg-1.0/xm2mg)
#define D2GDR2DT 0.5*R*(log(xm1co*xm2co/xm1mg/xm2mg))
#define D2GDR2DP (VR3) + 2.0*(VR3R3)*r[2] + \
(VR1R3)*r[0] + (VR2R3)*r[1] + (VR3R4)*r[3] +(VR3R5)*r[4] + \
(VR3S1)*s[0] + (VR3S2)*s[1] + (VR3S3)*s[2] +(VR3S4)*s[3]

#define D2GDR3R3 2.0*(GR4R4) + 0.25*R*t*( \
1.0/xm1mg + 1.0/xm2mg + 1.0/xm1ni + 1.0/xm2ni)
#define D2GDR3R4 (GR4R5) + 0.5*R*t*(1.0/xm2mg)
#define D2GDR3S0 (GR4S1) + 0.25*R*t*(1.0/xm2mg-1.0/xm1mg)
#define D2GDR3S1 (GR4S2) + 0.25*R*t*(1.0/xm2mg-1.0/xm1mg)
#define D2GDR3S2 (GR4S3) + 0.25*R*t*(1.0/xm1mg-1.0/xm2mg)
#define D2GDR3S3 (GR4S4) + 0.25*R*t*( \
1.0/xm1ni - 1.0/xm2ni + 1.0/xm1mg - 1.0/xm2mg)
#define D2GDR3DT 0.5*R*(log(xm1ni*xm2ni/xm1mg/xm2mg))
#define D2GDR3DP (VR4) + 2.0*(VR4R4)*r[3] + \
(VR1R4)*r[0] + (VR2R4)*r[1] + (VR3R4)*r[2] +(VR4R5)*r[4] + \
(VR4S1)*s[0] + (VR4S2)*s[1] + (VR4S3)*s[2] +(VR4S4)*s[3]

#define D2GDR4R4 2.0*(GR5R5) + R*t*(1.0/xm2mg + 1.0/xm2ca)
#define D2GDR4S0 (GR5S1) + 0.5*R*t*(1.0/xm2mg)
#define D2GDR4S1 (GR5S2) + 0.5*R*t*(1.0/xm2mg)
#define D2GDR4S2 (GR5S3) - 0.5*R*t*(1.0/xm2mg)
#define D2GDR4S3 (GR5S4) - 0.5*R*t*(1.0/xm2mg)
#define D2GDR4DT R*log(xm2ca/xm2mg)
#define D2GDR4DP (VR5) + 2.0*(VR5R5)*r[4] + \
(VR1R5)*r[0] + (VR2R5)*r[1] + (VR3R5)*r[2] + (VR4R5)*r[3] + \
(VR5S1)*s[0] + (VR5S2)*s[1] + (VR5S3)*s[2] + (VR5S4)*s[3]

#define D2GDS0S0 2.0*(GS1S1) + \
0.25*R*t*(1.0/xm1mn + 1.0/xm2mn + 1.0/xm1mg + 1.0/xm2mg)
#define D2GDS0S1 (GS1S2) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDS0S2 (GS1S3) - 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDS0S3 (GS1S4) - 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDS0DT 0.5*R*(log(xm2mn*xm1mg/xm1mn/xm2mg))
#define D2GDS0DP (VS1) + 2.0*(VS1S1)*s[0] + \
(VR1S1)*r[0] + (VR2S1)*r[1] + (VR3S1)*r[2] + (VR4S1)*r[3] + \
(VR5S1)*r[4] + (VS1S2)*s[1] + (VS1S3)*s[2] + (VS1S4)*s[3]

#define D2GDS1S1 2.0*(GS2S2) + \
0.25*R*t*(1.0/xm1fe + 1.0/xm2fe + 1.0/xm1mg +1.0/xm2mg)
#define D2GDS1S2 (GS2S3) - 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDS1S3 (GS2S4) - 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDS1DT 0.5*R*(log(xm2fe*xm1mg/xm1fe/xm2mg))
#define D2GDS1DP (VS2) + 2.0*(VS2S2)*s[1] + \
(VR1S2)*r[0] + (VR2S2)*r[1] + (VR3S2)*r[2] + (VR4S2)*r[3] + \
(VR5S2)*r[4] + (VS1S2)*s[0] + (VS2S3)*s[2] + (VS2S4)*s[3]

#define D2GDS2S2 2.0*(GS3S3) + \
0.25*R*t*(1.0/xm1co + 1.0/xm2co + 1.0/xm1mg +1.0/xm2mg)
#define D2GDS2S3 (GS3S4) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDS2DT 0.5*R*(log(xm1co*xm2mg/xm2co/xm1mg))
#define D2GDS2DP (VS3) + 2.0*(VS3S3)*s[2] + \
(VR1S3)*r[0] + (VR2S3)*r[1] + (VR3S3)*r[2] + (VR4S3)*r[3] + \
(VR5S3)*r[4] + (VS1S3)*s[0] + (VS2S3)*s[1] + (VS3S4)*s[3]

#define D2GDS3S3 2.0*(GS4S4) + \
0.25*R*t*(1.0/xm1ni + 1.0/xm2ni + 1.0/xm1mg +1.0/xm2mg)
#define D2GDS3DT 0.5*R*(log(xm1ni*xm2mg/xm2ni/xm1mg))
#define D2GDS3DP (VS4) + 2.0*(VS4S4)*s[3] + \
(VR1S4)*r[0] + (VR2S4)*r[1] + (VR3S4)*r[2] + (VR4S4)*r[3] + \
(VR5S4)*r[4] + (VS1S4)*s[0] + (VS2S4)*s[1] + (VS3S4)*s[2]

#define D2GDT2   0.0
#define D2GDTDP  0.0
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/
#define D3GDR0R0R0   0.125*R*t*(-1.0/SQUARE(xm1mn)+1.0/SQUARE(xm2mg) \
-1.0/SQUARE(xm2mn)+1.0/SQUARE(xm1mg))
#define D3GDR0R0R1   0.125*R*t*(1.0/SQUARE(xm1mg)+1.0/SQUARE(xm2mg))
#define D3GDR0R0S0   0.125*R*t*(1.0/SQUARE(xm1mn)+1.0/SQUARE(xm2mg) \
-1.0/SQUARE(xm2mn)-1.0/SQUARE(xm1mg))
#define D3GDR0R0S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R0S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R0S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R0DT   0.25*R*( \
1.0/xm1mn + 1.0/xm2mn + 1.0/xm1mg + 1.0/xm2mg)
#define D3GDR0R0DP   2.0*(VR1R1)

#define D3GDR0R1S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R1S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R1S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R1S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R1DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDR0R1DP   (VR1R2)

#define D3GDR0R2S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R2S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R2S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R2S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R2DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDR0R2DP   (VR1R3)

#define D3GDR0R3S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R3S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R3S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R3S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R3DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDR0R3DP   (VR1R4)

#define D3GDR0R4S0   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR0R4S1   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR0R4S2  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR0R4S3  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR0R4DT   0.5 *R*(1.0/xm2mg)
#define D3GDR0R4DP   (VR1R5)

#define D3GDR0S0S0   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg) \
-1.0/SQUARE(xm2mn)-1.0/SQUARE(xm1mn))
#define D3GDR0S0S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S0S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S0S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S0DT   0.25*R*( \
-1.0/xm1mn + 1.0/xm2mn - 1.0/xm1mg + 1.0/xm2mg)
#define D3GDR0S0DP   (VR1S1)

#define D3GDR0S1S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S1S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S1S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S1DT   0.25*R*(1.0/xm2mg-1.0/xm1mg)
#define D3GDR0S1DP   (VR1S2)

#define D3GDR0S2S2   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S2S3   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S2DT   0.25*R*(1.0/xm1mg-1.0/xm2mg)
#define D3GDR0S2DP   (VR1S3)

#define D3GDR0S3S3   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S3DT   0.25*R*(  1.0/xm1mg - 1.0/xm2mg)

#define D3GDR0S3DP   (VR1S4)

#define D3GDR1R1R1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg) \
-1.0/SQUARE(xm1fe)-1.0/SQUARE(xm2fe))
#define D3GDR1R1S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR1R1S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg) \
+1.0/SQUARE(xm1fe)-1.0/SQUARE(xm2fe))
#define D3GDR1R1S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR1R1S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR1R1DT   0.25*R*( \
1.0/xm1fe + 1.0/xm2fe + 1.0/xm1mg + 1.0/xm2mg)
#define D3GDR1R1DP   2.0*(VR2R2)

#define D3GDR1R2S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR1R2S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR1R2S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR1R2S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR1R2DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDR1R2DP   (VR2R3)

#define D3GDR1R3S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR1R3S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR1R3S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR1R3S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR1R3DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDR1R3DP   (VR2R4)

#define D3GDR1R4S0   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR1R4S1   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR1R4S2  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR1R4S3  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR1R4DT   0.5*R* (1.0/xm2mg)
#define D3GDR1R4DP   (VR2R5)

#define D3GDR1S0S0   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S0S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S0S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S0S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S0DT   0.25*R*(1.0/xm2mg - 1.0/xm1mg)
#define D3GDR1S0DP   (VR2S1)

#define D3GDR1S1S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg) \
-1.0/SQUARE(xm2fe)-1.0/SQUARE(xm1fe))
#define D3GDR1S1S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S1S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S1DT   0.25*R*( \
-1.0/xm1fe + 1.0/xm2fe - 1.0/xm1mg + 1.0/xm2mg)
#define D3GDR1S1DP   (VR2S2)

#define D3GDR1S2S2   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S2S3   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S2DT   0.25*R*(  1.0/xm1mg - 1.0/xm2mg)
#define D3GDR1S2DP   (VR2S3)

#define D3GDR1S3S3   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S3DT   0.25*R*(  1.0/xm1mg - 1.0/xm2mg)
#define D3GDR1S3DP   (VR2S4)

#define D3GDR2R2R2   0.125*R*t*(1.0/SQUARE(xm1mg)+1.0/SQUARE(xm2mg) \
-1.0/SQUARE(xm2co)-1.0/SQUARE(xm1co))
#define D3GDR2R2S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR2R2S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR2R2S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg) \
+1.0/SQUARE(xm2co)-1.0/SQUARE(xm1co))
#define D3GDR2R2S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR2R2DT   0.25*R*( \
1.0/xm1co + 1.0/xm2co + 1.0/xm1mg + 1.0/xm2mg)
#define D3GDR2R2DP   2.0*(VR3R3)

#define D3GDR2R3S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR2R3S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR2R3S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR2R3S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR2R3DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDR2R3DP   (VR3R4)

#define D3GDR2R4S0   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR2R4S1   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR2R4S2  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR2R4S3  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR2R4DT   0.5*R* (1.0/xm2mg)
#define D3GDR2R4DP   (VR3R5)

#define D3GDR2S0S0   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S0S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S0S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S0S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S0DT   0.25*R*(1.0/xm2mg-1.0/xm1mg)
#define D3GDR2S0DP   (VR3S1)

#define D3GDR2S1S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S1S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S1S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S1DT   0.25*R*(1.0/xm2mg-1.0/xm1mg)
#define D3GDR2S1DP   (VR3S2)

#define D3GDR2S2S2   0.125*R*t*(1.0/SQUARE(xm1mg)+1.0/SQUARE(xm2mg) \
-1.0/SQUARE(xm1co)-1.0/SQUARE(xm2co))
#define D3GDR2S2S3   0.125*R*t*(1.0/SQUARE(xm1mg)+1.0/SQUARE(xm2mg))
#define D3GDR2S2DT   0.25*R*( \
1.0/xm1co - 1.0/xm2co + 1.0/xm1mg - 1.0/xm2mg)
#define D3GDR2S2DP   (VR3S3)

#define D3GDR2S3S2   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S3S3   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S3DT   0.25*R*(1.0/xm1mg-1.0/xm2mg)
#define D3GDR2S3DP   (VR3S4)

#define D3GDR3R3R3   0.125*R*t*(1.0/SQUARE(xm1mg)+1.0/SQUARE(xm2mg) \
-1.0/SQUARE(xm2ni)-1.0/SQUARE(xm1ni))
#define D3GDR3R3S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR3R3S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR3R3S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR3R3S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg) \
+1.0/SQUARE(xm2ni)-1.0/SQUARE(xm1ni))
#define D3GDR3R3DT   0.25*R*( \
1.0/xm1mg + 1.0/xm2mg + 1.0/xm1ni + 1.0/xm2ni)
#define D3GDR3R3DP   2.0*(VR4R4)

#define D3GDR3R4S0   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR3R4S1   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR3R4S2  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR3R4S3  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR3R4DT   0.5*R* (1.0/xm2mg)
#define D3GDR3R4DP   (VR4R5)

#define D3GDR3S0S0   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S0S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S0S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S0S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S0DT   0.25*R*(1.0/xm2mg-1.0/xm1mg)
#define D3GDR3S0DP   (VR4S1)

#define D3GDR3S1S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S1S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S1S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S1DT   0.25*R*(1.0/xm2mg-1.0/xm1mg)
#define D3GDR3S1DP   (VR4S2)

#define D3GDR3S2S2   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S2S3   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S2DT   0.25*R*(1.0/xm1mg-1.0/xm2mg)
#define D3GDR3S2DP   (VR4S3)


#define D3GDR3S3S3   0.125*R*t*(1.0/SQUARE(xm1mg)+1.0/SQUARE(xm2mg) \
-1.0/SQUARE(xm1ni)-1.0/SQUARE(xm2ni))
#define D3GDR3S3DT   0.25*R*( \
1.0/xm1ni - 1.0/xm2ni + 1.0/xm1mg - 1.0/xm2mg)
#define D3GDR3S3DP   (VR4S4)

#define D3GDR4R4R4   R*t*(1.0/SQUARE(xm2mg) - 1.0/SQUARE(xm2ca))
#define D3GDR4R4R0   0.50*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4R0R0   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4R4S0   0.50*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4R4S1   0.50*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4R4S2  -0.50*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4R4S3  -0.50*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4R4DT   R*(1.0/xm2mg + 1.0/xm2ca)
#define D3GDR4R4DP   2.0*(VR5R5)

#define D3GDR4S0S0   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S0S1   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S0S2  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S0S3  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S0DT   0.5*R*(1.0/xm2mg)
#define D3GDR4S0DP   (VR5S1)

#define D3GDR4S1S1   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S1S2  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S1S3  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S1DT   0.5*R*(1.0/xm2mg)
#define D3GDR4S1DP   (VR5S2)

#define D3GDR4S2S2   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S2S3   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S2DT  -0.5*R*(1.0/xm2mg)
#define D3GDR4S2DP   (VR5S3)

#define D3GDR4S3S3   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S3DT  -0.5*R*(1.0/xm2mg)
#define D3GDR4S3DP   (VR5S4)

#define D3GDS0S0S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg) \
+1.0/SQUARE(xm1mn)-1.0/SQUARE(xm2mn))
#define D3GDS0S0S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS0S0S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS0S0S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS0S0DT   0.25*R*(1.0/xm1mn + 1.0/xm2mn + 1.0/xm1mg +1.0/xm2mg)
#define D3GDS0S0DP   2.0*(VS1S1)

#define D3GDS0S1S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS0S1S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS0S1S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS0S1DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDS0S1DP   (VS1S2)

#define D3GDS0S2S2   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS0S2S3   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS0S2DT  -0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDS0S2DP   (VS1S3)

#define D3GDS0S3S3   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS0S3DT  -0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDS0S3DP   (VS1S4)
#define D3GDS0DT2    0.0
#define D3GDS0DTDP   0.0
#define D3GDS0DP2    0.0

#define D3GDS1S1S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg) \
+1.0/SQUARE(xm1fe)-1.0/SQUARE(xm2fe))
#define D3GDS1S1S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS1S1S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS1S1DT   0.25*R*(1.0/xm1fe + 1.0/xm2fe + 1.0/xm1mg +1.0/xm2mg)
#define D3GDS1S1DP   2.0*(VS2S2)

#define D3GDS1S2S2   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS1S2S3   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS1S2DT  -0.25*R*   (1.0/xm1mg + 1.0/xm2mg)
#define D3GDS1S2DP   (VS2S3)

#define D3GDS1S3S3   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS1S3DT  -0.25*R*   (1.0/xm1mg + 1.0/xm2mg)
#define D3GDS1S3DP   (VS2S4)
#define D3GDS1DT2    0.0
#define D3GDS1DTDP   0.0
#define D3GDS1DP2    0.0

#define D3GDS2S2S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg) \
+1.0/SQUARE(xm2co)-1.0/SQUARE(xm1co))
#define D3GDS2S2S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS2S2DT   0.25*R*(1.0/xm1co + 1.0/xm2co + 1.0/xm1mg +1.0/xm2mg)
#define D3GDS2S2DP   2.0*(VS3S3)

#define D3GDS2S3S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS2S3DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDS2S3DP   (VS3S4)
#define D3GDS2DT2    0.0
#define D3GDS2DTDP   0.0
#define D3GDS2DP2    0.0

#define D3GDS3S3S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg) \
+1.0/SQUARE(xm2ni)-1.0/SQUARE(xm1ni))
#define D3GDS3S3DT   0.25*R*(1.0/xm1ni + 1.0/xm2ni + 1.0/xm1mg +1.0/xm2mg)
#define D3GDS3S3DP   2.0*(VS4S4)
#define D3GDS3DT2    0.0
#define D3GDS3DTDP   0.0
#define D3GDS3DP2    0.0

#define D3GDT3       0.0
#define D3GDT2DP     0.0
#define D3GDTDP2     0.0
#define D3GDP3       0.0

#define D3GDR0DT2  0.0  /* imported from spinel.c v2.0-2 */
#define D3GDR0DTDP 0.0
#define D3GDR0DP2  0.0
#define D3GDR1DT2  0.0
#define D3GDR1DTDP 0.0
#define D3GDR1DP2  0.0
#define D3GDR2DT2  0.0
#define D3GDR2DTDP 0.0
#define D3GDR2DP2  0.0
#define D3GDR3DT2  0.0
#define D3GDR3DTDP 0.0
#define D3GDR3DP2  0.0
#define D3GDR4DT2  0.0
#define D3GDR4DTDP 0.0
#define D3GDR4DP2  0.0
/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
d2gdr2[0][0] = (D2GDR0R0);     d2gdr2[0][1] = (D2GDR0R1); \
d2gdr2[0][2] = (D2GDR0R2);     d2gdr2[0][3] = (D2GDR0R3); \
d2gdr2[0][4] = (D2GDR0R4); \
d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = (D2GDR1R1); \
d2gdr2[1][2] = (D2GDR1R2);     d2gdr2[1][3] = (D2GDR1R3); \
d2gdr2[1][4] = (D2GDR1R4); \
d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2]; \
d2gdr2[2][2] = (D2GDR2R2);     d2gdr2[2][3] = (D2GDR2R3); \
d2gdr2[2][4] = (D2GDR2R4); \
d2gdr2[3][0] = d2gdr2[0][3]; d2gdr2[3][1] = d2gdr2[1][3]; \
d2gdr2[3][2] = d2gdr2[2][3]; d2gdr2[3][3] = (D2GDR3R3); \
d2gdr2[3][4] = (D2GDR3R4); \
d2gdr2[4][0] = d2gdr2[0][4]; d2gdr2[4][1] = d2gdr2[1][4]; \
d2gdr2[4][2] = d2gdr2[2][4]; d2gdr2[4][3] = d2gdr2[3][4]; \
d2gdr2[4][4] = (D2GDR4R4);

#define fillD2GDRDS \
d2gdrds[0][0] = (D2GDR0S0); d2gdrds[0][1] = (D2GDR0S1); \
d2gdrds[0][2] = (D2GDR0S2); d2gdrds[0][3] = (D2GDR0S3); \
d2gdrds[1][0] = (D2GDR1S0); d2gdrds[1][1] = (D2GDR1S1); \
d2gdrds[1][2] = (D2GDR1S2); d2gdrds[1][3] = (D2GDR1S3); \
d2gdrds[2][0] = (D2GDR2S0); d2gdrds[2][1] = (D2GDR2S1); \
d2gdrds[2][2] = (D2GDR2S2); d2gdrds[2][3] = (D2GDR2S3); \
d2gdrds[3][0] = (D2GDR3S0); d2gdrds[3][1] = (D2GDR3S1); \
d2gdrds[3][2] = (D2GDR3S2); d2gdrds[3][3] = (D2GDR3S3); \
d2gdrds[4][0] = (D2GDR4S0); d2gdrds[4][1] = (D2GDR4S1); \
d2gdrds[4][2] = (D2GDR4S2); d2gdrds[4][3] = (D2GDR4S3);

#define fillD2GDRDT \
d2gdrdt[0] = D2GDR0DT; d2gdrdt[1] = D2GDR1DT; d2gdrdt[2] = D2GDR2DT; \
d2gdrdt[3] = D2GDR3DT; d2gdrdt[4] = D2GDR4DT;

#define fillD2GDRDP \
d2gdrdp[0] = D2GDR0DP; d2gdrdp[1] = D2GDR1DP; d2gdrdp[2] = D2GDR2DP; \
d2gdrdp[3] = D2GDR3DP; d2gdrdp[4] = D2GDR4DP;

#define fillD2GDS2 \
d2gds2[0][0] = (D2GDS0S0);     d2gds2[0][1] = (D2GDS0S1); \
d2gds2[0][2] = (D2GDS0S2);     d2gds2[0][3] = (D2GDS0S3); \
d2gds2[1][0] = d2gds2[0][1]; d2gds2[1][1] = (D2GDS1S1); \
d2gds2[1][2] = (D2GDS1S2);     d2gds2[1][3] = (D2GDS1S3); \
d2gds2[2][0] = d2gds2[0][2]; d2gds2[2][1] = d2gds2[1][2]; \
d2gds2[2][2] = (D2GDS2S2);     d2gds2[2][3] = (D2GDS2S3); \
d2gds2[3][0] = d2gds2[0][3]; d2gds2[3][1] = d2gds2[1][3]; \
d2gds2[3][2] = d2gds2[2][3]; d2gds2[3][3] = (D2GDS3S3);

#define fillD2GDSDT \
d2gdsdt[0] = D2GDS0DT;  d2gdsdt[1] = D2GDS1DT; \
d2gdsdt[2] = D2GDS2DT;  d2gdsdt[3] = D2GDS3DT;

#define fillD2GDSDP \
d2gdsdp[0] = D2GDS0DP;  d2gdsdp[1] = D2GDS1DP; \
d2gdsdp[2] = D2GDS2DP;  d2gdsdp[3] = D2GDS3DP;

#define fillD3GDR3 \
d3gdr3[0][0][0] = D3GDR0R0R0;          d3gdr3[0][0][1] = D3GDR0R0R1; \
d3gdr3[0][0][2] = d3gdr3[0][0][1];     d3gdr3[0][0][3] = d3gdr3[0][0][1]; \
d3gdr3[0][0][4] = D3GDR4R0R0;          d3gdr3[0][1][0] = d3gdr3[0][0][1]; \
d3gdr3[0][1][1] = d3gdr3[0][0][1];     d3gdr3[0][1][2] = d3gdr3[0][0][1]; \
d3gdr3[0][1][3] = d3gdr3[0][0][1];     d3gdr3[0][1][4] = d3gdr3[0][0][4]; \
d3gdr3[0][2][0] = d3gdr3[0][0][1];     d3gdr3[0][2][1] = d3gdr3[0][0][1]; \
d3gdr3[0][2][2] = d3gdr3[0][0][1];     d3gdr3[0][2][3] = d3gdr3[0][0][1]; \
d3gdr3[0][2][4] = d3gdr3[0][0][4];     d3gdr3[0][3][0] = d3gdr3[0][0][1]; \
d3gdr3[0][3][1] = d3gdr3[0][0][1];     d3gdr3[0][3][2] = d3gdr3[0][0][1]; \
d3gdr3[0][3][3] = d3gdr3[0][0][1];     d3gdr3[0][3][4] = d3gdr3[0][0][4]; \
d3gdr3[0][4][0] = d3gdr3[0][0][4];     d3gdr3[0][4][1] = d3gdr3[0][0][4]; \
d3gdr3[0][4][2] = d3gdr3[0][0][4];     d3gdr3[0][4][3] = d3gdr3[0][0][4]; \
d3gdr3[0][4][4] = D3GDR4R4R0;          d3gdr3[1][0][0] = d3gdr3[0][0][1]; \
d3gdr3[1][0][1] = d3gdr3[0][0][1];     d3gdr3[1][0][2] = d3gdr3[0][0][1]; \
d3gdr3[1][0][3] = d3gdr3[0][0][1];     d3gdr3[1][0][4] = d3gdr3[0][0][4]; \
d3gdr3[1][1][0] = d3gdr3[0][0][1];     d3gdr3[1][1][1] = D3GDR1R1R1; \
d3gdr3[1][1][2] = d3gdr3[0][0][1];     d3gdr3[1][1][3] = d3gdr3[0][0][1]; \
d3gdr3[1][1][4] = d3gdr3[0][0][4];     d3gdr3[1][2][0] = d3gdr3[0][0][1]; \
d3gdr3[1][2][1] = d3gdr3[0][0][1];     d3gdr3[1][2][2] = d3gdr3[0][0][1]; \
d3gdr3[1][2][3] = d3gdr3[0][0][1];     d3gdr3[1][2][4] = d3gdr3[0][0][4]; \
d3gdr3[1][3][0] = d3gdr3[0][0][1];     d3gdr3[1][3][1] = d3gdr3[0][0][1]; \
d3gdr3[1][3][2] = d3gdr3[0][0][1];     d3gdr3[1][3][3] = d3gdr3[0][0][1]; \
d3gdr3[1][3][4] = d3gdr3[0][0][4];     d3gdr3[1][4][0] = d3gdr3[0][0][4]; \
d3gdr3[1][4][1] = d3gdr3[0][0][4];     d3gdr3[1][4][2] = d3gdr3[0][0][4]; \
d3gdr3[1][4][3] = d3gdr3[0][0][4];     d3gdr3[1][4][4] = d3gdr3[0][4][4]; \
d3gdr3[2][0][0] = d3gdr3[0][0][1];     d3gdr3[2][0][1] = d3gdr3[0][0][1]; \
d3gdr3[2][0][2] = d3gdr3[0][0][1];     d3gdr3[2][0][3] = d3gdr3[0][0][1]; \
d3gdr3[2][0][4] = d3gdr3[0][0][4];     d3gdr3[2][1][0] = d3gdr3[0][0][1]; \
d3gdr3[2][1][1] = d3gdr3[0][0][1];     d3gdr3[2][1][2] = d3gdr3[0][0][1]; \
d3gdr3[2][1][3] = d3gdr3[0][0][1];     d3gdr3[2][1][4] = d3gdr3[0][0][4]; \
d3gdr3[2][2][0] = d3gdr3[0][0][1];     d3gdr3[2][2][1] = d3gdr3[0][0][1]; \
d3gdr3[2][2][2] = D3GDR2R2R2;          d3gdr3[2][2][3] = d3gdr3[0][0][1]; \
d3gdr3[2][2][4] = d3gdr3[0][0][4];     d3gdr3[2][3][0] = d3gdr3[0][0][1]; \
d3gdr3[2][3][1] = d3gdr3[0][0][1];     d3gdr3[2][3][2] = d3gdr3[0][0][1]; \
d3gdr3[2][3][3] = d3gdr3[0][0][1];     d3gdr3[2][3][4] = d3gdr3[0][0][4]; \
d3gdr3[2][4][0] = d3gdr3[0][0][4];     d3gdr3[2][4][1] = d3gdr3[0][0][4]; \
d3gdr3[2][4][2] = d3gdr3[0][0][4];     d3gdr3[2][4][3] = d3gdr3[0][0][4]; \
d3gdr3[2][4][4] = d3gdr3[0][4][4];     d3gdr3[3][0][0] = d3gdr3[0][0][1]; \
d3gdr3[3][0][1] = d3gdr3[0][0][1];     d3gdr3[3][0][2] = d3gdr3[0][0][1]; \
d3gdr3[3][0][3] = d3gdr3[0][0][1];     d3gdr3[3][0][4] = d3gdr3[0][0][4]; \
d3gdr3[3][1][0] = d3gdr3[0][0][1];     d3gdr3[3][1][1] = d3gdr3[0][0][1]; \
d3gdr3[3][1][2] = d3gdr3[0][0][1];     d3gdr3[3][1][3] = d3gdr3[0][0][1]; \
d3gdr3[3][1][4] = d3gdr3[0][0][4];     d3gdr3[3][2][0] = d3gdr3[0][0][1]; \
d3gdr3[3][2][1] = d3gdr3[0][0][1];     d3gdr3[3][2][2] = d3gdr3[0][0][1]; \
d3gdr3[3][2][3] = d3gdr3[0][0][1];     d3gdr3[3][2][4] = d3gdr3[0][0][4]; \
d3gdr3[3][3][0] = d3gdr3[0][0][1];     d3gdr3[3][3][1] = d3gdr3[0][0][1]; \
d3gdr3[3][3][2] = d3gdr3[0][0][1];     d3gdr3[3][3][3] = D3GDR3R3R3; \
d3gdr3[3][3][4] = d3gdr3[0][0][4];     d3gdr3[3][4][0] = d3gdr3[0][0][4]; \
d3gdr3[3][4][1] = d3gdr3[0][0][4];     d3gdr3[3][4][2] = d3gdr3[0][0][4]; \
d3gdr3[3][4][3] = d3gdr3[0][0][4];     d3gdr3[3][4][4] = d3gdr3[0][4][4]; \
d3gdr3[4][0][0] = d3gdr3[0][0][4];     d3gdr3[4][0][1] = d3gdr3[0][0][4]; \
d3gdr3[4][0][2] = d3gdr3[0][0][4];     d3gdr3[4][0][3] = d3gdr3[0][0][4]; \
d3gdr3[4][0][4] = d3gdr3[0][4][4];     d3gdr3[4][1][0] = d3gdr3[0][0][4]; \
d3gdr3[4][1][1] = d3gdr3[0][0][4];     d3gdr3[4][1][2] = d3gdr3[0][0][4]; \
d3gdr3[4][1][3] = d3gdr3[0][0][4];     d3gdr3[4][1][4] = d3gdr3[0][4][4]; \
d3gdr3[4][2][0] = d3gdr3[0][0][4];     d3gdr3[4][2][1] = d3gdr3[0][0][4]; \
d3gdr3[4][2][2] = d3gdr3[0][0][4];     d3gdr3[4][2][3] = d3gdr3[0][0][4]; \
d3gdr3[4][2][4] = d3gdr3[0][4][4];     d3gdr3[4][3][0] = d3gdr3[0][0][4]; \
d3gdr3[4][3][1] = d3gdr3[0][0][4];     d3gdr3[4][3][2] = d3gdr3[0][0][4]; \
d3gdr3[4][3][3] = d3gdr3[0][0][4];     d3gdr3[4][3][4] = d3gdr3[0][4][4]; \
d3gdr3[4][4][0] = d3gdr3[0][4][4];     d3gdr3[4][4][1] = d3gdr3[0][4][4]; \
d3gdr3[4][4][2] = d3gdr3[0][4][4];     d3gdr3[4][4][3] = d3gdr3[0][4][4]; \
d3gdr3[4][4][4] = D3GDR4R4R4;

#define fillD3GDR2DS \
d3gdr2ds[0][0][0] = D3GDR0R0S0;        d3gdr2ds[0][0][1] = D3GDR0R0S1; \
d3gdr2ds[0][0][2] = D3GDR0R0S2;        d3gdr2ds[0][0][3] = D3GDR0R0S3; \
d3gdr2ds[0][1][0] = D3GDR0R1S0;        d3gdr2ds[0][1][1] = D3GDR0R1S1; \
d3gdr2ds[0][1][2] = D3GDR0R1S2;        d3gdr2ds[0][1][3] = D3GDR0R1S3; \
d3gdr2ds[0][2][0] = D3GDR0R2S0;        d3gdr2ds[0][2][1] = D3GDR0R2S1; \
d3gdr2ds[0][2][2] = D3GDR0R2S2;        d3gdr2ds[0][2][3] = D3GDR0R2S3; \
d3gdr2ds[0][3][0] = D3GDR0R3S0;        d3gdr2ds[0][3][1] = D3GDR0R3S1; \
d3gdr2ds[0][3][2] = D3GDR0R3S2;        d3gdr2ds[0][3][3] = D3GDR0R3S3; \
d3gdr2ds[0][4][0] = D3GDR0R4S0;        d3gdr2ds[0][4][1] = D3GDR0R4S1; \
d3gdr2ds[0][4][2] = D3GDR0R4S2;        d3gdr2ds[0][4][3] = D3GDR0R4S3; \
d3gdr2ds[1][0][0] = d3gdr2ds[0][1][0]; d3gdr2ds[1][0][1] = d3gdr2ds[0][1][1]; \
d3gdr2ds[1][0][2] = d3gdr2ds[0][1][2]; d3gdr2ds[1][0][3] = d3gdr2ds[0][1][3]; \
d3gdr2ds[1][1][0] = D3GDR1R1S0;        d3gdr2ds[1][1][1] = D3GDR1R1S1; \
d3gdr2ds[1][1][2] = D3GDR1R1S2;        d3gdr2ds[1][1][3] = D3GDR1R1S3; \
d3gdr2ds[1][2][0] = D3GDR1R2S0;        d3gdr2ds[1][2][1] = D3GDR1R2S1; \
d3gdr2ds[1][2][2] = D3GDR1R2S2;        d3gdr2ds[1][2][3] = D3GDR1R2S3; \
d3gdr2ds[1][3][0] = D3GDR1R3S0;        d3gdr2ds[1][3][1] = D3GDR1R3S1; \
d3gdr2ds[1][3][2] = D3GDR1R3S2;        d3gdr2ds[1][3][3] = D3GDR1R3S3; \
d3gdr2ds[1][4][0] = D3GDR1R4S0;        d3gdr2ds[1][4][1] = D3GDR1R4S1; \
d3gdr2ds[1][4][2] = D3GDR1R4S2;        d3gdr2ds[1][4][3] = D3GDR1R4S3; \
d3gdr2ds[2][0][0] = d3gdr2ds[0][2][0]; d3gdr2ds[2][0][1] = d3gdr2ds[0][2][1]; \
d3gdr2ds[2][0][2] = d3gdr2ds[0][2][2]; d3gdr2ds[2][0][3] = d3gdr2ds[0][2][3]; \
d3gdr2ds[2][1][0] = d3gdr2ds[1][2][0]; d3gdr2ds[2][1][1] = d3gdr2ds[1][2][1]; \
d3gdr2ds[2][1][2] = d3gdr2ds[1][2][2]; d3gdr2ds[2][1][3] = d3gdr2ds[1][2][3]; \
d3gdr2ds[2][2][0] = D3GDR2R2S0;        d3gdr2ds[2][2][1] = D3GDR2R2S1; \
d3gdr2ds[2][2][2] = D3GDR2R2S2;        d3gdr2ds[2][2][3] = D3GDR2R2S3; \
d3gdr2ds[2][3][0] = D3GDR2R3S0;        d3gdr2ds[2][3][1] = D3GDR2R3S1; \
d3gdr2ds[2][3][2] = D3GDR2R3S2;        d3gdr2ds[2][3][3] = D3GDR2R3S3; \
d3gdr2ds[2][4][0] = D3GDR2R4S0;        d3gdr2ds[2][4][1] = D3GDR2R4S1; \
d3gdr2ds[2][4][2] = D3GDR2R4S2;        d3gdr2ds[2][4][3] = D3GDR2R4S3; \
d3gdr2ds[3][0][0] = d3gdr2ds[0][3][0]; d3gdr2ds[3][0][1] = d3gdr2ds[0][3][1]; \
d3gdr2ds[3][0][2] = d3gdr2ds[0][3][2]; d3gdr2ds[3][0][3] = d3gdr2ds[0][3][3]; \
d3gdr2ds[3][1][0] = d3gdr2ds[1][3][0]; d3gdr2ds[3][1][1] = d3gdr2ds[1][3][1]; \
d3gdr2ds[3][1][2] = d3gdr2ds[1][3][2]; d3gdr2ds[3][1][3] = d3gdr2ds[1][3][3]; \
d3gdr2ds[3][2][0] = d3gdr2ds[2][3][0]; d3gdr2ds[3][2][1] = d3gdr2ds[2][3][1]; \
d3gdr2ds[3][2][2] = d3gdr2ds[2][3][2]; d3gdr2ds[3][2][3] = d3gdr2ds[2][3][3]; \
d3gdr2ds[3][3][0] = D3GDR3R3S0;        d3gdr2ds[3][3][1] = D3GDR3R3S1; \
d3gdr2ds[3][3][2] = D3GDR3R3S2;        d3gdr2ds[3][3][3] = D3GDR3R3S3; \
d3gdr2ds[3][4][0] = D3GDR3R4S0;        d3gdr2ds[3][4][1] = D3GDR3R4S1; \
d3gdr2ds[3][4][2] = D3GDR3R4S2;        d3gdr2ds[3][4][3] = D3GDR3R4S3; \
d3gdr2ds[4][0][0] = d3gdr2ds[0][4][0]; d3gdr2ds[4][0][1] = d3gdr2ds[0][4][1]; \
d3gdr2ds[4][0][2] = d3gdr2ds[0][4][2]; d3gdr2ds[4][0][3] = d3gdr2ds[0][4][3]; \
d3gdr2ds[4][1][0] = d3gdr2ds[1][4][0]; d3gdr2ds[4][1][1] = d3gdr2ds[1][4][1]; \
d3gdr2ds[4][1][2] = d3gdr2ds[1][4][2]; d3gdr2ds[4][1][3] = d3gdr2ds[1][4][3]; \
d3gdr2ds[4][2][0] = d3gdr2ds[2][4][0]; d3gdr2ds[4][2][1] = d3gdr2ds[2][4][1]; \
d3gdr2ds[4][2][2] = d3gdr2ds[2][4][2]; d3gdr2ds[4][2][3] = d3gdr2ds[2][4][3]; \
d3gdr2ds[4][3][0] = d3gdr2ds[3][4][0]; d3gdr2ds[4][3][1] = d3gdr2ds[3][4][1]; \
d3gdr2ds[4][3][2] = d3gdr2ds[3][4][2]; d3gdr2ds[4][3][3] = d3gdr2ds[3][4][3]; \
d3gdr2ds[4][4][0] = D3GDR4R4S0;        d3gdr2ds[4][4][1] = D3GDR4R4S1; \
d3gdr2ds[4][4][2] = D3GDR4R4S2;        d3gdr2ds[4][4][3] = D3GDR4R4S3;

#define fillD3GDR2DT \
d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT; \
d3gdr2dt[0][2] = D3GDR0R2DT;     d3gdr2dt[0][3] = D3GDR0R3DT; \
d3gdr2dt[0][4] = D3GDR0R4DT; \
d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT; \
d3gdr2dt[1][2] = D3GDR1R2DT;     d3gdr2dt[1][3] = D3GDR1R3DT; \
d3gdr2dt[1][4] = D3GDR1R4DT; \
d3gdr2dt[2][0] = d3gdr2dt[0][2]; d3gdr2dt[2][1] = d3gdr2dt[1][2]; \
d3gdr2dt[2][2] = D3GDR2R2DT;     d3gdr2dt[2][3] = D3GDR2R3DT; \
d3gdr2dt[2][4] = D3GDR2R4DT; \
d3gdr2dt[3][0] = d3gdr2dt[0][3]; d3gdr2dt[3][1] = d3gdr2dt[1][3]; \
d3gdr2dt[3][2] = d3gdr2dt[2][3]; d3gdr2dt[3][3] = D3GDR3R3DT;     \
d3gdr2dt[3][4] = D3GDR3R4DT; \
d3gdr2dt[4][0] = d3gdr2dt[0][4]; d3gdr2dt[4][1] = d3gdr2dt[1][4]; \
d3gdr2dt[4][2] = d3gdr2dt[2][4]; d3gdr2dt[4][3] = d3gdr2dt[3][4]; \
d3gdr2dt[4][4] = D3GDR4R4DT;

#define fillD3GDR2DP \
d3gdr2dp[0][0] = (D3GDR0R0DP);     d3gdr2dp[0][1] = (D3GDR0R1DP); \
d3gdr2dp[0][2] = (D3GDR0R2DP);     d3gdr2dp[0][3] = (D3GDR0R3DP); \
d3gdr2dp[0][4] = (D3GDR0R4DP); \
d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = (D3GDR1R1DP); \
d3gdr2dp[1][2] = (D3GDR1R2DP);     d3gdr2dp[1][3] = (D3GDR1R3DP); \
d3gdr2dp[1][4] = (D3GDR1R4DP); \
d3gdr2dp[2][0] = d3gdr2dp[0][2]; d3gdr2dp[2][1] = d3gdr2dp[1][2]; \
d3gdr2dp[2][2] = (D3GDR2R2DP);     d3gdr2dp[2][3] = (D3GDR2R3DP); \
d3gdr2dp[2][4] = (D3GDR2R4DP); \
d3gdr2dp[3][0] = d3gdr2dp[0][3]; d3gdr2dp[3][1] = d3gdr2dp[1][3]; \
d3gdr2dp[3][2] = d3gdr2dp[2][3]; d3gdr2dp[3][3] = (D3GDR3R3DP); \
d3gdr2dp[3][4] = (D3GDR3R4DP); \
d3gdr2dp[4][0] = d3gdr2dp[0][4]; d3gdr2dp[4][1] = d3gdr2dp[1][4]; \
d3gdr2dp[4][2] = d3gdr2dp[2][4]; d3gdr2dp[4][3] = d3gdr2dp[3][4]; \
d3gdr2dp[4][4] = (D3GDR4R4DP);

#define fillD3GDRDS2 \
d3gdrds2[0][0][0] = D3GDR0S0S0;        d3gdrds2[0][0][1] = D3GDR0S0S1; \
d3gdrds2[0][0][2] = D3GDR0S0S2;        d3gdrds2[0][0][3] = D3GDR0S0S3; \
d3gdrds2[0][1][0] = d3gdrds2[0][0][1]; d3gdrds2[0][1][1] = D3GDR0S1S1; \
d3gdrds2[0][1][2] = D3GDR0S1S2;        d3gdrds2[0][1][3] = D3GDR0S1S3; \
d3gdrds2[0][2][0] = d3gdrds2[0][0][2]; d3gdrds2[0][2][1] = d3gdrds2[0][1][2]; \
d3gdrds2[0][2][2] = D3GDR0S2S2;        d3gdrds2[0][2][3] = D3GDR0S2S3; \
d3gdrds2[0][3][0] = d3gdrds2[0][0][3]; d3gdrds2[0][3][1] = d3gdrds2[0][1][3]; \
d3gdrds2[0][3][2] = d3gdrds2[0][2][3]; d3gdrds2[0][3][3] = D3GDR0S3S3; \
d3gdrds2[1][0][0] = D3GDR1S0S0;        d3gdrds2[1][0][1] = D3GDR1S0S1; \
d3gdrds2[1][0][2] = D3GDR1S0S2;        d3gdrds2[1][0][3] = D3GDR1S0S3; \
d3gdrds2[1][1][0] = d3gdrds2[1][0][1]; d3gdrds2[1][1][1] = D3GDR1S1S1; \
d3gdrds2[1][1][2] = D3GDR1S1S2;        d3gdrds2[1][1][3] = D3GDR1S1S3; \
d3gdrds2[1][2][0] = d3gdrds2[1][0][2]; d3gdrds2[1][2][1] = d3gdrds2[1][1][2]; \
d3gdrds2[1][2][2] = D3GDR1S2S2;        d3gdrds2[1][2][3] = D3GDR1S2S3; \
d3gdrds2[1][3][0] = d3gdrds2[1][0][3]; d3gdrds2[1][3][1] = d3gdrds2[1][1][3]; \
d3gdrds2[1][3][2] = d3gdrds2[1][2][3]; d3gdrds2[1][3][3] = D3GDR1S3S3; \
d3gdrds2[2][0][0] = D3GDR2S0S0;        d3gdrds2[2][0][1] = D3GDR2S0S1; \
d3gdrds2[2][0][2] = D3GDR2S0S2;        d3gdrds2[2][0][3] = D3GDR2S0S3; \
d3gdrds2[2][1][0] = d3gdrds2[2][0][1]; d3gdrds2[2][1][1] = D3GDR2S1S1; \
d3gdrds2[2][1][2] = D3GDR2S1S2;        d3gdrds2[2][1][3] = D3GDR2S1S3; \
d3gdrds2[2][2][0] = d3gdrds2[2][0][2]; d3gdrds2[2][2][1] = d3gdrds2[2][1][2]; \
d3gdrds2[2][2][2] = D3GDR2S2S2;        d3gdrds2[2][2][3] = D3GDR2S2S3; \
d3gdrds2[2][3][0] = d3gdrds2[2][0][3]; d3gdrds2[2][3][1] = d3gdrds2[2][1][3]; \
d3gdrds2[2][3][2] = d3gdrds2[2][2][3]; d3gdrds2[2][3][3] = D3GDR2S3S3; \
d3gdrds2[3][0][0] = D3GDR3S0S0;        d3gdrds2[3][0][1] = D3GDR3S0S1; \
d3gdrds2[3][0][2] = D3GDR3S0S2;        d3gdrds2[3][0][3] = D3GDR3S0S3; \
d3gdrds2[3][1][0] = d3gdrds2[3][0][1]; d3gdrds2[3][1][1] = D3GDR3S1S1; \
d3gdrds2[3][1][2] = D3GDR3S1S2;        d3gdrds2[3][1][3] = D3GDR3S1S3; \
d3gdrds2[3][2][0] = d3gdrds2[3][0][2]; d3gdrds2[3][2][1] = d3gdrds2[3][1][2]; \
d3gdrds2[3][2][2] = D3GDR3S2S2;        d3gdrds2[3][2][3] = D3GDR3S2S3; \
d3gdrds2[3][3][0] = d3gdrds2[3][0][3]; d3gdrds2[3][3][1] = d3gdrds2[3][1][3]; \
d3gdrds2[3][3][2] = d3gdrds2[3][2][3]; d3gdrds2[3][3][3] = D3GDR3S3S3; \
d3gdrds2[4][0][0] = D3GDR4S0S0;        d3gdrds2[4][0][1] = D3GDR4S0S1; \
d3gdrds2[4][0][2] = D3GDR4S0S2;        d3gdrds2[4][0][3] = D3GDR4S0S3; \
d3gdrds2[4][1][0] = d3gdrds2[4][0][1]; d3gdrds2[4][1][1] = D3GDR4S1S1; \
d3gdrds2[4][1][2] = D3GDR4S1S2;        d3gdrds2[4][1][3] = D3GDR4S1S3; \
d3gdrds2[4][2][0] = d3gdrds2[4][0][2]; d3gdrds2[4][2][1] = d3gdrds2[4][1][2]; \
d3gdrds2[4][2][2] = D3GDR4S2S2;        d3gdrds2[4][2][3] = D3GDR4S2S3; \
d3gdrds2[4][3][0] = d3gdrds2[4][0][3]; d3gdrds2[4][3][1] = d3gdrds2[4][1][3]; \
d3gdrds2[4][3][2] = d3gdrds2[4][2][3]; d3gdrds2[4][3][3] = D3GDR4S3S3; \

#define fillD3GDRDSDT \
d3gdrdsdt[0][0] = D3GDR0S0DT;  d3gdrdsdt[0][1] = D3GDR0S1DT; \
d3gdrdsdt[0][2] = D3GDR0S2DT;  d3gdrdsdt[0][3] = D3GDR0S3DT; \
d3gdrdsdt[1][0] = D3GDR1S0DT;  d3gdrdsdt[1][1] = D3GDR1S1DT; \
d3gdrdsdt[1][2] = D3GDR1S2DT;  d3gdrdsdt[1][3] = D3GDR1S3DT; \
d3gdrdsdt[2][0] = D3GDR2S0DT;  d3gdrdsdt[2][1] = D3GDR2S1DT; \
d3gdrdsdt[2][2] = D3GDR2S2DT;  d3gdrdsdt[2][3] = D3GDR2S3DT; \
d3gdrdsdt[3][0] = D3GDR3S0DT;  d3gdrdsdt[3][1] = D3GDR3S1DT; \
d3gdrdsdt[3][2] = D3GDR3S2DT;  d3gdrdsdt[3][3] = D3GDR3S3DT; \
d3gdrdsdt[4][0] = D3GDR4S0DT;  d3gdrdsdt[4][1] = D3GDR4S1DT; \
d3gdrdsdt[4][2] = D3GDR4S2DT;  d3gdrdsdt[4][3] = D3GDR4S3DT;

#define fillD3GDRDSDP \
d3gdrdsdp[0][0] = D3GDR0S0DP;  d3gdrdsdp[0][1] = D3GDR0S1DP; \
d3gdrdsdp[0][2] = D3GDR0S2DP;  d3gdrdsdp[0][3] = D3GDR0S3DP; \
d3gdrdsdp[1][0] = D3GDR1S0DP;  d3gdrdsdp[1][1] = D3GDR1S1DP; \
d3gdrdsdp[1][2] = D3GDR1S2DP;  d3gdrdsdp[1][3] = D3GDR1S3DP; \
d3gdrdsdp[2][0] = D3GDR2S0DP;  d3gdrdsdp[2][1] = D3GDR2S1DP; \
d3gdrdsdp[2][2] = D3GDR2S2DP;  d3gdrdsdp[2][3] = D3GDR2S3DP; \
d3gdrdsdp[3][0] = D3GDR3S0DP;  d3gdrdsdp[3][1] = D3GDR3S1DP; \
d3gdrdsdp[3][2] = D3GDR3S2DP;  d3gdrdsdp[3][3] = D3GDR3S3DP; \
d3gdrdsdp[4][0] = D3GDR4S0DP;  d3gdrdsdp[4][1] = D3GDR4S1DP; \
d3gdrdsdp[4][2] = D3GDR4S2DP;  d3gdrdsdp[4][3] = D3GDR4S3DP;

#define fillD3GDS3 \
d3gds3[0][0][0] = D3GDS0S0S0;      d3gds3[0][0][1] = D3GDS0S0S1; \
d3gds3[0][0][2] = D3GDS0S0S2;      d3gds3[0][0][3] = D3GDS0S0S3; \
d3gds3[0][1][0] = d3gds3[0][0][1]; d3gds3[0][1][1] = D3GDS0S1S1; \
d3gds3[0][1][2] = D3GDS0S1S2;      d3gds3[0][1][3] = D3GDS0S1S3; \
d3gds3[0][2][0] = d3gds3[0][0][2]; d3gds3[0][2][1] = d3gds3[0][1][2]; \
d3gds3[0][2][2] = D3GDS0S2S2;      d3gds3[0][2][3] = D3GDS0S2S3; \
d3gds3[0][3][0] = d3gds3[0][0][3]; d3gds3[0][3][1] = d3gds3[0][1][3]; \
d3gds3[0][3][2] = d3gds3[0][2][3]; d3gds3[0][3][3] = D3GDS0S3S3; \
d3gds3[1][0][0] = d3gds3[0][0][1]; d3gds3[1][0][1] = d3gds3[0][1][1]; \
d3gds3[1][0][2] = d3gds3[0][1][2]; d3gds3[1][0][3] = d3gds3[0][1][3]; \
d3gds3[1][1][0] = d3gds3[0][1][1]; d3gds3[1][1][1] = D3GDS1S1S1; \
d3gds3[1][1][2] = D3GDS1S1S2;      d3gds3[1][1][3] = D3GDS1S1S3; \
d3gds3[1][2][0] = d3gds3[0][1][2]; d3gds3[1][2][1] = d3gds3[1][1][2]; \
d3gds3[1][2][2] = D3GDS1S2S2;      d3gds3[1][2][3] = D3GDS1S2S3; \
d3gds3[1][3][0] = d3gds3[0][1][3]; d3gds3[1][3][1] = d3gds3[1][1][3]; \
d3gds3[1][3][2] = d3gds3[1][2][3]; d3gds3[1][3][3] = D3GDS1S3S3; \
d3gds3[2][0][0] = d3gds3[0][0][2]; d3gds3[2][0][1] = d3gds3[0][1][2]; \
d3gds3[2][0][2] = d3gds3[0][2][2]; d3gds3[2][0][3] = d3gds3[0][2][3]; \
d3gds3[2][1][0] = d3gds3[0][1][2]; d3gds3[2][1][1] = d3gds3[1][1][2]; \
d3gds3[2][1][2] = d3gds3[1][2][2]; d3gds3[2][1][3] = d3gds3[1][2][3]; \
d3gds3[2][2][0] = d3gds3[0][2][2]; d3gds3[2][2][1] = d3gds3[1][2][2]; \
d3gds3[2][2][2] = D3GDS2S2S2;      d3gds3[2][2][3] = D3GDS2S2S3; \
d3gds3[2][3][0] = d3gds3[0][2][3]; d3gds3[2][3][1] = d3gds3[1][2][3]; \
d3gds3[2][3][2] = d3gds3[2][2][3]; d3gds3[2][3][3] = D3GDS2S3S3; \
d3gds3[3][0][0] = d3gds3[0][0][3]; d3gds3[3][0][1] = d3gds3[0][1][3]; \
d3gds3[3][0][2] = d3gds3[0][2][3]; d3gds3[3][0][3] = d3gds3[0][3][3]; \
d3gds3[3][1][0] = d3gds3[0][1][3]; d3gds3[3][1][1] = d3gds3[1][1][3]; \
d3gds3[3][1][2] = d3gds3[1][2][3]; d3gds3[3][1][3] = d3gds3[1][3][3]; \
d3gds3[3][2][0] = d3gds3[0][2][3]; d3gds3[3][2][1] = d3gds3[1][2][3]; \
d3gds3[3][2][2] = d3gds3[2][2][3]; d3gds3[3][2][3] = d3gds3[2][3][3]; \
d3gds3[3][3][0] = d3gds3[0][3][3]; d3gds3[3][3][1] = d3gds3[1][3][3]; \
d3gds3[3][3][2] = d3gds3[2][3][3]; d3gds3[3][3][3] = D3GDS3S3S3;

#define fillD3GDS2DT \
d3gds2dt[0][0] = D3GDS0S0DT;     d3gds2dt[0][1] = D3GDS0S1DT; \
d3gds2dt[0][2] = D3GDS0S2DT;     d3gds2dt[0][3] = D3GDS0S3DT; \
d3gds2dt[1][0] = d3gds2dt[0][1]; d3gds2dt[1][1] = D3GDS1S1DT; \
d3gds2dt[1][2] = D3GDS1S2DT;     d3gds2dt[1][3] = D3GDS1S3DT; \
d3gds2dt[2][0] = d3gds2dt[0][2]; d3gds2dt[2][1] = d3gds2dt[1][2]; \
d3gds2dt[2][2] = D3GDS2S2DT;     d3gds2dt[2][3] = D3GDS2S3DT; \
d3gds2dt[3][0] = d3gds2dt[0][3]; d3gds2dt[3][1] = d3gds2dt[1][3]; \
d3gds2dt[3][2] = d3gds2dt[2][3]; d3gds2dt[3][3] = D3GDS3S3DT;

#define fillD3GDS2DP \
d3gds2dp[0][0] = D3GDS0S0DP;     d3gds2dp[0][1] = D3GDS0S1DP; \
d3gds2dp[0][2] = D3GDS0S2DP;     d3gds2dp[0][3] = D3GDS0S3DP; \
d3gds2dp[1][0] = d3gds2dp[0][1]; d3gds2dp[1][1] = D3GDS1S1DP; \
d3gds2dp[1][2] = D3GDS1S2DP;     d3gds2dp[1][3] = D3GDS1S3DP; \
d3gds2dp[2][0] = d3gds2dp[0][2]; d3gds2dp[2][1] = d3gds2dp[1][2]; \
d3gds2dp[2][2] = D3GDS2S2DP;     d3gds2dp[2][3] = D3GDS2S3DP; \
d3gds2dp[3][0] = d3gds2dp[0][3]; d3gds2dp[3][1] = d3gds2dp[1][3]; \
d3gds2dp[3][2] = d3gds2dp[2][3]; d3gds2dp[3][3] = D3GDS3S3DP;

#define fillD3GDSDT2 \
d3gdsdt2[0] = D3GDS0DT2; d3gdsdt2[1] = D3GDS1DT2; \
d3gdsdt2[2] = D3GDS2DT2; d3gdsdt2[3] = D3GDS3DT2;

#define fillD3GDSDTDP \
d3gdsdtdp[0] = D3GDS0DTDP; d3gdsdtdp[1] = D3GDS1DTDP; \
d3gdsdtdp[2] = D3GDS2DTDP; d3gdsdtdp[3] = D3GDS3DTDP;

#define fillD3GDSDP2 \
d3gdsdp2[0] = D3GDS0DP2; d3gdsdp2[1] = D3GDS1DP2; \
d3gdsdp2[2] = D3GDS2DP2; d3gdsdp2[3] = D3GDS3DP2;

#define fillD3GDRDT2 \
d3gdrdt2[0] = D3GDR0DT2; d3gdrdt2[1] = D3GDR1DT2; \
d3gdrdt2[2] = D3GDR2DT2; d3gdrdt2[3] = D3GDR3DT2; \
d3gdrdt2[4] = D3GDR4DT2;

#define fillD3GDRDTDP \
d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP; \
d3gdrdtdp[2] = D3GDR2DTDP; d3gdrdtdp[3] = D3GDR3DTDP; \
d3gdrdtdp[4] = D3GDR4DTDP;

#define fillD3GDRDP2 \
d3gdrdp2[0] = D3GDR0DP2; d3gdrdp2[1] = D3GDR1DP2; \
d3gdrdp2[2] = D3GDR2DP2; d3gdrdp2[3] = D3GDR3DP2; \
d3gdrdp2[4] = D3GDR4DP2;

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
	int i, j;

	/* look-up or compute the current ordering state */
	if ( (t != tOld)       || (p != pOld) ||
		(r[0] != rOld[0]) || (r[1] != rOld[1]) || (r[2] != rOld[2]) ||
		(r[3] != rOld[3]) || (r[4] != rOld[4]) ) {
		double dgds[NS], sNew[NS], sMax[NS], sMin[NS];

		for (i=0; i<NS; i++) sOld[i] = 2.0;

		/* calculate initial guesses for ordering variables  */

		sMax[0]  = MIN(1.0-r[0],1.0+r[0]);
		sMax[0]  = MIN(sMax[0],1.0-r[0]-2.0*r[4]);
		sMin[0]  = MAX(r[0]-1.0,-r[0]-1.0);

		sMax[1]  = MIN(1.0-r[1],1.0+r[1]);
		sMax[1]  = MIN(sMax[1],1.0-r[1]-2.0*r[4]);
		sMin[1]  = MAX(r[1]-1.0,-r[1]-1.0);

		sMax[2]  = MIN(1.0-r[2],1.0+r[2]);
		sMin[2]  = MAX(r[2]-1.0,-r[2]-1.0);
		sMin[2]  = MAX(sMin[2],r[2]-1.0-2.0*r[4]);

		sMax[3]  = MIN(1.0-r[3],1.0+r[3]);
		sMin[3]  = MAX(r[3]-1.0,-r[3]-1.0);
		sMin[3]  = MAX(sMin[3],r[3]-1.0-2.0*r[4]);


		sNew[0] = (r[0]>0.0) ? sMax[0]/3.0  :0.35*sMin[0]+0.65*sMax[0];
		sNew[1] = (r[1]>0.0) ? sMax[1]/15.0 :0.45*sMin[1]+0.55*sMax[1];
		sNew[2] = (r[2]>0.0) ? sMax[2]/3.0  :0.25*sMin[2]+0.75*sMax[2];
		sNew[3] = (r[3]>0.0) ? sMax[3]/2.0  :0.25*sMin[3]+0.75*sMax[3];

		/* END OF INITIAL GUESSES    */

		while ((fabs(sNew[0]-sOld[0]) > 100.0*DBL_EPSILON) ||
			   (fabs(sNew[1]-sOld[1]) > 100.0*DBL_EPSILON) ||
			   (fabs(sNew[2]-sOld[2]) > 100.0*DBL_EPSILON) ||
			   (fabs(sNew[3]-sOld[3]) > 100.0*DBL_EPSILON) ) {
			double s[NS];

			for (i=0; i<NS; i++) s[i] = sNew[i];

			if (r[0] == -1.0) s[0]=0.0;
			if (r[1] == -1.0) s[1]=0.0;
			if (r[2] == -1.0) s[2]=0.0;
			if (r[3] == -1.0) s[3]=0.0;

			xm1mn= (r[0]-s[0]+1.0)/2.0;
			xm2mn= (r[0]+s[0]+1.0)/2.0;
			xm1fe= (r[1]-s[1]+1.0)/2.0;
			xm2fe= (r[1]+s[1]+1.0)/2.0;
			xm1co= (r[2]+s[2]+1.0)/2.0;
			xm2co= (r[2]-s[2]+1.0)/2.0;
			xm1ni= (r[3]+s[3]+1.0)/2.0;
			xm2ni= (r[3]-s[3]+1.0)/2.0;
			xm2ca=  r[4];
			xm1mg = (1.0 - xm1mn - xm1fe - xm1co - xm1ni);
			xm2mg = (1.0 - xm2mn - xm2fe - xm2co - xm2ni - xm2ca);

			if (xm1mn   <= DBL_EPSILON) xm1mn   = DBL_EPSILON;
			if (xm1fe   <= DBL_EPSILON) xm1fe   = DBL_EPSILON;
			if (xm1co   <= DBL_EPSILON) xm1co   = DBL_EPSILON;
			if (xm1ni   <= DBL_EPSILON) xm1ni   = DBL_EPSILON;
			if (xm1mg   <= DBL_EPSILON) xm1mg   = DBL_EPSILON;
			if (xm2mn   <= DBL_EPSILON) xm2mn   = DBL_EPSILON;
			if (xm2fe   <= DBL_EPSILON) xm2fe   = DBL_EPSILON;
			if (xm2co   <= DBL_EPSILON) xm2co   = DBL_EPSILON;
			if (xm2ni   <= DBL_EPSILON) xm2ni   = DBL_EPSILON;
			if (xm2mg   <= DBL_EPSILON) xm2mg   = DBL_EPSILON;
			if (xm2ca   <= DBL_EPSILON) xm2ca   = DBL_EPSILON;

			if (fabs(r[0]+1.0)<10.0*DBL_EPSILON) {
				dgds[0] = 0.0;
				invd2gds2[0][0] = 1.0;
				xm1mn = DBL_EPSILON;
				xm2mn = xm1mn;
			}
			else {
				dgds[0] = DGDS0;
				invd2gds2[0][0] = D2GDS0S0;
			}
			if  (fabs(r[1]+1.0)<10.0*DBL_EPSILON) {
				dgds[1] = 0.0;
				invd2gds2[1][1]=1.0;
				xm1fe = DBL_EPSILON;
				xm2fe = xm1fe;
			}
			else {
				dgds[1] = DGDS1;
				invd2gds2[1][1] = D2GDS1S1;
			}
			if (fabs(r[2]+1.0)<10.0*DBL_EPSILON) {
				dgds[2] = 0.0;
				invd2gds2[2][2]=1.0;
				xm2co = DBL_EPSILON;
				xm1co = xm2co;
			}
			else {
				dgds[2]      = DGDS2;
				invd2gds2[2][2]  = D2GDS2S2;
			}
			if (fabs(r[3]+1.0)<10.0*DBL_EPSILON) {
				dgds[3]    =0.0;
				invd2gds2[3][3]=1.0;
				xm2ni = DBL_EPSILON;
				xm1ni = xm2ni;
			}
			else {
				dgds[3]      = DGDS3;
				invd2gds2[3][3]  = D2GDS3S3;
			}

			invd2gds2[0][1] =
			((fabs(r[0]+1.0)>10.0*DBL_EPSILON)&&(fabs(r[1]+1.0)>10.0*DBL_EPSILON))
            ? D2GDS0S1 : 0.0;
			invd2gds2[0][2] =
			((fabs(r[0]+1.0)>10.0*DBL_EPSILON)&&(fabs(r[2]+1.0)>10.0*DBL_EPSILON))
            ? D2GDS0S2 : 0.0;
			invd2gds2[0][3] =
			((fabs(r[0]+1.0)>10.0*DBL_EPSILON)&&(fabs(r[3]+1.0)>10.0*DBL_EPSILON))
            ? D2GDS0S3 : 0.0;
			invd2gds2[1][2] =
			((fabs(r[1]+1.0)>10.0*DBL_EPSILON)&&(fabs(r[2]+1.0)>10.0*DBL_EPSILON))
            ? D2GDS1S2 : 0.0;
			invd2gds2[1][3] =
			((fabs(r[1]+1.0)>10.0*DBL_EPSILON)&&(fabs(r[3]+1.0)>10.0*DBL_EPSILON))
            ? D2GDS1S3 : 0.0;
			invd2gds2[2][3] =
			((fabs(r[2]+1.0)>10.0*DBL_EPSILON)&&(fabs(r[3]+1.0)>10.0*DBL_EPSILON))
            ? D2GDS2S3 : 0.0;
			invd2gds2[1][0] = invd2gds2[0][1];
			invd2gds2[2][0] = invd2gds2[0][2];
			invd2gds2[3][0] = invd2gds2[0][3];
			invd2gds2[2][1] = invd2gds2[1][2];
			invd2gds2[3][1] = invd2gds2[1][3];
			invd2gds2[3][2] = invd2gds2[2][3];

			for (i=0; i<NS; i++) sOld[i] = s[i];

			[self gaussj:invd2gds2];

			if (fabs(r[0]+1.0)<10.0*DBL_EPSILON)invd2gds2[0][0]=0.0;
			if (fabs(r[1]+1.0)<10.0*DBL_EPSILON)invd2gds2[1][1]=0.0;
			if (fabs(r[2]+1.0)<10.0*DBL_EPSILON)invd2gds2[2][2]=0.0;
			if (fabs(r[3]+1.0)<10.0*DBL_EPSILON)invd2gds2[3][3]=0.0;


			for (i=0; i<NS; i++) {
				for(j=0; j<NS; j++) s[i] += - invd2gds2[i][j]*dgds[j];
			}

			if (fabs(r[0]+1.0)>10.0*DBL_EPSILON) {
				s[0] = MAX(s[0], sMin[0]+3.0*DBL_EPSILON);
				s[0] = MIN(s[0], sMax[0]-3.0*DBL_EPSILON);
			}
			if (fabs(r[1]+1.0)>10.0*DBL_EPSILON) {
				s[1] = MAX(s[1], sMin[1]+3.0*DBL_EPSILON);
				s[1] = MIN(s[1], sMax[1]-3.0*DBL_EPSILON);
			}
			if (fabs(r[2]+1.0)>10.0*DBL_EPSILON) {
				s[2] = MAX(s[2], sMin[2]+3.0*DBL_EPSILON);
				s[2] = MIN(s[2], sMax[2]-3.0*DBL_EPSILON);
			}
			if (fabs(r[3]+1.0)>10.0*DBL_EPSILON) {
				s[3] = MAX(s[3], sMin[3]+3.0*DBL_EPSILON);
				s[3] = MIN(s[3], sMax[3]-3.0*DBL_EPSILON);
			}

			for (i=0; i<NS; i++) sNew[i] = s[i];
		}
		tOld = t;
		pOld = p;
		for (i=0; i<NR; i++) rOld[i] = r[i];

		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
        if (debug) {
			for (i=0; i<NS; i++) {
				if (dgds[i] > sqrt(DBL_EPSILON) && ABS(sNew[i]) > sqrt(DBL_EPSILON)) {
					NSLog(@"ERROR in OLIVINE.C (function ORDER). Failed to converge!");
					NSLog(@"  r1    = %13.6g, r2    = %13.6g, r3    = %13.6g", r[0], r[1], r[2]);
					NSLog(@"  r4    = %13.6g, r5    = %13.6g", r[3], r[4]);
					NSLog(@"  s1    = %13.6g, s2    = %13.6g", sOld[0], sOld[1]);
					NSLog(@"  s3    = %13.6g, s4    = %13.6g", sOld[2], sOld[3]);
					NSLog(@"  dgds1 = %13.6g, dgds2 = %13.6g", dgds[0], dgds[1]);
					NSLog(@"  dgds3 = %13.6g, dgds4 = %13.6g", dgds[2], dgds[3]);

					NSLog(@"xm1 mn = %13.6g, xm2 mn = %13.6g", xm1mn, xm2mn);
					NSLog(@"xm1 fe = %13.6g, xm2 fe = %13.6g", xm1fe, xm2fe);
					NSLog(@"xm1 co = %13.6g, xm2 co = %13.6g", xm1co, xm2co);
					NSLog(@"xm1 ni = %13.6g, xm2 ni = %13.6g", xm1ni, xm2ni);
					NSLog(@"xm1 mg = %13.6g, xm2 mg = %13.6g", xm1mg, xm2mg);
					NSLog(@"xm2 ca = %13.6g",                  xm2ca);
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
		double *s = sOld;
		double d2gdsdp[NS];

		fillD2GDSDP

		for (i=0; i<NS; i++) {
			dp[i] = 0.0;
			for (j=0; j<NS; j++) dp[i] += - invd2gds2[i][j]*d2gdsdp[j];
		}
	}
	if (mask & FIFTH  ) {   /* compute d2s/dr2 */
		double d2gdrds[NR][NS];
		double d3gdr2ds[NR][NR][NS];
		double d3gdrds2[NR][NS][NS];
		double d3gds3[NS][NS][NS];
		double dsdr[NS][NR], temp[NS];
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
		double d2gdrds[NR][NS];
		double d2gdsdt[NS];
		double d3gdrds2[NR][NS][NS];
		double d3gdrdsdt[NR][NS];
		double d3gds2dt[NS][NS];
		double d3gds3[NS][NS][NS];
		double dsdr[NS][NR], dsdt[NS], temp[NS];
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
		double d2gdsdt[NS];
		double d3gds3[NS][NS][NS];
		double d3gds2dt[NS][NS];
		double d3gdsdt2[NS];
		double dsdt[NS], temp[NS];
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
		double d2gdsdt[NS];
		double d2gdsdp[NS];
		double d3gds3[NS][NS][NS];
		double d3gds2dt[NS][NS];
		double d3gds2dp[NS][NS];
		double d3gdsdtdp[NS];
		double dsdt[NS], dsdp[NS], temp[NS];
		int k, l;

		fillD2GDSDT
		fillD2GDSDP
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDSDTDP
		fillD3GDS2DP

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
	const char *NAMES[NA]    = { "tephroite", "fayalite", "co-olivine", "ni-olivine", "monticellite","forsterite" };
	const char *FORMULAS[NA] = { "Mn2SiO4", "Fe2SiO4", "Co2SiO4", "Ni2SiO4", "CaMgSiO4", "Mg2SiO4" };
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
		result = result && (r[0] >= -1.0) && (r[0] <= 1.0-r[4]);
		result = result && (r[1] >= -1.0) && (r[1] <= 1.0-r[4]);
		result = result && (r[2] >= -1.0) && (r[2] <= 1.0-r[4]);
		result = result && (r[3] >= -1.0) && (r[3] <= 1.0-r[4]);
		result = result && (r[4] >= 0.0) && (r[4] <= 1.0);
		result = result && (r[0]+r[1]+r[2]+r[3]+2.0*r[4] <= -2.0-r[4]);
		result = result && (r[0]+r[1]+r[2]+r[3]+2.0*r[4] >= -4.0);
	}
	/* Check bounds on moles of endmember components */
	if (mask & SIXTH) {
		for (i=0, sum=0.0; i<NA; i++) sum += m[i];
		result = result && (sum >= 0.0);

		if (sum > 0.0) {
			result = result && (m[0]/sum >= 0.0) && (m[0]/sum <=  1.0);
			result = result && (m[1]/sum >= 0.0) && (m[1]/sum <=  1.0);
			result = result && (m[2]/sum >= 0.0) && (m[2]/sum <=  1.0);
			result = result && (m[3]/sum >= 0.0) && (m[3]/sum <=  1.0);
			result = result && (m[4]/sum >= 0.0) && (m[4]/sum <=  1.0);
			result = result && (m[5]/sum >= -0.5*m[4]/sum) && (m[5]/sum <= 1.0);
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
	 (2)  SECOND           THIRD | FOURTH | FIFTH | SIXTH | EIGHTH
	 (3)  THIRD            FOURTH | SEVENTH

	 (1) converts a vector of moles of elements into a vector of moles of
	 endmember olivine components.
	 (2) calculates from a vector of moles of endmember components, one or
	 all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
	 (3) calculates from a vector of independent compositional variables
	 mole fractions of endmember components and/or the Jacobian matrix
	 dx[]/dr[]

	 In this routine it is assumed that the elements are in the order of atomic
	 numbers and that the order of olivine components has been verified as:
	 m[0] = tephroite    (Mn2SiO4) ,
	 m[1] = fayalite     (Fe2SiO4) ,
	 m[2] = Co-olivine   (Co2SiO4) ,
	 m[3] = Ni-olivine   (Ni2SiO4) ,
	 m[4] = monticellite (CaMgSiO4) and
	 m[5] = forsterite   (Mg2SiO4)

	 ----------------------------------------------------------------------------*/

	int i, j, k;

	if (inpMask == FIRST && outMask == SECOND) {
		/* Converts a vector of moles of elements into a vector of moles of
		 end-member components.                                                 */
		static const int Mg = 12;
		static const int Ca = 20;
		static const int Mn = 25;
		static const int Fe = 26;
		static const int Co = 27;
		static const int Ni = 28;

		m[0] =  e[Mn]/2.0;        /* moles of Mn2SiO4                        */
		m[1] =  e[Fe]/2.0;        /* moles of Fe2SiO4                        */
		m[2] =  e[Co]/2.0;        /* moles of Co2SiO4                        */
		m[3] =  e[Ni]/2.0;        /* moles of Ni2SiO4                        */
		m[4] =  e[Ca];            /* Moles of CaMgSiO4                       */
		m[5] = (e[Mg]-e[Ca])/2.0; /* Moles of Mg2SiO4                        */

	} else if (inpMask == SECOND) {
		double sum;

		if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
			NSLog(@"Illegal call to conOlv with inpMask = %o and outMask = %o", inpMask, outMask);

		for (i=0, sum=0.0; i<NA; i++) sum += m[i];

		if (outMask & THIRD) {
			/* Converts a vector of moles of end-member components (m) into a vector
			 of independent compositional variables                               */
			r[0] = (sum != 0.0) ? 2.0*m[0]/sum - 1.0 : 0.0;
			r[1] = (sum != 0.0) ? 2.0*m[1]/sum - 1.0 : 0.0;
			r[2] = (sum != 0.0) ? 2.0*m[2]/sum - 1.0 : 0.0;
			r[3] = (sum != 0.0) ? 2.0*m[3]/sum - 1.0 : 0.0;
			r[4] = (sum != 0.0) ? m[4]/sum           : 0.0;
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
					dm[0][j] = (j == 0) ? 2.0*(1.0-m[0]/sum)/sum : -2.0*m[0]/SQUARE(sum);
					dm[1][j] = (j == 1) ? 2.0*(1.0-m[1]/sum)/sum : -2.0*m[1]/SQUARE(sum);
					dm[2][j] = (j == 2) ? 2.0*(1.0-m[2]/sum)/sum : -2.0*m[2]/SQUARE(sum);
					dm[3][j] = (j == 3) ? 2.0*(1.0-m[3]/sum)/sum : -2.0*m[3]/SQUARE(sum);
					dm[4][j] = (j == 4) ?     (1.0-m[4]/sum)/sum : -1.0*m[4]/SQUARE(sum);
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
						d2m[0][j][k]  = 4.0*m[0]/CUBE(sum);
						d2m[0][j][k] -= (j == 0) ? 2.0/SQUARE(sum) : 0.0;
						d2m[0][j][k] -= (k == 0) ? 2.0/SQUARE(sum) : 0.0;
						d2m[1][j][k]  = 4.0*m[1]/CUBE(sum);
						d2m[1][j][k] -= (j == 1) ? 2.0/SQUARE(sum) : 0.0;
						d2m[1][j][k] -= (k == 1) ? 2.0/SQUARE(sum) : 0.0;
						d2m[2][j][k]  = 4.0*m[2]/CUBE(sum);
						d2m[2][j][k] -= (j == 2) ? 2.0/SQUARE(sum) : 0.0;
						d2m[2][j][k] -= (k == 2) ? 2.0/SQUARE(sum) : 0.0;
						d2m[3][j][k]  = 4.0*m[3]/CUBE(sum);
						d2m[3][j][k] -= (j == 3) ? 2.0/SQUARE(sum) : 0.0;
						d2m[3][j][k] -= (k == 3) ? 2.0/SQUARE(sum) : 0.0;
						d2m[4][j][k]  = 2.0*m[4]/CUBE(sum);
						d2m[4][j][k] -= (j == 4) ? 1.0/SQUARE(sum) : 0.0;
						d2m[4][j][k] -= (k == 4) ? 1.0/SQUARE(sum) : 0.0;
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
								d3m[i][j][k][l] = -12.0*m[i]/QUARTIC(sum);
								d3m[i][j][k][l] += (i == j) ? 4.0/CUBE(sum) : 0.0;
								d3m[i][j][k][l] += (i == k) ? 4.0/CUBE(sum) : 0.0;
								d3m[i][j][k][l] += (i == l) ? 4.0/CUBE(sum) : 0.0;
								if (i == 4) d3m[i][j][k][l] /= 2.0;
							}
						}
					}
				}
			}
		}

	} else if (inpMask == THIRD) {

		if (outMask & ~(FOURTH | SEVENTH))
			NSLog(@"Illegal call to conOlv with inpMask = %o and outMask = %o", inpMask, outMask);

		if (outMask & FOURTH) {
			/* Converts a vector of independent compositional variables (r) into
			 a vector of mole fractions of endmember components (x).              */

			x[0] = (1.0+r[0])/2.0;
			x[1] = (1.0+r[1])/2.0;
			x[2] = (1.0+r[2])/2.0;
			x[3] = (1.0+r[3])/2.0;
			x[4] = r[4];
			x[5] = -1.0 - (r[0]+r[1]+r[2]+r[3])/2.0 - r[4];

			for (i=0; i<NA; i++) if (fabs(x[i]) < DBL_EPSILON) x[i] = 0.0;
		}

		if (outMask & SEVENTH) {
			/* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j]                 */
			for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
			dr[0][0] =  0.5;
			dr[1][1] =  0.5;
			dr[2][2] =  0.5;
			dr[3][3] =  0.5;
			dr[4][4] =  1.0;

			dr[5][0] = -0.5;
			dr[5][1] = -0.5;
			dr[5][2] = -0.5;
			dr[5][3] = -0.5;
			dr[5][4] = -1.0;
		}

	} else  {
		NSLog(@"Illegal call to conOlv with inpMask = %o and outMask = %o", inpMask, outMask);
	}

}
-(NSString *)displayFormula:(double)t
						  p:(double)p
						  r:(double [NA])r
{
	double totCa, totFe2, totMg,totMn,totCo,totNi;

	totCa  = 0.5*r[4];
	totFe2 = (1.0+r[1])/2.0;
	totMn  = (1.0+r[0])/2.0;
	totCo  = (1.0+r[2])/2.0;
	totNi  = (1.0+r[3])/2.0;

	totMg  = (1.0-totFe2-totMn-totCo-totNi);

	return [NSString stringWithFormat:@"(Ca%4.2fMg%4.2fFe''%4.2fMn%4.2fCo%4.2fNi%4.2f)2SiO4", totCa, totMg, totFe2, totMn, totCo, totNi];
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
		fr[i][0] = FR0(i);
		fr[i][1] = FR1(i);
		fr[i][2] = FR2(i);
		fr[i][3] = FR3(i);
		fr[i][4] = FR4(i);

	}

	[self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

	g       = G;
	dgdr[0] = DGDR0;
	dgdr[1] = DGDR1;
	dgdr[2] = DGDR2;
	dgdr[3] = DGDR3;
	dgdr[4] = DGDR4;

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
		double d2gdr2[NR][NR];
		double d2gdrds[NR][NS];
		double d2gds2[NS][NS];
		double dsdr[NS][NR], dfrdr[NA][NR], gs[NA][NS], dgsds[NA][NS], sum;
		int k, l;

		fillD2GDR2
		fillD2GDRDS
		fillD2GDS2

		for(i=0; i<NA; i++) {
			gs[i][0] = GGS0(i); /* s1 */
			gs[i][1] = GGS1(i); /* s2 */
			gs[i][2] = GGS2(i); /* s3 */
			gs[i][3] = GGS3(i); /* s4 */
			dfrdr[i][0] = DFR0DR0(i); /* r1 */
			dfrdr[i][1] = DFR1DR1(i); /* r2 */
			dfrdr[i][2] = DFR2DR2(i); /* r3 */
			dfrdr[i][3] = DFR3DR3(i); /* r4 */
			dfrdr[i][4] = DFR4DR4(i); /* r5 */

			dgsds[i][0] = DGS0DS0(i); /* s1 */
			dgsds[i][1] = DGS1DS1(i); /* s2 */
			dgsds[i][2] = DGS2DS2(i); /* s3 */
			dgsds[i][3] = DGS3DS3(i); /* s4 */
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
		dx[3] = DGDR3;
		dx[4] = DGDR4;
	}

	if(mask & THIRD) {
		double d2gdr2[NR][NR];
		double d2gdrds[NR][NS];
		double d2gds2[NS][NS];
		double dsdr[NS][NR];
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
		*smix = (S);
	}

	if(mask & SECOND) {
		double d2gdrds[NR][NS];
		double d2gdrdt[NR];
		double d2gds2[NS][NS];
		double d2gdsdt[NS];
		double dsdr[NS][NR], dsdt[NS];
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

		double d2gdrds[NR][NS];
		double d2gds2[NS][NS];
		double d2gdsdt[NS];
		double d3gdr2ds[NR][NR][NS];
		double d3gdr2dt[NR][NR];
		double d3gdrds2[NR][NS][NS];
		double d3gdrdsdt[NR][NS];
		double d3gds3[NS][NS][NS];
		double d3gds2dt[NS][NS];

		double dsdr[NS][NR], dsdt[NS], d2sdr2[NS][NR][NR], d2sdrdt[NS][NR];
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
		double d3gds3[NS][NS][NS];
		double d3gds2dt[NS][NS];
		double d3gdsdt2[NS];
		double d3gdt3 = D3GDT3;
		double d2sdt2[NS], temp;
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
		*vmix = (DGDP);

	}

	if(mask & SECOND) {
		double d2gdrds[NR][NS];
		double d2gdrdp[NR];
		double d2gds2[NS][NS];
		double d2gdsdp[NS];
		double dsdr[NS][NR], dsdp[NS];
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
		double d2gdrds[NR][NS];
		double d2gds2[NS][NS];
		double d2gdsdp[NS];
		double d3gdr2ds[NR][NR][NS];
		double d3gdr2dp[NR][NR];
		double d3gdrds2[NR][NS][NS];
		double d3gdrdsdp[NR][NS];
		double d3gds3[NS][NS][NS];
		double d3gds2dp[NS][NS];
		double dsdr[NS][NR],dsdp[NS],d2sdr2[NS][NR][NR],d2sdrdp[NS][NR];
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
		double d2gds2[NS][NS];
		double d2gdsdt[NS];
		double d2gdsdp[NS];
		double d2gdtdp = D2GDTDP;
		double dsdt[NS], dsdp[NS];
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
		double d2gds2[NS][NS];
		double d2gdsdp[NS];
		double d2gdp2 = D2GDP2;
		double dsdp[NS];
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
		double d2gds2[NS][NS];
		double d2gdsdt[NS];
		double d2gdsdp[NS];
		double d3gds3[NS][NS][NS];
		double d3gds2dt[NS][NS];
		double d3gdsdt2[NS];
		double d3gds2dp[NS][NS];
		double d3gdsdtdp[NS];
		double d3gdt2dp        = D3GDT2DP;
		double dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS];
		int i, j, k;

		fillD2GDS2
		fillD2GDSDT
		fillD2GDSDP
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDSDT2
		fillD3GDS2DP
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
		double d2gds2[NS][NS];
		double d2gdsdt[NS];
		double d2gdsdp[NS];
		double d3gds3[NS][NS][NS];
		double d3gds2dt[NS][NS];
		double d3gds2dp[NS][NS];
		double d3gdsdtdp[NS];
		double d3gdsdp2[NS];
		double d3gdtdp2             = D3GDTDP2;
		double dsdt[NS], dsdp[NS], d2sdtdp[NS], d2sdp2[NS];
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
		double d2gds2[NS][NS];
		double d2gdsdp[NS];
		double d3gds3[NS][NS][NS];
		double d3gds2dp[NS][NS];
		double d3gdsdp2[NS];
		double d3gdp3          = D3GDP3;
		double dsdp[NS], d2sdp2[NS];
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

-(NSUInteger)numberOfSolutionSpecies {
	return NA;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
	switch (index) {
		case 0:
			return @"tephroite";
			break;
		case 1:
			return @"fayalite";
			break;
		case 2:
			return @"co-olivine";
			break;
		case 3:
			return @"ni-olivine";
			break;
		case 4:
			return @"monticellite";
			break;
		case 5:
			return @"forsterite";
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
			for (i=0; i<nz; i++) if (xNz[i] < DBL_EPSILON) xNz[i] = DBL_EPSILON;

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
