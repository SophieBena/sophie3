//
//  LiquidpMelts.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/16/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import "LiquidpMelts.h"
#import "BermanProperties.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"
#import "LiquidpMeltsGenericEM.h"
#import "LiquidpMeltsH2O.h"


@implementation LiquidpMelts

static NSArray *endmembers;

static double modelParameters[] = {
       15094.7,           // W(TiO2      ,Si4O8     )
     -296975.2,           // W(Al4O6     ,Si4O8     )
     -164027.4,           // W(Fe2O3     ,Si4O8     )
       37459.2,           // W(MgCr2O4   ,Si4O8     )
      -18841.4,           // W(Fe2SiO4   ,Si4O8     )
           0.0,           // W(MnSi0.5O2 ,Si4O8     )
      -33833.5,           // W(Mg2SiO4   ,Si4O8     )
           0.0,           // W(NiSi0.5O2 ,Si4O8     )
           0.0,           // W(CoSi0.5O2 ,Si4O8     )
      -34232.9,           // W(Ca2Si2O6  ,Si4O8     )
      -59822.7,           // W(NaSi0.5O1.,Si4O8     )
     -102706.5,           // W(KAlSiO4   ,Si4O8     )
       37519.9,           // W(Ca3(PO4)2 ,Si4O8     )
      -45181.6,           // W(H2O       ,Si4O8     )

     -144804.9,           // W(Al4O6     ,TiO2      )
     -212292.3,           // W(Fe2O3     ,TiO2      )
      -22455.8,           // W(MgCr2O4   ,TiO2      )
        9324.2,           // W(Fe2SiO4   ,TiO2      )
           0.0,           // W(MnSi0.5O2 ,TiO2      )
       16335.6,           // W(Mg2SiO4   ,TiO2      )
           0.0,           // W(NiSi0.5O2 ,TiO2      )
           0.0,           // W(CoSi0.5O2 ,TiO2      )
       -9471.5,           // W(Ca2Si2O6  ,TiO2      )
       22194.2,           // W(NaSi0.5O1.,TiO2      )
       -3744.0,           // W(KAlSiO4   ,TiO2      )
       65544.0,           // W(Ca3(PO4)2 ,TiO2      )
       70663.0,           // W(H2O       ,TiO2      )

     -393566.0,           // W(Fe2O3     ,Al4O6     )
     -269339.7,           // W(MgCr2O4   ,Al4O6     )
     -200788.1,           // W(Fe2SiO4   ,Al4O6     )
           0.0,           // W(MnSi0.5O2 ,Al4O6     )
     -192709.0,           // W(Mg2SiO4   ,Al4O6     )
           0.0,           // W(NiSi0.5O2 ,Al4O6     )
           0.0,           // W(CoSi0.5O2 ,Al4O6     )
     -270700.8,           // W(Ca2Si2O6  ,Al4O6     )
     -205068.6,           // W(NaSi0.5O1.,Al4O6     )
     -114506.5,           // W(KAlSiO4   ,Al4O6     )
     -176584.1,           // W(Ca3(PO4)2 ,Al4O6     )
     -161944.4,           // W(H2O       ,Al4O6     )

      201536.3,           // W(MgCr2O4   ,Fe2O3     )
     -211493.4,           // W(Fe2SiO4   ,Fe2O3     )
           0.0,           // W(MnSi0.5O2 ,Fe2O3     )
     -196914.9,           // W(Mg2SiO4   ,Fe2O3     )
           0.0,           // W(NiSi0.5O2 ,Fe2O3     )
           0.0,           // W(CoSi0.5O2 ,Fe2O3     )
     -146008.1,           // W(Ca2Si2O6  ,Fe2O3     )
     -123728.7,           // W(NaSi0.5O1.,Fe2O3     )
     -130847.5,           // W(KAlSiO4   ,Fe2O3     )
     -126339.8,           // W(Ca3(PO4)2 ,Fe2O3     )
     -114508.6,           // W(H2O       ,Fe2O3     )

      -74759.0,           // W(Fe2SiO4   ,MgCr2O4   )
           0.0,           // W(MnSi0.5O2 ,MgCr2O4   )
       -3638.5,           // W(Mg2SiO4   ,MgCr2O4   )
           0.0,           // W(NiSi0.5O2 ,MgCr2O4   )
           0.0,           // W(CoSi0.5O2 ,MgCr2O4   )
       48337.5,           // W(Ca2Si2O6  ,MgCr2O4   )
      -43302.5,           // W(NaSi0.5O1.,MgCr2O4   )
      124517.4,           // W(KAlSiO4   ,MgCr2O4   )
       13004.3,           // W(Ca3(PO4)2 ,MgCr2O4   )
         -18.9,           // W(H2O       ,MgCr2O4   )

           0.0,           // W(MnSi0.5O2 ,Fe2SiO4   )
      -28736.4,           // W(Mg2SiO4   ,Fe2SiO4   )
           0.0,           // W(NiSi0.5O2 ,Fe2SiO4   )
           0.0,           // W(CoSi0.5O2 ,Fe2SiO4   )
      -28573.8,           // W(Ca2Si2O6  ,Fe2SiO4   )
       -4723.9,           // W(NaSi0.5O1.,Fe2SiO4   )
       22245.0,           // W(KAlSiO4   ,Fe2SiO4   )
        4909.8,           // W(Ca3(PO4)2 ,Fe2SiO4   )
        9769.4,           // W(H2O       ,Fe2SiO4   )

           0.0,           // W(Mg2SiO4   ,MnSi0.5O2 )
           0.0,           // W(NiSi0.5O2 ,MnSi0.5O2 )
           0.0,           // W(CoSi0.5O2 ,MnSi0.5O2 )
           0.0,           // W(Ca2Si2O6  ,MnSi0.5O2 )
           0.0,           // W(NaSi0.5O1.,MnSi0.5O2 )
           0.0,           // W(KAlSiO4   ,MnSi0.5O2 )
           0.0,           // W(Ca3(PO4)2 ,MnSi0.5O2 )
           0.0,           // W(H2O       ,MnSi0.5O2 )

           0.0,           // W(NiSi0.5O2 ,Mg2SiO4   )
           0.0,           // W(CoSi0.5O2 ,Mg2SiO4   )
         574.1,           // W(Ca2Si2O6  ,Mg2SiO4   )
        9272.3,           // W(NaSi0.5O1.,Mg2SiO4   )
       36512.7,           // W(KAlSiO4   ,Mg2SiO4   )
       -6766.8,           // W(Ca3(PO4)2 ,Mg2SiO4   )
       24630.1,           // W(H2O       ,Mg2SiO4   )

           0.0,           // W(CoSi0.5O2 ,NiSi0.5O2 )
           0.0,           // W(Ca2Si2O6  ,NiSi0.5O2 )
           0.0,           // W(NaSi0.5O1.,NiSi0.5O2 )
           0.0,           // W(KAlSiO4   ,NiSi0.5O2 )
           0.0,           // W(Ca3(PO4)2 ,NiSi0.5O2 )
           0.0,           // W(H2O       ,NiSi0.5O2 )

           0.0,           // W(Ca2Si2O6  ,CoSi0.5O2 )
           0.0,           // W(NaSi0.5O1.,CoSi0.5O2 )
           0.0,           // W(KAlSiO4   ,CoSi0.5O2 )
           0.0,           // W(Ca3(PO4)2 ,CoSi0.5O2 )
           0.0,           // W(H2O       ,CoSi0.5O2 )

        7430.3,           // W(NaSi0.5O1.,Ca2Si2O6  )
       19927.4,           // W(KAlSiO4   ,Ca2Si2O6  )
       88993.1,           // W(Ca3(PO4)2 ,Ca2Si2O6  )
       -1583.7,           // W(H2O       ,Ca2Si2O6  )

       -1102.3,           // W(KAlSiO4   ,NaSi0.5O1.)
      -13062.6,           // W(Ca3(PO4)2 ,NaSi0.5O1.)
       13043.1,           // W(H2O       ,NaSi0.5O1.)

       85064.0,           // W(Ca3(PO4)2 ,KAlSiO4   )
       35572.8,           // W(H2O       ,KAlSiO4   )

       53448.7            // W(H2O       ,Ca3(PO4)2 )
};

#define NR      14    // Number of independent mole fraction variables
#define NA      15    // Number of liquid components
#define NATOMS 5.0    // Average number of atoms in the formula unit (guess)
#define SMX  SHRT_MAX

#define corrH01  69752.597
#define corrH02 -16940.793
#define corrH03  43463.013
#define corrH04 163652.907
#define corrH05   7238.581
#define corrH06  13739.014
#define corrH07      0.000
#define corrH08   7911.937
#define corrH09      0.000
#define corrH10      0.000
#define corrH11 -21862.658
#define corrH12 -20028.750
#define corrH13 -16015.286
#define corrH14  13197.286
#define corrH15      0.000
#define corrH16      0.000
#define corrH17      0.000
#define corrH18      0.000

#define corrS01     45.580
#define corrS02      0.000
#define corrS03    -25.166
#define corrS04     25.413
#define corrS05      0.000
#define corrS06      6.424
#define corrS07      0.000
#define corrS08     12.777
#define corrS09      0.000
#define corrS10      0.000
#define corrS11     -2.266
#define corrS12      0.000
#define corrS13     -7.311
#define corrS14      0.000
#define corrS15      0.000
#define corrS16      0.000
#define corrS17      0.000
#define corrS18      0.000

#define corrV01      0.000
#define corrV02      0.000
#define corrV03      0.000
#define corrV04      0.000
#define corrV05      0.000
#define corrV06      0.000
#define corrV07      0.000
#define corrV08      0.000
#define corrV09      0.000
#define corrV10      0.000
#define corrV11      0.000
#define corrV12      0.000
#define corrV13      0.000
#define corrV14      0.000
#define corrV15      0.000
#define corrV16      0.000
#define corrV17      0.000
#define corrV18      0.000

+(void)initialize {
	if (self == [LiquidpMelts class]) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"Initialize(LiquidpMelts) - entry ...");
		NSMutableArray *mutableEndmembers = [NSMutableArray arrayWithCapacity:NA];

#define SIO2_MULT 4.0
		LiquidpMeltsGenericEM *si4o8 = [[LiquidpMeltsGenericEM alloc] initWithH:-906377.0*SIO2_MULT+corrH01
																			S:46.029*SIO2_MULT+corrS01
																		   k0:83.51*SIO2_MULT
																		   k1:-3.747e2*SIO2_MULT
																		   k2:-24.554e5*SIO2_MULT
																		   k3:28.007e7*SIO2_MULT
																		   l1:0.0
																		   l2:0.0
																		   Tt:0.0
																	   deltaH:0.0
																		 vLiq:2.690*SIO2_MULT+corrV01
																	  dvdtLiq:0.0
																	  dvdpLiq:-1.89e-5*SIO2_MULT
																   d2vdtdpLiq:1.3e-8*SIO2_MULT
																	d2vdp2Liq:3.6e-10*SIO2_MULT
																	  tFusion:1999.0
																	  sFusion:4.46*SIO2_MULT
																		cpLiq:82.6*SIO2_MULT];
#undef SIO2_MULT
		[si4o8 setPhaseName:@"Si4O8"];
		[si4o8 setPhaseFormula:@"Si4O8"];
		[mutableEndmembers addObject:si4o8];
		if (debug) NSLog(@"... allocated si4o8 ...");

		LiquidpMeltsGenericEM *tio2 = [[LiquidpMeltsGenericEM alloc] initWithH:-944750.0+corrH02
																		   S:50.460+corrS02
																		  k0:77.84
																		  k1:0.0
																		  k2:-33.678e5
																		  k3:40.294e7
																		  l1:0.0
																		  l2:0.0
																		  Tt:0.0
																	  deltaH:0.0
																		vLiq:2.316+corrV02
																	 dvdtLiq:7.246e-4
																	 dvdpLiq:-2.310e-5
																  d2vdtdpLiq:0.0
																   d2vdp2Liq:5.0e-10
																	 tFusion:1870.0
																	 sFusion:35.824
																	   cpLiq:109.2];
		[tio2 setPhaseName:@"TiO2"];
		[tio2 setPhaseFormula:@"TiO2"];
		[mutableEndmembers addObject:tio2];
		if (debug) NSLog(@"... allocated tio2 ...");

#define AL2O3_MULT 2.0
		LiquidpMeltsGenericEM *al4o6 = [[LiquidpMeltsGenericEM alloc] initWithH:-1675700.0*AL2O3_MULT+corrH03
																			S:50.82*AL2O3_MULT+corrS03
																		   k0:155.02*AL2O3_MULT
																		   k1:-8.284E2*AL2O3_MULT
																		   k2:-38.614E5*AL2O3_MULT
																		   k3:40.908E7*AL2O3_MULT
																		   l1:0.0
																		   l2:0.0
																		   Tt:0.0
																	   deltaH:0.0
																		 vLiq:3.711*AL2O3_MULT+corrV03
																	  dvdtLiq:2.62e-4*AL2O3_MULT
																	  dvdpLiq:-2.26e-5*AL2O3_MULT
																   d2vdtdpLiq:2.7e-8*AL2O3_MULT
																	d2vdp2Liq:4.0e-10*AL2O3_MULT
																	  tFusion:2319.65
																	  sFusion:48.61*AL2O3_MULT
																		cpLiq:170.3*AL2O3_MULT];
#undef AL2O3_MULT
		[al4o6 setPhaseName:@"Al4O6"];
		[al4o6 setPhaseFormula:@"Al4O6"];
		[mutableEndmembers addObject:al4o6];
		if (debug) NSLog(@"... allocated al4o6 ...");

		LiquidpMeltsGenericEM *fe2o3 = [[LiquidpMeltsGenericEM alloc] initWithH:-822000.00+corrH04
																			S:87.40+corrS04
																		   k0:146.86
																		   k1:0.0
																		   k2:-55.768E5
																		   k3:52.563E7
																		   l1:-7.403E-2
																		   l2:27.921E-5
																		   Tt:955.0
																	   deltaH:1287.0
																		 vLiq:4.213+corrV04
																	  dvdtLiq:9.09E-4
																	  dvdpLiq:-2.53E-5
																   d2vdtdpLiq:3.1E-8
																	d2vdp2Liq:4.4e-10
																	  tFusion:1895.0
																	  sFusion:60.41
																		cpLiq:240.9];
		[fe2o3 setPhaseName:@"Fe2O3"];
		[fe2o3 setPhaseFormula:@"Fe2O3"];
		[mutableEndmembers addObject:fe2o3];
		if (debug) NSLog(@"... allocated fe2o3 ...");

		LiquidpMeltsGenericEM *mgcr2o4 = [[LiquidpMeltsGenericEM alloc] initWithH:-1783640.0+corrH05
																			  S:106.02+corrS05
																			 k0:201.981
																			 k1:-5.519E2
																			 k2:-57.844E5
																			 k3:57.729E7
																			 l1:0.0
																			 l2:0.0
																			 Tt:0.0
																		 deltaH:0.0
																		   vLiq:5.358+corrV05
																		dvdtLiq:11.71e-4
																		dvdpLiq:-2.26e-5
																	 d2vdtdpLiq:1.8e-8
																	  d2vdp2Liq:4.67e-10
																		tFusion:2673.15
																		sFusion:73.22
																		  cpLiq:335.1];
		[mgcr2o4 setPhaseName:@"MgCr2O4"];
		[mgcr2o4 setPhaseFormula:@"MgCr2O4"];
		[mutableEndmembers addObject:mgcr2o4];
		if (debug) NSLog(@"... allocated mgcr2o4 ...");

		LiquidpMeltsGenericEM *fe2sio4 = [[LiquidpMeltsGenericEM alloc] initWithH:-1479360.0+corrH06
																			  S:150.930+corrS06
																			 k0:248.93
																			 k1:-19.239E2
																			 k2:0.0
																			 k3:-13.910E7
																			 l1:0.0
																			 l2:0.0
																			 Tt:0.0
																		 deltaH:0.0
																		   vLiq:5.420+corrV06
																		dvdtLiq:5.84e-4
																		dvdpLiq:-2.79e-5
																	 d2vdtdpLiq:-2.3e-8
																	  d2vdp2Liq:14.6e-10
																		tFusion:1490.0
																		sFusion:59.9
																		  cpLiq:240.2];
		[fe2sio4 setPhaseName:@"Fe2SiO4"];
		[fe2sio4 setPhaseFormula:@"Fe2SiO4"];
		[mutableEndmembers addObject:fe2sio4];
		if (debug) NSLog(@"... allocated fe2sio4 ...");

		LiquidpMeltsGenericEM *mn2sio4 = [[LiquidpMeltsGenericEM alloc] initWithH:-1732000.0/2.0+corrH07
																			  S:155.9/2.0+corrS07
																			 k0:219.89/2.0
																			 k1:-12.710E2/2.0
																			 k2:-20.496E5/2.0
																			 k3:17.652E7/2.0
																			 l1:0.0
																			 l2:0.0
																			 Tt:0.0
																		 deltaH:0.0
																		   vLiq:2.84+corrV07
																		dvdtLiq:5.84e-4/2.0
																		dvdpLiq:-2.79e-5/2.0
																	 d2vdtdpLiq:-2.3e-8/2.0
																	  d2vdp2Liq:14.6e-10/2.0
																		tFusion:1620.0
																		sFusion:27.6
																		  cpLiq:121.6];
		[mn2sio4 setPhaseName:@"MnSi0.5O2"];
		[mn2sio4 setPhaseFormula:@"MnSi0.5O2"];
		[mutableEndmembers addObject:mn2sio4];
		if (debug) NSLog(@"... allocated mn2sio4 ...");

		LiquidpMeltsGenericEM *mg2sio4 = [[LiquidpMeltsGenericEM alloc] initWithH:-2174420.0+corrH08
																			  S:94.010+corrS08
																			 k0:238.64
																			 k1:-20.013E2
																			 k2:0.0
																			 k3:-11.624E7
																			 l1:0.0
																			 l2:0.0
																			 Tt:0.0
																		 deltaH:0.0
																		   vLiq:4.980+corrV08
																		dvdtLiq:5.24e-4
																		dvdpLiq:-1.35e-5
																	 d2vdtdpLiq:-1.3e-8
																	  d2vdp2Liq:4.14e-10
																		tFusion:2163.0
																		sFusion:57.2
																		  cpLiq:271.0];
		[mg2sio4 setPhaseName:@"Mg2SiO4"];
		[mg2sio4 setPhaseFormula:@"Mg2SiO4"];
		[mutableEndmembers addObject:mg2sio4];
		if (debug) NSLog(@"... allocated mgsio4 ...");

		LiquidpMeltsGenericEM *ni2sio4 = [[LiquidpMeltsGenericEM alloc] initWithH:-1395300.0/2.0+corrH09
																			  S:128.1/2.0+corrS09
																			 k0:214.997/2.0
																			 k1:-10.3075E2/2.0
																			 k2:-49.4453E5/2.0
																			 k3:62.375E7/2.0
																			 l1:0.0
																			 l2:0.0
																			 Tt:0.0
																		 deltaH:0.0
																		   vLiq:2.48+corrV09
																		dvdtLiq:5.84e-4/2.0
																		dvdpLiq:-2.79e-5/2.0
																	 d2vdtdpLiq:-2.3e-8/2.0
																	  d2vdp2Liq:14.6e-10/2.0
																		tFusion:1923.0
																		sFusion:29.0
																		  cpLiq:119.3];
		[ni2sio4 setPhaseName:@"NiSi0.5O2"];
		[ni2sio4 setPhaseFormula:@"NiSi0.5O2"];
		[mutableEndmembers addObject:ni2sio4];
		if (debug) NSLog(@"... allocated ni2sio4 ...");

		LiquidpMeltsGenericEM *co2sio4 = [[LiquidpMeltsGenericEM alloc] initWithH:-1414100.0/2.0+corrH10
																			  S:142.6/2.0+corrS10
																			 k0:201.048/2.0
																			 k1:-0.369E2/2.0
																			 k2:-71.81E5/2.0
																			 k3:90.05E7/2.0
																			 l1:0.0
																			 l2:0.0
																			 Tt:0.0
																		 deltaH:0.0
																		   vLiq:2.30+corrV10
																		dvdtLiq:5.84e-4/2.0
																		dvdpLiq:-2.79e-5/2.0
																	 d2vdtdpLiq:-2.3e-8/2.0
																	  d2vdp2Liq:14.6e-10/2.0
																		tFusion:1688.0
																		sFusion:29.0
																		  cpLiq:125.3];
		[co2sio4 setPhaseName:@"CoSi0.5O2"];
		[co2sio4 setPhaseFormula:@"CoSi0.5O2"];
		[mutableEndmembers addObject:co2sio4];
		if (debug) NSLog(@"... allocated co2sio4 ...");

#define CASIO3_MULT 2.0
		LiquidpMeltsGenericEM *casio3 = [[LiquidpMeltsGenericEM alloc] initWithH:-1627427.0*CASIO3_MULT+corrH11
																			 S:85.279*CASIO3_MULT+corrS11
																			k0:141.16*CASIO3_MULT
																			k1:-4.172e2*CASIO3_MULT
																			k2:-58.576e5*CASIO3_MULT
																			k3:94.074e7*CASIO3_MULT
																			l1:0.0
																			l2:0.0
																			Tt:0.0
																		deltaH:0.0
																		  vLiq:4.347*CASIO3_MULT+corrV11
																	   dvdtLiq:2.92e-4*CASIO3_MULT
																	   dvdpLiq:-1.55e-5*CASIO3_MULT
																	d2vdtdpLiq:-1.6e-8*CASIO3_MULT
																	 d2vdp2Liq:3.89e-10*CASIO3_MULT
																	   tFusion:1817.0
																	   sFusion:31.5*CASIO3_MULT
																		 cpLiq:172.4*CASIO3_MULT];
#undef CASIO3_MULT
		[casio3 setPhaseName:@"Ca2Si2O6"];
		[casio3 setPhaseFormula:@"Ca2Si2O6"];
		[mutableEndmembers addObject:casio3];
		if (debug) NSLog(@"... allocated ca2si2o6 ...");

#define NA2SIO3_MULT 0.5
		LiquidpMeltsGenericEM *na2sio3 = [[LiquidpMeltsGenericEM alloc] initWithH:-373190.0*4.184*NA2SIO3_MULT+corrH12
																			  S:27.21*4.184*NA2SIO3_MULT+corrS12
																			 k0:234.77*NA2SIO3_MULT
																			 k1:-22.189E2*NA2SIO3_MULT
																			 k2:0.0*NA2SIO3_MULT
																			 k3:13.530E7*NA2SIO3_MULT
																			 l1:0.0
																			 l2:0.0
																			 Tt:0.0
																		 deltaH:0.0
																		   vLiq:5.568*NA2SIO3_MULT+corrV12
																		dvdtLiq:7.41e-4*NA2SIO3_MULT
																		dvdpLiq:-4.29e-5*NA2SIO3_MULT
																	 d2vdtdpLiq:-5.3e-8*NA2SIO3_MULT
																	  d2vdp2Liq:8.4e-10*NA2SIO3_MULT
																		tFusion:1361.0
																		sFusion:38.34*NA2SIO3_MULT
																		  cpLiq:180.2*NA2SIO3_MULT];
#undef NA2SIO3_MULT
		[na2sio3 setPhaseName:@"NaSi0.5O1.5"];
		[na2sio3 setPhaseFormula:@"NaSi0.5O1.5"];
		[mutableEndmembers addObject:na2sio3];
		if (debug) NSLog(@"... allocated nasi0.5o1.5 ...");

		LiquidpMeltsGenericEM *kalsio4 = [[LiquidpMeltsGenericEM alloc] initWithH:-2111813.55+corrH13
																			  S:133.9653+corrS13
																			 k0:186.0
																			 k1:0.0
																			 k2:-131.067E5
																			 k3:213.893E7
																			 l1:-7.096454E-2
																			 l2:21.682E-5
																			 Tt:800.15
																		 deltaH:1154.0
																		   vLiq:6.8375+corrV13
																		dvdtLiq:7.265e-4
																		dvdpLiq:-6.395e-5
																	 d2vdtdpLiq:-4.6e-8
																	  d2vdp2Liq:12.1e-10
																		tFusion:2023.15
																		sFusion:24.5
																		  cpLiq:217.0];
		[kalsio4 setPhaseName:@"KAlSiO4"];
		[kalsio4 setPhaseFormula:@"KAlSiO4"];
		[mutableEndmembers addObject:kalsio4];
		if (debug) NSLog(@"... allocated kalsio4 ...");

		LiquidpMeltsGenericEM *ca3po4 = [[LiquidpMeltsGenericEM alloc] initWithH:-4097169.0+corrH14
																			 S:235.978+corrS14
																			k0:402.997
																			k1:-28.0835E2
																			k2:0.0
																			k3:-32.6230E7
																			l1:2.5427E-2
																			l2:19.255E-5
																			Tt:1373.0
																		deltaH:14059.0
																		  vLiq:10.7382+corrV14
																	   dvdtLiq:0.0
																	   dvdpLiq:0.0
																	d2vdtdpLiq:0.0
																	 d2vdp2Liq:0.0
																	   tFusion:1943.15
																	   sFusion:35.690
																		 cpLiq:574.67];
		[ca3po4 setPhaseName:@"Ca3(PO4)2"];
		[ca3po4 setPhaseFormula:@"Ca3(PO4)2"];
		[mutableEndmembers addObject:ca3po4];
		if (debug) NSLog(@"... allocated ca3(po4)2 ...");

		LiquidpMeltsH2O *h2o = [[LiquidpMeltsH2O alloc] init];
		[h2o setPhaseName:@"H2O"];     // must be set here, because not initialized in the class instance
		[h2o setPhaseFormula:@"H2O"];  // (above)
		[mutableEndmembers addObject:h2o];
		if (debug) NSLog(@"... allocated h2o and exiting.");

		endmembers = [NSArray arrayWithArray:mutableEndmembers];
	}
}

#undef corrH01
#undef corrH02
#undef corrH03
#undef corrH04
#undef corrH05
#undef corrH06
#undef corrH07
#undef corrH08
#undef corrH09
#undef corrH10
#undef corrH11
#undef corrH12
#undef corrH13
#undef corrH14
#undef corrH15
#undef corrH16
#undef corrH17
#undef corrH18

#undef corrS01
#undef corrS02
#undef corrS03
#undef corrS04
#undef corrS05
#undef corrS06
#undef corrS07
#undef corrS08
#undef corrS09
#undef corrS10
#undef corrS11
#undef corrS12
#undef corrS13
#undef corrS14
#undef corrS15
#undef corrS16
#undef corrS17
#undef corrS18

#undef corrV01
#undef corrV02
#undef corrV03
#undef corrV04
#undef corrV05
#undef corrV06
#undef corrV07
#undef corrV08
#undef corrV09
#undef corrV10
#undef corrV11
#undef corrV12
#undef corrV13
#undef corrV14
#undef corrV15
#undef corrV16
#undef corrV17
#undef corrV18

-(id)init {
	if ((self = [super init])) {
		BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
		if (debug) NSLog(@"init(LiquidpMelts) ... entry ...");
		computeMixingQuantities = NO;
		[self setPhaseName:@"Liquid"];
		if (debug) NSLog(@"... exiting.");
	}
	return self;
}

/**
 "C" code from MELTS
 */

/*
 * Array to convert W(i,j) indexes to entries in the array of structures
 * modelParameters[k], where k is:
 */

static const short int Index[NA][NA] = {
	{ SMX,   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13 },
	{   0, SMX,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26 },
	{   1,  14, SMX,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38 },
	{   2,  15,  27, SMX,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49 },
	{   3,  16,  28,  39, SMX,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59 },
	{   4,  17,  29,  40,  50, SMX,  60,  61,  62,  63,  64,  65,  66,  67,  68 },
	{   5,  18,  30,  41,  51,  60, SMX,  69,  70,  71,  72,  73,  74,  75,  76 },
	{   6,  19,  31,  42,  52,  61,  69, SMX,  77,  78,  79,  80,  81,  82,  83 },
	{   7,  20,  32,  43,  53,  62,  70,  77, SMX,  84,  85,  86,  87,  88,  89 },
	{   8,  21,  33,  44,  54,  63,  71,  78,  84, SMX,  90,  91,  92,  93,  94 },
	{   9,  22,  34,  45,  55,  64,  72,  79,  85,  90, SMX,  95,  96,  97,  98 },
	{  10,  23,  35,  46,  56,  65,  73,  80,  86,  91,  95, SMX,  99, 100, 101 },
	{  11,  24,  36,  47,  57,  66,  74,  81,  87,  92,  96,  99, SMX, 102, 103 },
	{  12,  25,  37,  48,  58,  67,  75,  82,  88,  93,  97, 100, 102, SMX, 104 },
	{  13,  26,  38,  49,  59,  68,  76,  83,  89,  94,  98, 101, 103, 104, SMX }
};

#undef SMX
#define WH(i,j) modelParameters[Index[i][j]]
#define WS(i,j) 0.0
#define WV(i,j) 0.0

#define R 8.3143

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
	const char *NAMES[NA] = { "SiO2", "TiO2", "Al2O3", "Fe2O3", "Cr2O3", "FeO", "MnO", "MgO",
		"NiO", "CoO", "CaO", "Na2O", "K2O", "P2O5", "H2O" };
	const char *FORMULAS[NA] = { "Si4O8", "TiO2", "Al4O6", "Fe2O3", "MgCr2O4", "Fe2SiO4", "MnSi0.5O2", "Mg2SiO4",
		"NiSi0.5O2", "CoSi0.5O2", "Ca2Si2O6", "NaSi0.5O1.5", "KAlSiO4", "Ca3(PO4)2", "H2O" };
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
	 endmember olivine components.
	 (2) calculates from a vector of moles of endmember components, one or
	 all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
	 (3) calculates from a vector of independent compositional variables
	 mole fractions of endmember components and/or the Jacobian matrix
	 dx[]/dr[]
	 ----------------------------------------------------------------------------*/

	int i, j, k;

	if (inpMask == FIRST && outMask == SECOND) {
		/* Converts a vector of moles of elements into a vector of moles of
		 end-member components. */
		static const int  H =  1;
		static const int  O =  8;
		static const int Na = 11;
		static const int Mg = 12;
		static const int Al = 13;
		static const int Si = 14;
		static const int  P = 15;
		static const int  K = 19;
		static const int Ca = 20;
		static const int Ti = 22;
		static const int Cr = 24;
		static const int Mn = 25;
		static const int Fe = 26;
		static const int Co = 27;
		static const int Ni = 28;
		double ox;


		m[ 0] = (e[Si] - e[Mn]/2.0 - (e[Mg]-e[Cr]/2.0)/2.0 - e[Ni]/2.0 - e[Co]/2.0 - (e[Ca]-3.0*e[P]/2.0) - e[Na]/2.0 - e[K])/4.0; // Si4O8
		m[ 1] = e[Ti];                    // TiO2
		m[ 2] = (e[Al]-e[K])/4.0;         // Al4O6
		m[ 4] = e[Cr]/2.0;                // MgCr2O4
		m[ 6] = e[Mn];                    // MnSi0.5O2
		m[ 7] = (e[Mg]-e[Cr]/2.0)/2.0;    // Mg2SiO4
		m[ 8] = e[Ni];                    // NiSi0.5O2
		m[ 9] = e[Co];                    // CoSi0.5O2
		m[10] = (e[Ca]-3.0*e[P]/2.0)/2.0; // Ca2Si2O6
		m[11] = e[Na];                    // NaSi0.5O1.5
		m[12] = e[K];                     // KAlSiO4
		m[13] = e[P]/2.0;                 // Ca3(PO4)2
		m[14] = e[H]/2.0;                 // H2O

		ox = 8.0*m[0] + 2.0*m[1] + 6.0*m[2] + 4.0*m[4] + 2.0*m[6] + 4.0*m[7] + 2.0*m[8]
		   + 2.0*m[9] + 6.0*m[10] + 3.0*m[11]/2.0 + 4.0*m[12] + 8.0*m[13] + m[14];

		if (e[Fe] > 10.0*DBL_EPSILON) {
			m[3] = e[O] - ox - e[Fe];       // Fe2O3
			if (fabs(m[3]) < 100.0*DBL_EPSILON) m[3] = 0.0;
			m[5] = (e[Fe] - 2.0*m[3])/2.0;  // Fe2SiO4
			if (fabs(m[5]) < 100.0*DBL_EPSILON) m[5] = 0.0;
			m[0] -= m[5]/4.0;               // readjust Si4O8 moles for the Fe2SiO4 just made
		} else {
			m[3] = 0.0;
			m[5] = 0.0;
		}


	} else if (inpMask == SECOND) {
		double sum;

		if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
			NSLog(@"Illegal call to conLiq_v34 with inpMask = %o and outMask = %o", inpMask, outMask);

		for (i=0, sum=0.0; i<NA; i++) sum += m[i];

		if (outMask & THIRD) {
			/* Converts a vector of moles of end-member components (m) into a vector
			 of independent compositional variables (r) required as input for the
			 remaining public functions.
			 The dependent variable is taken to be SiO2 (1st component), as this
			 component will never have a mole fraction of zero.                   */

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
			NSLog(@"Illegal call to convert with inpMask = %o and outMask = %o", inpMask, outMask);

		if (outMask & FOURTH) {
			/* Converts a vector of independent compositional variables (r)
			 into a vector of mole fractions of end-member components (x)            */

			for (i=0, x[0] = 1.0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }
		}

		if (outMask & SEVENTH) {
			/* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
			for (j=0; j<NR; j++) dr[0][j] = -1.0;
			for (i=1; i<NA; i++) {
				for (j=0; j<NR; j++) dr[i][j] = 0.0;
				dr[i][i-1] = 1.0;
			}
		}

	} else {
		NSLog(@"Illegal call to conLiq_v34 with inpMask = %o and outMask = %o", inpMask, outMask);
	}

}

-(NSString *)displayFormula:(double)t
						  p:(double)p
						  r:(double [NA])r
{
	double sio2 = (1.0-r[0]-r[1]-r[2]-r[3]-r[4]-r[5]-r[6]-r[7]-r[8]-r[9]-r[10]-r[11]-r[12]-r[13])*4.0*60.0843
	            + (r[4] + r[5]/2.0 + r[6] + r[7]/2.0 + r[8]/2.0 + 2.0*r[9] + r[10]/2.0 + r[11])*60.0843;
	double tio2  = r[0]*79.8658;
	double al2o3 = (r[1]*2.0 + r[11]/2.0)*101.9613;
	double fe2o3 = r[2]*159.6882;
	double cr2o3 = r[3]*151.9904;
	double feo   = 2.0*r[4]*71.8444;
	double mno   = r[5]*70.93745;
	double mgo   = (r[3] + 2.0*r[6])*40.3044;
	double nio   = r[7]*74.6928;
	double coo   = r[8]*74.9326;
	double cao   = (r[9]*2.0 + 3.0*r[12])*56.0774;
	double na2o  = r[10]*61.97894/2.0;
	double k2o   = r[11]*94.196/2.0;
	double p2o5  = r[12]*141.9445;
	double h2o   = r[13]*18.01528;
	double sum = sio2 + tio2 + al2o3 + fe2o3 + cr2o3 + feo + mno + mgo + nio + coo + cao + na2o + k2o + p2o5 + h2o;
	if (sum == 0.0) sum = 1.0;

	return [NSString stringWithFormat:@"wt%%:SiO2 %5.2f TiO2 %5.2f Al2O3 %5.2f Fe2O3 %5.2f Cr2O3 %5.2f FeO %5.2f MnO %5.2f MgO %5.2f NiO %5.2f CoO %5.2f CaO %5.2f Na2O %5.2f K2O %5.2f P2O5 %5.2f H2O %5.2f",
			sio2*100.0/sum, tio2*100.0/sum, al2o3*100.0/sum, fe2o3*100.0/sum, cr2o3*100.0/sum, feo*100.0/sum, mno*100.0/sum,
			mgo*100.0/sum, nio*100.0/sum, coo*100.0/sum, cao*100.0/sum, na2o*100.0/sum, k2o*100.0/sum, p2o5*100.0/sum, h2o*100.0/sum];
}

-(void)activity:(int)mask
			  t:(double)t
			  p:(double)p
			  r:(double [NA])r
			  a:(double [NA])a      // (pointer to a[]) activities              BINARY MASK: 0001
			 mu:(double [NA])mu     // (pointer to mu[]) chemical potentials    BINARY MASK: 0010
			 dx:(double [NA][NR])dr // (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100
{
	double x[NA], gex;
	int i, j;

	/* x[0] --> x[NA-1] is an array of mole fractions of liquid components      */
	for (x[0] = 1.0, i=0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }
	for (gex=0.0, i=0; i<NA; i++) {
		for (j=i+1; j<NA; j++) gex += x[i]*x[j]*(WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j));
	}

	if (mask & FIRST) {
		for (i=0; i<NA; i++) {
			a[i] = - gex;
			for (j=0;   j<i;  j++) a[i] += (WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j))*x[j];
			for (j=i+1; j<NA; j++) a[i] += (WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j))*x[j];
			a[i] = (x[i] != 0.0) ? x[i]*exp(a[i]/(R*t)) : 0.0;
			a[i] = (i != NA-1)   ? (1.0 - x[NA-1])*a[i] : x[NA-1]*a[NA-1];
		}
	}

	if (mask & SECOND) {
		for (i=0; i<NA; i++) {
			mu[i] = - gex;
			for (j=0;   j<i;  j++) mu[i] += (WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j))*x[j];
			for (j=i+1; j<NA; j++) mu[i] += (WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j))*x[j];
			mu[i] = (x[i] != 0.0) ? mu[i] + R*t*log(x[i]) : 0.0;
			if (i != NA-1)           mu[i]    += R*t*log(1.0-x[NA-1]);
			else if (x[NA-1] != 0.0) mu[NA-1] += R*t*log(x[NA-1]);
		}
	}

	if (mask & THIRD) {
		double a[NA], dgexdr[NR];

		for (i=0; i<NA; i++) {
			a[i] = - gex;
			for (j=0;   j<i;  j++) a[i] += (WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j))*x[j];
			for (j=i+1; j<NA; j++) a[i] += (WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j))*x[j];
			a[i] = (x[i] != 0.0) ? x[i]*exp(a[i]/(R*t)) : 0.0;
			a[i] = (i != NA-1)   ? (1.0 - x[NA-1])*a[i] : x[NA-1]*a[NA-1];
		}
		for (i=0; i<NR; i++) {
			for (dgexdr[i] = x[0]*(WH(0,i+1)-t*WS(0,i+1)+(p-1.0)*WV(0,i+1)), j=0;
				 j<NR; j++) dgexdr[i] += (i != j) ? r[j]*((WH(i+1,j+1)-t*WS(i+1,j+1)
														   +(p-1.0)*WV(i+1,j+1)) - (WH(0,j+1)-t*WS(0,j+1)+(p-1.0)*WV(0,j+1)))
				: - r[j]*(WH(0,j+1)-t*WS(0,j+1)+(p-1.0)*WV(0,j+1));
		}

		/* Special case for component 0 (SiO2)                                    */
		for (j=0; j<NR; j++) {
			dr[0][j] = - R*t/x[0] + (WH(0,j+1)-t*WS(0,j+1)+(p-1.0)*WV(0,j+1))
			- dgexdr[j];
			if (j == NR-1) dr[0][j] += - R*t/(1.0-x[NA-1]);
			dr[0][j] *= a[0]/(R*t);
		}

		/* All other cases                                                        */
		for (i=1; i<NA; i++) {
			for (j=0; j<NR; j++) {
				dr[i][j] = (i == (j+1)) ? ((r[j] != 0.0) ? R*t/r[j] : 0.0) :
				WH(i,j+1)-t*WS(i,j+1)+(p-1.0)*WV(i,j+1);
				dr[i][j] += - (WH(0,i)-t*WS(0,i)+(p-1.0)*WV(0,i)) - dgexdr[j];
				if (i != NA-1 && j == NR-1) dr[i][NR-1] += - R*t/(1.0-x[NA-1]);
				else if (i == NA-1 && j == NR-1 && x[NA-1] != 0.0)
					dr[NA-1][NR-1] += R*t/x[NA-1];
				dr[i][j] *= a[i]/(R*t);
			}
		}

	}

}

-(void)gmix:(int)mask
		  t:(double)t
		  p:(double)p
		  r:(double [NR])r
	   gmix:(double *)gmix           // Gibbs energy of mixing             BINARY MASK: 0001
		 dx:(double [NR])dr          // (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010
		dx2:(double [NR][NR])dr2     // (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100
		dx3:(double [NR][NR][NR])dx3 // (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000
{
	double x[NA];
	int i, j;

	/* x[0] --> x[NA-1] is an array of mole fractions of liquid components      */
	for (x[0] = 1.0, i=0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }

	if (mask & FIRST) {
		for (*gmix = 0.0, i=0; i<NA; i++) {
			for (j=i+1; j<NA; j++) *gmix += x[i]*x[j]*(WH(i,j) - t*WS(i,j)
													   + (p-1.0)*WV(i,j));
			*gmix += (x[i] != 0.0) ? R*t*x[i]*log(x[i]) : 0.0;
		}
		*gmix += (x[NA-1] != 0.0) ?
		R*t*(x[NA-1]*log(x[NA-1]) + (1.0-x[NA-1])*log(1.0-x[NA-1])) : 0.0;
	}

	if(mask & SECOND) {
		for (i=0; i<NR; i++) {
			dr[i] = (r[i] != 0.0) ? R*t*(log(r[i]) - log(x[0])) : 0.0;
			dr[i] += x[0]*(WH(0,i+1)-t*WS(0,i+1)+(p-1.0)*WV(0,i+1));
			for (j=0; j<NR; j++)
				dr[i] += (i != j) ? r[j]*((WH(i+1,j+1) - t*WS(i+1,j+1)
										   + (p-1.0)*WV(i+1,j+1))
										  - (WH(0,j+1)-t*WS(0,j+1)+(p-1.0)*WV(0,j+1)))
				: - r[j]*(WH(0,j+1)-t*WS(0,j+1)+(p-1.0)*WV(0,j+1));
		}
		dr[NR-1] += (x[NA-1] != 0.0) ? R*t*(log(x[NA-1])-log(1.0-x[NA-1])) : 0.0;
	}

	if(mask & THIRD) {
		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) {
				dr2[i][j] = R*t/x[0] - (WH(0,i+1)-t*WS(0,i+1)+(p-1.0)*WV(0,i+1))
				- (WH(0,j+1)-t*WS(0,j+1)+(p-1.0)*WV(0,j+1));
				dr2[i][j] += (i != j) ? (WH(i+1,j+1)-t*WS(i+1,j+1)+(p-1.0)*WV(i+1,j+1))
				: 0.0;
			}
			dr2[i][i] += (r[i] != 0.0) ? R*t/r[i] : 0.0;
		}
		dr2[NR-1][NR-1] += (x[NA-1] != 0.0) ? R*t*(1.0/x[NA-1] + 1.0/(1.0-x[NA-1])) : 0.0;
	}

	// To be completed
	if(mask & FOURTH) {
		int k;
		for (i=0; i<NR; i++) for (j=0; j<NR; j++) for (k=0; k<NR; k++) dx3[i][j][k] = 0.0;
	}
}

-(void)hmix:(int)mask
		  t:(double)t
		  p:(double)p
		  r:(double [NR])r
	   hmix:(double *)hmix // Enthalpy of mixing BINARY MASK: 1
{
	double x[NA];
	int i, j;

	/* x[0] --> x[NA-1] is an array of mole fractions of liquid components      */
	for (x[0] = 1.0, i=0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }

	for (*hmix = 0.0, i=0; i<NA; i++) {
		for (j=i+1; j<NA; j++) *hmix += x[i]*x[j]*(WH(i,j)+(p-1.0)*WV(i,j));
	}
}

-(void)smix:(int)mask
		  t:(double)t
		  p:(double)p
		  r:(double [NR])r
	   smix:(double *)smix        // Entropy of mixing                  BINARY MASK: 001
		 dx:(double [NR])dr       // (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010
		dx2:(double [NR][NR])dr2  // (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100
{
	double x[NA];
	int i, j;

	/* x[0] --> x[NA-1] is an array of mole fractions of liquid components      */
	for (x[0] = 1.0, i=0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }

	if (mask & FIRST) {
		for (*smix = 0.0, i=0; i<NA; i++) {
			for (j=i+1; j<NA; j++) *smix += x[i]*x[j]*WS(i,j);
			*smix += (x[i] != 0.0) ? - R*x[i]*log(x[i]) : 0.0;
		}
		*smix += (x[NA-1] != 0.0) ?
		-R*(x[NA-1]*log(x[NA-1]) + (1.0-x[NA-1])*log(1.0-x[NA-1])) : 0.0;
	}

	if(mask & SECOND) {
		for (i=0; i<NR; i++) {
			dr[i] = (r[i] != 0.0) ? R*(log(x[0]) - log(r[i])) : 0.0;
			dr[i] += x[0]*WS(0,i+1);
			for (j=0; j<NR; j++)
				dr[i] += (i != j) ? r[j]*(WS(i+1,j+1) - WS(0,j+1)) : - r[j]*WS(0,j+1);
		}
		dr[NR-1] += (x[NA-1] != 0.0) ? -R*(log(x[NA-1])-log(1.0-x[NA-1])) : 0.0;
	}

	if(mask & THIRD) {
		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) {
				dr2[i][j] = -R/x[0] - WS(0,i+1) - WS(0,j+1);
				dr2[i][j] += (i != j) ? WS(i+1,j+1) : 0.0;
			}
			dr2[i][i] += (r[i] != 0.0) ? - R/r[i] : 0.0;
		}
		dr2[NR-1][NR-1] +=
		(x[NA-1] != 0.0) ? -R*(1.0/x[NA-1] + 1.0/(1.0-x[NA-1])) : 0.0;
	}
}

-(void)cpmix:(int)mask
		   t:(double)t
		   p:(double)p
		   r:(double [NR])r
	   cpmix:(double *)cpmix // Heat capacity of mixing         BINARY MASK: 001
		  dt:(double *)dt    // d(cp)/d(t)                      BINARY MASK: 010
		  dx:(double [NR])dr // d(cp)/d(x[])d(t)                BINARY MASK: 100
{
	double x[NA];
	int i;

	/* x[0] --> x[NA-1] is an array of mole fractions of liquid components      */
	for (x[0] = 1.0, i=0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }

	if (mask & FIRST) {
		*cpmix = 0.0;
	}

	if(mask & SECOND) {
		*dt = 0.0;
	}

	if(mask & THIRD) {
		for(i=0; i<NR; i++) dr[i] = 0.0;
	}
}

-(void)vmix:(int)mask
		  t:(double)t
		  p:(double)p
		  r:(double [NR])r
	   vmix:(double *)vmix       // Volume of mixing                BINARY MASK: 0000000001
		 dx:(double [NR])dr      // pointer to dx[]) d(v)/d(x[])    BINARY MASK: 0000000010
		dx2:(double [NR][NR])dr2 // pointer to dx2[][]) d(v)/d(x[])2BINARY MASK: 0000000100
		 dt:(double *)dt         // d(v)/d(t)                       BINARY MASK: 0000001000
		 dp:(double *)dp         // d(v)/d(p)                       BINARY MASK: 0000010000
		dt2:(double *)dt2        // d2(v)/d(t)2                     BINARY MASK: 0000100000
	   dtdp:(double *)dtdp       // d2(v)/d(t)d(p)                  BINARY MASK: 0001000000
		dp2:(double *)dp2        // d2(v)/d(p)2                     BINARY MASK: 0010000000
	   dxdt:(double [NR])drdt    // d2(v)/d(x[])d(t)                BINARY MASK: 0100000000
	   dxdp:(double [NR])drdp    // d2(v)/d(x[])d(p)                BINARY MASK: 1000000000
{
	double x[NA];
	int i, j;

	/* x[0] --> x[NA-1] is an array of mole fractions of liquid components      */
	for (x[0] = 1.0, i=0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }

	if (mask & FIRST) {
		for (*vmix = 0.0, i=0; i<NA; i++)
			for (j=i+1; j<NA; j++) *vmix += x[i]*x[j]*WV(i,j);
	}

	if(mask & SECOND) {
		for (i=0; i<NR; i++) {
			dr[i] = x[0]*WV(0,i+1);
			for (j=0; j<NR; j++)
				dr[i] += (i != j) ? r[j]*(WV(i+1,j+1) - WV(0,j+1)) : - r[j]*WV(0,j+1);
		}
	}

	if(mask & THIRD) {
		for (i=0; i<NR; i++) {
			for (j=0; j<NR; j++) {
				dr2[i][j]  = - WV(0,i+1) - WV(0,j+1);
				dr2[i][j] += (i != j) ? WV(i+1,j+1) : 0.0;
			}
		}
	}

	if(mask & FOURTH) {
		*dt = 0.0;
	}

	if(mask & FIFTH) {
		*dp = 0.0;
	}

	if(mask & SIXTH) {
		*dt2 = 0.0;
	}

	if(mask & SEVENTH) {
		*dtdp = 0.0;
	}

	if(mask & EIGHTH) {
		*dp2 = 0.0;
	}

	if(mask & NINTH) {
		for (i=0; i<NR; i++) drdt[i] = 0.0;
	}

	if(mask & TENTH) {
		for (i=0; i<NR; i++) drdp[i] = 0.0;
	}
}

-(NSUInteger)numberOfSolutionSpecies {
	return NA;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
	switch (index) {
		case 0:
			return @"Si4O8";
			break;
		case 1:
			return @"TiO2";
			break;
		case 2:
			return @"Al4O6";
			break;
		case 3:
			return @"Fe2O3";
			break;
		case 4:
			return @"MgCr2O4";
			break;
		case 5:
			return @"Fe2SiO4";
			break;
		case 6:
			return @"MnSi0.5O2";
			break;
		case 7:
			return @"Mg2SiO4";
			break;
		case 8:
			return @"NiSi0.5O2";
			break;
		case 9:
			return @"CoSi0.5O2";
			break;
		case 10:
			return @"Ca2Si2O6";
			break;
		case 11:
			return @"NaSi0.5O1.5";
			break;
		case 12:
			return @"KAlSiO4";
			break;
		case 13:
			return @"Ca3(PO4)2";
			break;
		case 14:
			return @"H2O";
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

-(void)correctActivityCoefficients:(double [NA])gamma forComposition:(double [NA])x {
	for (NSUInteger i=0, j=0; i<NA; i++) if (x[i] != 0.0) {
		if      (gamma[j] > 10000.0) gamma[j] = 10000.0;
		else if (gamma[j] < 0.0001)  gamma[j] = 0.0001;
		j++;
	}
}

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

#include "SolutionPhase.h"

@end
