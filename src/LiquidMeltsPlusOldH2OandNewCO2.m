//
//  LiquidMeltsPlusOldH2OandNewCO2.m
//  PhaseObjC
//
//  Created by Mark Ghiorso on 7/6/17.
//  Copyright © 2017 Mark Ghiorso. All rights reserved.
//

#import "LiquidMeltsPlusOldH2OandNewCO2.h"
#import "BermanProperties.h"
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "DoubleTensor.h"
#import "LiquidMeltsGenericEM.h"
#import "LiquidMeltsSiO2.h"
#import "LiquidMeltsH2O.h"
#import "LiquidMeltsCO2.h"

@implementation LiquidMeltsPlusOldH2OandNewCO2

@synthesize endmembers;

static double referenceValuesOfModelParameters[] = {
    26266.7,  //   0 W(TiO2      ,SiO2      )
    -39120.0,  //   1 W(Al2O3     ,SiO2      )
    8110.3,  //   2 W(Fe2O3     ,SiO2      )
    27886.3,  //   3 W(MgCr2O4   ,SiO2      )
    23660.9,  //   4 W(Fe2SiO4   ,SiO2      )
    18393.9,  //   5 W(MnSi0.5O2 ,SiO2      )
    3421.0,  //   6 W(Mg2SiO4   ,SiO2      )
    25197.4,  //   7 W(NiSi0.5O2 ,SiO2      )
    14802.8,  //   8 W(CoSi0.5O2 ,SiO2      )
    -863.7,  //   9 W(CaSiO3    ,SiO2      )
    -99039.0,  //  10 W(Na2SiO3   ,SiO2      )
    -33921.7,  //  11 W(KAlSiO4   ,SiO2      )
    61891.6,  //  12 W(Ca3(PO4)2 ,SiO2      )
    30967.3,  //  13 W(H2O       ,SiO2      )
    0.0,  //  14 W(CO2       ,SiO2      )
    6.3281206517489e+01*1000.0,  //  15 W(CaCO3     ,SiO2      )

    -29449.8,  //  16 W(Al2O3     ,TiO2      )
    -84756.9,  //  17 W(Fe2O3     ,TiO2      )
    -72303.4,  //  18 W(MgCr2O4   ,TiO2      )
    5209.1,  //  19 W(Fe2SiO4   ,TiO2      )
    -16123.5,  //  20 W(MnSi0.5O2 ,TiO2      )
    -4178.3,  //  21 W(Mg2SiO4   ,TiO2      )
    3614.8,  //  22 W(NiSi0.5O2 ,TiO2      )
    -1640.0,  //  23 W(CoSi0.5O2 ,TiO2      )
    -35372.5,  //  24 W(CaSiO3    ,TiO2      )
    -15415.6,  //  25 W(Na2SiO3   ,TiO2      )
    -48094.6,  //  26 W(KAlSiO4   ,TiO2      )
    25938.8,  //  27 W(Ca3(PO4)2 ,TiO2      )
    81879.1,  //  28 W(H2O       ,TiO2      )
    -1.9265668647972e+01*1000.0,  //  29 W(CO2       ,TiO2      )
    -7.9202707514741e+01*1000.0,  //  30 W(CaCO3     ,TiO2      )

    -17089.4,  //  31 W(Fe2O3     ,Al2O3     )
    -31770.3,  //  32 W(MgCr2O4   ,Al2O3     )
    -30509.0,  //  33 W(Fe2SiO4   ,Al2O3     )
    -53874.9,  //  34 W(MnSi0.5O2 ,Al2O3     )
    -32880.3,  //  35 W(Mg2SiO4   ,Al2O3     )
    2985.2,  //  36 W(NiSi0.5O2 ,Al2O3     )
    -2677.4,  //  37 W(CoSi0.5O2 ,Al2O3     )
    -57917.9,  //  38 7(CaSiO3    ,Al2O3     )
    -130785.0,  //  39 W(Na2SiO3   ,Al2O3     )
    -25859.2,  //  40 W(KAlSiO4   ,Al2O3     )
    52220.8,  //  41 W(Ca3(PO4)2 ,Al2O3     )
    -16098.1,  //  42 W(H2O       ,Al2O3     )
    0.0,  //  43 W(CO2       ,Al2O3     )
    4.6716121152125e+01*1000.0,  //  44 W(CaCO3     ,Al2O3     )

    21605.9,  //  45 W(MgCr2O4   ,Fe2O3     )
    -179064.9,  //  46 W(Fe2SiO4   ,Fe2O3     )
    3907.9,  //  47 W(MnSi0.5O2 ,Fe2O3     )
    -71518.6,  //  48 W(Mg2SiO4   ,Fe2O3     )
    408.7,  //  49 W(NiSi0.5O2 ,Fe2O3     )
    -223.7,  //  50 W(CoSi0.5O2 ,Fe2O3     )
    12076.6,  //  51 W(CaSiO3    ,Fe2O3     )
    -149662.2,  //  52 W(Na2SiO3   ,Fe2O3     )
    57555.9,  //  53 W(KAlSiO4   ,Fe2O3     )
    -4213.9,  //  54 W(Ca3(PO4)2 ,Fe2O3     )
    31405.5,  //  55 W(H2O       ,Fe2O3     )
    -3.1868817914962e+00*1000.0,  //  56 W(CO2       ,Fe2O3     )
    6.5508729974818e+01*1000.0,  //  57 W(CaCO3     ,Fe2O3     )

    -82971.8,  //  58 W(Fe2SiO4   ,MgCr2O4   )
    182.4,  //  59 W(MnSi0.5O2 ,MgCr2O4   )
    46049.2,  //  60 W(Mg2SiO4   ,MgCr2O4   )
    -266.0,  //  61 W(NiSi0.5O2 ,MgCr2O4   )
    -384.0,  //  62 W(CoSi0.5O2 ,MgCr2O4   )
    30704.7,  //  63 W(CaSiO3    ,MgCr2O4   )
    113646.0,  //  64 W(Na2SiO3   ,MgCr2O4   )
    75709.1,  //  65 W(KAlSiO4   ,MgCr2O4   )
    5341.8,  //  66 W(Ca3(PO4)2 ,MgCr2O4   )
    0.0,  //  67 W(H2O       ,MgCr2O4   )
    0.0,  //  68 W(CO2       ,MgCr2O4   )
    0.0,  //  69 W(CaCO3     ,MgCr2O4   )

    -6823.9,  //  70 W(MnSi0.5O2 ,Fe2SiO4   )
    -37256.7,  //  71 W(Mg2SiO4   ,Fe2SiO4   )
    -17019.8,  //  72 W(NiSi0.5O2 ,Fe2SiO4   )
    -11746.3,  //  73 W(CoSi0.5O2 ,Fe2SiO4   )
    -12970.8,  //  74 W(CaSiO3    ,Fe2SiO4   )
    -90533.8,  //  75 W(Na2SiO3   ,Fe2SiO4   )
    23649.4,  //  76 W(KAlSiO4   ,Fe2SiO4   )
    87410.3,  //  77 W(Ca3(PO4)2 ,Fe2SiO4   )
    28873.6,  //  78 W(H2O       ,Fe2SiO4   )
    -3.2464545980393e+01*1000.0,  //  79 W(CO2       ,Fe2SiO4   )
    -7.2996815311226e+01*1000.0,  //  80 W(CaCO3     ,Fe2SiO4   )

    -13040.1,  //  81 W(Mg2SiO4   ,MnSi0.5O2 )
    785.8,  //  82 W(NiSi0.5O2 ,MnSi0.5O2 )
    -50.6,  //  83 W(CoSi0.5O2 ,MnSi0.5O2 )
    2934.6,  //  84 W(CaSiO3    ,MnSi0.5O2 )
    -15780.8,  //  85 W(Na2SiO3   ,MnSi0.5O2 )
    23727.4,  //  86 W(KAlSiO4   ,MnSi0.5O2 )
    0.0,  //  87 W(Ca3(PO4)2 ,MnSi0.5O2 )
    0.0,  //  88 W(H2O       ,MnSi0.5O2 )
    0.0,  //  89 W(CO2       ,MnSi0.5O2 )
    0.0,  //  90 W(CaCO3     ,MnSi0.5O2 )

    -21175.5,  //  91 W(NiSi0.5O2 ,Mg2SiO4   )
    -14994.9,  //  92 W(CoSi0.5O2 ,Mg2SiO4   )
    -31731.9,  //  93 W(CaSiO3    ,Mg2SiO4   )
    -41876.9,  //  94 W(Na2SiO3   ,Mg2SiO4   )
    22323.1,  //  95 W(KAlSiO4   ,Mg2SiO4   )
    -23208.8,  //  96 W(Ca3(PO4)2 ,Mg2SiO4   )
    35633.7,  //  97 W(H2O       ,Mg2SiO4   )
    -4.0853595138597e+01*1000.0,  //  98 W(CO2       ,Mg2SiO4   )
    -2.4872562549739e+01*1000.0,  //  99 W(CaCO3     ,Mg2SiO4   )

    258.9,  // 100 W(CoSi0.5O2 ,NiSi0.5O2 )
    7027.5,  // 101 W(CaSiO3    ,NiSi0.5O2 )
    -3647.8,  // 102 W(Na2SiO3   ,NiSi0.5O2 )
    4261.4,  // 103 W(KAlSiO4   ,NiSi0.5O2 )
    0.0,  // 104 W(Ca3(PO4)2 ,NiSi0.5O2 )
    0.0,  // 105 W(H2O       ,NiSi0.5O2 )
    0.0,  // 106 W(CO2       ,NiSi0.5O2 )
    0.0,  // 107 W(CaCO3     ,NiSi0.5O2 )

    -26685.7,  // 108 W(CaSiO3    ,CoSi0.5O2 )
    531.2,  // 109 W(Na2SiO3   ,CoSi0.5O2 )
    265.7,  // 110 W(KAlSiO4   ,CoSi0.5O2 )
    0.0,  // 111 W(Ca3(PO4)2 ,CoSi0.5O2 )
    0.0,  // 112 W(H2O       ,CoSi0.5O2 )
    0.0,  // 113 W(CO2       ,CoSi0.5O2 )
    0.0,  // 114 W(CaCO3     ,CoSi0.5O2 )

    -13247.1,  // 115 W(Na2SiO3   ,CaSiO3    )
    17111.1,  // 116 W(KAlSiO4   ,CaSiO3    )
    37070.3,  // 117 W(Ca3(PO4)2 ,CaSiO3    )
    20374.6,  // 118 W(H2O       ,CaSiO3    )
    3.0012481472120e+01*1000.0,  // 119 W(CO2       ,CaSiO3    )
    3.7534127781795e+01*1000.0,  // 120 W(CaCO3     ,CaSiO3    )

    6522.8,  // 121 W(KAlSiO4   ,Na2SiO3   )
    15571.9,  // 122 W(Ca3(PO4)2 ,Na2SiO3   )
    -96937.6,  // 123 W(H2O       ,Na2SiO3   )
    0.0,  // 124 W(CO2       ,Na2SiO3   )
    -3.1101129932991e+02*1000.0,  // 126 W(CaCO3     ,Na2SiO3   )

    17100.6,  // 127 W(Ca3(PO4)2 ,KAlSiO4   )
    10374.2,  // 128 W(H2O       ,KAlSiO4   )
    0.0,  // 129 W(CO2       ,KAlSiO4   )
    -2.7864967599360e+01*1000.0,  // 130 W(CaCO3     ,KAlSiO4   )

    43451.3,  // 131 W(H2O       ,Ca3(PO4)2 )
    -3.4728165916175e+00*1000.0,  // 132 W(CO2       ,Ca3(PO4)2 )
    2.0119124588637e+00*1000.0,  // 133 W(CaCO3     ,Ca3(PO4)2 )

    2.3255499914663e+01*1000.0,  // 134 W(CO2       ,H2O       )
    7.8729683594722e+00*1000.0,  // 135 W(CaCO3     ,H2O       )

    0.0,  // 180 W(CaCO3     ,CO2       )
};

#ifdef DEBUG
#undef DEBUG
#endif

#define NR      15      // Number of independent mole fraction variables
#define NS       1      // Three ordering parameters
#define NA      16      // Number of liquid components
#define NE      (NA+NS) // Total number of species
#define NATOMS  5.0     // Average number of atoms in the formula unit (guess)
#define mIdxH2O NA-2    // This is the index number of H2O in the endmember array
#define rIdxH2O NR-2    // This is the index number of H2O in the independent variable array
#define mIdxCO2 NA-1    // This is the index number of CO2 in the endmember array
#define rIdxCO2 NR-1    // This is the index number of CO2 in the independent variable array
#define rIdxCaSiO3    9
#define SMX     SHRT_MAX

#pragma mark -
#pragma mark class methods

+(void)initialize {
    if (self == [LiquidMeltsPlusOldH2OandNewCO2 class]) {
        BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
        if (debug) NSLog(@"Initialize(LiquidMelts) - entry ...");
    }
}

#pragma mark -
#pragma mark instance methods

-(id)init {
    if ((self = [super init])) {
        BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.VERBOSE"];
        if (debug) NSLog(@"init(LiquidMeltsPlusCO2) ... entry ...");

        NSMutableArray *mutableEndmembers = [NSMutableArray arrayWithCapacity:NA];

        LiquidMeltsSiO2 *sio2 = [[LiquidMeltsSiO2 alloc] init];
        [mutableEndmembers addObject:sio2];
        if (debug) NSLog(@"... allocated sio2 ...");

        LiquidMeltsGenericEM *tio2 = [[LiquidMeltsGenericEM alloc] initWithH:-944750.0
                                                                           S:50.460
                                                                          k0:77.84
                                                                          k1:0.0
                                                                          k2:-33.678e5
                                                                          k3:40.294e7
                                                                          l1:0.0
                                                                          l2:0.0
                                                                          Tt:0.0
                                                                      deltaH:0.0
                                                                        vLiq:2.316
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

        LiquidMeltsGenericEM *al2o3 = [[LiquidMeltsGenericEM alloc] initWithH:-1675700.0
                                                                            S:50.82
                                                                           k0:155.02
                                                                           k1:-8.284E2
                                                                           k2:-38.614E5
                                                                           k3:40.908E7
                                                                           l1:0.0
                                                                           l2:0.0
                                                                           Tt:0.0
                                                                       deltaH:0.0
                                                                         vLiq:3.711
                                                                      dvdtLiq:2.62e-4
                                                                      dvdpLiq:-2.26e-5
                                                                   d2vdtdpLiq:2.7e-8
                                                                    d2vdp2Liq:4.0e-10
                                                                      tFusion:2319.65
                                                                      sFusion:48.61
                                                                        cpLiq:170.3];
        [al2o3 setPhaseName:@"Al2O3"];
        [al2o3 setPhaseFormula:@"Al2O3"];
        [mutableEndmembers addObject:al2o3];
        if (debug) NSLog(@"... allocated al2o3 ...");

        LiquidMeltsGenericEM *fe2o3 = [[LiquidMeltsGenericEM alloc] initWithH:-822000.00
                                                                            S:87.40
                                                                           k0:146.86
                                                                           k1:0.0
                                                                           k2:-55.768E5
                                                                           k3:52.563E7
                                                                           l1:-7.403E-2
                                                                           l2:27.921E-5
                                                                           Tt:955.0
                                                                       deltaH:1287.0
                                                                         vLiq:4.213
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

        LiquidMeltsGenericEM *mgcr2o4 = [[LiquidMeltsGenericEM alloc] initWithH:-1783640.0
                                                                              S:106.02
                                                                             k0:201.981
                                                                             k1:-5.519E2
                                                                             k2:-57.844E5
                                                                             k3:57.729E7
                                                                             l1:0.0
                                                                             l2:0.0
                                                                             Tt:0.0
                                                                         deltaH:0.0
                                                                           vLiq:5.358
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

        LiquidMeltsGenericEM *fe2sio4 = [[LiquidMeltsGenericEM alloc] initWithH:-1479360.0
                                                                              S:150.930
                                                                             k0:248.93
                                                                             k1:-19.239E2
                                                                             k2:0.0
                                                                             k3:-13.910E7
                                                                             l1:0.0
                                                                             l2:0.0
                                                                             Tt:0.0
                                                                         deltaH:0.0
                                                                           vLiq:5.420
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

        LiquidMeltsGenericEM *mn2sio4 = [[LiquidMeltsGenericEM alloc] initWithH:-1732000.0/2.0
                                                                              S:155.9/2.0
                                                                             k0:219.89/2.0
                                                                             k1:-12.710E2/2.0
                                                                             k2:-20.496E5/2.0
                                                                             k3:17.652E7/2.0
                                                                             l1:0.0
                                                                             l2:0.0
                                                                             Tt:0.0
                                                                         deltaH:0.0
                                                                           vLiq:2.84
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

        LiquidMeltsGenericEM *mg2sio4 = [[LiquidMeltsGenericEM alloc] initWithH:-2174420.0
                                                                              S:94.010
                                                                             k0:238.64
                                                                             k1:-20.013E2
                                                                             k2:0.0
                                                                             k3:-11.624E7
                                                                             l1:0.0
                                                                             l2:0.0
                                                                             Tt:0.0
                                                                         deltaH:0.0
                                                                           vLiq:4.980
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

        LiquidMeltsGenericEM *ni2sio4 = [[LiquidMeltsGenericEM alloc] initWithH:-1395300.0/2.0
                                                                              S:128.1/2.0
                                                                             k0:214.997/2.0
                                                                             k1:-10.3075E2/2.0
                                                                             k2:-49.4453E5/2.0
                                                                             k3:62.375E7/2.0
                                                                             l1:0.0
                                                                             l2:0.0
                                                                             Tt:0.0
                                                                         deltaH:0.0
                                                                           vLiq:2.48
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

        LiquidMeltsGenericEM *co2sio4 = [[LiquidMeltsGenericEM alloc] initWithH:-1414100.0/2.0
                                                                              S:142.6/2.0
                                                                             k0:201.048/2.0
                                                                             k1:-0.369E2/2.0
                                                                             k2:-71.81E5/2.0
                                                                             k3:90.05E7/2.0
                                                                             l1:0.0
                                                                             l2:0.0
                                                                             Tt:0.0
                                                                         deltaH:0.0
                                                                           vLiq:2.30
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

        LiquidMeltsGenericEM *casio3 = [[LiquidMeltsGenericEM alloc] initWithH:-1627427.0
                                                                             S:85.279
                                                                            k0:141.16
                                                                            k1:-4.172e2
                                                                            k2:-58.576e5
                                                                            k3:94.074e7
                                                                            l1:0.0
                                                                            l2:0.0
                                                                            Tt:0.0
                                                                        deltaH:0.0
                                                                          vLiq:4.347
                                                                       dvdtLiq:2.92e-4
                                                                       dvdpLiq:-1.55e-5
                                                                    d2vdtdpLiq:-1.6e-8
                                                                     d2vdp2Liq:3.89e-10
                                                                       tFusion:1817.0
                                                                       sFusion:31.5
                                                                         cpLiq:172.4];
        [casio3 setPhaseName:@"CaSiO3"];
        [casio3 setPhaseFormula:@"CaSiO3"];
        [mutableEndmembers addObject:casio3];
        if (debug) NSLog(@"... allocated casio3 ...");

        LiquidMeltsGenericEM *na2sio3 = [[LiquidMeltsGenericEM alloc] initWithH:-373190.0*4.184
                                                                              S:27.21*4.184
                                                                             k0:234.77
                                                                             k1:-22.189E2
                                                                             k2:0.0
                                                                             k3:13.530E7
                                                                             l1:0.0
                                                                             l2:0.0
                                                                             Tt:0.0
                                                                         deltaH:0.0
                                                                           vLiq:5.568
                                                                        dvdtLiq:7.41e-4
                                                                        dvdpLiq:-4.29e-5
                                                                     d2vdtdpLiq:-5.3e-8
                                                                      d2vdp2Liq:8.4e-10
                                                                        tFusion:1361.0
                                                                        sFusion:38.34
                                                                          cpLiq:180.2];
        [na2sio3 setPhaseName:@"Na2SiO3"];
        [na2sio3 setPhaseFormula:@"Na2SiO3"];
        [mutableEndmembers addObject:na2sio3];
        if (debug) NSLog(@"... allocated na2sio3 ...");

        LiquidMeltsGenericEM *kalsio4 = [[LiquidMeltsGenericEM alloc] initWithH:-2111813.55
                                                                              S:133.9653
                                                                             k0:186.0
                                                                             k1:0.0
                                                                             k2:-131.067E5
                                                                             k3:213.893E7
                                                                             l1:-7.096454E-2
                                                                             l2:21.682E-5
                                                                             Tt:800.15
                                                                         deltaH:1154.0
                                                                           vLiq:6.8375
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

        LiquidMeltsGenericEM *ca3po4 = [[LiquidMeltsGenericEM alloc] initWithH:-4097169.0
                                                                             S:235.978
                                                                            k0:402.997
                                                                            k1:-28.0835E2
                                                                            k2:0.0
                                                                            k3:-32.6230E7
                                                                            l1:2.5427E-2
                                                                            l2:19.255E-5
                                                                            Tt:1373.0
                                                                        deltaH:14059.0
                                                                          vLiq:10.7382
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

        LiquidMeltsH2O *h2o = [[LiquidMeltsH2O alloc] init];
        [h2o setPhaseName:@"H2O"];
        [h2o setPhaseFormula:@"H2O"];
        [mutableEndmembers addObject:h2o];
        if (debug) NSLog(@"... allocated h2o ...");

        LiquidMeltsCO2 *co2 = [[LiquidMeltsCO2 alloc] init];
        [mutableEndmembers addObject:co2];
        if (debug) NSLog(@"... allocated co2 and exiting.");

        endmembers = [NSArray arrayWithArray:mutableEndmembers];

        NSUInteger lengthOfModelParametersVector = NE*(NE-1)/2;
        _modelParametersWrapper = [[DoubleVector alloc] initWithSize:lengthOfModelParametersVector];
        double *modelParameters = [_modelParametersWrapper pointerToDouble];
        for (NSUInteger i=0; i<lengthOfModelParametersVector; i++) modelParameters[i] = referenceValuesOfModelParameters[i];
        modelParametersHaveBeenAlteredSinceLastOrderingDetermination = NO;

        computeMixingQuantities = NO;
        tOld = -9999.0;
        pOld = -9999.0;
        for (NSUInteger i=0; i<NR; i++) rOld[i] = -9999.0;
        for (NSUInteger i=0; i<NS; i++) sOld[i] = 2.0;
        [self setPhaseName:@"Liquid"];
        deltaHsp[0] = -1.7574497522747e+01*1000.0;
        deltaSsp[0] = 0.0;
        deltaVsp[0] = -1.9034060173857e+00; // CaCO3
        if (debug) NSLog(@"... exiting.");
    }
    return self;
}

#pragma mark -
#pragma mark NSSecureCoding protocol methods

static NSString *kmodelParametersWrapper = @"modelParametersWrapper";
static NSString *kendmembers = @"endmembers";
static NSString *kdeltaHsp_0 = @"deltaHsp_0";
static NSString *kdeltaSsp_0 = @"deltaSsp_0";
static NSString *kdeltaVsp_0 = @"deltaVsp_0";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
#ifdef __APPLE__
        _modelParametersWrapper = (DoubleVector *) [aDecoder decodeObjectOfClass:[DoubleVector class] forKey:kmodelParametersWrapper];
        endmembers = (NSArray *) [aDecoder decodeObjectOfClass:[NSArray class] forKey:kendmembers];
#else
        _modelParametersWrapper = (DoubleVector *) [aDecoder decodeObjectForKey:kmodelParametersWrapper];
        endmembers = (NSArray *) [aDecoder decodeObjectForKey:kendmembers];
#endif
        deltaHsp[0] = [aDecoder decodeDoubleForKey:kdeltaHsp_0];
        deltaSsp[0] = [aDecoder decodeDoubleForKey:kdeltaSsp_0];
        deltaVsp[0] = [aDecoder decodeDoubleForKey:kdeltaVsp_0];
        computeMixingQuantities = NO;
        tOld = -9999.0;
        pOld = -9999.0;
        for (NSUInteger i=0; i<NR; i++) rOld[i] = -9999.0;
        for (NSUInteger i=0; i<NS; i++) sOld[i] = 2.0;
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeObject:_modelParametersWrapper forKey:kmodelParametersWrapper];
        [aCoder encodeObject:endmembers forKey:kendmembers];
        [aCoder encodeDouble:deltaHsp[0] forKey:kdeltaHsp_0];
        [aCoder encodeDouble:deltaSsp[0] forKey:kdeltaSsp_0];
        [aCoder encodeDouble:deltaVsp[0] forKey:kdeltaVsp_0];
    } // else [NSException raise:NSInvalidArchiveOperationException format:@"Class %@ only supports NSKeyedArchiver coders.", [self className]];
}

#pragma mark -
#pragma mark original C code from MELTS

/**
 "C" code from MELTS
 */

/*
 * Array to convert W(i,j) indexes to entries in the array of structures
 * modelParameters[k], where k is:
 */

static const short int Index[NE][NE] = {
    { SMX,   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15 },
    {   0, SMX,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30 },
    {   1,  16, SMX,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44 },
    {   2,  17,  31, SMX,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57 },
    {   3,  18,  32,  45, SMX,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69 },
    {   4,  19,  33,  46,  58, SMX,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80 },
    {   5,  20,  34,  47,  59,  70, SMX,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90 },
    {   6,  21,  35,  48,  60,  71,  81, SMX,  91,  92,  93,  94,  95,  96,  97,  98,  99 },
    {   7,  22,  36,  49,  61,  72,  82,  91, SMX, 100, 101, 102, 103, 104, 105, 106, 107 },
    {   8,  23,  37,  50,  62,  73,  83,  92, 100, SMX, 108, 109, 110, 111, 112, 113, 114 },
    {   9,  24,  38,  51,  63,  74,  84,  93, 101, 108, SMX, 115, 116, 117, 118, 119, 120 },
    {  10,  25,  39,  52,  64,  75,  85,  94, 102, 109, 115, SMX, 121, 122, 123, 124, 125 },
    {  11,  26,  40,  53,  65,  76,  86,  95, 103, 110, 116, 121, SMX, 126, 127, 128, 129 },
    {  12,  27,  41,  54,  66,  77,  87,  96, 104, 111, 117, 122, 126, SMX, 130, 131, 132 },
    {  13,  28,  42,  55,  67,  78,  88,  97, 105, 112, 118, 123, 127, 130, SMX, 133, 134 },
    {  14,  29,  43,  56,  68,  79,  89,  98, 106, 113, 119, 124, 128, 131, 133, SMX, 135 },
    {  15,  30,  44,  57,  69,  80,  90,  99, 107, 114, 120, 125, 129, 132, 134, 135, SMX }
};

#undef SMX
#define WH(i,j) modelParameters[Index[i][j]]

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

#define MAX_ITER 200    /* Maximum number of iterations allowed in order */

#define SQUARE(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

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

#define fillXSPECIES \
double xSpecies[NE]; \
xSpecies[ 0] = 1.0 - rTotal + s[0]; \
xSpecies[ 1] = r[ 0];               \
xSpecies[ 2] = r[ 1];               \
xSpecies[ 3] = r[ 2];               \
xSpecies[ 4] = r[ 3];               \
xSpecies[ 5] = r[ 4];               \
xSpecies[ 6] = r[ 5];               \
xSpecies[ 7] = r[ 6];               \
xSpecies[ 8] = r[ 7];               \
xSpecies[ 9] = r[ 8];               \
xSpecies[10] = r[ 9] - s[0];        \
xSpecies[11] = r[10];               \
xSpecies[12] = r[11];               \
xSpecies[13] = r[12];               \
xSpecies[14] = r[13];               \
xSpecies[15] = r[14] - s[0];        \
xSpecies[16] = s[0];                                           \
for (i=0; i<NE; i++) if (xSpecies[i] <= 0.0) xSpecies[i] = DBL_EPSILON;

static const double dxSpeciesds[NE][NS] = {
    {  1.0},
    {  0.0},
    {  0.0},
    {  0.0},
    {  0.0},
    {  0.0},
    {  0.0},
    {  0.0},
    {  0.0},
    {  0.0},
    { -1.0},
    {  0.0},
    {  0.0},
    {  0.0},
    {  0.0},
    { -1.0},
    {  1.0},
};

static const double dxSpeciesdr[NE][NR] = {
    { -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0},
    {  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    {  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    {  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    {  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    {  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    {  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0},
    {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0},
    {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0},
    {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0},
    {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0},
    {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
};

#define fillG \
g = 0.0; \
for (i=0; i<NS; i++)                          g += (deltaHsp[i] - t*deltaSsp[i] + (p-1.0)*deltaVsp[i])*xSpecies[NA+i]; \
for (i=0; i<NE-1; i++) for (j=i+1; j<NE; j++) g += (WH(i, j))*xSpecies[i]*xSpecies[j]; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)   g += R*t*xSpecies[i]*log(xSpecies[i]); \
if (xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0) g += R*t*(xSpecies[mIdxH2O]*log(xSpecies[mIdxH2O]) + (1.0-xSpecies[mIdxH2O])*log(1.0-xSpecies[mIdxH2O]));

#define fillS \
entropy = 0.0; \
for (i=0; i<NS; i++)                          entropy += deltaSsp[i]*xSpecies[NA+i]; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)   entropy += -R*xSpecies[i]*log(xSpecies[i]); \
if (xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0) entropy += -R*(xSpecies[mIdxH2O]*log(xSpecies[mIdxH2O]) + (1.0-xSpecies[mIdxH2O])*log(1.0-xSpecies[mIdxH2O]));

#define fillV \
v = 0.0; \
for (i=0; i<NS; i++)                          v += deltaVsp[i]*xSpecies[NA+i];

#define fillDGDR \
for (k=0; k<NR; k++) { \
dgdr[k] = 0.0; \
for (i=0; i<NS; i++)                          dgdr[k] += (deltaHsp[i] - t*deltaSsp[i] + (p-1.0)*deltaVsp[i])*dxSpeciesdr[NA+i][k]; \
for (i=0; i<NE-1; i++) for (j=i+1; j<NE; j++) dgdr[k] += (WH(i, j))*(dxSpeciesdr[i][k]*xSpecies[j] + xSpecies[i]*dxSpeciesdr[j][k]); \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][k] != 0.0)                   dgdr[k] += R*t*dxSpeciesdr[i][k]*(1.0 + log(xSpecies[i])); \
if (xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0 && dxSpeciesdr[mIdxH2O][k] != 0.0) dgdr[k] += R*t*dxSpeciesdr[mIdxH2O][k]*(log(xSpecies[mIdxH2O]) - log(1.0-xSpecies[mIdxH2O])); \
}

#define fillDGDS \
for (i=0; i<NS; i++) { \
dgds[i] = 0.0; \
for (j=0; j<NS; j++)                          dgds[i] += (deltaHsp[j] - t*deltaSsp[j] + (p-1.0)*deltaVsp[j])*dxSpeciesds[NA+j][i]; \
for (j=0; j<NE-1; j++) for (k=j+1; k<NE; k++) dgds[i] += (WH(j, k))*(xSpecies[j]*dxSpeciesds[k][i] + dxSpeciesds[j][i]*xSpecies[k]); \
for (j=0; j<NE; j++) if (xSpecies[j] > 0.0)   dgds[i] += R*t*dxSpeciesds[j][i]*(1.0 + log(xSpecies[j])); \
if (xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0) dgds[i] += R*t*dxSpeciesds[mIdxH2O][i]*(log(xSpecies[mIdxH2O]) - log(1.0-xSpecies[mIdxH2O])); \
}

#define fillD2GDR2 \
for (k=0; k<NR; k++) { \
for (l=k; l<NR; l++) { \
d2gdr2[k][l]  = 0.0; \
for (i=0; i<NE-1; i++) for (j=i+1; j<NE; j++) d2gdr2[k][l] += (WH(i, j))*(dxSpeciesdr[i][k]*dxSpeciesdr[j][l] + dxSpeciesdr[i][l]*dxSpeciesdr[j][k]); \
} \
} \
for (j=0; j<NR; j++) { \
for (k=j; k<NR; k++) { \
double config = 0.0; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0 && dxSpeciesdr[i][k] != 0.0) \
config += dxSpeciesdr[i][j]*dxSpeciesdr[i][k]/xSpecies[i]; \
if (xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0 && dxSpeciesdr[mIdxH2O][j] != 0.0 && dxSpeciesdr[mIdxH2O][k] != 0.0) \
config += dxSpeciesdr[mIdxH2O][j]*dxSpeciesdr[mIdxH2O][k]*(1.0/xSpecies[mIdxH2O] + 1.0/(1.0-xSpecies[mIdxH2O])); \
d2gdr2[j][k] += R*t*config; \
d2gdr2[k][j]  = d2gdr2[j][k]; \
} \
}

#define fillD2GDRDS \
for (k=0; k<NR; k++) { \
for (l=0; l<NS; l++) { \
d2gdrds[k][l]  = 0.0; \
for (i=0; i<NE-1; i++) for (j=i+1; j<NE; j++) d2gdrds[k][l] += (WH(i, j))*(dxSpeciesdr[i][k]*dxSpeciesds[j][l] + dxSpeciesds[i][l]*dxSpeciesdr[j][k]); \
} \
} \
for (j=0; j<NR; j++) { \
for (k=0; k<NS; k++) { \
double config = 0.0; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)             config += (dxSpeciesdr[i][j]*dxSpeciesds[i][k]/xSpecies[i]); \
if (xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0) config += dxSpeciesds[mIdxH2O][k]*dxSpeciesdr[mIdxH2O][j]*(1.0/xSpecies[mIdxH2O] + 1.0/(1.0-xSpecies[mIdxH2O])); \
d2gdrds[j][k] += R*t*config; \
} \
}

#define fillD2GDRDT \
for (j=0; j<NR; j++) { \
double config = 0.0; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0)                   config += dxSpeciesdr[i][j]*(1.0 + log(xSpecies[i])); \
if (xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0 && dxSpeciesdr[mIdxH2O][j] != 0.0) config += dxSpeciesdr[mIdxH2O][j]*(log(xSpecies[mIdxH2O]) - log(1.0-xSpecies[mIdxH2O])); \
d2gdrdt[j] = R*config; \
}

#define fillD2GDRDP \
for (k=0; k<NR; k++) d2gdrdp[k] = 0.0;

#define fillD2GDS2 \
for (i=0; i<NS; i++) { \
for (l=0; l<NS; l++) { \
d2gds2[i][l]  = 0.0; \
for (j=0; j<NE-1; j++) for (k=j+1; k<NE; k++) d2gds2[i][l] += (WH(j, k))*(dxSpeciesds[j][l]*dxSpeciesds[k][i] + dxSpeciesds[j][i]*dxSpeciesds[k][l]); \
} \
} \
for (j=0; j<NS; j++) { \
for (k=j; k<NS; k++) { \
double config = 0.0; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)             config += dxSpeciesds[i][j]*dxSpeciesds[i][k]/xSpecies[i]; \
if (xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0) config += dxSpeciesds[mIdxH2O][j]*dxSpeciesds[mIdxH2O][k]*(1.0/xSpecies[mIdxH2O] + 1.0/(1.0-xSpecies[mIdxH2O])); \
d2gds2[j][k] += R*t*config; \
d2gds2[k][j]  = d2gds2[j][k]; \
} \
}

#define fillD2GDSDT \
for (i=0; i<NS; i++) { \
d2gdsdt[i] = 0.0; \
for (j=0; j<NS; j++)  d2gdsdt[i] += (-deltaSsp[j])*dxSpeciesds[NA+j][i]; \
} \
for (j=0; j<NS; j++) { \
double config = 0.0; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0) \
config += dxSpeciesds[i][j]*(1.0 + log(xSpecies[i])); \
if (xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0) config += dxSpeciesds[mIdxH2O][j]*(log(xSpecies[mIdxH2O]) - log(1.0-xSpecies[mIdxH2O])); \
d2gdsdt[j] += R*config; \
}

#define fillD2GDSDP \
for (i=0; i<NS; i++) { \
d2gdsdp[i] = 0.0; \
for (j=0; j<NS; j++) d2gdsdp[i] += (deltaVsp[j])*dxSpeciesds[NA+j][i]; \
}

#define D2GDT2 0.0

#define D2GDTDP 0.0

#define D2GDP2 0.0

#define fillD3GDR2DS \
for (j=0; j<NR; j++) { \
for (k=j; k<NR; k++) { \
for (l=0; l<NS; l++) { \
double config = 0.0; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0) \
config += - dxSpeciesdr[i][j]*dxSpeciesdr[i][k]*dxSpeciesds[i][l]/(xSpecies[i]*xSpecies[i]); \
if (mIdxH2O != -1 && xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0) \
config += (- dxSpeciesdr[mIdxH2O][j]*dxSpeciesdr[mIdxH2O][k]*dxSpeciesds[mIdxH2O][l] \
*(1.0/(xSpecies[mIdxH2O]*xSpecies[mIdxH2O]) - 1.0/((1.0-xSpecies[mIdxH2O])*(1.0-xSpecies[mIdxH2O]))) ); \
d3gdr2ds[j][k][l] = R*t*config; \
d3gdr2ds[k][j][l] = d3gdr2ds[j][k][l]; \
} \
} \
}

#define fillD3GDR2DT \
for (j=0; j<NR; j++) { \
for (k=j; k<NR; k++) { \
double config = 0.0; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0 && dxSpeciesdr[i][k] != 0.0) \
config += dxSpeciesdr[i][j]*dxSpeciesdr[i][k]/xSpecies[i]; \
if (mIdxH2O != -1 && xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0 && dxSpeciesdr[mIdxH2O][j] != 0.0 && dxSpeciesdr[mIdxH2O][k] != 0.0) \
config += dxSpeciesdr[mIdxH2O][j]*dxSpeciesdr[mIdxH2O][k]*(1.0/xSpecies[mIdxH2O] + 1.0/(1.0-xSpecies[mIdxH2O])); \
d3gdr2dt[j][k] = R*config; \
d3gdr2dt[k][j] = d3gdr2dt[j][k]; \
} \
}


#define fillD3GDR2DP \
for (k=0; k<NR; k++) for (l=0; l<NR; l++) d3gdr2dp[k][l] = 0.0;

#define fillD3GDRDS2 \
for (j=0; j<NR; j++) { \
for (k=0; k<NS; k++) { \
for (l=k; l<NS; l++) { \
double config = 0.0; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0) \
config += - dxSpeciesdr[i][j]*dxSpeciesds[i][k]*dxSpeciesds[i][l]/(xSpecies[i]*xSpecies[i]); \
if (mIdxH2O != -1 && xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0) \
config += ( - dxSpeciesds[mIdxH2O][k]*dxSpeciesdr[mIdxH2O][j]*dxSpeciesds[mIdxH2O][l] \
*(1.0/(xSpecies[mIdxH2O]*xSpecies[mIdxH2O]) - 1.0/((1.0-xSpecies[mIdxH2O])*(1.0-xSpecies[mIdxH2O]))) ); \
d3gdrds2[j][k][l] = R*t*config; \
d3gdrds2[j][l][k] = d3gdrds2[j][k][l]; \
} \
} \
}

#define fillD3GDRDSDT \
for (j=0; j<NR; j++) { \
for (k=0; k<NS; k++) { \
double config = 0.0; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0) \
config += dxSpeciesdr[i][j]*dxSpeciesds[i][k]/xSpecies[i]; \
if (mIdxH2O != -1 && xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0) \
config += dxSpeciesds[mIdxH2O][k]*dxSpeciesdr[mIdxH2O][j]*(1.0/xSpecies[mIdxH2O] + 1.0/(1.0-xSpecies[mIdxH2O])); \
d3gdrdsdt[j][k] = R*config; \
} \
}

#define fillD3GDRDSDP \
for (k=0; k<NR; k++) for (l=0; l<NS; l++) d3gdrdsdp[k][l] = 0.0;

#define fillD3GDRDT2 \
for (i=0; i<NR; i++) d3gdrdt2[i] = 0.0;

#define fillD3GDRDTDP \
for (i=0; i<NR; i++) d3gdrdtdp[i] = 0.0;

#define fillD3GDRDP2 \
for (i=0; i<NR; i++) d3gdrdp2[i] = 0.0;

#define fillD3GDS3 \
for (j=0; j<NS; j++) { \
for (k=j; k<NS; k++) { \
for (l=k; l<NS; l++) { \
double config = 0.0; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0) config += -dxSpeciesds[i][j]*dxSpeciesds[i][k]*dxSpeciesds[i][l]/(xSpecies[i]*xSpecies[i]); \
if (mIdxH2O != -1 && xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0) \
config += - dxSpeciesds[mIdxH2O][j]*dxSpeciesds[mIdxH2O][k]*dxSpeciesds[mIdxH2O][l] \
*(1.0/(xSpecies[mIdxH2O]*xSpecies[mIdxH2O]) - 1.0/((1.0-xSpecies[mIdxH2O])*(1.0-xSpecies[mIdxH2O]))); \
d3gds3[j][k][l] = R*t*config; \
d3gds3[k][j][l] = d3gds3[j][k][l]; \
d3gds3[l][j][k] = d3gds3[j][k][l]; \
d3gds3[l][k][j] = d3gds3[j][k][l]; \
d3gds3[j][l][k] = d3gds3[j][k][l]; \
d3gds3[k][l][j] = d3gds3[j][k][l]; \
} \
} \
}

#define fillD3GDS2DT \
for (j=0; j<NS; j++) { \
for (k=j; k<NS; k++) { \
double config = 0.0; \
for (i=0; i<NE; i++) if (xSpecies[i] > 0.0) config += dxSpeciesds[i][j]*dxSpeciesds[i][k]/xSpecies[i]; \
if (mIdxH2O != -1 && xSpecies[mIdxH2O] > 0.0 && xSpecies[mIdxH2O] < 1.0) \
config += dxSpeciesds[mIdxH2O][j]*dxSpeciesds[mIdxH2O][k]*(1.0/xSpecies[mIdxH2O] + 1.0/(1.0-xSpecies[mIdxH2O])); \
d3gds2dt[j][k] = R*config; \
d3gds2dt[k][j] = d3gds2dt[j][k]; \
} \
}

#define fillD3GDS2DP \
for (i=0; i<NS; i++) { \
for (l=0; l<NS; l++) d3gds2dp[i][l] = 0.0; \
}

#define fillD3GDSDT2 \
for (i=0; i<NS; i++) d3gdsdt2[i] = 0.0;

#define fillD3GDSDTDP \
for (i=0; i<NS; i++) d3gdsdtdp[i] = 0.0;

#define fillD3GDSDP2 \
for (i=0; i<NS; i++) d3gdsdp2[i] = 0.0;

#define D3GDT3 0.0

#define D3GDT2DP 0.0

#define D3GDTDP2 0.0

#define D3GDP3 0.0

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
    int i, j, k, l, iter = 0;
    double *modelParameters = [[self modelParametersWrapper] pointerToDouble];

    /* look-up or compute the current ordering state */
    if (modelParametersHaveBeenAlteredSinceLastOrderingDetermination || (t != tOld) || (p != pOld) ||
        (r[ 0] != rOld[ 0]) || (r[ 1] != rOld[ 1]) || (r[ 2] != rOld[ 2]) || (r[ 3] != rOld[ 3]) || (r[ 4] != rOld[ 4]) ||
        (r[ 5] != rOld[ 5]) || (r[ 6] != rOld[ 6]) || (r[ 7] != rOld[ 7]) || (r[ 8] != rOld[ 8]) || (r[ 9] != rOld[ 9]) ||
        (r[10] != rOld[10]) || (r[11] != rOld[11]) || (r[12] != rOld[12]) || (r[13] != rOld[13]) || (r[14] != rOld[14]) ) {
        if (r[rIdxCO2] > 0.0) {
            double dgds[NS], sNew[NS], rTotal = 0.0;
            BOOL hasCaO  = ((r[rIdxCaSiO3]  > 0.0) && (r[rIdxCO2] > 0.0)) ? YES : NO;

            for (i=0; i<NS; i++) sOld[i] = 2.0;
            for (i=0; i<NR; i++) rTotal += r[i];

            double denom = 1.0 + ((hasCaO) ? 1.0 : 0.0);
            sNew[0] = (hasCaO ) ? r[rIdxCO2]/denom : 0.0; // if CaO  is present 1/denom of all CO2 is CaCO3

            if (sNew[0] > r[rIdxCaSiO3]) sNew[0] = r[rIdxCaSiO3]/2.0;

            {
                if( (1.0 - rTotal + sNew[0]) < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]",  0);
                if( (r[ 0])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]",  1);
                if( (r[ 1])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]",  2);
                if( (r[ 2])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]",  3);
                if( (r[ 3])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]",  4);
                if( (r[ 4])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]",  5);
                if( (r[ 5])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]",  6);
                if( (r[ 6])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]",  7);
                if( (r[ 7])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]",  8);
                if( (r[ 8])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]",  9);
                if( (r[ 9] - sNew[0])        < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]", 10);
                if( (r[10])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]", 11);
                if( (r[11])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]", 12);
                if( (r[12])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]", 13);
                if( (r[13])                  < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]", 14);
                if( (r[14] - sNew[0])        < 0.0) NSLog(@"Initial guess in odering invalid for xSpecies[%d]", 15);
            }

            while ( ((fabs(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON)) && (iter < MAX_ITER)) {
                double s[NS], deltaS[NS], lambda, d2gds2[NS][NS];

                for (i=0; i<NS; i++) s[i] = sNew[i];

                fillXSPECIES

                for (i=0; i<NE; i++) {
                    if      (xSpecies[i] <= 0.0) xSpecies[i] = DBL_EPSILON;
                    else if (xSpecies[i] >= 1.0) xSpecies[i] = 1.0 - DBL_EPSILON;
                }

                fillDGDS
                fillD2GDS2

                if (!hasCaO ) { dgds[0] = 0.0; } //d2gds2[0][0] = 1.0; }

                for (i=0; i<NS; i++) for (l=0; l<NS; l++) invd2gds2[i][l] = d2gds2[i][l];
                for (i=0; i<NS; i++) sOld[i] = s[i];

                if (NS == 1) invd2gds2[0][0] = (invd2gds2[0][0] != 0.0) ? 1.0/invd2gds2[0][0] : DBL_MAX;
                else         [self gaussj:invd2gds2];

                for (i=0; i<NS; i++) {
                    for(j=0; j<NS; j++) s[i] += - invd2gds2[i][j]*dgds[j];
                    deltaS[i] = s[i] - sOld[i];
                }

                lambda = 1.0;
                BOOL guessIsOkay = NO;
                while (!guessIsOkay && (lambda > DBL_EPSILON)) {
                    guessIsOkay = YES;
                    for (i=0; i<NS; i++) s[i] = sOld[i] + lambda*deltaS[i];
                    if( (1.0 - rTotal + s[0]) < 0.0) guessIsOkay = NO;
                    if( (r[ 0])               < 0.0) guessIsOkay = NO;
                    if( (r[ 1])               < 0.0) guessIsOkay = NO;
                    if( (r[ 2])               < 0.0) guessIsOkay = NO;
                    if( (r[ 3])               < 0.0) guessIsOkay = NO;
                    if( (r[ 4])               < 0.0) guessIsOkay = NO;
                    if( (r[ 5])               < 0.0) guessIsOkay = NO;
                    if( (r[ 6])               < 0.0) guessIsOkay = NO;
                    if( (r[ 7])               < 0.0) guessIsOkay = NO;
                    if( (r[ 8])               < 0.0) guessIsOkay = NO;
                    if( (r[ 9] - s[0])        < 0.0) guessIsOkay = NO;
                    if( (r[10])               < 0.0) guessIsOkay = NO;
                    if( (r[11])               < 0.0) guessIsOkay = NO;
                    if( (r[12])               < 0.0) guessIsOkay = NO;
                    if( (r[13])               < 0.0) guessIsOkay = NO;
                    if( (r[14] - s[0])        < 0.0) guessIsOkay = NO;
                    if( (s[0])                < 0.0) guessIsOkay = NO;
                    lambda = lambda/2.0;
                }

                for (i=0; i<NS; i++) sNew[i] = s[i];
                iter++;
            }
            tOld = t;
            pOld = p;
            for (i=0; i<NR; i++) rOld[i] = r[i];
            modelParametersHaveBeenAlteredSinceLastOrderingDetermination = NO;

            BOOL debug = false; // [[NSUserDefaults standardUserDefaults] boolForKey:@"DEBUG.SATLOOP.SIMPLE"];
            if (debug || (iter >= 200)) {
                for (i=0; i<NS; i++) {
                    if (dgds[i] > sqrt(DBL_EPSILON) && ABS(sOld[i]) > DBL_EPSILON) {
                        NSLog(@"ERROR in Liquid Class (function ORDER). Failed to converge!\n");
                        if (iter >= MAX_ITER) NSLog(@"  Iteration limit (%4d) exceeded.\n", iter);
                        NSLog(@"  s0    = %13.6g\n", sOld[0]);
                        NSLog(@"  dgds0 = %13.6g\n", dgds[0]);
                        break;
                    }
                }
            }

        } else {
            for (i=0; i<NS; i++) {
                sOld[i] = 0.0;
                for (j=0; j<NS; j++) invd2gds2[i][j] = 0.0;
                invd2gds2[i][i] = 0.0;
            }

        }

    }

    if (mask & FIRST  ) {   /* return s        */
        for (i=0; i<NS; i++) s[i] = sOld[i];
    }
    if (mask & SECOND ) {   /* compute ds/dr:  */
        double d2gdrds[NR][NS], rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        fillXSPECIES
        fillD2GDRDS

        for (i=0; i<NS; i++) {
            for (j=0; j<NR; j++) {
                dr[i][j] = 0.0;
                for (k=0; k<NS; k++) dr[i][j] += (s[i] > 0.0) ? - invd2gds2[i][k]*d2gdrds[j][k] : 0.0;
            }
        }
    }
    if (mask & THIRD  ) {   /* compute ds/dt:  */
        double d2gdsdt[NS], rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        fillXSPECIES
        fillD2GDSDT

        for (i=0; i<NS; i++) {
            dt[i] = 0.0;
            for (j=0; j<NS; j++) dt[i] += - invd2gds2[i][j]*d2gdsdt[j];
        }
    }
    if (mask & FOURTH ) {   /* compute ds/dp:  */
        double d2gdsdp[NS], rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        fillXSPECIES
        fillD2GDSDP

        for (i=0; i<NS; i++) {
            dp[i] = 0.0;
            for (j=0; j<NS; j++) dp[i] += - invd2gds2[i][j]*d2gdsdp[j];
        }
    }
    if (mask & FIFTH  ) {   /* compute d2s/dr2 */
        double d2gdrds[NR][NS], d3gdr2ds[NR][NR][NS], d3gdrds2[NR][NS][NS],
        d3gds3[NS][NS][NS], dsdr[NS][NR], temp[NS];
        int m, n;
        double rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        fillXSPECIES
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
        double rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        fillXSPECIES
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
        double rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        fillXSPECIES
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
        double rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        fillXSPECIES
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
        double rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        fillXSPECIES
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
        double rTotal = 0.0;
        for (i=0; i<NR; i++) rTotal += r[i];
        double s[NS];
        for (i=0; i<NS; i++) s[i] = sOld[i];

        fillXSPECIES
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
    const char *NAMES[NA] = { "SiO2", "TiO2", "Al2O3", "Fe2O3", "Cr2O3", "FeO", "MnO", "MgO",
        "NiO", "CoO", "CaO", "Na2O", "K2O", "P2O5", "H2O", "CO2" };
    const char *FORMULAS[NA] = { "SiO2", "TiO2", "Al2O3", "Fe2O3", "MgCr2O4", "Fe2SiO4", "MnSi0.5O2", "Mg2SiO4",
        "NiSi0.5O2", "CoSi0.5O2", "CaSiO3", "Na2SiO3", "KAlSiO4", "Ca3(PO4)2", "H2O", "CO2" };
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
        static const int  C =  6;
        static const int  O =  8;
        static const int  F =  9;
        static const int Na = 11;
        static const int Mg = 12;
        static const int Al = 13;
        static const int Si = 14;
        static const int  P = 15;
        //        static const int  S = 16;
        static const int Cl = 17;
        static const int  K = 19;
        static const int Ca = 20;
        static const int Ti = 22;
        static const int Cr = 24;
        static const int Mn = 25;
        static const int Fe = 26;
        static const int Co = 27;
        static const int Ni = 28;
        double ox;


        m[ 0] = e[Si] - e[Mn]/2.0 - (e[Mg]/2.0-e[Cr]/4.0) - e[Ni]/2.0 - e[Co]/2.0 - (e[Ca]-3.0*e[P]/2.0) - e[Na]/2.0 - e[K];     // SiO2
        m[ 1] = e[Ti];                 // TiO2
        m[ 2] = e[Al]/2.0-e[K]/2.0;    // Al2O3
        m[ 4] = e[Cr]/2.0;             // MgCr2O4
        m[ 6] = e[Mn];                 // MnSi0.5O2
        m[ 7] = e[Mg]/2.0-e[Cr]/4.0;   // Mg2SiO4
        m[ 8] = e[Ni];                 // NiSi0.5O2
        m[ 9] = e[Co];                 // CoSi0.5O2
        m[10] = e[Ca]-3.0*e[P]/2.0;    // CaSiO3
        m[11] = e[Na]/2.0;             // Na2SiO3
        m[12] = e[K];                  // KAlSiO4
        m[13] = e[P]/2.0;              // Ca3(PO4)2
        m[14] = (e[H]-e[F]-e[Cl])/2.0; // H2O
        m[15] = e[C];                  // CO2

        //        m[16] = e[F];                  // HF
        //        m[17] = e[Cl];                 // HCl
        //        m[18] = e[S];                  // SO2

        if (fabs(m[14]) < 100.0*DBL_EPSILON) m[14] = 0.0;

        ox = 2.0*m[0] + 2.0*m[1] + 3.0*m[2] + 4.0*m[4] + 2.0*m[6] + 4.0*m[7] + 2.0*m[8]
        + 2.0*m[9] + 3.0*m[10] + 3.0*m[11] + 4.0*m[12] + 8.0*m[13] + m[14] + m[15]*2.0;

        // e[O] - (ox - 2 m[5]) = 3 m[3] + 4 m[5]  O balance, note production of m[5] reduces the amount of m[0], hence the correction to ox
        // e[Fe] = 2 m[3] + 2 m[5]                 Fe balance

        if (e[Fe] > 10.0*DBL_EPSILON) {
            m[3] = e[O] - ox - e[Fe];       // Fe2O3
            if (fabs(m[3]) < 100.0*DBL_EPSILON) m[3] = 0.0;
            m[5] = (e[Fe] - 2.0*m[3])/2.0;  // Fe2SiO4
            if (fabs(m[5]) < 100.0*DBL_EPSILON) m[5] = 0.0;
            m[0] -= m[5];                   // readjust SiO2 moles for the Fe2SiO4 just made
        } else {
            m[3] = 0.0;
            m[5] = 0.0;
        }


    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
            NSLog(@"Illegal call to convert with inpMask = %o and outMask = %o", inpMask, outMask);

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
    double sio2 =  (1.0-r[0]-r[1]-r[2]-r[3]-r[5]/2.0-r[7]/2.0-r[8]/2.0-r[12]-r[13]-r[14])*60.0843;
    double tio2  = r[0]*79.8658;
    double al2o3 = (r[1] + r[11]/2.0)*101.9613;
    double fe2o3 = r[2]*159.6882;
    double cr2o3 = r[3]*151.9904;
    double feo   = 2.0*r[4]*71.8444;
    double mno   = r[5]*70.93745;
    double mgo   = (r[3] + 2.0*r[6])*40.3044;
    double nio   = r[7]*74.6928;
    double coo   = r[8]*74.9326;
    double cao   = (r[9] + 3.0*r[12])*56.0774;
    double na2o  = r[10]*61.97894;
    double k2o   = r[11]*94.196/2.0;
    double p2o5  = r[12]*141.9445;
    double h2o   = r[13]*18.01528;
    double co2   = r[14]*44.0095;
    double sum = sio2 + tio2 + al2o3 + fe2o3 + cr2o3 + feo + mno + mgo + nio + coo + cao + na2o + k2o + p2o5 + h2o + co2;
    if (sum == 0.0) sum = 1.0;

    return [NSString stringWithFormat:@"wt%%:SiO2 %5.2f TiO2 %5.2f Al2O3 %5.2f Fe2O3 %5.2f Cr2O3 %5.2f FeO %5.2f MnO %5.2f MgO %5.2f NiO %5.2f CoO %5.2f CaO %5.2f Na2O %5.2f K2O %5.2f P2O5 %5.2f H2O %5.2f CO2 %5.2f",
            sio2*100.0/sum, tio2*100.0/sum, al2o3*100.0/sum, fe2o3*100.0/sum, cr2o3*100.0/sum, feo*100.0/sum, mno*100.0/sum,
            mgo*100.0/sum, nio*100.0/sum, coo*100.0/sum, cao*100.0/sum, na2o*100.0/sum, k2o*100.0/sum, p2o5*100.0/sum,
            h2o*100.0/sum, co2*100.0/sum];
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
    int i, j, k;
    double *modelParameters = [[self modelParametersWrapper] pointerToDouble];
    double rTotal = 0.0;
    for (i=0; i<NR; i++) rTotal += r[i];

    for (i=0; i<NR; i++) fr[0][i] = -r[i];
    for (j=1; j<NA; j++) for (i=0; i<NR; i++) fr[j][i] = (i+1 == j) ? 1.0 - r[i] : - r[i];

    [self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

    fillXSPECIES
    fillG
    fillDGDR

    if (mask & FIRST) {
        for (i=0; i<NA; i++) {
            for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
            a[i] = exp(a[i]/(R*t));
        }
    }

    if (mask & SECOND) {
        for (i=0; i<NA; i++) {
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

        for (j=0; j<NA; j++) for (i=0; i<NR; i++) dfrdr[j][i] = -1.0;
        for (j=0; j<NA; j++) for (i=0; i<NS; i++) gs[j][i]    = -s[i];
        for (j=0; j<NA; j++) for (i=0; i<NS; i++) dgsds[j][i] = -1.0;

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
         dx:(double [NR])dr          // (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010
        dx2:(double [NR][NR])dr2     // (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100
        dx3:(double [NR][NR][NR])dx3 // (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000
{
    double s[NS];
    double rTotal = 0.0;
    int i;
    for (i=0; i<NR; i++) rTotal += r[i];
    double *modelParameters = [[self modelParametersWrapper] pointerToDouble];

    [self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    fillXSPECIES

    if (mask & FIRST) {
        double g;
        int j;

        fillG
        *gmix = g;
    }

    if(mask & SECOND) {
        double dgdr[NR];
        int j, k;

        fillDGDR
        for (i=0; i<NR; i++) dr[i] = dgdr[i];
    }

    if(mask & THIRD) {
        double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR];
        int j, k, l;

        fillD2GDR2
        fillD2GDRDS
        fillD2GDS2

        [self order:SECOND t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dr2[i][j] = d2gdr2[i][j];
                for (k=0; k<NS; k++) {
                    dr2[i][j] += d2gdrds[i][k]*dsdr[k][j] + d2gdrds[j][k]*dsdr[k][i];
                    for (l=0; l<NS; l++) dr2[i][j] += d2gds2[k][l]*dsdr[k][i]*dsdr[l][j];
                }
            }
        }
    }

    // To be completed
    if(mask & FOURTH) {
        double d3gdr3[NR][NR][NR], d3gdr2ds[NR][NR][NS], d3gdrds2[NR][NS][NS];
        double d3gds3[NS][NS][NS], dsdr[NS][NR];
        int j, k, l, m, n;

        for (i=0; i<NR; i++) for (j=0; j<NR; j++) for (k=0; k<NR; k++) d3gdr3[i][j][k] = 0.0;
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
{
    double s[NS], g, entropy;
    double rTotal = 0.0;
    int i, j;
    for (i=0; i<NR; i++) rTotal += r[i];
    double *modelParameters = [[self modelParametersWrapper] pointerToDouble];

    [self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    fillXSPECIES

    fillG
    fillS
    *hmix = g + t*entropy;
}

-(void)smix:(int)mask
          t:(double)t
          p:(double)p
          r:(double [NR])r
       smix:(double *)smix        // Entropy of mixing                  BINARY MASK: 001
         dx:(double [NR])dr       // (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010
        dx2:(double [NR][NR])dr2  // (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100
{
    double s[NS];
    double rTotal = 0.0;
    int i;
    for (i=0; i<NR; i++) rTotal += r[i];
    double *modelParameters = [[self modelParametersWrapper] pointerToDouble];

    [self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    fillXSPECIES

    if (mask & FIRST) {
        double entropy;

        fillS
        *smix = entropy;
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NS], d2gdrdt[NR], d2gds2[NS][NS], d2gdsdt[NS], dsdr[NS][NR], dsdt[NS];
        int j, k, l;

        fillD2GDRDS
        fillD2GDRDT
        fillD2GDS2
        fillD2GDSDT

        [self order:SECOND | THIRD t:t p:p r:r s:NULL dr:dsdr dt:dsdt dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        for (i=0; i<NR; i++) {
            dr[i] = d2gdrdt[i];
            for (k=0; k<NS; k++) {
                dr[i] += d2gdrds[i][k]*dsdt[k] + d2gdsdt[k]*dsdr[k][i];
                for (l=0; l<NS; l++) dr[i] += d2gds2[k][l]*dsdt[k]*dsdr[l][i] ;
            }
            dr[i] *= -1.0;
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
                dr2[i][j] = d3gdr2dt[i][j];
                for (k=0; k<NS; k++) {
                    dr2[i][j] += d3gdr2ds[i][j][k]*dsdt[k]
                    + d3gdrdsdt[i][k]*dsdr[k][j]
                    + d3gdrdsdt[j][k]*dsdr[k][i]
                    + d2gdsdt[k]*d2sdr2[k][i][j]
                    + d2gdrds[i][k]*d2sdrdt[k][j]
                    + d2gdrds[j][k]*d2sdrdt[k][i];
                    for (l=0; l<NS; l++) {
                        dr2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdt[l]
                        + d3gdrds2[j][k][l]*dsdr[k][i]*dsdt[l]
                        + d2gds2[k][l]*d2sdr2[k][i][j]*dsdt[l]
                        + d3gds2dt[k][l]*dsdr[k][i]*dsdr[l][j]
                        + d2gds2[k][l]*dsdr[k][i]*d2sdrdt[l][j]
                        + d2gds2[k][l]*dsdr[k][j]*d2sdrdt[l][i];
                        for (m=0; m<NS; m++)
                            dr2[i][j] += d3gds3[k][l][m]*dsdr[k][i]*dsdr[l][j]*dsdt[m];
                    }
                }
                dr2[i][j] *= -1.0;
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
          dx:(double [NR])dr // d(cp)/d(x[])d(t)                BINARY MASK: 100
{
    double s[NS], dsdt[NS], d2gdsdt[NS], d2gds2[NS][NS], d2gdt2;
    double rTotal = 0.0;
    int i, j, k, l;
    for (i=0; i<NR; i++) rTotal += r[i];
    double *modelParameters = [[self modelParametersWrapper] pointerToDouble];

    [self order:FIRST | THIRD t:t p:p r:r s:s dr:NULL dt:dsdt dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    fillXSPECIES
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
        double d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gdsdt2[NS], d2sdt2[NS], temp;
        double d3gdt3 = D3GDT3;

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

        fillD2GDRDS
        fillD3GDRDS2
        fillD3GDRDSDT
        fillD3GDRDT2
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2

        [self order:SECOND | SIXTH | EIGHTH t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:NULL dr2:NULL drt:d2sdrdt drp:NULL dt2:d2sdt2 dtp:NULL dp2:NULL];

        for (i=0; i<NR; i++) {
            for (j=0,dr[i]=d3gdrdt2[i]; j<NS; j++) {
                dr[i] += d3gdsdt2[j]*dsdr[j][i] + 2.0*d2gdsdt[j]*d2sdrdt[j][i] +
                2.0*d3gdrdsdt[i][j]*dsdt[j] + d2gdrds[i][j]*d2sdt2[j];
                for (k=0; k<NS; k++) {
                    dr[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdt[k] +
                    2.0*d2gds2[j][k]*dsdt[j]*d2sdrdt[k][i] +
                    2.0*d3gds2dt[j][k]*dsdr[j][i]*dsdt[k] +
                    d2gds2[j][k]*dsdr[j][i]*d2sdt2[k];
                    for (l=0; l<NS; l++)
                        dr[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdt[k]*dsdt[l];
                }
            }
            dr[i] *= -t;
        }
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
    double s[NS];
    double rTotal = 0.0;
    int i, j;
    for (i=0; i<NR; i++) rTotal += r[i];
    double *modelParameters = [[self modelParametersWrapper] pointerToDouble];

    [self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];
    fillXSPECIES

    if (mask & FIRST) {
        double v;

        fillV
        *vmix = v;
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NS], d2gdrdp[NR], d2gds2[NS][NS], d2gdsdp[NS], dsdr[NS][NR], dsdp[NS];
        int k, l;

        fillD2GDRDS
        fillD2GDRDP
        fillD2GDS2
        fillD2GDSDP

        [self order:SECOND | FOURTH t:t p:p r:r s:NULL dr:dsdr dt:NULL dp:dsdp dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

        for (i=0; i<NR; i++) {
            dr[i] = d2gdrdp[i];
            for (j=0; j<NS; j++) {
                dr[i] += d2gdrds[i][j]*dsdp[j] + d2gdsdp[j]*dsdr[j][i];
                for (k=0; k<NS; k++) dr[i] += d2gds2[j][k]*dsdp[j]*dsdr[k][i];
            }
        }
    }

    if(mask & THIRD) {
        double d2gdrds[NR][NS], d2gds2[NS][NS], d2gdsdp[NS], d3gdr2ds[NR][NR][NS],
        d3gdr2dp[NR][NR], d3gdrds2[NR][NS][NS], d3gdrdsdp[NR][NS],
        d3gds3[NS][NS][NS], d3gds2dp[NS][NS], dsdr[NS][NR], dsdp[NS],
        d2sdr2[NS][NR][NR], d2sdrdp[NS][NR];
        int k, l, m;

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
                dr2[i][j] = d3gdr2dp[i][j];
                for (k=0; k<NS; k++) {
                    dr2[i][j] += d3gdr2ds[i][j][k]*dsdp[k]
                    + d3gdrdsdp[i][k]*dsdr[k][j]
                    + d3gdrdsdp[j][k]*dsdr[k][i]
                    + d2gdsdp[k]*d2sdr2[k][i][j]
                    + d2gdrds[i][k]*d2sdrdp[k][j]
                    + d2gdrds[j][k]*d2sdrdp[k][i];
                    for (l=0; l<NS; l++) {
                        dr2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdp[l]
                        + d3gdrds2[j][k][l]*dsdr[k][i]*dsdp[l]
                        + d2gds2[k][l]*d2sdr2[k][i][j]*dsdp[l]
                        + d3gds2dp[k][l]*dsdr[k][i]*dsdr[l][j]
                        + d2gds2[k][l]*dsdr[k][i]*d2sdrdp[l][j]
                        + d2gds2[k][l]*dsdr[k][j]*d2sdrdp[l][i];
                        for (m=0; m<NS; m++)
                            dr2[i][j] += d3gds3[k][l][m]*dsdr[k][i]*dsdr[l][j]*dsdp[m];
                    }
                }
            }
        }
    }

    if(mask & FOURTH) {
        double d2gdtdp = D2GDTDP;
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], dsdt[NS], dsdp[NS];
        int k, l;

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
        int k, l;

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
        int k, l;

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
        int k, l;

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
        int k, l;

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
        int k, l;

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
            for (j=0,drdt[i]=d3gdrdtdp[i]; j<NS; j++) {
                drdt[i] += d3gdsdtdp[j]*dsdr[j][i] + d2gdsdt[j]*d2sdrdp[j][i] +
                d3gdrdsdt[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdtdp[j] +
                d3gdrdsdp[i][j]*dsdt[j] + d2gdsdp[j]*d2sdrdt[j][i];
                for (k=0; k<NS; k++) {
                    drdt[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdp[k] +
                    d2gds2[j][k]*dsdt[j]*d2sdrdp[k][i] +
                    d2gds2[j][k]*dsdp[j]*d2sdrdt[k][i] +
                    d3gds2dt[j][k]*dsdr[j][i]*dsdp[k] +
                    d3gds2dp[j][k]*dsdr[j][i]*dsdt[k] +
                    d2gds2[j][k]*dsdr[j][i]*d2sdtdp[k];
                    for (l=0; l<NS; l++)
                        drdt[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdt[k]*dsdp[l];
                }
            }
        }
    }

    if(mask & TENTH) {
        double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gds2dp[NS][NS],
        d2gdrds[NR][NS], dsdp[NS], dsdr[NS][NR], d2sdrdp[NS][NR], d2gds2[NS][NS],
        d3gdrdsdp[NR][NS], d3gdrdp2[NR], d3gdsdp2[NS], d2gdsdp[NS], d2sdp2[NS];
        int k, l;

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
            for (j=0,drdp[i]=d3gdrdp2[i]; j<NS; j++) {
                drdp[i] += d3gdsdp2[j]*dsdr[j][i] + d2gdsdp[j]*d2sdrdp[j][i] +
                2.0*d3gdrdsdp[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdp2[j] +
                d2gdsdp[j]*d2sdrdp[j][i];
                for (k=0; k<NS; k++) {
                    drdp[i] += d3gdrds2[i][j][k]*dsdp[j]*dsdp[k] +
                    2.0*d2gds2[j][k]*dsdp[j]*d2sdrdp[k][i] +
                    2.0*d3gds2dp[j][k]*dsdr[j][i]*dsdp[k] +
                    d2gds2[j][k]*dsdr[j][i]*d2sdp2[k];
                    for (l=0; l<NS; l++)
                        drdp[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdp[k]*dsdp[l];
                }
            }
        }
    }
}

#define NAS NE

-(NSUInteger)numberOfSolutionSpecies {
    return NE;
}

-(NSString *)nameOfSolutionSpeciesAtIndex:(NSUInteger)index {
    switch (index) {
        case 0:
            return @"SiO2";
            break;
        case 1:
            return @"TiO2";
            break;
        case 2:
            return @"Al2O3";
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
            return @"CaSiO3";
            break;
        case 11:
            return @"Na2SiO3";
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
        case 15:
            return @"CO2";
            break;
        case 16:
            return @"CaCO3";
            break;
        default:
            return @"";
            break;
    }
}

-(NSDictionary *)getSpeciesMoleFractionsForBulkComposition:(double *)m aT:(double)t andP:(double)p {
    double r[NR], s[NS];
    double rTotal = 0.0;
    int i;

    [self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:m r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
    for (i=0; i<NR; i++) rTotal += r[i];
    [self order:FIRST t:t p:p r:r s:s dr:NULL dt:NULL dp:NULL dr2:NULL drt:NULL drp:NULL dt2:NULL dtp:NULL dp2:NULL];

    fillXSPECIES

    NSMutableDictionary *result = [NSMutableDictionary dictionaryWithCapacity:NE];
    for (i=0; i<NE; i++) [result setObject:[NSNumber numberWithDouble:xSpecies[i]] forKey:[self nameOfSolutionSpeciesAtIndex:i]];

    return [NSDictionary dictionaryWithDictionary:result];
}

-(DoubleVector *)convertMolesOfSpeciesToMolesOfComponents:(double *)mSpecies {
    DoubleVector *mComponentWrapper =[[DoubleVector alloc] initWithSize:NA];
    double *mComponents = [mComponentWrapper pointerToDouble];
    for (NSUInteger i=0; i<NA; i++) mComponents[i] = mSpecies[i];
    mComponents[10] += mSpecies[16]; // CaCO3  -> CaSiO3
    mComponents[ 0] -= mSpecies[16]; // SiO2 must be rediuced
    mComponents[15] += mSpecies[16]; // CO2 must be increased
    return mComponentWrapper;
}

-(DoubleVector *)chemicalPotentialsOfSpeciesFromMolesOfComponents:(double *)m andT:(double)t andP:(double)p {
    DoubleVector *muSpeciesWrapper = [[DoubleVector alloc] initWithSize:NE];
    double *muSpecies = [muSpeciesWrapper pointerToDouble];
    double *muComponents = [[self getChemicalPotentialFromMolesOfComponents:m andT:t andP:p] pointerToDouble];
    for (NSUInteger i=0; i<NA; i++) muSpecies[i] = muComponents[i];
    // CaCO3
    muSpecies[16] = muSpecies[10] + muSpecies[15] - muSpecies[0];
    return muSpeciesWrapper;
}

-(DoubleVector *)elementalCompositionOfSpeciesAtIndex:(NSUInteger)index {
    DoubleVector *speciesElementArrayWrapper = nil;

    if (index < NA) speciesElementArrayWrapper = [[endmembers objectAtIndex:index] formulaAsElementArray];
    else if (index == 16) {
        speciesElementArrayWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
        double *speciesElementArray = [speciesElementArrayWrapper pointerToDouble];
        double *component = [[[endmembers objectAtIndex:10] formulaAsElementArray] pointerToDouble];
        for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] += component[i];
        component = [[[endmembers objectAtIndex:15] formulaAsElementArray] pointerToDouble];
        for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] += component[i];
        component = [[[endmembers objectAtIndex:0] formulaAsElementArray] pointerToDouble];
        for (NSUInteger i=1; i<107; i++) if (component[i] != 0.0) speciesElementArray[i] -= component[i];
    }
    return speciesElementArrayWrapper;
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
    double mu0[NAS], deltaMu[NAS], xNz[NAS], x[NAS], gamma[NAS], xLast[NAS], affinity = 0.0;
    NSUInteger i, j, nz = 0, index[NAS];

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
            nz++;
        }
        x[i] = 0.0;
        xLast[i] = 0.0;
    }

    // CaCO3
    if ((chemicalPotentials[0] != 0.0) && (chemicalPotentials[10] != 0.0) && (chemicalPotentials[15] != 0.0)) {
        mu0[16] = mu0[10] + mu0[15] - mu0[0] + (deltaHsp[0] - t*deltaSsp[0] + (p-1.0)*deltaVsp[0]);
        deltaMu[nz] = chemicalPotentials[10] + chemicalPotentials[15] - chemicalPotentials[0] - mu0[16];
        index[nz] = 16;
        gamma[nz] = 1.0;
        nz++;
    } else mu0[16] = 0.0;
    x[16] = 0.0;
    xLast[16] = 0.0;

    // There are no non-zero chemical potentials, so the phase can never form.  Return an affinity of zero and zero the mole fractions
    if (nz == 0) {
        [results addObject:[NSNumber numberWithDouble:affinity]];
        for (i=0; i<NA; i++) [results addObject:[NSNumber numberWithDouble:x[i]]];
        [results addObject:[NSNumber numberWithBool:YES]];          // convergence flag
        [results addObject:[NSNumber numberWithUnsignedInteger:0]]; // iteration count
        [results addObject:[NSNumber numberWithDouble:NATOMS]];     // number of atoms used to scale affinity
        [results addObject:[NSNumber numberWithDouble:0.0]];        // likely error in affinity

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
        double a[NAS], mu[NA], r[NR];
        for (i=0; i<NA; i++) xReduced[i] = x[i];
        xReduced[ 0] -= x[16];
        xReduced[10] += x[16];
        xReduced[15] += x[16];
        [self convert:SECOND outMask:THIRD t:0.0 p:0.0 e:NULL m:xReduced r:r x:NULL dm:NULL d2m:NULL dr:NULL d3m:NULL];
        [self activity:FIRST | SECOND t:t p:p r:r a:a mu:mu dx:NULL];
        if (![self testPermissibleValuesOfComponents:xReduced]) {
            NSLog(@"Composition estimate is infeasible.");
#if TARGET_OS_IPHONE || TARGET_IPHONE_SIMULATOR
            for (i=0; i<NAS; i++)
                NSLog(@"species X%lu %@ = %g", (unsigned long)i, [[self nameOfSolutionSpeciesAtIndex:i] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[i]);
#else
            for (i=0; i<NAS; i++)
                NSLog(@"species X%lu %@ = %g", i, [[self nameOfSolutionSpeciesAtIndex:i] stringByPaddingToLength:15 withString:@" " startingAtIndex:0], x[i]);
#endif
        }

        for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) gamma[j++] = a[i]/x[i];

        if ((chemicalPotentials[0] != 0.0) && (chemicalPotentials[10] != 0.0) && (chemicalPotentials[15] != 0.0)) {
            a[16] = exp((mu[10] + mu[15] - mu[0] + mu0[10] + mu0[15] - mu0[0] - mu0[16])/(R*t));
            gamma[j++] = a[16]/x[16];
        }

        [self correctActivityCoefficients:gamma forComposition:x];
        converged = (fabs(affinity-affinityLast) < 0.1);
        count++;

    } while (count < 50 && !converged);

    if (debugS) {
#if TARGET_OS_IPHONE || TARGET_IPHONE_SIMULATOR
        NSLog(@"... Terminated (converged %@) for phase %@ in %lu iterations with delta affinity %f J for %f atoms.",
              converged ? @"YES" : @"NO", [self phaseName], (unsigned long)count, fabs(affinity-affinityLast), NATOMS);
#else
        NSLog(@"... Terminated (converged %@) for phase %@ in %lu iterations with delta affinity %f J for %f atoms.",
              converged ? @"YES" : @"NO", [self phaseName], count, fabs(affinity-affinityLast), NATOMS);
#endif
        for (i=0, j=0; i<NA; i++) if (x[i] != 0.0) NSLog(@"... ... Activity coefficient of component %@ is %f with mole fraction %f",
                                                         [[endmembers objectAtIndex:i] phaseName], gamma[j++], x[i]);
        NSLog(@"... ... Activity coefficient of component %@ is %f with mole fraction %f", @"CaCO3",  gamma[j++], x[16]);
    }
    if (debugV) NSLog(@"Exiting [... affinityAndCompositionFromLiquidChemicalPotentialSum].");


    [results addObject:[NSNumber numberWithDouble:affinity]];                         // affinity in J
    for (i=0; i<NA; i++) [results addObject:[NSNumber numberWithDouble:xReduced[i]]]; // composition in mole fraction of endmembers
    [results addObject:[NSNumber numberWithBool:converged]];                          // convergence flag
    [results addObject:[NSNumber numberWithUnsignedInteger:count]];                   // iteration count
    [results addObject:[NSNumber numberWithDouble:NATOMS]];                           // number of atoms used to scale affinity
    [results addObject:[NSNumber numberWithDouble:fabs(affinity-affinityLast)]];      // likely error in affinity

    return [NSArray arrayWithArray:results];
}

#include "SolutionPhase.h"

@end
