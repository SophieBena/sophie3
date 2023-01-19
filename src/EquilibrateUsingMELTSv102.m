//
//  EquilibrateUsingMELTSv102.m
//  PhaseMELTSobjC
//
//  Created by Mark Ghiorso on 3/1/17.
//  Copyright © 2017 Mark Ghiorso. All rights reserved.
//

#import "EquilibrateUsingMELTSv102.h"

// This provides the class support to translate chemical formulas into molecular weights and moles of elements
#import "PhaseBase.h"

// Support the DoubleVector class for method override
#import "DoubleVector.h"

// This is the Equilibrate state that is eventually passed back to the main thread
#import "EquilibrateState.h"

// These are the phases known to MELTS
#import "BermanStoichiometricPhases.h"
#import "LiquidMelts.h"
#import "OlivineBerman.h"
#import "FeldsparBerman.h"
#import "SpinelBerman.h"
#import "CpxBerman.h"
#import "OpxBerman.h"
#import "LiquidMeltsH2O.h"
#import "GarnetBerman.h"
#import "RhombohedralBerman.h"
#import "LeuciteBerman.h"
#import "NephelineSSBerman.h"
#import "KalsiliteSSBerman.h"
#import "BiotiteBerman.h"
#import "HornblendeBerman.h"
#import "AlloyLiquid.h"
#import "AlloySolid.h"
#import "MeliliteBerman.h"
#import "OrthoOxideBerman.h"
#import "ClinoamphiboleBerman.h"
#import "OrthoamphiboleBerman.h"

@implementation EquilibrateUsingMELTSv102

-(EquilibrateStatePhase *)allocatePhaseWithName:(NSString *)name phaseClass:(Class)phaseClass {
    EquilibrateStatePhase *phase = [[EquilibrateStatePhase alloc] init];
    if ([phaseClass instancesRespondToSelector:@selector(initWithCompositionConstraint:)])
        [phase setPhaseClassInstance:[[phaseClass alloc] initWithCompositionConstraint:name]];
    else [phase setPhaseClassInstance:[[phaseClass alloc] init]];
    if ([phaseClass instancesRespondToSelector:@selector(setOperationParent:)]) [[phase phaseClassInstance] setOperationParent:@"Python"];
    return phase;
}

- (instancetype)init {
    [[NSUserDefaults standardUserDefaults] setInteger:1 forKey:@"PPCalculationDatabase"];
    self = [super initWithInformationalDebug:NO
                             andVerboseDebug:NO
                             andSystemOxides:[NSArray arrayWithObjects:@"SiO2", @"TiO2", @"Al2O3", @"Fe2O3",
                                              @"Cr2O3", @"FeO", @"MnO", @"MgO", @"NiO", @"CoO", @"CaO", @"Na2O",
                                              @"K2O", @"P2O5", @"H2O", nil]
               andDictionaryOfPhasesInSystem:[NSDictionary dictionaryWithObjectsAndKeys:
                    [self allocatePhaseWithName:@"Liquid"         phaseClass:[LiquidMelts class]         ], @"Liquid",
                    [self allocatePhaseWithName:@"Olivine"        phaseClass:[OlivineBerman class]       ], @"Olivine",
                    [self allocatePhaseWithName:@"Plagioclase"    phaseClass:[FeldsparBerman class]      ], @"Plagioclase",
                    [self allocatePhaseWithName:@"Sanidine"       phaseClass:[FeldsparBerman class]      ], @"Sanidine",
                    [self allocatePhaseWithName:@"Augite"         phaseClass:[CpxBerman class]           ], @"Augite",
                    [self allocatePhaseWithName:@"Pigeonite"      phaseClass:[CpxBerman class]           ], @"Pigeonite",
                    [self allocatePhaseWithName:@"Titanaugite"    phaseClass:[CpxBerman class]           ], @"Titanaugite",
                    [self allocatePhaseWithName:@"Spinel"         phaseClass:[SpinelBerman class]        ], @"Spinel",
                    [self allocatePhaseWithName:@"Orthopyroxene"  phaseClass:[OpxBerman class]           ], @"Orthopyroxene",
                    [self allocatePhaseWithName:@"Quartz"         phaseClass:[QuartzBerman class]        ], @"Quartz",
                    [self allocatePhaseWithName:@"Garnet"         phaseClass:[GarnetBerman class]        ], @"Garnet",
                    [self allocatePhaseWithName:@"Ilmenite ss"    phaseClass:[RhombohedralBerman class]  ], @"Ilmenite ss",
                    [self allocatePhaseWithName:@"Leucite"        phaseClass:[LeuciteBerman class]       ], @"Leucite",
                    [self allocatePhaseWithName:@"Nepheline ss"   phaseClass:[NephelineSSBerman class]   ], @"Nepheline ss",
                    [self allocatePhaseWithName:@"Panunzite"      phaseClass:[NephelineSSBerman class]   ], @"Panunzite",
                    [self allocatePhaseWithName:@"Kalsilite ss"   phaseClass:[KalsiliteSSBerman class]   ], @"Kalsilite ss",
                    [self allocatePhaseWithName:@"Melilite"       phaseClass:[MeliliteBerman class]      ], @"Melilite",
                    [self allocatePhaseWithName:@"Biotite"        phaseClass:[BiotiteBerman class]       ], @"Biotite",
                    [self allocatePhaseWithName:@"Hornblende"     phaseClass:[HornblendeBerman class]    ], @"Hornblende",
                    [self allocatePhaseWithName:@"OrthoOxide"     phaseClass:[OrthoOxideBerman class]    ], @"OrthoOxide",
                    [self allocatePhaseWithName:@"Actinolite"     phaseClass:[ClinoamphiboleBerman class]], @"Actinolite",
                    [self allocatePhaseWithName:@"Cummingtonite"  phaseClass:[ClinoamphiboleBerman class]], @"Cummingtonite",
                    [self allocatePhaseWithName:@"Anthophyllite"  phaseClass:[OrthoamphiboleBerman class]], @"Anthophyllite",
                    [self allocatePhaseWithName:@"Liquid Alloy"   phaseClass:[AlloyLiquid class]         ], @"Liquid Alloy",
                    [self allocatePhaseWithName:@"Solid Alloy"    phaseClass:[AlloySolid class]          ], @"Solid Alloy",
                    [self allocatePhaseWithName:@"Water"          phaseClass:[WaterMelts class]          ], @"Water",
                    [self allocatePhaseWithName:@"Aegirine"       phaseClass:[AegirineBerman class]      ], @"Aegirine",
                    [self allocatePhaseWithName:@"Aenigmatite"    phaseClass:[AenigmatiteBerman class]   ], @"Aenigmatite",
                    [self allocatePhaseWithName:@"Akermanite"     phaseClass:[AkermaniteBerman class]    ], @"Akermanite",
                    [self allocatePhaseWithName:@"Andalusite"     phaseClass:[AndalusiteBerman class]    ], @"Andalusite",
                    [self allocatePhaseWithName:@"Apatite"        phaseClass:[ApatiteBerman class]       ], @"Apatite",
                    [self allocatePhaseWithName:@"Chromite"       phaseClass:[ChromiteBerman class]      ], @"Chromite",
                    [self allocatePhaseWithName:@"Coesite"        phaseClass:[CoesiteBerman class]       ], @"Coesite",
                    [self allocatePhaseWithName:@"Corundum"       phaseClass:[CorundumBerman class]      ], @"Corundum",
                    [self allocatePhaseWithName:@"Cristobalite"   phaseClass:[CristobaliteBerman class]  ], @"Cristobalite",
                    [self allocatePhaseWithName:@"Fayalite"       phaseClass:[FayaliteBerman class]      ], @"Fayalite",
                    [self allocatePhaseWithName:@"Forsterite"     phaseClass:[ForsteriteBerman class]    ], @"Forsterite",
                    [self allocatePhaseWithName:@"Gehlenite"      phaseClass:[GehleniteBerman class]     ], @"Gehlenite",
                    [self allocatePhaseWithName:@"Hematite"       phaseClass:[HematiteBerman class]      ], @"Hematite",
                    [self allocatePhaseWithName:@"Ilmenite"       phaseClass:[IlmeniteBerman class]      ], @"Ilmenite",
                    [self allocatePhaseWithName:@"Kalsilite"      phaseClass:[KalsiliteBerman class]     ], @"Kalsilite",
                    [self allocatePhaseWithName:@"Kyanite"        phaseClass:[KyaniteBerman class]       ], @"Kyanite",
                    [self allocatePhaseWithName:@"Lime"           phaseClass:[LimeBerman class]          ], @"Lime",
                    [self allocatePhaseWithName:@"Magnetite"      phaseClass:[MagnetiteBerman class]     ], @"Magnetite",
                    [self allocatePhaseWithName:@"Muscovite"      phaseClass:[MuscoviteBerman class]     ], @"Muscovite",
                    [self allocatePhaseWithName:@"Nepheline"      phaseClass:[NephelineBerman class]     ], @"Nepheline",
                    [self allocatePhaseWithName:@"Periclase"      phaseClass:[PericlaseBerman class]     ], @"Periclase",
                    [self allocatePhaseWithName:@"Perovskite"     phaseClass:[PerovskiteBerman class]    ], @"Perovskite",
                    [self allocatePhaseWithName:@"Phlogopite"     phaseClass:[PhlogopiteBerman class]    ], @"Phlogopite",
                    [self allocatePhaseWithName:@"Rutile"         phaseClass:[RutileBerman class]        ], @"Rutile",
                    [self allocatePhaseWithName:@"Sillimanite"    phaseClass:[SillimaniteBerman class]   ], @"Sillimanite",
                    [self allocatePhaseWithName:@"Sphene"         phaseClass:[SpheneBerman class]        ], @"Sphene",
                    [self allocatePhaseWithName:@"Tridymite"      phaseClass:[TridymiteBerman class]     ], @"Tridymite",
                    [self allocatePhaseWithName:@"Whitlockite"    phaseClass:[WhitlockiteBerman class]   ], @"Whitlockite",
                    nil]];
    if (self) {
    }
    return self;
}

-(DoubleVector *)molesOfOxidesFromMolesOfElements:(double *)e {
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
    DoubleVector *mWrapper = [[DoubleVector alloc] initWithSize:15];
    double *m = [mWrapper pointerToDouble];


    m[ 0] = e[Si];     // SiO2
    m[ 1] = e[Ti];     // TiO2
    m[ 2] = e[Al]/2.0; // Al2O3
    m[ 4] = e[Cr]/2.0; // Cr2O3
    m[ 6] = e[Mn];     // MnO
    m[ 7] = e[Mg];     // MgO
    m[ 8] = e[Ni];     // NiO
    m[ 9] = e[Co];     // CoO
    m[10] = e[Ca];     // CaO
    m[11] = e[Na]/2.0; // Na2O
    m[12] = e[K]/2.0;  // K2O
    m[13] = e[P]/2.0;  // P2O5
    m[14] = e[H]/2.0;  // H2O

    ox = 2.0*m[0] + 2.0*m[1] + 3.0*m[2] + 3.0*m[4] + m[6] + m[7] + m[8] + m[9] + m[10] + m[11] + m[12] + 5.0*m[13] + m[14];

    // e[O] - ox = 3 m[3] + m[5]  O balance
    // e[Fe]     = 2 m[3] + m[5] Fe balance

    m[ 3] = e[O] - ox - e[Fe];  // Fe2O3
    m[ 5] = e[Fe] - 2.0*m[3];   // FeO

    return mWrapper;
}

@end
