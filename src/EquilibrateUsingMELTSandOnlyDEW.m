//
//  EquilibrateUsingMELTSandOnlyDEW.m
//  PhaseMELTSobjC
//
//  Created by Mark Ghiorso on 2/21/18.
//  Copyright © 2018 Mark Ghiorso. All rights reserved.
//

#import "EquilibrateUsingMELTSandOnlyDEW.h"

// This provides the class support to translate chemical formulas into molecular weights and moles of elements
#import "PhaseBase.h"
#import "BermanProperties.h"

// Support the DoubleVector class for method override
#import "DoubleVector.h"

// This is the Equilibrate state that is eventually passed back to the main thread
#import "EquilibrateState.h"

// These are the phases known to MELTS
#import "BermanStoichiometricPhases.h"
#import "LiquidMeltsPlusCO2.h"
#import "DEWFluid.h"

#import "LiquidMeltsSiO2.h"
#import "LiquidMeltsH2ORevised.h"
#import "LiquidMeltsCO2.h"

@implementation EquilibrateUsingMELTSandOnlyDEW

-(EquilibrateStatePhase *)allocatePhaseWithName:(NSString *)name phaseClass:(Class)phaseClass {
    EquilibrateStatePhase *phase = [[EquilibrateStatePhase alloc] init];
    [phase setPhaseClassInstance:[[phaseClass alloc] init]];
    return phase;
}

- (instancetype)init {
    self = [super initWithInformationalDebug:YES
                             andVerboseDebug:YES
                             andSystemOxides:[NSArray arrayWithObjects:@"SiO2", @"TiO2", @"Al2O3", @"Fe2O3",
                                              @"Cr2O3", @"FeO", @"MnO", @"MgO", @"NiO", @"CoO", @"CaO", @"Na2O",
                                              @"K2O", @"P2O5", @"H2O", @"CO2", nil]
               andDictionaryOfPhasesInSystem:[NSDictionary dictionaryWithObjectsAndKeys:
                    [self allocatePhaseWithName:@"Liquid"         phaseClass:[LiquidMeltsPlusCO2 class]  ], @"Liquid",
                    [self allocatePhaseWithName:@"DEWFluid"       phaseClass:[DEWFluid class]            ], @"DEWFluid",
                    nil]];
    if (self) {
        // All instances of subclasses of the BermanProperties class have properties calculated using the
        // Gibbs free energy reference option
        [BermanProperties enableGibbsFreeEnergyReferenceStateUsed];
        EquilibrateStatePhase *phase = (EquilibrateStatePhase *) [self.dictionaryOfPhasesInSystem objectForKey:@"Liquid"];
        LiquidMeltsPlusCO2 *liquid = (LiquidMeltsPlusCO2 *) phase.phaseClassInstance;
        [(LiquidMeltsSiO2 *)      [liquid.endmembers objectAtIndex:0 ] setGibbsFreeEnergyReferenceStateUsed:YES];
        [(LiquidMeltsH2ORevised *)[liquid.endmembers objectAtIndex:14] setGibbsFreeEnergyReferenceStateUsed:YES];
        [(LiquidMeltsCO2 *)       [liquid.endmembers objectAtIndex:15] setGibbsFreeEnergyReferenceStateUsed:YES];
        EquilibrateStatePhase *dewPhase = (EquilibrateStatePhase *) [self.dictionaryOfPhasesInSystem objectForKey:@"DEWFluid"];
        DEWFluid *dew = (DEWFluid *) dewPhase.phaseClassInstance;
        dew.debugS = NO;
        dew.debugV = NO;
    }
    return self;
}

-(DoubleVector *)molesOfOxidesFromMolesOfElements:(double *)e {
    static const int  H =  1;
    static const int  C =  6;
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
    DoubleVector *mWrapper = [[DoubleVector alloc] initWithSize:16];
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
    m[15] = e[C];      // CO2

    ox = 2.0*m[0] + 2.0*m[1] + 3.0*m[2] + 3.0*m[4] + m[6] + m[7] + m[8] + m[9] + m[10] + m[11] + m[12] + 5.0*m[13] + m[14] + 2.0*m[15];

    // e[O] - ox = 3 m[3] + m[5]  O balance
    // e[Fe]     = 2 m[3] + m[5] Fe balance

    m[ 3] = e[O] - ox - e[Fe];  // Fe2O3
    m[ 5] = e[Fe] - 2.0*m[3];   // FeO

    return mWrapper;
}

@end
