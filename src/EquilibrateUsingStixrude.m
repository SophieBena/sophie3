//
//  EquilibrateUsingStixrude.m
//  PhaseMELTSobjC
//
//  Created by Mark Ghiorso on 3/1/17.
//  Copyright © 2017 Mark Ghiorso. All rights reserved.
//

#import "EquilibrateUsingStixrude.h"

// This provides the class support to translate chemical formulas into molecular weights and moles of elements
#import "PhaseBase.h"

// This is the Equilibrate state that is eventually passed back to the main thread
#import "EquilibrateState.h"

// these are phases known to Stixrude
#import "StixrudeEndmembers.h"
#import "StixrudeSolutions.h"
#import "StixrudeSolutionPhase.h"
#import "Stixrude.h"
#import "StixrudeProperties.h"
#import "StixrudeStoichiometricPhases.h"

@implementation EquilibrateUsingStixrude

-(EquilibrateStatePhase *)allocatePhaseWithName:(NSString *)name phaseClass:(Class)phaseClass {
    EquilibrateStatePhase *phase = [[EquilibrateStatePhase alloc] init];
    [phase setPhaseClassInstance:[[phaseClass alloc] init]];
    return phase;
}

- (instancetype)init {
    [[NSUserDefaults standardUserDefaults] setInteger:3 forKey:@"PPCalculationDatabase"];
    self = [super initWithInformationalDebug:NO
                             andVerboseDebug:NO
                             andSystemOxides:[NSArray arrayWithObjects:@"SiO2", @"Al2O3", @"FeO", @"MgO", @"CaO", @"Na2O", nil]
               andDictionaryOfPhasesInSystem:[NSMutableDictionary dictionaryWithObjectsAndKeys:
                    [self allocatePhaseWithName:@"Olivine"          phaseClass:[OlivineStixrude          class]], @"Olivine",
                    [self allocatePhaseWithName:@"Wadsleyite"       phaseClass:[WadsleyiteStixrude       class]], @"Wadsleyite",
                    [self allocatePhaseWithName:@"Ringwoodite"      phaseClass:[RingwooditeStixrude      class]], @"Ringwoodite",
                    [self allocatePhaseWithName:@"Perovskite"       phaseClass:[PerovskiteStixrude       class]], @"Perovskite",
                    [self allocatePhaseWithName:@"Feldspar"         phaseClass:[FeldsparStixrude         class]], @"Feldspar",
                    [self allocatePhaseWithName:@"Spinel"           phaseClass:[SpinelStixrude           class]], @"Spinel",
                    [self allocatePhaseWithName:@"Orthopyroxene"    phaseClass:[OrthopyroxeneStixrude    class]], @"Orthopyroxene",
                    [self allocatePhaseWithName:@"Clinopyroxene"    phaseClass:[ClinopyroxeneStixrude    class]], @"Clinopyroxene",
                    [self allocatePhaseWithName:@"HP-Clinopyroxene" phaseClass:[HPClinopyroxeneStixrude  class]], @"HP-Clinopyroxene",
                    [self allocatePhaseWithName:@"Garnet"           phaseClass:[GarnetStixrude           class]], @"Garnet",
                    [self allocatePhaseWithName:@"Akimotoite"       phaseClass:[AkimotoiteStixrude       class]], @"Akimotoite",
                    [self allocatePhaseWithName:@"Post-Perovskite"  phaseClass:[PostPerovskiteStixrude   class]], @"Post-Perovskite",
                    [self allocatePhaseWithName:@"Ferropericlase"   phaseClass:[FerropericlaseStixrude   class]], @"Ferropericlase",
                    [self allocatePhaseWithName:@"CaFerritePhase"   phaseClass:[CaFerritePhaseStixrude   class]], @"CaFerritePhase",
                    // pure phases
                    [self allocatePhaseWithName:@"Quartz"           phaseClass:[QuartzStixrude           class]], @"Quartz",
                    [self allocatePhaseWithName:@"Coesite"          phaseClass:[CoesiteStixrude          class]], @"Coesite",
                    [self allocatePhaseWithName:@"Stishovite"       phaseClass:[StishoviteStixrude       class]], @"Stishovite",
                    [self allocatePhaseWithName:@"Seifertite"       phaseClass:[SeifertiteStixrude       class]], @"Seifertite",
                    [self allocatePhaseWithName:@"CaPerovskite"     phaseClass:[CaPerovskiteStixrude     class]], @"CaPerovskite",
                    [self allocatePhaseWithName:@"Kyanite"          phaseClass:[KyaniteStixrude          class]], @"Kyanite",
                    [self allocatePhaseWithName:@"Nepheline"        phaseClass:[NephelineStixrude        class]], @"Nepheline",
                    nil]];
    if (self) {
    }
    return self;
}

@end
