//
//  Equilibrate.h
//  ENKI PhaseObjC package (originally part of the PhasePlot package)
//
//  Originally created by Mark Ghiorso on 5/29/10.
//  Modiied for ENKI on 2/14/17
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"

@class PhaseBase;
@class EquilibrateState;
@class DoubleVector;

@interface Equilibrate : NSObject {
	EquilibrateState *equilibrateState;
	NSMutableDictionary *permissablePhasesState;
	NSDictionary *calculationOptions;
    BOOL initialPhaseAssemblageRequiresEquilibrationStep;
    NSMutableDictionary *pseudoPhaseDictionary;
    NSMutableSet *phasesThatShouldBeIncludedInEvaluateSaturationState;
    BOOL convergenceByMinimalEnergy;
}

@property (strong, readonly) NSArray <NSString *> *systemOxides;
@property (strong, readonly) NSArray <PhaseBase *> *bulkComposition;
@property (strong, readonly) NSArray <NSNumber *> *elementsKnownToSystem;
@property (strong, readwrite) NSDictionary <NSString *, PhaseBase *> *dictionaryOfPhasesInSystem;

@property (assign, readwrite) BOOL debugV;
@property (assign, readwrite) BOOL debugS;

@property (strong, readwrite) NSNumber *PPOPTIONS_DISTRIBUTION;
@property (strong, readwrite) NSNumber *PPOPTIONS_CHEM_POTENTIAL;
@property (strong, readwrite) NSNumber *PPOPTIONS_RANK;
@property (strong, readwrite) NSNumber *PPOPTIONS_ZERO_LINEAR;
@property (strong, readwrite) NSNumber *PPOPTIONS_QUAD_ITERS;
@property (strong, readwrite) NSNumber *PPOPTIONS_FORCE_COMPONENTS;
@property (strong, readwrite) NSNumber *PPOPTIONS_DETECT_MINIMAL;

@property (strong, readwrite) NSNumber *PPPARAMETERS_CHEM_POTENTIAL;
@property (strong, readwrite) NSNumber *PPPARAMETERS_RANK;
@property (strong, readwrite) NSNumber *PPPARAMETERS_QUAD_OPTIMAL;
@property (strong, readwrite) NSNumber *PPPARAMETERS_QUAD_SUBOPTIMAL;
@property (strong, readwrite) NSNumber *PPPARAMETERS_QUAD_MAX_ITERS;
@property (strong, readwrite) NSNumber *PPPARAMETERS_LINEAR_MAX_ITERS;
@property (strong, readwrite) NSNumber *PPPARAMETERS_LINEAR_RELTEST;
@property (strong, readwrite) NSNumber *PPPARAMETERS_LINEAR_MIN_STEPLENGTH;
@property (strong, readwrite) NSNumber *PPPARAMETERS_COMPONENT_MINIMUM;
@property (strong, readwrite) NSNumber *PPPARAMETERS_MASSOUT;
@property (strong, readwrite) NSNumber *PPPARAMETERS_MOLESIN;
@property (strong, readwrite) NSNumber *PPPARAMETERS_LINEAR_MINIMAL;
@property (strong, readwrite) NSNumber *PPPARAMETERS_CYCLES_MAX;

-(NSMutableArray *) oxideNames;
-(NSArray *) phaseNames;
-(DoubleVector *)molesOfOxidesFromMolesOfElements:(double *)e;
-(DoubleVector *)gramsOfOxidesFromMolesOfElements:(double *)e;
-(BOOL)compositionIsFeasible:(NSArray *)compositionInWtPercentOxides forSolution:(id <SolutionPhaseProtocol>)omniComponentPhase;

-(instancetype)init;
-(instancetype)initWithSystemOxides:(NSArray <NSString *> *)systemOxides
      andDictionaryOfPhasesInSystem:(NSDictionary <NSString *, PhaseBase *> *)dictionaryOfPhasesInSystem;
-(instancetype)initWithInformationalDebug:(BOOL)debugS
                          andVerboseDebug:(BOOL)debugV
                          andSystemOxides:(NSArray <NSString *> *)systemOxides
            andDictionaryOfPhasesInSystem:(NSDictionary <NSString *, PhaseBase *> *)dictionaryOfPhasesInSystem;

-(BOOL)evaluateSaturationState:(NSString *)lastPhaseRemovedFromSystem;
-(void)setTemperature:(double)t;
-(void)incrementTemperature:(double)t;
-(void)setPressure:(double)p;
-(void)incrementPressure:(double)p;
-(void)setEntropy:(double)s;
-(void)incrementEntropy:(double)s;
-(void)setVolume:(double)v;
-(void)incrementVolume:(double)v;
-(void)setComposition:(double *)wtOxides;
-(void)setPermissablePhasesState:(NSArray *)arrayOfStateValues;
-(NSDictionary *)getPermissablePhasesState;
-(void)setCalculationOptions:(NSDictionary *)dictionaryOfCalculationOptions;
-(NSDictionary *)execute;

-(double)chemicalPotentialsOfElements:(double *)elementMu;

-(NSString *)equilibrateResultsAsXML;

@end
