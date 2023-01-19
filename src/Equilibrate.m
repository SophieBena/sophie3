//
//  Equilibrate.m
//  ENKI PhaseObjC package
//
//  Created by Mark Ghiorso on 5/29/10.
//  Modified by Mark S. Ghiorso on 2/14/17.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "Equilibrate.h"

// These are support functions
#import "DoubleVector.h"
#import "DoubleMatrix.h"
#import "IntegerVector.h"
#import "MathSupport.h"

// This is a prototype to a fictious pseudo-liquid like phase
#import "PseudoPhase.h"

// This provides the class support to translate chemical formulas into molecular weights and moles of elements
#import "PhaseBase.h"

// This is the Equilibrate state that is eventually passed back to the main thread
#import "EquilibrateState.h"

// The accelerate framework (BLAS and LAPACK)
#ifdef __APPLE__
#import <Accelerate/Accelerate.h>
#else
#include "clapack.h"
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#pragma mark -
#pragma mark Equilibrate class

@implementation Equilibrate

#define KEY_ORDINATE      @"ordinate"
#define KEY_ABSCISSA      @"abscissa"
#define KEY_IMPOSE_BUFFER @"imposeBuffer"
#define KEY_BUFFER        @"buffer"
#define KEY_BUFFER_OFFSET @"bufferOffset"

#pragma mark -
#pragma mark public instance methods

-(instancetype)init {
    return [self initWithInformationalDebug:NO
                            andVerboseDebug:NO
                            andSystemOxides:nil
              andDictionaryOfPhasesInSystem:nil];
}

-(instancetype)initWithSystemOxides:(NSArray <NSString *> *)systemOxides
      andDictionaryOfPhasesInSystem:(NSDictionary <NSString *, PhaseBase *> *)dictionaryOfPhasesInSystem {
    return [self initWithInformationalDebug:NO
                            andVerboseDebug:NO
                            andSystemOxides:systemOxides
              andDictionaryOfPhasesInSystem:dictionaryOfPhasesInSystem];
}

-(instancetype)initWithInformationalDebug:(BOOL)debugS
                          andVerboseDebug:(BOOL)debugV
                          andSystemOxides:(NSArray <NSString *> *)systemOxides
            andDictionaryOfPhasesInSystem:(NSDictionary <NSString *, PhaseBase *> *)dictionaryOfPhasesInSystem {
	if ((self = [super init])) {
        _debugS = debugS;
        _debugV = debugV;
        [[NSUserDefaults standardUserDefaults] setBool:NO forKey:@"DEBUG.SATLOOP.VERBOSE"];
        [[NSUserDefaults standardUserDefaults] setBool:NO forKey:@"DEBUG.SATLOOP.SIMPLE"];
		if (_debugV) NSLog(@"[Equilibrate init] Entering method ...");
        _PPOPTIONS_DISTRIBUTION     = @NO;
        _PPOPTIONS_CHEM_POTENTIAL   = @NO;
        _PPOPTIONS_RANK             = @NO;
        _PPOPTIONS_ZERO_LINEAR      = @NO;
        _PPOPTIONS_QUAD_ITERS       = @NO;
        _PPOPTIONS_FORCE_COMPONENTS = @NO;
        _PPOPTIONS_DETECT_MINIMAL   = @NO;

        _PPPARAMETERS_CHEM_POTENTIAL        = [NSNumber numberWithDouble:pow(DBL_EPSILON, 0.75)];
        _PPPARAMETERS_RANK                  = [NSNumber numberWithDouble:10.0*DBL_EPSILON];
        _PPPARAMETERS_QUAD_OPTIMAL          = [NSNumber numberWithDouble:pow(DBL_EPSILON, 0.75)];
        _PPPARAMETERS_QUAD_SUBOPTIMAL       = [NSNumber numberWithDouble:sqrt(DBL_EPSILON)];
        _PPPARAMETERS_QUAD_MAX_ITERS        = [NSNumber numberWithUnsignedInteger:100];
        _PPPARAMETERS_LINEAR_MAX_ITERS      = [NSNumber numberWithUnsignedInteger:100];
        _PPPARAMETERS_LINEAR_RELTEST        = [NSNumber numberWithDouble:10.0*DBL_EPSILON];
        _PPPARAMETERS_LINEAR_MIN_STEPLENGTH = [NSNumber numberWithDouble:DBL_EPSILON];
        _PPPARAMETERS_COMPONENT_MINIMUM     = [NSNumber numberWithDouble:DBL_EPSILON];
        _PPPARAMETERS_MASSOUT               = [NSNumber numberWithDouble:1.0e-5];
        _PPPARAMETERS_MOLESIN               = [NSNumber numberWithDouble:1.0e-6];
        _PPPARAMETERS_LINEAR_MINIMAL        = [NSNumber numberWithUnsignedInteger:5];
        _PPPARAMETERS_CYCLES_MAX            = [NSNumber numberWithUnsignedInteger:3];

        _systemOxides = [NSArray arrayWithArray:systemOxides];
        NSMutableArray *mutableBulkComposition = [NSMutableArray arrayWithCapacity:[systemOxides count]];
        BOOL foundElement[107]; for (NSUInteger i=0; i<107; i++) foundElement[i] = NO;

        for (NSString *entry in systemOxides) {
            PhaseBase *oxide = [[PhaseBase alloc] init];
            [oxide setPhaseName:entry];
            [oxide setPhaseFormula:entry];
            double *elements = [[oxide formulaAsElementArray] pointerToDouble];
            for (NSUInteger i=1; i<107; i++) if (elements[i] != 0.0) foundElement[i] = YES;
            [mutableBulkComposition addObject:oxide];
        }
        _bulkComposition = [NSArray arrayWithArray:mutableBulkComposition];
        NSMutableArray *mutableElementsKnownToSystem = [NSMutableArray arrayWithCapacity:20];
        for (NSUInteger i=1; i<107; i++) if (foundElement[i]) [mutableElementsKnownToSystem addObject:[NSNumber numberWithUnsignedInteger:i]];
        _elementsKnownToSystem = [NSArray arrayWithArray:mutableElementsKnownToSystem];

        if (_debugV) NSLog(@"... creating instance of EquilibrateState class.");
        equilibrateState = [[EquilibrateState alloc] initWithComponentNumber:[_systemOxides count]];

        if (_debugV) NSLog(@"... allocating instances of potential phases.");
        [equilibrateState setPotentialPhaseList: [NSMutableDictionary dictionaryWithDictionary:dictionaryOfPhasesInSystem]];
        _dictionaryOfPhasesInSystem = [NSDictionary dictionaryWithDictionary:dictionaryOfPhasesInSystem];

        permissablePhasesState = [NSMutableDictionary dictionaryWithCapacity:1];
        for (NSString *key in dictionaryOfPhasesInSystem.allKeys) {
            [permissablePhasesState setObject:[NSNumber numberWithBool:YES] forKey:key];
        }

		calculationOptions = [NSDictionary dictionaryWithObjectsAndKeys:@"Temperature (°C)",KEY_ORDINATE,
							  @"Pressure (MPa)", KEY_ABSCISSA,
							  [NSNumber numberWithBool:NO], KEY_IMPOSE_BUFFER,
							  @"QFM", KEY_BUFFER,
							  [NSNumber numberWithDouble:0.0], KEY_BUFFER_OFFSET,
							  nil];
        initialPhaseAssemblageRequiresEquilibrationStep = NO;
        pseudoPhaseDictionary = [NSMutableDictionary dictionaryWithCapacity:0];
        phasesThatShouldBeIncludedInEvaluateSaturationState = [NSMutableSet setWithCapacity:0];
        convergenceByMinimalEnergy = NO;
		if (_debugV) NSLog(@"... exiting [Equilibrate init].");
    }
    return self;
}

/**
 Provides an array of NSString objects that correspond to the names of oxide composition variables.

 @return
 An array of Cocoa strings containing names of known system oxides.
 */
-(NSArray *) oxideNames {
    NSArray *array = [NSArray arrayWithArray:self.systemOxides];
    return array;
}

/**
 Provides a test on the feasibility of an input composition

 @param
 compositionInWtPercentOxides array of NSNumber objects containing wt% values of system oxides

 @return
 YES = feasible, NO = unfeasible
 */
-(BOOL)compositionIsFeasible:(NSArray *)compositionInWtPercentOxides forSolution:(id <SolutionPhaseProtocol>)omniComponentPhase {
    BOOL result = YES;

    NSUInteger nc = [compositionInWtPercentOxides count];
    DoubleVector *molesOxidesWrapper = [[DoubleVector alloc] initWithSize:nc];
    double *molesOxides = [molesOxidesWrapper pointerToDouble];
    for (NSUInteger i=0; i<nc; i++) {
        molesOxides[i] = [[compositionInWtPercentOxides objectAtIndex:i] doubleValue]/[[self.bulkComposition objectAtIndex:i] mw];
        if (molesOxides[i] < 0.0) result = NO;
    }

    if (result) {
        DoubleVector *molesElementsWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
        double *molesElements = [molesElementsWrapper pointerToDouble];
        for (NSUInteger i=0; i<nc; i++) {
            [EquilibrateState accumulateIntoVectorOfElements:molesElements
                                                     aScaler:molesOxides[i]
                                 timesElementReferenceVector:[[[self.bulkComposition objectAtIndex:i] formulaAsElementArray] pointerToDouble]];
        }

        result = [omniComponentPhase testPermissibleValuesOfComponents:
                  [[omniComponentPhase convertElementsToMoles:molesElements] pointerToDouble]];
    }

    return result;
}

/**
 Provides an array of NSString objects that correspond to the names of phases known to this class.

 @return
 An array of Cocoa strings.
 */
-(NSArray *) phaseNames {
    return [[self.dictionaryOfPhasesInSystem allKeys] sortedArrayUsingSelector:@selector(localizedCaseInsensitiveCompare:)];
}

/**
 Converts moles of elements into moles of oxides

 This is a generic implementation that should be overridden by a subclass.  This version assumes that there are no multiple
 oxidation states of a given element.  Oxygen balance is ignored, and multiple instances, like FeO and Fe2O3, are
 assigned values without regard to mass and oxygen mole balance

 @param
 e vector of length 107 containing moles of elements

 @return
 Double vector of length [self.systemOxides count] containing moles of oxides
 */
-(DoubleVector *)molesOfOxidesFromMolesOfElements:(double *)e {
    DoubleVector *mWrapper = [[DoubleVector alloc] initWithSize:[self.systemOxides count]];
    double *m = [mWrapper pointerToDouble];

    NSUInteger nOx = 0;
    for (PhaseBase *oxide in self.bulkComposition) {
        DoubleVector *formulaAsElementArrayWrapper = oxide.formulaAsElementArray;
        double *formulaAsElementArray = [formulaAsElementArrayWrapper pointerToDouble];
        for (NSUInteger i=1; i<107; i++) {
            if ( (i != 8) && (formulaAsElementArray[i] > 0.0) && (e[i] > 0.0)) {
                m[nOx] = e[i]/formulaAsElementArray[i];
                break;
            }
        }
        nOx++;
    }
    return mWrapper;
}

/**
 Converts moles of oxides into moles of elements

 @param
 e vector of length [self.systemOxides count] containing moles of oxides

 @return
 Double vector of length 107 containing moles of elements
 */
-(DoubleVector *)gramsOfOxidesFromMolesOfElements:(double *)e {
    DoubleVector *mWrapper = [self molesOfOxidesFromMolesOfElements:e];
    double *m = [mWrapper pointerToDouble];
    DoubleVector *gWrapper = [[DoubleVector alloc] initWithSize:[self.systemOxides count]];
    double *g = [gWrapper pointerToDouble];
    for (NSUInteger i=0; i<[self.systemOxides count]; i++)
        g[i] = m[i]*[(PhaseBase *)[self.bulkComposition objectAtIndex:i] mw];
    return gWrapper;
}

-(NSArray *)equilibrateResultsHelper:(NSArray *)molesOfElements {
    double moles[107];
    moles[0] = 0.0;
    for (NSUInteger i=0; i<106; i++) moles[i+1] = [[molesOfElements objectAtIndex:i] doubleValue];
    DoubleVector *gramsOfOxidesWrapper = [self gramsOfOxidesFromMolesOfElements:moles];
    double *gramsOfOxides = [gramsOfOxidesWrapper pointerToDouble];
    NSMutableArray *result = [NSMutableArray arrayWithCapacity:0];

    NSArray *oxideNames = [self oxideNames];
    NSUInteger n = [oxideNames count];
    for (NSUInteger i=0; i<n; i++) [result addObject:[NSArray arrayWithObjects:[oxideNames objectAtIndex:i], [NSNumber numberWithDouble:gramsOfOxides[i]], nil]];

    return result;
}

-(NSString *)equilibrateResultsAsXML {
    NSXMLElement *root = (NSXMLElement *)[NSXMLNode elementWithName:@"EquilibrateState"];
    NSXMLDocument *xmlDoc = [[NSXMLDocument alloc] initWithRootElement:root];
    [xmlDoc setVersion:@"1.0"];
    [xmlDoc setCharacterEncoding:@"UTF-8"];
    [root addChild:[NSXMLNode elementWithName:@"Date" stringValue:[[NSDate date] descriptionWithLocale:nil]]];
    [root addChild:[equilibrateState getEquilibrateStateAsXMLElement:self selector:@selector(equilibrateResultsHelper:)]];
    return [[NSString alloc] initWithData:[xmlDoc XMLData] encoding:NSUTF8StringEncoding];
}

/**
 Sets the temperature

 @param
   t in K
 */
-(void)setTemperature:(double)t {
	[equilibrateState setT:t];
    [equilibrateState setReferenceEntropyOfSystem:0.0];
    [equilibrateState setCorrectionToReferenceEntropyOfSystem:0.0];
}

/**
 Increments the temperature (assumes a previously successful execute step)

 @param
 t in K
 */
-(void)incrementTemperature:(double)t {
    [self setTemperature:t];
    initialPhaseAssemblageRequiresEquilibrationStep = YES;
}

/**
 Sets the pressure

 @param
   p in bars
 */
-(void)setPressure:(double)p {
	[equilibrateState setP:p];
    [equilibrateState setReferenceVolumeOfSystem:0.0];
    [equilibrateState setCorrectionToReferenceVolumeOfSystem:0.0];
}

/**
 Increments the pressure (assumes a previously successful execute step)

 @param
 p in bars
 */
-(void)incrementPressure:(double)p {
    [self setPressure:p];
    initialPhaseAssemblageRequiresEquilibrationStep = YES;
}

/**
 Sets the entropy

 @param
 s in J/K
 */
-(void)setEntropy:(double)s {
    //[equilibrateState setT:0.0];
    [equilibrateState setReferenceEntropyOfSystem:s];
    [equilibrateState setCorrectionToReferenceEntropyOfSystem:0.0];
}

/**
 Increments the entropy (assumes a previously successful execute step)

 @param
 s in J/K
 */
-(void)incrementEntropy:(double)s {
    [self setEntropy:s];
    initialPhaseAssemblageRequiresEquilibrationStep = YES;
}

/**
 Sets the volume

 @param
 v in J/bar
 */
-(void)setVolume:(double)v {
    //[equilibrateState setP:0.0];
    [equilibrateState setReferenceVolumeOfSystem:v];
    [equilibrateState setCorrectionToReferenceVolumeOfSystem:0.0];
}

/**
 Increments the volume (assumes a previously successful execute step)

 @param
 v in J/bar
 */
-(void)incrementVolume:(double)v {
    [self setVolume:v];
    initialPhaseAssemblageRequiresEquilibrationStep = YES;
}

enum modeForEstimateLiquidAbsentEquilibriumAssemblage {
    kDoNotForceInclusionOfAnyEndmemberOfAnyPhaseInSystem,
    kForceInclusionOfSolidSolutionEndmembersForViablePhasesInSystem,
    kForceInclusionOfEveryViablePhaseInSystem
};

-(void)estimateLiquidAbsentEquilibriumAssemblage:(NSUInteger)mode {
    NSUInteger numberOfElementsKnownToSystem;
    NSMutableArray *constraints;
    NSSet *refElementSet;
    numberOfElementsKnownToSystem = [self.elementsKnownToSystem count];
    constraints = [NSMutableArray arrayWithCapacity:numberOfElementsKnownToSystem];
    refElementSet = [NSSet setWithArray:self.elementsKnownToSystem];

    NSUInteger numberOfInequalityConstraints = 0;
	NSUInteger numberOfEqualityConstraints = 0;

    for (EquilibrateStatePhase *phaseWrapper in [[equilibrateState potentialPhaseList] allValues]) {
        id phase = [phaseWrapper phaseClassInstance];

        BOOL includePhase = YES;
        if (permissablePhasesState) includePhase = [[permissablePhasesState valueForKey:[phase phaseName]] boolValue];

        if (includePhase) {
            if (![phase conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
                numberOfEqualityConstraints++;
                [constraints addObject:[NSArray arrayWithObjects:[phase phaseName],
                                        [NSNumber numberWithUnsignedInteger:1],
                                        @"",
                                        [NSValue valueWithPointer:[[phase formulaAsElementArray] pointerToDouble]],
                                        [NSNumber numberWithDouble:[phase getGibbsFreeEnergyFromT:[equilibrateState T] andP:[equilibrateState P]]], nil]];
            } else {
                NSUInteger na = [phase numberOfSolutionComponents];
                NSMutableArray *phaseEntryForConstraintArray = [NSMutableArray arrayWithCapacity:2+4*na];
                [phaseEntryForConstraintArray addObject:[phase phaseName]];
                [phaseEntryForConstraintArray addObject:[NSNumber numberWithUnsignedInteger:na]];

                for (NSUInteger i=0; i<na; i++) {
                    [phaseEntryForConstraintArray addObject:[phase nameOfSolutionSpeciesAtIndex:i]];
                    [phaseEntryForConstraintArray addObject:[NSValue valueWithPointer:[[phase elementalCompositionOfSpeciesAtIndex:i] pointerToDouble]]];
                    [phaseEntryForConstraintArray addObject:[NSNumber numberWithDouble:[[phase componentAtIndex:i] getGibbsFreeEnergyFromT:[equilibrateState T] andP:[equilibrateState P]]]];

                    NSMutableSet *elementsInFormulaSet = [NSMutableSet setWithCapacity:numberOfElementsKnownToSystem];
                    double *elementArray = [[phase elementalCompositionOfSpeciesAtIndex:i] pointerToDouble];
                    for (NSUInteger j=1; j<107; j++) if (elementArray[j] != 0.0) [elementsInFormulaSet addObject:[NSNumber numberWithUnsignedInteger:j]];
                    BOOL activeInequalityConstraint = [elementsInFormulaSet isSubsetOfSet:refElementSet];
                    [phaseEntryForConstraintArray addObject:[NSNumber numberWithBool:activeInequalityConstraint]];
                    if (activeInequalityConstraint && !(mode == kDoNotForceInclusionOfAnyEndmemberOfAnyPhaseInSystem)) numberOfInequalityConstraints++;
                }

                numberOfEqualityConstraints += na;
                [constraints addObject:phaseEntryForConstraintArray];
            }
        }
    }

    DoubleMatrix *aSMPXWrapper = [[DoubleMatrix alloc]initWithRowSize:numberOfInequalityConstraints+numberOfElementsKnownToSystem+3
                                                    andWithColumnSize:numberOfEqualityConstraints+2
                                                      andInitialValue:0.0];
	double **aSMPX = [aSMPXWrapper pointerToPointerToDouble];

    // Inequality constraints
	NSUInteger colIndex = 2;
	NSUInteger rowIndex = 2;

    if (!(mode == kDoNotForceInclusionOfAnyEndmemberOfAnyPhaseInSystem)) {
        for (NSArray *object in constraints) {
            NSUInteger na = [[object objectAtIndex:1] unsignedIntValue];
            if (self.debugV) NSLog(@"Inequality %@ na:%lu colIndex = %lu rowIndex = %lu", [object objectAtIndex:0], na, colIndex, rowIndex);
            if (na > 1) {
                NSUInteger ii=0;
                for (NSUInteger i=0; i<na; i++) {
                    if ([[object objectAtIndex:2+4*i+3] boolValue]) {
                        double molesIn = self.PPPARAMETERS_MOLESIN.doubleValue;
                        NSUInteger jj=0;
                        for (NSUInteger j=0; j<na; j++) {
                            if ([[object objectAtIndex:2+4*j+3] boolValue]) {
                                if (mode == kForceInclusionOfEveryViablePhaseInSystem) {
                                    if (ii == jj) {
                                        aSMPX[rowIndex+ii][colIndex+j] = -1.0;
                                        aSMPX[rowIndex+ii][1] = molesIn;
                                        if (self.debugV) NSLog(@"... a(rowIndex+%lu, colIndex+%lu) = %lf", ii, j, -1.0);
                                        if (self.debugV) NSLog(@"... a(rowIndex+%lu, 1)            = %lf", ii, molesIn);
                                    }
                                } else {
                                    double factor = -molesIn;
                                    if (ii == jj) factor = 1.0 - molesIn;
                                    aSMPX[rowIndex+ii][colIndex+j] = -factor;
                                    if (self.debugV) NSLog(@"... a(rowIndex+%lu, colIndex+%lu) = %lf", ii, j, -factor);
                                }
                                jj++;
                            }
                        }
                        ii++;
                    }
                }
                rowIndex += ii;
            }
            colIndex += na;
            if (self.debugV) NSLog(@"... exit colIndex = %lu rowIndex = %lu", colIndex, rowIndex);
        }
        NSAssert((rowIndex-2) == numberOfInequalityConstraints, @"Improper construction of inequality constraints (rows) for SIMPLX tableau.");
        NSAssert((colIndex-2) == numberOfEqualityConstraints,   @"Improper construction of inequality constraints (columns) for SIMPLX tableau.");
        colIndex = 2;
    }

    // objective function and equality constraints
	for (NSArray *object in constraints) {
		NSUInteger na = [[object objectAtIndex:1] unsignedIntValue];
        if (self.debugV) NSLog(@"Equality %@ na:%lu colIndex = %lu rowIndex = %lu", [object objectAtIndex:0], na, colIndex, rowIndex);
		for (NSUInteger i=0; i<na; i++) {
            NSUInteger ii = rowIndex;
            double *elementArray = [[object objectAtIndex:2+4*i+1] pointerValue];
            for (NSNumber *entry in self.elementsKnownToSystem) {
                NSUInteger elmNo = [entry unsignedIntValue];
                aSMPX[ii][colIndex] = -elementArray[elmNo];
                if (self.debugV) NSLog(@"... a(%lu, %lu) = %lf", ii, colIndex, -elementArray[elmNo]);
                ii++;
            }
            aSMPX[1][colIndex] = -[[object objectAtIndex:2+4*i+2] doubleValue]; // objective functionnegative because function is maximized
            colIndex++;
        }
        if (self.debugV) NSLog(@"... exit colIndex = %lu rowIndex = %lu", colIndex, rowIndex);
	}
    NSAssert((colIndex-2) == numberOfEqualityConstraints,   @"Improper construction of equality constraints (columns) for SIMPLX tableau.");

    // rhs vector, objective function, inequality constraints, equality constraints
	aSMPX[1][1] = 0.0;
	double *bulkCompositionArray = [[equilibrateState bulkCompositionInElements] pointerToDouble];
    for (NSNumber *entry in self.elementsKnownToSystem) {
        NSUInteger elmNo = [entry unsignedIntValue];
        aSMPX[rowIndex++][1] = bulkCompositionArray[elmNo];
    }
    NSAssert((rowIndex-2) == (numberOfElementsKnownToSystem+numberOfInequalityConstraints), @"Improper number of rows for SIMPLX tableau");

	NSInteger icaseSMPX;
    IntegerVector *izrovSMPXWrapper = [[IntegerVector alloc] initWithSize:numberOfEqualityConstraints+1 andInitialValue:0];
    IntegerVector *iposvSMPXWrapper = [[IntegerVector alloc] initWithSize:numberOfElementsKnownToSystem+numberOfInequalityConstraints+1 andInitialValue:0];
	NSInteger *izrovSMPX = [izrovSMPXWrapper pointerToInteger];
	NSInteger *iposvSMPX = [iposvSMPXWrapper pointerToInteger];

	[MathSupport simplx:aSMPX
					  m:numberOfElementsKnownToSystem+numberOfInequalityConstraints
					  n:numberOfEqualityConstraints
					 m1:0
					 m2:numberOfInequalityConstraints
					 m3:numberOfElementsKnownToSystem
				  icase:&icaseSMPX
				  izrov:izrovSMPX
				  iposv:iposvSMPX];

    DoubleVector *xSMPXWrapper = [[DoubleVector alloc] initWithSize:numberOfEqualityConstraints];
	double *xSMPX = [xSMPXWrapper pointerToDouble];
	if (self.debugS) NSLog(@"...return from simplex method with icase = %ld. Objective function = %20.13g. Solution:", icaseSMPX, aSMPX[1][1]);
	for (NSUInteger j=1, i; j<=numberOfElementsKnownToSystem+numberOfInequalityConstraints; j++) if ((i = iposvSMPX[j]) <= numberOfEqualityConstraints) xSMPX[i-1] = aSMPX[j+1][1];
	for (NSUInteger j=1, i; j<=numberOfEqualityConstraints; j++) if ((i = izrovSMPX[j]) <= numberOfEqualityConstraints) xSMPX[i-1] = 0.0;

	colIndex = 0;
    DoubleVector *totalElementsInAllPhasesWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
	double *totalElementsInAllPhases = [totalElementsInAllPhasesWrapper pointerToDouble];
	for (NSArray *object in constraints) {
		NSString *key = [object objectAtIndex:0];
		NSUInteger na = [[object objectAtIndex:1] unsignedIntValue];
        DoubleVector *elementsInPhaseWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
		double *elementsInPhase = [elementsInPhaseWrapper pointerToDouble];

		for (NSUInteger i=0; i<na; i++) {
			if (xSMPX[colIndex] > sqrt(DBL_EPSILON)) {
				if (self.debugS) NSLog(@"...%@[%@] = %g", [[object objectAtIndex:0] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
					  [[object objectAtIndex:2+4*i] stringByPaddingToLength:20 withString:@" " startingAtIndex:0], xSMPX[colIndex]);
				[EquilibrateState accumulateIntoVectorOfElements:elementsInPhase aScaler:xSMPX[colIndex] timesElementReferenceVector:[[object objectAtIndex:2+4*i+1] pointerValue]];
			}
			colIndex++;
		}

		double sumOfElementsInPhase = 0.0;
		for (NSUInteger i=1; i<107; i++) sumOfElementsInPhase += elementsInPhase[i];

		if (sumOfElementsInPhase > sqrt(DBL_EPSILON)) {
			EquilibrateStatePhase *phaseWrapper = [[equilibrateState potentialPhaseList] objectForKey:key];
			[[equilibrateState potentialPhaseList] removeObjectForKey:key];

			double *phaseCompositionInElements  = [[phaseWrapper bulkCompositionInElements] pointerToDouble];
			for (NSUInteger i=1; i<107; i++) phaseCompositionInElements[i] = elementsInPhase[i];
			// [EquilibrateState accumulateIntoVectorOfElements:phaseCompositionInElements aScaler:1.0 timesElementReferenceVector:elementsInPhase];
			[EquilibrateState accumulateIntoVectorOfElements:totalElementsInAllPhases aScaler:1.0 timesElementReferenceVector:elementsInPhase];

			double massOfPhase = 0.0;
			if ([[phaseWrapper phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
				massOfPhase = [[phaseWrapper phaseClassInstance] convertElementsToTotalMass:phaseCompositionInElements];
				double *m = [[[phaseWrapper phaseClassInstance] convertElementsToMoles:phaseCompositionInElements] pointerToDouble];
				if (self.debugS) NSLog(@"...Formula: %@", [[phaseWrapper phaseClassInstance] getFormulaFromMolesOfComponents:m andT:[equilibrateState T] andP:[equilibrateState P]]);
			} else {
				massOfPhase = [[phaseWrapper phaseClassInstance] convertElementsToMassOfPhase:phaseCompositionInElements];
			}
			[phaseWrapper setMass:massOfPhase];

			[[equilibrateState phasesInSystem] setObject:phaseWrapper forKey:key];
			if ([[phaseWrapper phaseClassInstance] respondsToSelector:@selector(incrementInstanceCountOfPhase)]) [[phaseWrapper phaseClassInstance] incrementInstanceCountOfPhase];
		}
	}

	if (self.debugS) {
		double *totalElementsInSystem = [[equilibrateState bulkCompositionInElements] pointerToDouble];
		for (NSUInteger i=1; i<107; i++) if (totalElementsInSystem[i] != 0.0) {
			NSLog(@"...element[%2@]: System: %13.6g  NNLS: %13.6g  delta: %13.6g",
				  [[PhaseBase elementNameFromAtomicNumber:i] stringByPaddingToLength:2 withString:@" " startingAtIndex:0], totalElementsInSystem[i], totalElementsInAllPhases[i],
				  totalElementsInSystem[i] - totalElementsInAllPhases[i]);
		}
	}

    // recast computed assemblage into chemical potentials of elements in the system
    DoubleVector *elementMuWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
    double *elementMu = [elementMuWrapper pointerToDouble];
    double residualNorm = [self chemicalPotentialsOfElements:elementMu];
    if (self.debugS) NSLog(@"Residual norm of recasted system composition: %lf", residualNorm);

    // remove all phases from the system and add them back to the potential phase list
    [phasesThatShouldBeIncludedInEvaluateSaturationState removeAllObjects];
    for (NSString *key in [[equilibrateState phasesInSystem] allKeys]) {
        EquilibrateStatePhase *phaseWrapper = [[equilibrateState phasesInSystem] objectForKey:key];
        [EquilibrateState scaleVectorOfElements:[[phaseWrapper bulkCompositionInElements] pointerToDouble] aScaler:0.0];
        [phaseWrapper setMass:0.0];
        [[equilibrateState potentialPhaseList] setObject:phaseWrapper forKey:key];
        [phasesThatShouldBeIncludedInEvaluateSaturationState addObject:key];
    }
    [[equilibrateState phasesInSystem] removeAllObjects];

    // Now make element pseudophases with the chemical potentials calculated above
    for (NSNumber *entry in self.elementsKnownToSystem) {
        NSUInteger elmNo = [entry unsignedIntValue];
        if (totalElementsInAllPhases[elmNo] > 0.0) {
            NSString *name = [PhaseBase elementNameFromAtomicNumber:elmNo];

            PseudoPhase *pseudoPhase = [[PseudoPhase alloc] init];
            [pseudoPhase setPhaseName:[name stringByAppendingString:@"PseudoPhase"]];
            [pseudoPhase setPhaseFormula:name];
            [pseudoPhase setG:elementMu[elmNo]];
            [pseudoPhase setGLast:elementMu[elmNo]-5000.0];

            EquilibrateStatePhase *phaseWrapper = [[EquilibrateStatePhase alloc] init];
            double *bc = [[phaseWrapper bulkCompositionInElements] pointerToDouble];
            bc[elmNo] = totalElementsInAllPhases[elmNo];
            [phaseWrapper setPhaseClassInstance:pseudoPhase];
            [phaseWrapper setMass:totalElementsInAllPhases[elmNo]*[pseudoPhase mw]];

            [[equilibrateState phasesInSystem] setObject:phaseWrapper forKey:[name stringByAppendingString:@"PseudoPhase"]];
            [pseudoPhaseDictionary setObject:pseudoPhase forKey:name];
        }
    }

    initialPhaseAssemblageRequiresEquilibrationStep = NO;
}

-(void)estimateLiquidAbsentEquilibriumAssemblageUsingLDPandMinimumConcentration:(double)minimumConcentration {
	NSMutableArray *constraints = [NSMutableArray arrayWithCapacity:[[equilibrateState potentialPhaseList] count]];
    BOOL elementFoundInPotentialPhaseList[107]; for (NSUInteger i=0; i<107; i++) elementFoundInPotentialPhaseList[i] = NO;

    NSUInteger numberOfPhaseComponentEndmembers = 0;
    for (EquilibrateStatePhase *phaseWrapper in [[equilibrateState potentialPhaseList] allValues]) {
        id phase = [phaseWrapper phaseClassInstance];

        BOOL includePhase = YES;
        if (permissablePhasesState) includePhase = [[permissablePhasesState valueForKey:[phase phaseName]] boolValue];

        if (includePhase) {
            if (![phase conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
                numberOfPhaseComponentEndmembers++;
                [constraints addObject:[NSArray arrayWithObjects:[phase phaseName],
                                        [NSNumber numberWithUnsignedInteger:1],
                                        @"",
                                        [NSValue valueWithPointer:[[phase formulaAsElementArray] pointerToDouble]], nil]];
                double *formulaOfPhase = [[phase formulaAsElementArray] pointerToDouble];
                for (NSUInteger i=1; i<107; i++) if (formulaOfPhase[i] != 0.0) elementFoundInPotentialPhaseList[i] = YES;
            } else {
                NSUInteger na = [phase numberOfSolutionComponents];
                NSMutableArray *phaseEntryForConstraintArray = [NSMutableArray arrayWithCapacity:2+2*na];
                [phaseEntryForConstraintArray addObject:[phase phaseName]];
                [phaseEntryForConstraintArray addObject:[NSNumber numberWithUnsignedInteger:na]];

                for (NSUInteger i=0; i<na; i++) {
                    [phaseEntryForConstraintArray addObject:[phase nameOfSolutionSpeciesAtIndex:i]];
                    [phaseEntryForConstraintArray addObject:[NSValue valueWithPointer:[[phase elementalCompositionOfSpeciesAtIndex:i] pointerToDouble]]];
                    double *formulaOfPhase = [[phase elementalCompositionOfSpeciesAtIndex:i] pointerToDouble];
                    for (NSUInteger i=1; i<107; i++) if (formulaOfPhase[i] != 0.0) elementFoundInPotentialPhaseList[i] = YES;
                }

                numberOfPhaseComponentEndmembers += na;
                [constraints addObject:phaseEntryForConstraintArray];
            }
        }
    }

	NSUInteger numberOfElementsKnownToSystem = 0;
    for (NSUInteger i=1; i<107; i++) if (elementFoundInPotentialPhaseList[i]) numberOfElementsKnownToSystem++;

    NSUInteger numberOfInequalityConstraints = numberOfPhaseComponentEndmembers + 2*numberOfElementsKnownToSystem;
    DoubleMatrix *gLDPWrapper = [[DoubleMatrix alloc] initWithRowSize:numberOfInequalityConstraints
                                                    andWithColumnSize:numberOfPhaseComponentEndmembers
                                                      andInitialValue:0.0];
	double **gLDP = [gLDPWrapper pointerToPointerToDouble];
	DoubleVector *hLDPWrapper = [[DoubleVector alloc] initWithSize:numberOfInequalityConstraints andInitialValue:0.0];
    double *hLDP = [hLDPWrapper pointerToDouble];

    for (NSUInteger i=0; i<numberOfPhaseComponentEndmembers; i++) {
        gLDP[i][i] = 1.0;
        hLDP[i]    = minimumConcentration;
    }

	NSUInteger rowIndex = numberOfPhaseComponentEndmembers;

    double *bulkCompositionArray = [[equilibrateState bulkCompositionInElements] pointerToDouble];
	for (NSUInteger elmNo=1; elmNo<107; elmNo++) if (elementFoundInPotentialPhaseList[elmNo]) {
        NSUInteger colIndex = 0;

        for (NSArray *object in constraints) {
            NSUInteger na = [[object objectAtIndex:1] unsignedIntValue];
            for (NSUInteger i=0; i<na; i++) {
                double *bc = [[object objectAtIndex:3+2*i] pointerValue];
                gLDP[rowIndex  ][colIndex] =  bc[elmNo];
                gLDP[rowIndex+1][colIndex] = -bc[elmNo];
                colIndex++;
            }
        }
        NSAssert(colIndex == numberOfPhaseComponentEndmembers, @"Improper construction of inequality constraints (columns) for LDP matrix.");

		hLDP[rowIndex  ] =  bulkCompositionArray[elmNo]*(1.0-sqrt(DBL_EPSILON));
        hLDP[rowIndex+1] = -bulkCompositionArray[elmNo]*(1.0+sqrt(DBL_EPSILON));
        NSLog(@"element %lu moles = %20.13e", elmNo, bulkCompositionArray[elmNo]);

        rowIndex += 2;
	}
    NSAssert(rowIndex == numberOfInequalityConstraints, @"Improper construction of inequality constraints (rows) for LDP tableau.");

    if (self.debugS) {
        for (NSUInteger i=0; i<numberOfPhaseComponentEndmembers; i++) {
            NSMutableString *temp = [NSMutableString stringWithCapacity:0];
            for (NSUInteger j=0; j<20; j++) if (gLDP[i][j] != 0.0) [temp appendFormat:@" %5.1lf", gLDP[i][j]]; else [temp appendString:@"     *"];
            NSLog(@"h: %20.13e %@", hLDP[i], temp);
            temp = [NSMutableString stringWithCapacity:0];
            for (NSUInteger j=20; j<40; j++) if (gLDP[i][j] != 0.0) [temp appendFormat:@" %5.1lf", gLDP[i][j]]; else [temp appendString:@"     *"];
            NSLog(@"                        %@", temp);
            temp = [NSMutableString stringWithCapacity:0];
            for (NSUInteger j=40; j<numberOfPhaseComponentEndmembers; j++) if (gLDP[i][j] != 0.0) [temp appendFormat:@" %5.1lf", gLDP[i][j]]; else [temp appendString:@"     *"];
            NSLog(@"                        %@", temp);
        }

        for (NSUInteger i=numberOfPhaseComponentEndmembers; i<numberOfInequalityConstraints; i++) {
            NSMutableString *temp = [NSMutableString stringWithCapacity:0];
            for (NSUInteger j=0; j<20; j++) if (gLDP[i][j] != 0.0) [temp appendFormat:@" %5.1lf", gLDP[i][j]]; else [temp appendString:@"     *"];
            NSLog(@"h: %20.13e %@", hLDP[i], temp);
            temp = [NSMutableString stringWithCapacity:0];
            for (NSUInteger j=20; j<40; j++) if (gLDP[i][j] != 0.0) [temp appendFormat:@" %5.1lf", gLDP[i][j]]; else [temp appendString:@"     *"];
            NSLog(@"                        %@", temp);
            temp = [NSMutableString stringWithCapacity:0];
            for (NSUInteger j=40; j<numberOfPhaseComponentEndmembers; j++) if (gLDP[i][j] != 0.0) [temp appendFormat:@" %5.1lf", gLDP[i][j]]; else [temp appendString:@"     *"];
            NSLog(@"                        %@", temp);
        }
    }

    DoubleVector *xLDPWrapper = [[DoubleVector alloc] initWithSize:numberOfPhaseComponentEndmembers];
	double *xLDP  = [xLDPWrapper pointerToDouble];
    DoubleVector *wLDPWrapper = [[DoubleVector alloc] initWithSize:2*numberOfPhaseComponentEndmembers + 2*numberOfInequalityConstraints + 2];
    double *wLDP  = [wLDPWrapper pointerToDouble];
    DoubleMatrix *eLDPWrapper = [[DoubleMatrix alloc] initWithRowSize:numberOfPhaseComponentEndmembers+1
                                             andWithColumnSize:numberOfInequalityConstraints];
    double **eLDP = [eLDPWrapper pointerToPointerToDouble];
    IntegerVector *indexLDPWrapper = [[IntegerVector alloc] initWithSize:numberOfInequalityConstraints andInitialValue:0];
	NSInteger *indexLDP = [indexLDPWrapper pointerToInteger];
    double xNorm;

    if (self.debugS) NSLog(@"Call to LDP with %lu rows and %lu columns.", numberOfInequalityConstraints, numberOfPhaseComponentEndmembers);
	if (![MathSupport ldpWithConstraintMatrix:gLDP
                             andRowDimension:numberOfInequalityConstraints
                          andColumnDimension:numberOfPhaseComponentEndmembers
                         andConstraintVector:hLDP
                           andSolutionVector:xLDP
                       andSolutionVectorNorm:&xNorm
                      andDoubleWorkingVector:wLDP
                      andDoubleWorkingMatrix:eLDP
                     andIntegerWorkingVector:indexLDP]) {
        NSLog(@"Failure in LDP (estimateLiquidAbsentEquilibriumAssemblageUsingLDP)");
    }

	NSUInteger colIndex = 0;
	for (NSArray *object in constraints) {
		NSString *key = [object objectAtIndex:0];
		NSUInteger na = [[object objectAtIndex:1] unsignedIntValue];

        EquilibrateStatePhase *phaseWrapper = [[equilibrateState potentialPhaseList] objectForKey:key];
        [[equilibrateState potentialPhaseList] removeObjectForKey:key];
        double *phaseCompositionInElements  = [[phaseWrapper bulkCompositionInElements] pointerToDouble];
        for (NSUInteger i=1; i<107; i++) phaseCompositionInElements[i] = 0.0;

        for (NSUInteger i=0; i<na; i++, colIndex++) {
            double *elementsInPhase = [[object objectAtIndex:3+2*i] pointerValue];
			[EquilibrateState accumulateIntoVectorOfElements:phaseCompositionInElements aScaler:xLDP[colIndex] timesElementReferenceVector:elementsInPhase];
            if (self.debugS) NSLog(@"...solution moles for %@ [%@] is %20.13e", [key stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
                              [[object objectAtIndex:2+2*i] stringByPaddingToLength:20 withString:@" " startingAtIndex:0], xLDP[colIndex]);
        }

        double massOfPhase = 0.0;
        if ([[phaseWrapper phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
            massOfPhase = [[phaseWrapper phaseClassInstance] convertElementsToTotalMass:phaseCompositionInElements];
            double *m = [[[phaseWrapper phaseClassInstance] convertElementsToMoles:phaseCompositionInElements] pointerToDouble];
            if (self.debugS) NSLog(@"...Formula: %@", [[phaseWrapper phaseClassInstance] getFormulaFromMolesOfComponents:m andT:[equilibrateState T] andP:[equilibrateState P]]);
        } else {
            massOfPhase = [[phaseWrapper phaseClassInstance] convertElementsToMassOfPhase:phaseCompositionInElements];
        }
        if (self.debugS) NSLog(@"...Mass of phase: %20.13e", massOfPhase);
        [phaseWrapper setMass:massOfPhase];

        [[equilibrateState phasesInSystem] setObject:phaseWrapper forKey:key];
        if ([[phaseWrapper phaseClassInstance] respondsToSelector:@selector(incrementInstanceCountOfPhase)]) [[phaseWrapper phaseClassInstance] incrementInstanceCountOfPhase];

	}
    NSAssert(colIndex == numberOfPhaseComponentEndmembers, @"Improper solution vector reallocation (columns) for LDP matrix.");

    initialPhaseAssemblageRequiresEquilibrationStep = YES;
}

/**
 Sets the bulk composition of the system.

 Wipes the system clean of all phases and allocates an instance of a liquid phase (which may be metastable).

 @param
   wtOxides NSArray of length [self.systemOxides count] containing objects of class NSNumber.
   Each object holds a doubleValue of grams of an oxide.
 */
-(void)setComposition:(double *)wtOxides {
	double *bulkCompOxides = [[equilibrateState bulkCompositionInOxides] pointerToDouble];
	double *bulkCompElements = [[equilibrateState bulkCompositionInElements] pointerToDouble];
    double mass = 0.0;
    [EquilibrateState scaleVectorOfElements:bulkCompElements aScaler: 0.0];
    for (NSUInteger i=0; i<[self.systemOxides count]; i++) {
        PhaseBase *oxide = [self.bulkComposition objectAtIndex:i];

        double value = wtOxides[i];
        double moles = value/[oxide mw];

        bulkCompOxides[i] = moles;
        mass += value;
        [EquilibrateState accumulateIntoVectorOfElements:bulkCompElements
                                                 aScaler:moles
                             timesElementReferenceVector:[[oxide formulaAsElementArray] pointerToDouble]];
    }
    [equilibrateState setMass:mass];

	if ([[equilibrateState phasesInSystem] count] > 0) {
		if (self.debugV) NSLog(@"The phasesInSystem dictionary has a non-zero length.  Removing all objects.");
		//
		// To do:  eliminate more than one entry for a given phase by looking at key name
		[[equilibrateState potentialPhaseList] addEntriesFromDictionary:[equilibrateState phasesInSystem]]; // transfer all phases back to potential list
		[[equilibrateState phasesInSystem] removeAllObjects];                                         // wipe all the phases from the system
	}

	BOOL liquidAllowed = YES;
	if (permissablePhasesState) {
		NSNumber *liquidEntry = [permissablePhasesState valueForKey:@"Liquid"];
		if (liquidEntry) liquidAllowed = [liquidEntry boolValue];
		else             liquidAllowed = NO;
	}

	if (self.debugV) NSLog(@"The phasesInSystem dictionary has length zero or has been zeroed.");

	if (liquidAllowed) {
		if (self.debugS) NSLog(@"Adding a liquid phase to the system.");
		EquilibrateStatePhase *phase = [[equilibrateState potentialPhaseList] valueForKey:@"Liquid"];
		[[equilibrateState potentialPhaseList] removeObjectForKey:@"Liquid"];

		[phase setMass:[equilibrateState mass]];                                                       // its mass is the system mass
		double *bulkCompositionInElements = [[phase bulkCompositionInElements] pointerToDouble]; // its composition is the bulk composition of the system
        [EquilibrateState scaleVectorOfElements:bulkCompositionInElements aScaler:0.0];
		[EquilibrateState accumulateIntoVectorOfElements:bulkCompositionInElements
										   aScaler:1.0
					   timesElementReferenceVector:[[equilibrateState bulkCompositionInElements] pointerToDouble]];

		[[equilibrateState phasesInSystem] setObject:phase forKey:@"Liquid"];
		if ([[phase phaseClassInstance] respondsToSelector:@selector(incrementInstanceCountOfPhase)]) [[phase phaseClassInstance] incrementInstanceCountOfPhase];
		if (self.debugV) NSLog(@"... Liquid phase added to initial system.");

	} else {
		if (self.debugS) NSLog(@"Performing an equilibration step using stoichiometric phases (liquid absent case).");
        //[self estimateLiquidAbsentEquilibriumAssemblageUsingLDPandMinimumConcentration:1.0e-4];
		[self estimateLiquidAbsentEquilibriumAssemblage:kDoNotForceInclusionOfAnyEndmemberOfAnyPhaseInSystem];
		if (self.debugV) NSLog(@"...Exiting set composition (liquid absent branch).");
	}
}

/**
 Sets the permissable phase state.

 Turns on or off entries in the permissablePhasesState array that determine inclusion of phases in the system.

 @param
 arrayOfStateValues NSArray of containing objects of class NSNumber.  Each object holds a BOOL value. The order is the same as returned by the method phaseNames.
 */
-(void)setPermissablePhasesState:(NSArray *)arrayOfStateValues {
	if (self.debugV) NSLog(@"Entering [Equilibrate setPermissablePhasesState ...");
	NSArray *permissablePhasesName = nil;
	permissablePhasesName = self.phaseNames;
	NSUInteger nPhases = [permissablePhasesName count];
	for (NSUInteger i=0; i<nPhases; i++) {
		BOOL value = [[arrayOfStateValues objectAtIndex:i] boolValue];
		NSString *key = [permissablePhasesName objectAtIndex:i];
		[permissablePhasesState setObject:[NSNumber numberWithBool:value] forKey:key];
	}
}
-(NSDictionary *)getPermissablePhasesState {
    return [NSDictionary dictionaryWithDictionary:permissablePhasesState];
}

/**
 Set the options for the Phase equilibration.

 @param
   dictionaryOfCalculationOptions NSDictionary containing options corresponding to known keys:
   - the value of the key @"ordinate" is either @"Temperature", @"Enthalpy", or @"Entropy"
   - the value of the key @"abscissa" is either @"Pressure" or @"Volume"
   - the value of the key @"buffer: is either @"HM", @"IW", @"NNO", or @"QFM"
   - the value of the key @"imposeBuffer" is either the boolean YES or NO contained in an NSNumber object
   - the value of the key @"bufferOffset" is a double contained in an NSNumber object

 */
-(void)setCalculationOptions:(NSDictionary *)dictionaryOfCalculationOptions {
	if (self.debugV) NSLog(@"Entering [Equilibrate setCalculationOptions ...");
	calculationOptions = [NSDictionary dictionaryWithDictionary:dictionaryOfCalculationOptions];
    if (self.debugV) {
        for (NSString *key in [dictionaryOfCalculationOptions allKeys]) {
            NSLog(@"%@ %@", key, dictionaryOfCalculationOptions[key]);
        }
    }
}

/**
 Re-distributes a phase composition (in moles of elements) amongst all other phases in the system.

 The algorithm uses HFTI least squares (Lawson and Hanson, 1974) to compute a correction vector to the system phase
 composition.  This is a speed improvement over Asimow and Ghiorso (1995) who utilized Singular Value Decomposition for the
 same purpose.  Alternatively, the algorithm can be made to use NNLS (non-negative least squares, Lawson and Hansen, 1974)
 to insure that the redistribution always generates positive additions to the base phase assemblage.

 @param
 elementVector of moles of elements on length 107.

 @param
 scaleFactor of adding this mass to the system phases (+) or subtracting this mass from the system phases (-).

 @param
 useNNLS BOOL value.  YES forces use of NNLS.  NO forces use of HFTI.

 @return
 BOOL value for successful (YES) or unsuccessful (NO) redistribution of mass.

 */
-(BOOL)distributeMolesOfElementsAmongstSystemPhases:(double *)elementVector withMultiplier:(double)scaleFactor withExcludedPhases:(NSSet *)exclude withNonNegativeSolution:(BOOL)useNNLS {
	NSArray *nonZeroElements = [equilibrateState hashTableOfNonZeroElementEntries];
	NSUInteger numberOfRowConstraints = [nonZeroElements count];
	NSUInteger numberOfColumnConstraints = 0;
	NSMutableDictionary *phaseList = [NSMutableDictionary dictionaryWithCapacity:0];
	NSDictionary *phasesAvailable = [equilibrateState phasesInSystem];
	for (NSString *key in [phasesAvailable allKeys]) {
		if (![exclude containsObject:key]) {
			id phaseClassInstance = [[phasesAvailable valueForKey:key] phaseClassInstance];

			NSMutableArray *nonZero = [NSMutableArray arrayWithCapacity:0];
			if ([phaseClassInstance conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
				NSUInteger na = [phaseClassInstance numberOfSolutionComponents];
                DoubleVector *mWrapper = [phaseClassInstance convertElementsToMoles:[[[phasesAvailable valueForKey:key] bulkCompositionInElements] pointerToDouble]];
				double *m = [mWrapper pointerToDouble];
				for (NSUInteger i=0; i<na; i++) if (m[i] != 0.0) {
					double *bc = [[[phaseClassInstance componentAtIndex:i] formulaAsElementArray] pointerToDouble];
					[nonZero addObject:[NSValue valueWithPointer:bc]];
				}
			} else {
				double *ePure = [[(PhaseBase *)phaseClassInstance formulaAsElementArray] pointerToDouble];
				[nonZero addObject:[NSValue valueWithPointer:ePure]];
			}

			numberOfColumnConstraints += [nonZero count];
			if ([nonZero count] > 0) [phaseList setObject:nonZero forKey:key];
		}
	}

	NSUInteger maxNumberOfConstraints = (numberOfRowConstraints < numberOfColumnConstraints) ? numberOfColumnConstraints : numberOfRowConstraints;
    DoubleMatrix *aHFTIWrapper = [[DoubleMatrix alloc] initWithRowSize:maxNumberOfConstraints
                                                     andWithColumnSize:numberOfColumnConstraints
                                                       andInitialValue:0.0];
	double **aHFTI = [aHFTIWrapper pointerToPointerToDouble];
    DoubleMatrix *bHFTIWrapper = [[DoubleMatrix alloc] initWithRowSize:maxNumberOfConstraints
                                                     andWithColumnSize:1];
	double **bHFTI = [bHFTIWrapper pointerToPointerToDouble];
    DoubleVector *hHFTIWrapper = [[DoubleVector alloc] initWithSize:maxNumberOfConstraints];
	double  *hHFTI = [hHFTIWrapper pointerToDouble];
    DoubleVector *gHFTIWrapper = [[DoubleVector alloc] initWithSize:maxNumberOfConstraints];
	double  *gHFTI = [gHFTIWrapper pointerToDouble];
    IntegerVector *pHFTIWrapper = [[IntegerVector alloc] initWithSize:maxNumberOfConstraints andInitialValue:0];
	NSInteger     *pHFTI = [pHFTIWrapper pointerToInteger];

	NSUInteger rowIndex = 0;
	for (NSNumber *entry in nonZeroElements) {
		NSUInteger elmNo = [entry unsignedIntValue];
		bHFTI[rowIndex][0] = elementVector[elmNo];

		NSUInteger colIndex = 0;
		for (NSString *key in [phaseList allKeys]) {
			NSArray *nonZero = [phaseList objectForKey:key];
			for (NSValue *object in nonZero) {
				double *elementsStoichiometry = [object pointerValue];
				aHFTI[rowIndex][colIndex] = elementsStoichiometry[elmNo];
				colIndex++;
			}
		}

		rowIndex++;
	}

	double tolerance = 100.0*DBL_EPSILON, rNormHFTI;
	NSInteger pseudoRank;

    if (useNNLS) {
        DoubleVector *bNNLSWrapper = [[DoubleVector alloc] initWithSize:maxNumberOfConstraints andInitialValue:0.0];
        double *bNNLS = [bNNLSWrapper pointerToDouble];
        DoubleVector *xNNLSWrapper = [[DoubleVector alloc] initWithSize:maxNumberOfConstraints andInitialValue:0.0];
        double *xNNLS = [xNNLSWrapper pointerToDouble];
        for (NSUInteger i=0; i<numberOfColumnConstraints; i++) bNNLS[i] = bHFTI[i][0];
        [MathSupport nnlsWithConstraintMatrix:aHFTI
                              andRowDimension:numberOfRowConstraints
                           andColumnDimension:numberOfColumnConstraints
                          andConstraintVector:bNNLS
                            andSolutionVector:xNNLS
                              andResidualNorm:&rNormHFTI
                        andDualSolutionVector:hHFTI
                       andDoubleWorkingVector:gHFTI
                     andUIntegerWorkingVector:pHFTI];
        for (NSUInteger i=0; i<numberOfColumnConstraints; i++) bHFTI[i][0] = xNNLS[i];
        pseudoRank = numberOfColumnConstraints;
    } else {
        [MathSupport hfti:aHFTI
                        m:numberOfRowConstraints
                        n:numberOfColumnConstraints
                        b:bHFTI
                       nb:1
                      tau:tolerance
                        k:&pseudoRank
                    rnorm:&rNormHFTI
                        h:hHFTI
                        g:gHFTI
                        p:pHFTI];
    }

	if (self.debugS) {
		NSLog(@"%@ phase re-distribution algorithm.  Results (scaleFactor = %g):", (useNNLS) ? @"NNLS" : @"HFTI", scaleFactor);
		NSLog(@"...rNorm = %g, pseudorank = %ld, col = %lu, row = %lu", rNormHFTI, pseudoRank, numberOfColumnConstraints, numberOfRowConstraints);
	}
	if (self.debugV) for (NSUInteger i=0; i<numberOfColumnConstraints; i++) NSLog(@"...x[%2.2lu] = %13.6g", i, bHFTI[i][0]);
	BOOL failedToBalanceElementTotals = (rNormHFTI > tolerance);

	if (failedToBalanceElementTotals || (useNNLS && self.debugS)) {
		NSUInteger rowIndex = 0;
		for (NSNumber *entry in nonZeroElements) {
			NSUInteger elmNo = [entry unsignedIntValue];
			double residual = elementVector[elmNo];

			NSUInteger colIndex = 0;
			for (NSString *key in [phaseList allKeys]) {
				NSArray *nonZero = [phaseList objectForKey:key];
				for (NSValue *object in nonZero) {
					double *elementsStoichiometry = [object pointerValue];
					residual -= elementsStoichiometry[elmNo]*bHFTI[colIndex][0];
					colIndex++;
				}
			}

			NSLog(@"Error for element %@ is %g", [[PhaseBase elementNameFromAtomicNumber:elmNo] stringByPaddingToLength:2 withString:@" " startingAtIndex:0], residual);
            if (fabs(residual) > tolerance) {
                NSString *pp = [[PhaseBase elementNameFromAtomicNumber:elmNo] stringByAppendingString:@"PseudoPhase"];
                id object = [[equilibrateState phasesInSystem] objectForKey:pp];
                if (object) {
                    NSLog(@"... PseudoPhase %@ is in system phase list.", pp);
                } else if (![[equilibrateState phasesInSystem] objectForKey:@"Liquid"]) {
                    NSLog(@"... Liquid is not present so the PseudoPhase %@ should be addd to balance the elemental abundances.  Terminating", pp);
                }
            }
			rowIndex++;
		}
		if (failedToBalanceElementTotals) return NO;
	}

	NSMutableSet *phasesWhoseCompositionsAreNoLongerValid = [NSMutableSet setWithCapacity:1];
	NSUInteger colIndex = 0;
	for (NSString *key in [phaseList allKeys]) {
		NSArray *nonZero = [phaseList objectForKey:key];
		BOOL phaseNeedsUpdating = NO;
        DoubleVector *deltaPhaseCompositionInElementsWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
		double *deltaPhaseCompositionInElements = [deltaPhaseCompositionInElementsWrapper pointerToDouble];
		for (NSValue *object in nonZero) {
			if (bHFTI[colIndex][0] != 0.0) {
				double *elementsStoichiometry = [object pointerValue];
				[EquilibrateState accumulateIntoVectorOfElements:deltaPhaseCompositionInElements
												   aScaler:bHFTI[colIndex][0]
							   timesElementReferenceVector:elementsStoichiometry];
				phaseNeedsUpdating = YES;
				if (self.debugS) NSLog(@"...x[%@] = %13.6g", [key stringByPaddingToLength:10 withString:@" " startingAtIndex:0], bHFTI[colIndex][0]);
			}
			colIndex++;
		}
		if (phaseNeedsUpdating) {
			EquilibrateStatePhase *phaseWrapper = [[equilibrateState phasesInSystem] valueForKey:key];
			double *phaseCompositionInElements  = [[phaseWrapper bulkCompositionInElements] pointerToDouble];
			[EquilibrateState accumulateIntoVectorOfElements:phaseCompositionInElements
											   aScaler:scaleFactor
						   timesElementReferenceVector:deltaPhaseCompositionInElements];
			double massOfPhase = 0.0;
			if ([[phaseWrapper phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
				massOfPhase = [[phaseWrapper phaseClassInstance] convertElementsToTotalMass:phaseCompositionInElements];
				if (![[phaseWrapper phaseClassInstance] testPermissibleValuesOfComponents:[[[phaseWrapper phaseClassInstance] convertElementsToMoles:phaseCompositionInElements] pointerToDouble]] ) {
					if (self.debugS) NSLog(@"... in distributeMolesOfElementsAmongstSystemPhases. Phase %@ has illegal composition.", key);
					[phasesWhoseCompositionsAreNoLongerValid addObject:key];
				}
			} else {
				massOfPhase = [[phaseWrapper phaseClassInstance] convertElementsToMassOfPhase:phaseCompositionInElements];
			}
			[phaseWrapper setMass:massOfPhase];
		}
	}

	if ([phasesWhoseCompositionsAreNoLongerValid count] > 0) {
		if (self.debugS) NSLog(@"... in distributeMolesOfElementsAmongstSystemPhases. Calling method recursively to reverse composition correction.");
		BOOL recursiveCall = [self distributeMolesOfElementsAmongstSystemPhases:elementVector
                                                                 withMultiplier:-scaleFactor
                                                             withExcludedPhases:exclude
                                                        withNonNegativeSolution:useNNLS];
		if (!recursiveCall) {
			if (self.debugS) NSLog(@"... in distributeMolesOfElementsAmongstSystemPhases. Cannot resurrect unaltered phase assemblage. Exiting.");
			return NO;
		}
		if (self.debugS) NSLog(@"... in distributeMolesOfElementsAmongstSystemPhases. Calling method recursively with illegal phase excluded.");
		[phasesWhoseCompositionsAreNoLongerValid unionSet:exclude];
		recursiveCall = [self distributeMolesOfElementsAmongstSystemPhases:elementVector
                                                            withMultiplier:scaleFactor
                                                        withExcludedPhases:phasesWhoseCompositionsAreNoLongerValid
                                                   withNonNegativeSolution:useNNLS];
		if (!recursiveCall) {
			if (self.debugS) NSLog(@"... in distributeMolesOfElementsAmongstSystemPhases. Failure with supplemented excluded set.");
			return NO;
		}
	}

	return YES;
}

/**
 Computes chemical potentials of elements from current system phase assemblage

 The solution will be unique if the system is in equilibrium and the norm of the residuals will be zero.
 A non-zero residual norm indicates a disequilibrium state.

 @param
 elementMu of moles of elements on length 107. Solution is returned in this vector

 @return
 Residual norm.
 */
-(double)chemicalPotentialsOfElements:(double *)elementMu {
	NSArray *nonZeroElements = [equilibrateState hashTableOfNonZeroElementEntries];

	NSUInteger numberOfRowConstraints    = 0;
	NSUInteger numberOfColumnConstraints = [nonZeroElements count];

	NSMutableDictionary *phaseList = [NSMutableDictionary dictionaryWithCapacity:0];
	for (NSString *key in [[equilibrateState phasesInSystem] allKeys]) if (![key hasSuffix:@"PseudoPhase"])  {
        EquilibrateStatePhase *object = [[equilibrateState phasesInSystem] valueForKey:key];
        NSArray	*nonZero = [object nonZeroMolesOfEndmemberComponents:(independentT | independentP) andT:[equilibrateState T] andP:[equilibrateState P]];
        numberOfRowConstraints += [nonZero count];
        [phaseList setObject:nonZero forKey:key];
	}

    DoubleMatrix *cMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:(numberOfRowConstraints < numberOfColumnConstraints) ? numberOfColumnConstraints : numberOfRowConstraints
                                                       andWithColumnSize:numberOfColumnConstraints
                                                         andInitialValue:0.0];
	double **cMatrix = [cMatrixWrapper pointerToPointerToDouble];
    DoubleMatrix *dVectorWrapper = [[DoubleMatrix alloc] initWithRowSize:(numberOfRowConstraints < numberOfColumnConstraints) ? numberOfColumnConstraints : numberOfRowConstraints
                                                       andWithColumnSize:1];
	double **dVector = [dVectorWrapper pointerToPointerToDouble];

	NSUInteger rowIndex = 0;
	double fNorm = 0.0;
	for (NSString *key in [phaseList allKeys]) {
		NSArray *nonZero = [phaseList objectForKey:key];
		for (NSArray *object in nonZero) {
			double *elementsStoichiometry = [[object objectAtIndex:1] pointerValue];

			NSInteger colIndex = 0;
			for (NSNumber *entry in nonZeroElements) {
				NSUInteger elmNo = [entry unsignedIntValue];
				cMatrix[rowIndex][colIndex] = elementsStoichiometry[elmNo];
				colIndex++;
			}

			dVector[rowIndex][0] = [[object objectAtIndex:2] doubleValue];
			fNorm += pow(dVector[rowIndex][0], 2.0);
			rowIndex++;
		}
	}
	fNorm = sqrt(fNorm);

    DoubleVector *hVectorWrapper = [[DoubleVector alloc] initWithSize:numberOfColumnConstraints];
	double *hVector = [hVectorWrapper pointerToDouble];
    DoubleVector *gVectorWrapper = [[DoubleVector alloc] initWithSize:numberOfColumnConstraints];
	double *gVector = [gVectorWrapper pointerToDouble];
    IntegerVector *pVectorWrapper = [[IntegerVector alloc] initWithSize:numberOfColumnConstraints andInitialValue:0];
	NSInteger    *pVector = [pVectorWrapper pointerToInteger];

	double tolerance = 10.0*DBL_EPSILON, rNorm;
	NSInteger pseudoRank;
	[MathSupport hfti:cMatrix
					m:numberOfRowConstraints
					n:numberOfColumnConstraints
					b:dVector
				   nb:1
				  tau:tolerance
					k:&pseudoRank
				rnorm:&rNorm
					h:hVector
					g:gVector
					p:pVector];
	if (self.debugS) NSLog(@"chemicalPotentialsOfElements: HFTI completed. Dimension = %lu, %lu, pseudo rank = %ld, rNorm = %g, (%g)",
					  numberOfRowConstraints, numberOfColumnConstraints, pseudoRank, rNorm, (fNorm > 0.0) ? rNorm/fNorm : 0.0);

	if (self.debugS) {
		rowIndex = 0;
		for (NSString *key in [phaseList allKeys]) {
			NSArray *nonZero = [phaseList objectForKey:key];
			for (NSArray *object in nonZero) {
				double *elementsStoichiometry = [[object objectAtIndex:1] pointerValue];

				double sum = 0.0;
				NSUInteger colIndex = 0;
				for (NSNumber *entry in nonZeroElements) {
					NSUInteger elmNo = [entry unsignedIntValue];
					sum += elementsStoichiometry[elmNo]*dVector[colIndex][0];
					colIndex++;
				}

				NSLog(@"[%@][%2.2lu] = %13.6g, est = %13.6g, delta = %13.6g", [key stringByPaddingToLength:25 withString:@" " startingAtIndex:0],
					  rowIndex, [[object objectAtIndex:2] doubleValue], sum, [[object objectAtIndex:2] doubleValue] - sum);
				rowIndex++;
			}
		}
	}

	NSInteger colIndex = 0;
	for (NSNumber *entry in nonZeroElements) {
		NSUInteger elmNo = [entry unsignedIntValue];
		elementMu[elmNo] = dVector[colIndex][0];

        if (elementMu[elmNo] == 0.0) {
            NSString *name = [PhaseBase elementNameFromAtomicNumber:elmNo];
            if ([pseudoPhaseDictionary objectForKey:name] != nil) {
                elementMu[elmNo] = [(PseudoPhase *)[pseudoPhaseDictionary objectForKey:name] g];
                if (self.debugS) NSLog(@"...Chemical potential for %@ defaults to pseudoPhase value.", name);
            }
        }

		if (self.debugV) NSLog(@"... HFTI solution ElmNo[%@] = %g", [[PhaseBase elementNameFromAtomicNumber:elmNo] stringByPaddingToLength:2 withString:@" " startingAtIndex:0],
						  dVector[colIndex][0]);
		colIndex++;
	}

	return ((fNorm > 0.0) ? rNorm/fNorm : 0.0);
}

/**
 Determines if a phase is supersaturated in the system realative to a liquid phase or a collection of solid phases.

 The most supersaturated phase is transferred from the potentialPhaseList dictionary to the phasesInSystem dictionary.

 @return
   A YES if supersaturation is detected.  A NO if the current system is stable.
 */
-(BOOL)evaluateSaturationState:(NSString *)lastPhaseRemovedFromSystem {
	double t = [equilibrateState T];
	double p = [equilibrateState P];
	BOOL hasLiquid = YES;

	EquilibrateStatePhase *liquidWrapper;
	id <SolutionPhaseProtocol> liquid;
    DoubleVector *elementMuWrapper = [[DoubleVector alloc] initWithSize:107 andInitialValue:0.0];
    DoubleVector *liqMolesWrapper, *liqMuWrapper;
	double *elementMu = [elementMuWrapper pointerToDouble];
	NSUInteger nl=0;

	if (![[equilibrateState phasesInSystem] objectForKey:@"Liquid"]) {
		if (self.debugS) NSLog(@"... in evaluateSaturationState (liquid not present in phasesInSystem.");
		hasLiquid = NO;
        double residualNorm = [self chemicalPotentialsOfElements:elementMu];
		if (!self.PPOPTIONS_CHEM_POTENTIAL.boolValue) {
            BOOL hasPseudoPhase = NO;
            for (NSString *pn in [[equilibrateState phasesInSystem] allKeys]) if ([pn hasSuffix:@"PseudoPhase"]) { hasPseudoPhase = YES; break; }
            if (!hasPseudoPhase) {
                double defaultTolerance = self.PPPARAMETERS_CHEM_POTENTIAL.doubleValue;
                if (!convergenceByMinimalEnergy && (residualNorm > defaultTolerance)) return NO;  // bipass calculation and force exit from calling program
                convergenceByMinimalEnergy = NO;
            }
		}
        for (NSUInteger i=1; i<107; i++) if (elementMu[i] != 0.0) {
            PseudoPhase *pseudoPhase = [pseudoPhaseDictionary objectForKey:[PhaseBase elementNameFromAtomicNumber:i]];
            if (pseudoPhase) [pseudoPhase setG:elementMu[i]];
        }
	} else {
		liquidWrapper = [[equilibrateState phasesInSystem] valueForKey:@"Liquid"];
		liquid          = [liquidWrapper phaseClassInstance];
		liqMolesWrapper = [liquid convertElementsToMoles:[[liquidWrapper bulkCompositionInElements] pointerToDouble]];
		liqMuWrapper    = [liquid getChemicalPotentialFromMolesOfComponents:[liqMolesWrapper pointerToDouble] andT:t andP:p];
		nl              = [liquid numberOfSolutionComponents];
	}

	double lowestAffinity = 0.0;
	NSString *phaseWithLowestAffinity = @"";
    NSMutableSet *allPhasesWithNegativeAffinities = [NSMutableSet setWithCapacity:0];

	for (EquilibrateStatePhase *phaseWrapper in [[equilibrateState potentialPhaseList] allValues]) {
		id phase = [phaseWrapper phaseClassInstance];

		BOOL includePhase = YES;
		if (permissablePhasesState) includePhase = [[permissablePhasesState valueForKey:[phase phaseName]] boolValue];
		if ([lastPhaseRemovedFromSystem isEqualToString:[phase phaseName]]) includePhase = NO;

		if ([phase conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
			if (self.debugV) NSLog(@"Potential phase %@ is a solution phase.", [phase phaseName]);
			NSUInteger na = [phase numberOfSolutionComponents];
            DoubleVector *chemicalPotentialSumWrapper = [[DoubleVector alloc] initWithSize:na andInitialValue:0.0];
			double *chemicalPotentialSum = [chemicalPotentialSumWrapper pointerToDouble];
			for (NSUInteger i=0; i<na; i++) {
				chemicalPotentialSum[i] = 0.0;
				if (hasLiquid) {
					DoubleVector *coeffWrapper = [liquid convertElementsToMoles:[[[phase componentAtIndex:i] formulaAsElementArray] pointerToDouble]];
                    double *coeff    = [coeffWrapper pointerToDouble];
                    double *liqMoles = [liqMolesWrapper pointerToDouble];
                    double *liqMu    = [liqMuWrapper pointerToDouble];
					for (NSUInteger j=0; j<nl; j++) if (coeff[j] != 0.0) {
						if (liqMoles[j] != 0.0) {
                            if (self.debugV) NSLog(@"... ... i = %lu j = %lu coeff = %g", i, j, coeff[j]);
							chemicalPotentialSum[i] += coeff[j]*liqMu[j];
						} else {
							chemicalPotentialSum[i] = 0.0;
							break;
						}
					}
				} else {
                    DoubleVector *coeffWrapper = [[phase componentAtIndex:i] formulaAsElementArray];
					double *coeff = [coeffWrapper pointerToDouble];
					for (NSUInteger j=1; j<107; j++) if (coeff[j] != 0.0)  {
						if (elementMu[j] != 0.0) chemicalPotentialSum[i] += coeff[j]*elementMu[j];
						else {
							chemicalPotentialSum[i] = 0.0;
							break;
						}
					}
				}
				if (self.debugV) NSLog(@"... whose implied liquid chemical potential sum for component %@ is %f",
					  [[phase componentAtIndex:i] phaseName], chemicalPotentialSum[i]);

			}
			NSArray *affinityAndComposition = [phase affinityAndCompositionFromLiquidChemicalPotentialSum:chemicalPotentialSum andT:t andP:p];
			if (self.debugS) {
				NSLog(@"Converged: %@ - Phase %@ has the affinity is %f",
                      ([[affinityAndComposition objectAtIndex:na+1] boolValue]) ? @"YES" : @"NO",
                      [[phase phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
					  [[affinityAndComposition objectAtIndex:0] doubleValue]);
                DoubleVector *mWrapper = [[DoubleVector alloc] initWithSize:na];
				double *m = [mWrapper pointerToDouble];
				for (NSUInteger i=0; i<na; i++) {
					NSLog(@"... estimated composition of %@ is %g", [[[phase componentAtIndex:i] phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0],
						  [[affinityAndComposition objectAtIndex:1+i] doubleValue]);
					m[i] = [[affinityAndComposition objectAtIndex:1+i] doubleValue];
				}
				NSLog(@"... Formula: %@", [phase getFormulaFromMolesOfComponents:m andT:t andP:p]);
				if ([lastPhaseRemovedFromSystem isEqualToString:[phase phaseName]])
					NSLog(@"... Phase will be suppressed because it was recently removed from the system.");
			}

			// converged, so copy and save affinity and composition info (scale the composition by MOLESIN)
			if ([[affinityAndComposition objectAtIndex:na+1] boolValue]) {
				double affinity = [[affinityAndComposition objectAtIndex:0] doubleValue];
				[phaseWrapper setAffinity:affinity];
                DoubleVector *xWrapper = [[DoubleVector alloc] initWithSize:na];
				double *x = [xWrapper pointerToDouble];
				double molesIn = self.PPPARAMETERS_MOLESIN.doubleValue;
				for (NSUInteger i=0; i<na; i++) x[i] = [[affinityAndComposition objectAtIndex:1+i] doubleValue]*molesIn;
				[phaseWrapper setBulkCompositionInElements:[phase convertMolesToElements:x]];

				affinity /= [[affinityAndComposition objectAtIndex:na+3] doubleValue];  // scale this affinity with # atoms
                if (affinity < 1.0) [allPhasesWithNegativeAffinities addObject:[phase phaseName]];
				if (includePhase && (affinity < lowestAffinity)) {
					lowestAffinity = affinity;
					phaseWithLowestAffinity = [NSString stringWithString:[phase phaseName]];
				}
			} else {
				[phaseWrapper setAffinity:0.0];
                DoubleVector *xWrapper = [[DoubleVector alloc] initWithSize:na];
				double *x = [xWrapper pointerToDouble];
				double molesIn = self.PPPARAMETERS_MOLESIN.doubleValue;
				for (NSUInteger i=0; i<na; i++) x[i] = [[affinityAndComposition objectAtIndex:1+i] doubleValue]*molesIn;
				[phaseWrapper setBulkCompositionInElements:[phase convertMolesToElements:x]];
			}


		} else {
			if (self.debugV) NSLog(@"Potential phase %@ is a stoichiometric phase", [phase phaseName]);
			double *refElements = [[phase formulaAsElementArray] pointerToDouble];
			double affinity = [phase getGibbsFreeEnergyFromT:t andP:p];
			if (hasLiquid) {
                DoubleVector *coeffWrapper = [liquid convertElementsToMoles:refElements];
				double *coeff    = [coeffWrapper pointerToDouble];
                double *liqMoles = [liqMolesWrapper pointerToDouble];
                double *liqMu    = [liqMuWrapper pointerToDouble];
				for (NSUInteger i=0; i<nl; i++) if ( ((coeff[i] != 0.0) && (liqMoles[i] != 0.0)) || (coeff[i] == 0.0) ) affinity -= coeff[i]*liqMu[i];
				else { affinity = 0.0; break; }
			} else {
				for (NSUInteger i=1; i<107; i++) if (refElements[i] != 0.0) {
					if (elementMu[i] != 0.0) affinity -= refElements[i]*elementMu[i];
					else {
						affinity = 0.0;
						break;
					}
				}
			}
			[phaseWrapper setAffinity:affinity];
			double *elements = [[phaseWrapper bulkCompositionInElements] pointerToDouble];
			double molesIn = self.PPPARAMETERS_MOLESIN.doubleValue;
			[EquilibrateState accumulateIntoVectorOfElements:elements aScaler:molesIn timesElementReferenceVector:refElements];
			if (self.debugS) NSLog(@"Phase %@ has chemical affinity is %f", [[phase phaseName] stringByPaddingToLength:20 withString:@" " startingAtIndex:0], affinity);

			double numberOfAtomsInFormulaUnit = 0.0;
			for (NSUInteger i=1; i<107; i++) numberOfAtomsInFormulaUnit += refElements[i];
			affinity /= numberOfAtomsInFormulaUnit;
            if (affinity < 1.0) [allPhasesWithNegativeAffinities addObject:[phase phaseName]];
			if (includePhase && (affinity < lowestAffinity)) {
				lowestAffinity = affinity;
				phaseWithLowestAffinity = [NSString stringWithString:[phase phaseName]];
			}
		}

	}

	if (lowestAffinity < 0.0) {
        NSMutableSet *phasesThatMustBeAdded = [NSMutableSet setWithSet:[allPhasesWithNegativeAffinities objectsPassingTest:^(id obj, BOOL *stop) {
            return [phasesThatShouldBeIncludedInEvaluateSaturationState containsObject:obj];
        }]];
        [phasesThatMustBeAdded addObject:phaseWithLowestAffinity];
        [phasesThatShouldBeIncludedInEvaluateSaturationState removeAllObjects];

		if (self.debugS) NSLog(@"Phase %@ will be added to the assemblage (scaled affinity = %f)", phaseWithLowestAffinity, lowestAffinity);
        if ([phasesThatMustBeAdded count] > 1) {
            if (self.debugS) {
                NSLog(@"(Special case) - All these phases will be added to the assemblage:");
                for (NSString *key in phasesThatMustBeAdded) NSLog(@"... %@", key);
            }
        }

        for (NSString *nameOfPhaseToAddToAssemblage in phasesThatMustBeAdded) {
            EquilibrateStatePhase *phaseWrapper = [[equilibrateState potentialPhaseList] objectForKey:nameOfPhaseToAddToAssemblage];
            [[equilibrateState potentialPhaseList] removeObjectForKey:nameOfPhaseToAddToAssemblage];

            double *phaseCompositionInElements = [[phaseWrapper bulkCompositionInElements] pointerToDouble];
            double massOfPhase = 0.0;
            if ([[phaseWrapper phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)]) {
                massOfPhase = [[phaseWrapper phaseClassInstance] convertElementsToTotalMass:phaseCompositionInElements];
            } else {
                massOfPhase = [[phaseWrapper phaseClassInstance] convertElementsToMassOfPhase:phaseCompositionInElements];
            }
            [phaseWrapper setMass:massOfPhase];

            if (hasLiquid) {
                double *liquidCompositionInElements = [[liquidWrapper bulkCompositionInElements] pointerToDouble];
                [EquilibrateState accumulateIntoVectorOfElements:liquidCompositionInElements
                                                   aScaler:-1.0
                               timesElementReferenceVector:phaseCompositionInElements];
                if (![[liquidWrapper phaseClassInstance] testPermissibleValuesOfComponents:[[[liquidWrapper phaseClassInstance] convertElementsToMoles:liquidCompositionInElements] pointerToDouble]]) {
                    [EquilibrateState accumulateIntoVectorOfElements:liquidCompositionInElements
                                                       aScaler:1.0
                                   timesElementReferenceVector:phaseCompositionInElements];
                    if (![self distributeMolesOfElementsAmongstSystemPhases:phaseCompositionInElements
                                                             withMultiplier:-1.0
                                                         withExcludedPhases:[NSSet setWithObject:@"Liquid"]
                                                    withNonNegativeSolution:NO]) {
                        if (!self.PPOPTIONS_DISTRIBUTION.boolValue) return NO;
                    }
                } else [liquidWrapper setMass:[liquidWrapper mass]-massOfPhase];
            } else {
                if (![self distributeMolesOfElementsAmongstSystemPhases:phaseCompositionInElements
                                                         withMultiplier:-1.0
                                                     withExcludedPhases:phasesThatMustBeAdded  // was withExcludedPhases:[NSSet setWithObject:@"none"]]
                                                withNonNegativeSolution:NO]) {
                    if (!self.PPOPTIONS_DISTRIBUTION.boolValue) return NO;
                }
            }

            [[equilibrateState phasesInSystem] setObject:phaseWrapper forKey:nameOfPhaseToAddToAssemblage];
            if ([[phaseWrapper phaseClassInstance] respondsToSelector:@selector(incrementInstanceCountOfPhase)]) [[phaseWrapper phaseClassInstance] incrementInstanceCountOfPhase];
        }
		return YES;
	} else {
        NSMutableSet *pseudoPhasesThatCanBeDestabilized = [NSMutableSet setWithCapacity:0];
        for (NSString *key in [[equilibrateState phasesInSystem] allKeys]) if ([key hasSuffix:@"PseudoPhase"]) [pseudoPhasesThatCanBeDestabilized addObject:key];
        if ([pseudoPhasesThatCanBeDestabilized count] > 0) {
            if (self.debugS) NSLog(@"The unconstrained PseudoPhases will now be destabilized...");
            BOOL forceAnotherIteration = NO;
            for (NSString *key in pseudoPhasesThatCanBeDestabilized) {
                PseudoPhase *pseudoPhaseInstance = (PseudoPhase *)[[[equilibrateState phasesInSystem] objectForKey:key] phaseClassInstance];
                double oldG     = [pseudoPhaseInstance gLast];
                double currentG = [pseudoPhaseInstance g];
                if (fabs(currentG-oldG) > 1.0*DBL_EPSILON) {
                    [pseudoPhaseInstance setG:currentG+5000.0];
                    [pseudoPhaseInstance setGLast:currentG];
                    forceAnotherIteration = YES;
                    if (self.debugS) NSLog(@"... %@ is destablized by 5000 J from %13.3lf => %13.3lf (cur-old %13.6e)",
                                      [key stringByPaddingToLength:13 withString:@" " startingAtIndex:0], currentG, currentG+5000.0, currentG-oldG);
                }
            }
            return forceAnotherIteration;
        }
    }

	return NO;
}

// Initial guess failure flags for non-linear constraints
static const NSUInteger kNlConstraintsSuccess                         = 001LU;
static const NSUInteger kNlConstraintsFailureForTemperatureLowerBound = 002LU;
static const NSUInteger kNlConstraintsFailureForTemperatureUpperBound = 004LU;
static const NSUInteger kNlConstraintsFailureForPressureLowerBound    = 010LU;
static const NSUInteger kNlConstraintsFailureForPressureUpperBound    = 020LU;
static const NSUInteger kNlConstraintsFailureGeneric                  = 040LU;

-(NSUInteger)computeTemperatureFromSystemEntropy {
    double t = [equilibrateState T];
    double p = [equilibrateState P];
    double referenceEntropy = [equilibrateState referenceEntropyOfSystem] + [equilibrateState correctionToReferenceEntropyOfSystem];
    double sumOfComputedEntropies = 0.0;
    NSUInteger iterations = 0;

    if (referenceEntropy == 0.0) return kNlConstraintsFailureGeneric;
    if (t == 0.0) t = 1000.0;

    if (self.debugS) NSLog(@"Entry to computeTemperatureFromSystemEntropy ... (S = %g, T = %g", referenceEntropy, t-273.15);
    do {
        sumOfComputedEntropies = 0.0;
        double sumOfComputedHeatCapacities = 0.0;
        for (NSString *key in [[equilibrateState phasesInSystem] allKeys]) {
            EquilibrateStatePhase *object = [[equilibrateState phasesInSystem] valueForKey:key];
            double entropy      = [object entropyFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
            double heatCapacity = [object heatCapacityFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
            if (self.debugV) NSLog(@"*.*.*.*. Phase:%@ S = %g Cp = %g", [[object phaseClassInstance] phaseName], entropy, heatCapacity);
            if (entropy      != compositionIsUnacceptable) sumOfComputedEntropies      += entropy;
            if (heatCapacity != compositionIsUnacceptable) sumOfComputedHeatCapacities += heatCapacity;
        }
        if (self.debugV) NSLog(@"*.*.*.*. Iter = %lu Entropy %20.13e, Cp = %20.13e", iterations, sumOfComputedEntropies, sumOfComputedHeatCapacities);

        t -= (sumOfComputedEntropies-referenceEntropy)*t/sumOfComputedHeatCapacities;
        if (self.debugV) NSLog(@"*.*.*.*. T = %g", t-273.15);
        if      (t <  300.0) t =  300.0;
        else if (t > 5000.0) t = 5000.0;

        iterations++;
    } while ((iterations < 100) && (fabs(sumOfComputedEntropies-referenceEntropy) > 1.0e-9*referenceEntropy));

    [equilibrateState setT:t];
    if (self.debugS) NSLog(@"...Exit from computeTemperatureFromSystemEntropy with T = %g", t-273.15);
    if (iterations < 100) return kNlConstraintsSuccess;
    else {
        NSUInteger result = 0;
        if      (t ==  300.0) result |= kNlConstraintsFailureForTemperatureLowerBound;
        else if (t == 5000.0) result |= kNlConstraintsFailureForTemperatureUpperBound;
        else                  result |= kNlConstraintsFailureGeneric;
        return result;
    }
};

-(NSUInteger)computePressureFromSystemVolume {
    double t = [equilibrateState T];
    double p = [equilibrateState P];
    double referenceVolume = [equilibrateState referenceVolumeOfSystem] + [equilibrateState correctionToReferenceVolumeOfSystem];
    double sumOfComputedVolumes = 0.0;
    NSUInteger iterations = 0;

    if (referenceVolume == 0.0) return kNlConstraintsFailureGeneric;
    if (t == 0.0) t = 1000.0;
    if (p == 0.0) p = 2000.0;

    if (self.debugS) NSLog(@"Entry to computePressureFromSystemVolume ... (V = %g, P = %g", referenceVolume, p);
    do {
        sumOfComputedVolumes = 0.0;
        double sumOfComputedDvdp = 0.0;
        for (NSString *key in [[equilibrateState phasesInSystem] allKeys]) {
            EquilibrateStatePhase *object = [[equilibrateState phasesInSystem] valueForKey:key];
            double volume = [object volumeFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
            double dvdp   = [object dvdpFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
            if (self.debugV) NSLog(@"*.*.*.*. Phase:%@ V = %g dvdp = %g", [[object phaseClassInstance] phaseName], volume, dvdp);
            if (volume != compositionIsUnacceptable) sumOfComputedVolumes += volume;
            if (dvdp   != compositionIsUnacceptable) sumOfComputedDvdp    += dvdp;
        }
        if (self.debugV) NSLog(@"*.*.*.*. Iter = %lu Volume %20.13e, dvdp = %20.13e", iterations, sumOfComputedVolumes, sumOfComputedDvdp);

        p -= (sumOfComputedVolumes-referenceVolume)/sumOfComputedDvdp;
        if (self.debugV) NSLog(@"*.*.*.*. P = %g", p);
        if      (p <      1.0) p =      1.0;
        else if (p > 500000.0) p = 500000.0;

        iterations++;
    } while ((iterations < 100) && (fabs(sumOfComputedVolumes-referenceVolume) > 1.0e-9*referenceVolume));

    [equilibrateState setP:p];
    if (self.debugS) NSLog(@"...Exit from computePressureFromSystemVolume with P = %g", p);
    if (iterations < 100) return kNlConstraintsSuccess;
    else {
        NSUInteger result = 0;
        if      (p ==      1.0) result |= kNlConstraintsFailureForPressureLowerBound;
        else if (p == 500000.0) result |= kNlConstraintsFailureForPressureUpperBound;
        else                    result |= kNlConstraintsFailureGeneric;
        return result;
    }
};

-(NSUInteger)computeTemperatureAndPressureFromSystemEntropyAndVolume {
    double t = [equilibrateState T];
    double p = [equilibrateState P];
    double referenceEntropy = [equilibrateState referenceEntropyOfSystem] + [equilibrateState correctionToReferenceEntropyOfSystem];
    double referenceVolume  = [equilibrateState referenceVolumeOfSystem]  + [equilibrateState correctionToReferenceVolumeOfSystem];
    double sumOfComputedEntropies = 0.0, sumOfComputedVolumes = 0.0;
    NSUInteger iterations = 0;

    if ((referenceEntropy == 0.0) || (referenceVolume == 0.0)) return kNlConstraintsFailureGeneric;
    if (t == 0.0) t = 1000.0;
    if (p == 0.0) p = 2000.0;

    if (self.debugS) NSLog(@"Entry to computeTemperatureAndPressureFromSystemEntropyAndVolume ... S = %g V = %g, T = %g P = %g", referenceEntropy, referenceVolume, t-273.15, p);
    do {
        sumOfComputedEntropies = 0.0;
        sumOfComputedVolumes   = 0.0;
        double sumOfComputedHeatCapacity = 0.0, sumOfComputedDvdp = 0.0, sumOfComputedDvdt = 0.0;
        for (NSString *key in [[equilibrateState phasesInSystem] allKeys]) {
            EquilibrateStatePhase *object = [[equilibrateState phasesInSystem] valueForKey:key];
            double entropy      = [object entropyFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
            double volume       = [object volumeFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
            double heatCapacity = [object heatCapacityFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
            double dvdp         = [object dvdpFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
            double dvdt         = [object dvdtFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
            if (self.debugV) NSLog(@"*.*.*.*. Phase:%@ S = %g V = %g Cp = %g dvdp = %g dvdt = %g",
                              [[object phaseClassInstance] phaseName], entropy, volume, heatCapacity, dvdp, dvdt);
            if (entropy      != compositionIsUnacceptable) sumOfComputedEntropies    += entropy;
            if (volume       != compositionIsUnacceptable) sumOfComputedVolumes      += volume;
            if (heatCapacity != compositionIsUnacceptable) sumOfComputedHeatCapacity += heatCapacity;
            if (dvdp         != compositionIsUnacceptable) sumOfComputedDvdp         += dvdp;
            if (dvdt         != compositionIsUnacceptable) sumOfComputedDvdt         += dvdt;
        }
        if (self.debugV) NSLog(@"*.*.*.*. Iter = %lu S %20.13e V %20.13e Cp %20.13e dvdp %20.13e dvdt %20.13e",
                          iterations, sumOfComputedEntropies, sumOfComputedVolumes, sumOfComputedHeatCapacity, sumOfComputedDvdp, sumOfComputedDvdt);
        double d11 = sumOfComputedHeatCapacity/t;
        double d12 = -sumOfComputedDvdt;
        double d21 = sumOfComputedDvdt;
        double d22 = sumOfComputedDvdp;
        double determ = d11*d22 - d12*d21;
        t -= (  (sumOfComputedEntropies-referenceEntropy)*d22 - (sumOfComputedVolumes-referenceVolume)*d12 )/determ;
        p -= ( -(sumOfComputedEntropies-referenceEntropy)*d21 + (sumOfComputedVolumes-referenceVolume)*d11 )/determ;
        if (self.debugV) NSLog(@"*.*.*.*. T = %g P = %g", t-273.15, p);
        if      (t <    300.0) t =    300.0;
        else if (t >   5000.0) t =   5000.0;
        if      (p <      1.0) p =      1.0;
        else if (p > 500000.0) p = 500000.0;

        iterations++;
    } while ((iterations < 100) &&
             (fabs(sumOfComputedVolumes-referenceVolume) > 1.0e-9*referenceVolume) &&
             (fabs(sumOfComputedEntropies-referenceEntropy) > 1.0e-9*referenceEntropy));
    if (self.debugS) {
        NSLog(@"...Iterations %lu", iterations);
        NSLog(@"...Vref = %g V error = %g %% = %g", referenceVolume, sumOfComputedVolumes-referenceVolume, 100.0*(sumOfComputedVolumes-referenceVolume)/referenceVolume);
        NSLog(@"...Sref = %g S error = %g %% = %g", referenceEntropy, sumOfComputedEntropies-referenceEntropy, 100.0*(sumOfComputedEntropies-referenceEntropy)/referenceEntropy);
    }

    [equilibrateState setT:t];
    [equilibrateState setP:p];
    if (self.debugS) NSLog(@"...Exit from computeTemperatureAndPressureFromSystemEntropyAndVolume with T = %g and P = %g", t-273.15, p);
    if (iterations < 100) return kNlConstraintsSuccess;
    else {
        NSUInteger result = 0;
        if (fabs(sumOfComputedVolumes-referenceVolume)    > 1.0e-9*referenceVolume) {
            if      (p ==      1.0) result |= kNlConstraintsFailureForPressureLowerBound;
            else if (p == 500000.0) result |= kNlConstraintsFailureForPressureUpperBound;
            else                    result |= kNlConstraintsFailureGeneric;
        } else if (fabs(sumOfComputedEntropies-referenceEntropy) > 1.0e-9*referenceEntropy) {
            if      (t ==  300.0) result |= kNlConstraintsFailureForTemperatureLowerBound;
            else if (t == 5000.0) result |= kNlConstraintsFailureForTemperatureUpperBound;
            else                  result |= kNlConstraintsFailureGeneric;
        }
        return result;
    }
};

-(void)computeAndFillUsefulOutputQuantitiesForEquilibrateStateClass {
    if ([equilibrateState referenceEntropyOfSystem] == 0.0) {
        double t = [equilibrateState T];
        double p = [equilibrateState P];
        double sumOfComputedEntropies = 0.0;
        for (NSString *key in [[equilibrateState phasesInSystem] allKeys]) {
            EquilibrateStatePhase *object = [[equilibrateState phasesInSystem] valueForKey:key];
            double entropy = [object entropyFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
            if (entropy != compositionIsUnacceptable) sumOfComputedEntropies += entropy;
        }
        [equilibrateState setReferenceEntropyOfSystem:sumOfComputedEntropies];
    }
    if ([equilibrateState referenceVolumeOfSystem] == 0.0) {
        double t = [equilibrateState T];
        double p = [equilibrateState P];
        double sumOfComputedVolumes = 0.0;
        for (NSString *key in [[equilibrateState phasesInSystem] allKeys]) {
            EquilibrateStatePhase *object = [[equilibrateState phasesInSystem] valueForKey:key];
            double volume = [object volumeFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
            if (volume != compositionIsUnacceptable) sumOfComputedVolumes += volume;
        }
        [equilibrateState setReferenceVolumeOfSystem:sumOfComputedVolumes];
    }
    if ([equilibrateState referenceEnthalpyOfSystem] == 0.0) {
        double t = [equilibrateState T];
        double p = [equilibrateState P];
        double sumOfComputedEnthalpies = 0.0;
        for (NSString *key in [[equilibrateState phasesInSystem] allKeys]) {
            EquilibrateStatePhase *object = [[equilibrateState phasesInSystem] valueForKey:key];
            double enthalpy = [object enthalpyFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
            if (enthalpy != compositionIsUnacceptable) sumOfComputedEnthalpies += enthalpy;
        }
        [equilibrateState setReferenceEnthalpyOfSystem:sumOfComputedEnthalpies];
    }
    if ([[equilibrateState phasesInSystem] objectForKey:@"Liquid"]) {
        const double t0 = 1673.15,        /* K       */
        a =    0.196,
        b =    1.1492e4,                  /* K       */
        c =   -6.675,
        e =   -3.364,
        f =   -7.01e-7  * 1.0e5,          /* K/bar   */
        g =   -1.54e-10 * 1.0e5,          /* 1/bar   */
        h =    3.85e-17 * 1.0e5 * 1.0e5;  /* K/bar^2 */
        const NSUInteger indexFeO = 5, indexFe2O3 = 3;

        DoubleVector *elementWrapper = [[[equilibrateState phasesInSystem] objectForKey:@"Liquid"] bulkCompositionInElements];
        double *element = [elementWrapper pointerToDouble];
        DoubleVector *moleWrapper = [self molesOfOxidesFromMolesOfElements:element];
        double *mole = [moleWrapper pointerToDouble];

        if ((mole[indexFeO] != 0.0) && (mole[indexFe2O3] != 0.0))  {
            double sum = 0.0;
            for (NSUInteger i=0; i<15; i++) sum += mole[i];
            sum += mole[indexFe2O3];
            if (sum != 0.0) {
                const NSUInteger indexAl2O3 = 2, indexCaO = 10, indexNa2O = 11, indexK2O = 12;
                double t = [equilibrateState T];
                double p = [equilibrateState P];
                double temp = b/t + c + e*(1.0-t0/t - log(t/t0)) + f*p/t + g*(t-t0)*p/t + h*p*p/t
                - 2.243*mole[indexAl2O3]/sum
                - 1.828*mole[indexFeO]/sum
                + 3.201*mole[indexCaO]/sum
                + 5.854*mole[indexNa2O]/sum
                + 6.215*mole[indexK2O]/sum
                - 1.828*2.0*mole[indexFe2O3]/sum;

                double logfo2 = (log(mole[indexFe2O3]/mole[indexFeO]) - temp)/(a*log(10.0));
                [equilibrateState setFo2:logfo2];
                [equilibrateState setFo2Path:FO2_NNO];
                [equilibrateState setFo2Delta:logfo2 - (-24930.0/t +  9.360)];
            }
        }
    }
}

#pragma mark -
#pragma mark public dependent instance methods

/**
 Executes an equilibration calculation.

 @return
   An NSDictionary of status and output items:
     - status key:  A flag indicating success or failure of the calculation
     - results key: A EquilibrateState object
 */
-(NSDictionary *)execute {
	NSString *exitStatus = @"success, Trivial case with no quadratic search.";

	double *bulkComp = [[equilibrateState bulkCompositionInOxides] pointerToDouble];
    if (self.debugV) {
        for (NSUInteger i=0; i<[self.systemOxides count]; i++) NSLog(@"... Moles of %@ is %g", [self.systemOxides objectAtIndex:i], bulkComp[i]);
        NSLog(@"There are %lu phases in the phasesInSystem dictionary", [[equilibrateState phasesInSystem] count]);
        NSLog(@"There are %lu phases in the potentialPhaseList dictionary", [[equilibrateState potentialPhaseList] count]);
    }
	// BOOL hasLiquid = ([[equilibrateState phasesInSystem] objectForKey:@"Liquid"]) ? YES : NO;
    NSUInteger independentVariablesOfSystemPotential = (independentT | independentP);
    if ([[calculationOptions objectForKey:KEY_ORDINATE] isEqualToString:@"Entropy (J/K-kg)"]) { independentVariablesOfSystemPotential &= ~independentT; independentVariablesOfSystemPotential |= independentS; }
    if ([[calculationOptions objectForKey:KEY_ABSCISSA] isEqualToString:@"Volume (cc/kg)"])   { independentVariablesOfSystemPotential &= ~independentP; independentVariablesOfSystemPotential |= independentV; }
    if (self.debugS) NSLog(@"The potential independent variables are T:%@ S:%@ P:%@ V:%@",
                      (independentVariablesOfSystemPotential & independentT) ? @"Yes" : @"No",
                      (independentVariablesOfSystemPotential & independentS) ? @"Yes" : @"No",
                      (independentVariablesOfSystemPotential & independentP) ? @"Yes" : @"No",
                      (independentVariablesOfSystemPotential & independentV) ? @"Yes" : @"No");

    NSUInteger nlInitialGuessIsFeasible = kNlConstraintsSuccess;
    if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentV)) nlInitialGuessIsFeasible = [self computeTemperatureAndPressureFromSystemEntropyAndVolume];
    else if (independentVariablesOfSystemPotential & independentS) nlInitialGuessIsFeasible = [self computeTemperatureFromSystemEntropy];
    else if (independentVariablesOfSystemPotential & independentV) nlInitialGuessIsFeasible = [self computePressureFromSystemVolume];

    //  Attempt via a Kobayashi Maru to construct a feasible guess to the non-linear constraints
    while (nlInitialGuessIsFeasible != kNlConstraintsSuccess) {
        if (nlInitialGuessIsFeasible & kNlConstraintsFailureGeneric) {
            exitStatus = @"failure, Cannot for a valid initial guess for T or P.";
            [self computeAndFillUsefulOutputQuantitiesForEquilibrateStateClass];
            NSDictionary *dictionary = [NSDictionary dictionaryWithObjectsAndKeys:exitStatus, @"status", equilibrateState, @"results", nil];
            return dictionary;

        } else if ( (nlInitialGuessIsFeasible & kNlConstraintsFailureForPressureLowerBound) || (nlInitialGuessIsFeasible & kNlConstraintsFailureForPressureUpperBound)) {
            double signOfCorrectionTerm = (nlInitialGuessIsFeasible & kNlConstraintsFailureForPressureLowerBound) ? -1.0 : 1.0;
            NSUInteger iterations = 0;
            while ((nlInitialGuessIsFeasible & kNlConstraintsFailureForPressureLowerBound) || (nlInitialGuessIsFeasible & kNlConstraintsFailureForPressureUpperBound) || (iterations > 100)) {
                [equilibrateState setCorrectionToReferenceVolumeOfSystem:[equilibrateState correctionToReferenceVolumeOfSystem] + signOfCorrectionTerm*0.1];

                if (self.debugS) NSLog(@"Make an adjustment to the system reference volume to attempt a feasibile initial guess: %g J/bar.", [equilibrateState correctionToReferenceVolumeOfSystem]);
                if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentV)) nlInitialGuessIsFeasible = [self computeTemperatureAndPressureFromSystemEntropyAndVolume];
                else if (independentVariablesOfSystemPotential & independentV) nlInitialGuessIsFeasible = [self computePressureFromSystemVolume];

                iterations++;
            }
            if (iterations > 100) nlInitialGuessIsFeasible &= kNlConstraintsFailureGeneric;

        } else if ( (nlInitialGuessIsFeasible & kNlConstraintsFailureForTemperatureLowerBound) || (nlInitialGuessIsFeasible & kNlConstraintsFailureForTemperatureUpperBound)) {
            double signOfCorrectionTerm = (nlInitialGuessIsFeasible & kNlConstraintsFailureForTemperatureLowerBound) ? 1.0 : -1.0;
            NSUInteger iterations = 0;
            while ((nlInitialGuessIsFeasible & kNlConstraintsFailureForTemperatureLowerBound) || (nlInitialGuessIsFeasible & kNlConstraintsFailureForTemperatureUpperBound) || (iterations > 100)) {
                [equilibrateState setCorrectionToReferenceEntropyOfSystem:[equilibrateState correctionToReferenceEntropyOfSystem] + signOfCorrectionTerm*10.0];

                if (self.debugS) NSLog(@"Make an adjustment to the system reference entropy to attempt a feasibile initial guess: %g J/K.", [equilibrateState correctionToReferenceEntropyOfSystem]);
                if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentV)) nlInitialGuessIsFeasible = [self computeTemperatureAndPressureFromSystemEntropyAndVolume];
                else if (independentVariablesOfSystemPotential & independentS) nlInitialGuessIsFeasible = [self computeTemperatureFromSystemEntropy];

                iterations++;
            }
            if (iterations > 100) nlInitialGuessIsFeasible &= kNlConstraintsFailureGeneric;

        }
    }

    BOOL converged = NO;
	NSString *lastPhaseRemovedFromSystem = @"";
	NSCountedSet *setOfPhasesRemovedFromSystem = [[NSCountedSet alloc] initWithCapacity:1];
	if (!initialPhaseAssemblageRequiresEquilibrationStep) converged = ![self evaluateSaturationState:lastPhaseRemovedFromSystem];

	NSUInteger quadraticIterations = 0;
	double lastPotentialCalculatedInLinearSearch = 0.0;
	double averageLastPotentialCalculatedInLinearSearch = 0.0;
	NSMutableArray *historyPotentialCalculatedInLinearSearch = [NSMutableArray arrayWithCapacity:100];

	while (!converged) {
		if (self.debugS) {
			NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><");
			NSLog(@"Begin Quadratic Iteration Loop %3.3lu (T = %8.2f °C, P = %8.2f Mb)", quadraticIterations, [equilibrateState T]-273.15, [equilibrateState P]/10.0);
			NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><");
		}
		if (quadraticIterations == 0) [historyPotentialCalculatedInLinearSearch removeAllObjects];

        if ([equilibrateState correctionToReferenceVolumeOfSystem] != 0.0) {
            double sumOfComputedDvdp = 0.0;
            double t = [equilibrateState T];
            double p = [equilibrateState P];
            for (NSString *key in [[equilibrateState phasesInSystem] allKeys]) {
                EquilibrateStatePhase *object = [[equilibrateState phasesInSystem] valueForKey:key];
                double dvdp   = [object dvdpFunctionAtT:t andP:p andDeltaMolesOfNonZeroComponents:NULL andScalerMultiplierForCorrection:0.0];
                if (dvdp   != compositionIsUnacceptable) sumOfComputedDvdp += dvdp;
            }
            double pCorrection = -[equilibrateState correctionToReferenceVolumeOfSystem]/sumOfComputedDvdp;
            if (self.debugS) NSLog(@"Pressure in the absence of corection to the reference volume is = %8.2f Mb +corr = %8.2f Mb", p/10.0, (p+pCorrection)/10.0);
            if      (p+pCorrection     > 1.0) [equilibrateState setCorrectionToReferenceVolumeOfSystem:0.0];
            else if (p+pCorrection/2.0 > 1.0) [equilibrateState setCorrectionToReferenceVolumeOfSystem:[equilibrateState correctionToReferenceVolumeOfSystem]/2.0];
            else if (p+pCorrection/4.0 > 1.0) [equilibrateState setCorrectionToReferenceVolumeOfSystem:[equilibrateState correctionToReferenceVolumeOfSystem]*3.0/4.0];
        }

		// create the constraint matrix
		NSArray *nonZeroElements = [equilibrateState hashTableOfNonZeroElementEntries];
		NSUInteger numberOfRowConstraints = [nonZeroElements count];
        if (independentVariablesOfSystemPotential & independentV) numberOfRowConstraints++;
        if (independentVariablesOfSystemPotential & independentS) numberOfRowConstraints++;
		if (self.debugV) NSLog(@"Constraint rows = %lu", numberOfRowConstraints);

		NSUInteger numberOfColumnConstraints = 0;
		NSMutableDictionary *phaseList = [NSMutableDictionary dictionaryWithCapacity:0];
		NSDictionary *phasesAvailable = [equilibrateState phasesInSystem];
		for (NSString *key in [phasesAvailable allKeys]) {
			EquilibrateStatePhase *object = [phasesAvailable valueForKey:key];
			NSArray	*nonZero = [object nonZeroMolesOfEndmemberComponents:independentVariablesOfSystemPotential andT:[equilibrateState T] andP:[equilibrateState P]];
			numberOfColumnConstraints += [nonZero count];
			[phaseList setObject:nonZero forKey:key];
		}
        if (independentVariablesOfSystemPotential & independentV) numberOfColumnConstraints++;
        if (independentVariablesOfSystemPotential & independentS) numberOfColumnConstraints++;
		if (self.debugV) NSLog(@"Constraint columns = %lu", numberOfColumnConstraints);

        DoubleMatrix *cMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:numberOfRowConstraints
                                                           andWithColumnSize:numberOfColumnConstraints
                                                             andInitialValue:0.0];
		double **cMatrix = [cMatrixWrapper pointerToPointerToDouble];
        DoubleVector *dVectorWrapper = [[DoubleVector alloc] initWithSize:numberOfRowConstraints andInitialValue:0.0];
		double  *dVector = [dVectorWrapper pointerToDouble];

		NSUInteger rowIndex = 0;
		double *bulkCompositionInElements = [[equilibrateState bulkCompositionInElements] pointerToDouble];
		for (NSNumber *entry in nonZeroElements) {
			NSUInteger elmNo = [entry unsignedIntValue];
            if ((independentVariablesOfSystemPotential & independentT) && (independentVariablesOfSystemPotential & independentP)) dVector[rowIndex] = bulkCompositionInElements[elmNo];

			NSUInteger colIndex = 0;
			for (NSString *key in [phaseList allKeys]) {
				NSArray *nonZero = [phaseList objectForKey:key];
				for (NSArray *object in nonZero) {
					double *elementsStoichiometry = [[object objectAtIndex:1] pointerValue];
					cMatrix[rowIndex][colIndex] = elementsStoichiometry[elmNo];
					colIndex++;
				}
			}

			rowIndex++;
		}

        // Add non-linear constraint row
        if ( ((independentVariablesOfSystemPotential & independentT) && (independentVariablesOfSystemPotential & independentV)) ||
             ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentP)) ) {
            NSUInteger colIndex = 0;
            for (NSString *key in [phaseList allKeys]) {
                NSArray *nonZero = [phaseList objectForKey:key];
                for (NSArray *object in nonZero) {
                    cMatrix[rowIndex][colIndex] = [[object objectAtIndex:5] doubleValue];
                    colIndex++;
                }
            }
            for (NSString *key in [phaseList allKeys]) {
                NSArray *object = [[phaseList objectForKey:key] objectAtIndex:0];
                cMatrix[rowIndex][colIndex] += [[object objectAtIndex:9] doubleValue];
            }
        } else if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentV)) {
            NSUInteger colIndex = 0;
            for (NSString *key in [phaseList allKeys]) {
                NSArray *nonZero = [phaseList objectForKey:key];
                for (NSArray *object in nonZero) {
                    cMatrix[rowIndex  ][colIndex] = [[[object objectAtIndex:5] objectAtIndex:0] doubleValue];
                    cMatrix[rowIndex+1][colIndex] = [[[object objectAtIndex:5] objectAtIndex:1] doubleValue];
                    colIndex++;
                }
            }
            for (NSString *key in [phaseList allKeys]) {
                NSArray *object = [[phaseList objectForKey:key] objectAtIndex:0];
                cMatrix[rowIndex  ][colIndex  ] += [[[object objectAtIndex:9] objectAtIndex:0] doubleValue];
                cMatrix[rowIndex  ][colIndex+1] += [[[object objectAtIndex:9] objectAtIndex:1] doubleValue];
                cMatrix[rowIndex+1][colIndex  ] += [[[object objectAtIndex:9] objectAtIndex:2] doubleValue];
                cMatrix[rowIndex+1][colIndex+1] += [[[object objectAtIndex:9] objectAtIndex:3] doubleValue];
            }
        }

        // Store matrix transpose for Lagrange multiplier extraction
        DoubleMatrix *cTMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:numberOfColumnConstraints
                                                            andWithColumnSize:numberOfRowConstraints
                                                              andInitialValue:0.0];
        double **cTMatrix = [cTMatrixWrapper pointerToPointerToDouble];
        if ((independentVariablesOfSystemPotential & independentV) || (independentVariablesOfSystemPotential & independentS)) {
            for (NSUInteger i=0; i<numberOfRowConstraints; i++) {
                for (NSUInteger j=0; j<numberOfColumnConstraints; j++) cTMatrix[j][i] = cMatrix[i][j];
            }
        }

        // Perform an orthogonal decomposition of the constraint matrix
        DoubleVector *hVectorWrapper = [[DoubleVector alloc] initWithSize:numberOfRowConstraints];
		double *hVector = [hVectorWrapper pointerToDouble];
        DoubleVector *yVectorWrapper = [[DoubleVector alloc] initWithSize:numberOfRowConstraints];
		double *yVector = [yVectorWrapper pointerToDouble];
		// Householder decompose the constraint matrix
		for (NSUInteger i=0; i<numberOfRowConstraints; i++) [MathSupport householderRowRow:HOUSEHOLDER_CALC_MODE_H1
																						 p:i
																						 l:i+1
																						 m:numberOfColumnConstraints-1
																						 v:cMatrix
																					  vRow:i
																						 h:&hVector[i]
																						 c:cMatrix
																				 cRowStart:i+1
																				   cRowEnd:numberOfRowConstraints-1];
		// back-solve the lower-triangular system and store the soln in yVector
		for (NSUInteger i=0; i<numberOfRowConstraints; i++) {
			yVector[i] = dVector[i];
			for (NSUInteger j=0; j<i; j++) yVector[i] -= cMatrix[i][j]*yVector[j];
			if (fabs(cMatrix[i][i]) > 10.0*DBL_EPSILON) {  // was cMatrix[i][i] != 0.0
				yVector[i] /= cMatrix[i][i];
                if (self.debugV) NSLog(@"... yVector[%lu] = %20.13e cMatrix[%lu][%lu] = %20.13e", i, yVector[i], i, i, cMatrix[i][i]);
			} else {
				numberOfRowConstraints--;
				if (self.debugS) NSLog(@"...Constraint matrix row %lu is redundant.  numberOfRowConstraints reset to %lu", i, numberOfRowConstraints);
			}
		}
		if (self.debugV) NSLog(@"Constraint phase completed.");

        DoubleMatrix *eMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:numberOfColumnConstraints
                                                           andWithColumnSize:numberOfColumnConstraints
                                                             andInitialValue:0.0];
		double **eMatrix = [eMatrixWrapper pointerToPointerToDouble];
        DoubleMatrix *bMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:numberOfColumnConstraints
                                                           andWithColumnSize:1
                                                             andInitialValue:0.0];
		double **bMatrix = [bMatrixWrapper pointerToPointerToDouble];
        DoubleVector *referenceMolesWrapper = [[DoubleVector alloc] initWithSize:numberOfColumnConstraints];
		double *referenceMoles = [referenceMolesWrapper pointerToDouble];
        DoubleMatrix *eLambdaMatrixWrapper, *eLambdaSMatrixWrapper, *eLambdaVMatrixWrapper;
        double **eLambdaMatrix=NULL, **eLambdaSMatrix=NULL, **eLambdaVMatrix=NULL;

        if ( ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentP)) ||
             ((independentVariablesOfSystemPotential & independentT) && (independentVariablesOfSystemPotential & independentV)) ) {
            eLambdaMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:numberOfColumnConstraints
                                                 andWithColumnSize:numberOfColumnConstraints
                                                   andInitialValue:0.0];
            eLambdaMatrix = [eLambdaMatrixWrapper pointerToPointerToDouble];
        } else if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentV)) {
            eLambdaSMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:numberOfColumnConstraints
                                                       andWithColumnSize:numberOfColumnConstraints
                                                         andInitialValue:0.0];
            eLambdaSMatrix = [eLambdaSMatrixWrapper pointerToPointerToDouble];
            eLambdaVMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:numberOfColumnConstraints
                                                        andWithColumnSize:numberOfColumnConstraints
                                                          andInitialValue:0.0];
            eLambdaVMatrix = [eLambdaVMatrixWrapper pointerToPointerToDouble];
        }

		NSUInteger colIndex = 0;
		rowIndex = 0;
		for (NSString *key in [phaseList allKeys]) {
			NSArray *nonZero = [phaseList objectForKey:key];
			NSUInteger colBase = colIndex;

			for (NSArray *object in nonZero) {
				BOOL hessianCompositionalEntry = [[[phasesAvailable objectForKey:key] phaseClassInstance] conformsToProtocol:@protocol(SolutionPhaseProtocol)];
				referenceMoles[rowIndex] = [[object objectAtIndex:0] doubleValue]; // moles of component/phase
				bMatrix[rowIndex][0] = -[[object objectAtIndex:2] doubleValue];    // - dp

				if (hessianCompositionalEntry) {
					NSArray *hessian = [object objectAtIndex:3];
					colIndex = colBase;
					for (NSNumber *entry in hessian) {
						eMatrix[rowIndex][colIndex] = [entry doubleValue];
						colIndex++;
					}
                    if ( ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentP)) ||
                         ((independentVariablesOfSystemPotential & independentT) && (independentVariablesOfSystemPotential & independentV)) ) {
                        NSArray *hessian = [object objectAtIndex:6];
                        colIndex = colBase;
                        for (NSNumber *entry in hessian) {
                            eLambdaMatrix[rowIndex][colIndex] = [entry doubleValue];
                            colIndex++;
                        }
                    } else if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentV)) {
                        NSArray *hessianS = [[object objectAtIndex:6] objectAtIndex:0];
                        colIndex = colBase;
                        for (NSNumber *entry in hessianS) {
                            eLambdaSMatrix[rowIndex][colIndex] = [entry doubleValue];
                            colIndex++;
                        }
                        NSArray *hessianV = [[object objectAtIndex:6] objectAtIndex:1];
                        colIndex = colBase;
                        for (NSNumber *entry in hessianV) {
                            eLambdaVMatrix[rowIndex][colIndex] = [entry doubleValue];
                            colIndex++;
                        }
                    }
				} else colIndex++;

                if ( ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentP)) ||
                     ((independentVariablesOfSystemPotential & independentT) && (independentVariablesOfSystemPotential & independentV)) ) {
                    eMatrix[rowIndex][numberOfColumnConstraints-1] = [[object objectAtIndex:4] doubleValue];
                    eMatrix[numberOfColumnConstraints-1][rowIndex] = eMatrix[rowIndex][numberOfColumnConstraints-1];

                    eLambdaMatrix[rowIndex][numberOfColumnConstraints-1] = [[object objectAtIndex:7] doubleValue];
                    eLambdaMatrix[numberOfColumnConstraints-1][rowIndex] = eLambdaMatrix[rowIndex][numberOfColumnConstraints-1];
                } else if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentV)) {
                    eMatrix[rowIndex][numberOfColumnConstraints-2] = [[[object objectAtIndex:4] objectAtIndex:0] doubleValue];
                    eMatrix[rowIndex][numberOfColumnConstraints-1] = [[[object objectAtIndex:4] objectAtIndex:1] doubleValue];
                    eMatrix[numberOfColumnConstraints-2][rowIndex] = eMatrix[rowIndex][numberOfColumnConstraints-2];
                    eMatrix[numberOfColumnConstraints-1][rowIndex] = eMatrix[rowIndex][numberOfColumnConstraints-1];

                    eLambdaSMatrix[rowIndex][numberOfColumnConstraints-2] = [[[object objectAtIndex:7] objectAtIndex:0] doubleValue];
                    eLambdaSMatrix[rowIndex][numberOfColumnConstraints-1] = [[[object objectAtIndex:7] objectAtIndex:1] doubleValue];
                    eLambdaVMatrix[rowIndex][numberOfColumnConstraints-2] = [[[object objectAtIndex:7] objectAtIndex:2] doubleValue];
                    eLambdaVMatrix[rowIndex][numberOfColumnConstraints-1] = [[[object objectAtIndex:7] objectAtIndex:3] doubleValue];

                    eLambdaSMatrix[numberOfColumnConstraints-2][rowIndex] = eLambdaSMatrix[rowIndex][numberOfColumnConstraints-2];
                    eLambdaSMatrix[numberOfColumnConstraints-1][rowIndex] = eLambdaSMatrix[rowIndex][numberOfColumnConstraints-1];
                    eLambdaVMatrix[numberOfColumnConstraints-2][rowIndex] = eLambdaVMatrix[rowIndex][numberOfColumnConstraints-2];
                    eLambdaVMatrix[numberOfColumnConstraints-1][rowIndex] = eLambdaVMatrix[rowIndex][numberOfColumnConstraints-1];
                }

				rowIndex++;
			}

            if ( ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentP)) ||
                 ((independentVariablesOfSystemPotential & independentT) && (independentVariablesOfSystemPotential & independentV)) ) {
                NSArray *object = [nonZero objectAtIndex:0];
                bMatrix[numberOfColumnConstraints-1][0] -= [[object objectAtIndex:8] doubleValue];  // -dp
                eMatrix[numberOfColumnConstraints-1][numberOfColumnConstraints-1] += [[object objectAtIndex:10] doubleValue];
                eLambdaMatrix[numberOfColumnConstraints-1][numberOfColumnConstraints-1] += [[object objectAtIndex:11] doubleValue];

            } else if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentV)) {
                NSArray *object = [nonZero objectAtIndex:0];

                bMatrix[numberOfColumnConstraints-2][0] -= [[[object objectAtIndex:8] objectAtIndex:0] doubleValue];  // -dp
                bMatrix[numberOfColumnConstraints-1][0] -= [[[object objectAtIndex:8] objectAtIndex:1] doubleValue];  // -dp

                eMatrix[numberOfColumnConstraints-2][numberOfColumnConstraints-2] += [[[object objectAtIndex:10] objectAtIndex:0] doubleValue];
                eMatrix[numberOfColumnConstraints-2][numberOfColumnConstraints-1] += [[[object objectAtIndex:10] objectAtIndex:1] doubleValue];
                eMatrix[numberOfColumnConstraints-1][numberOfColumnConstraints-2] += [[[object objectAtIndex:10] objectAtIndex:1] doubleValue];
                eMatrix[numberOfColumnConstraints-1][numberOfColumnConstraints-1] += [[[object objectAtIndex:10] objectAtIndex:2] doubleValue];

                eLambdaSMatrix[numberOfColumnConstraints-2][numberOfColumnConstraints-2] += [[[object objectAtIndex:11] objectAtIndex:0] doubleValue];
                eLambdaSMatrix[numberOfColumnConstraints-2][numberOfColumnConstraints-1] += [[[object objectAtIndex:11] objectAtIndex:1] doubleValue];
                eLambdaSMatrix[numberOfColumnConstraints-1][numberOfColumnConstraints-2] += [[[object objectAtIndex:11] objectAtIndex:1] doubleValue];
                eLambdaSMatrix[numberOfColumnConstraints-1][numberOfColumnConstraints-1] += [[[object objectAtIndex:11] objectAtIndex:2] doubleValue];

                eLambdaVMatrix[numberOfColumnConstraints-2][numberOfColumnConstraints-2] += [[[object objectAtIndex:11] objectAtIndex:3] doubleValue];
                eLambdaVMatrix[numberOfColumnConstraints-2][numberOfColumnConstraints-1] += [[[object objectAtIndex:11] objectAtIndex:4] doubleValue];
                eLambdaVMatrix[numberOfColumnConstraints-1][numberOfColumnConstraints-2] += [[[object objectAtIndex:11] objectAtIndex:4] doubleValue];
                eLambdaVMatrix[numberOfColumnConstraints-1][numberOfColumnConstraints-1] += [[[object objectAtIndex:11] objectAtIndex:5] doubleValue];
            }

		}
        if (independentVariablesOfSystemPotential & independentS) { rowIndex++; colIndex++; }
        if (independentVariablesOfSystemPotential & independentV) { rowIndex++; colIndex++; }
		if (self.debugS) {
			NSLog(@"Quadratic search direction matrix constructed dimension check: %lu = %lu = %lu", rowIndex, numberOfColumnConstraints, colIndex);
			NSLog(@"Structure of eMatrix, bMatrix[][0]:");
		}
		for (NSUInteger i=0; i<numberOfColumnConstraints; i++) {
			NSMutableString *rowStructure = [NSMutableString stringWithCapacity:0];
			for (NSUInteger j=0; j<numberOfColumnConstraints; j++) {
				if (eMatrix[i][j] != 0.0) [rowStructure appendString:@"x"]; else [rowStructure appendString:@"0"];
			}
			if (self.debugS) NSLog(@"eMatrix[%2.2lu][] = %@ bMatrix[%2.2lu][0] = %f", i, rowStructure, i, bMatrix[i][0]);
		}

        //Compute Lagrange multipliers if the constraints are non-linear
        double lambdaS = 0.0, lambdaV = 0.0;
        if ((independentVariablesOfSystemPotential & independentV) || (independentVariablesOfSystemPotential & independentS)) {
            DoubleMatrix *gMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:numberOfColumnConstraints
                                                                 andWithColumnSize:1
                                                                   andInitialValue:0.0];
            double **gMatrix = [gMatrixWrapper pointerToPointerToDouble];
            for (NSUInteger i=0; i<numberOfColumnConstraints; i++) gMatrix[i][0] = -bMatrix[i][0];

            DoubleVector *hVectorWrapper = [[DoubleVector alloc] initWithSize:numberOfRowConstraints];
			double *hVector = [hVectorWrapper pointerToDouble];
            DoubleVector *gVectorWrapper = [[DoubleVector alloc] initWithSize:numberOfRowConstraints];
			double *gVector = [gVectorWrapper pointerToDouble];
            IntegerVector *pVectorWrapper = [[IntegerVector alloc] initWithSize:numberOfRowConstraints andInitialValue:0];
			NSInteger *pVector = [pVectorWrapper pointerToInteger];
            double tolerance = DBL_EPSILON, rNorm;
			NSInteger pseudoRank;
			[MathSupport hfti:cTMatrix
							m:numberOfColumnConstraints
							n:numberOfRowConstraints
							b:gMatrix
						   nb:1
						  tau:tolerance
							k:&pseudoRank
						rnorm:&rNorm
							h:hVector
							g:gVector
							p:pVector];
            if (self.debugS) {
				NSLog(@"Lagrange HFTI step completed. Dimension = %lu, pseudo rank = %ld, rNorm = %g", numberOfRowConstraints, pseudoRank, rNorm);
                NSUInteger i = 0;
                for (NSNumber *entry in nonZeroElements) {
                    NSUInteger elmNo = [entry unsignedIntValue];
                    NSLog(@"Lagrange multiplier[%@] = %13.6g", [[PhaseBase elementNameFromAtomicNumber:elmNo] stringByPaddingToLength:2 withString:@" " startingAtIndex:0], gMatrix[i][0]);
                    i++;
                }
			}


            NSUInteger i = numberOfRowConstraints - 1;
            if (independentVariablesOfSystemPotential & independentV)  lambdaV  = gMatrix[i--][0];
            if (independentVariablesOfSystemPotential & independentS)  lambdaS  = gMatrix[i--][0];
            if (self.debugS) NSLog(@"Lagrange multiplier[*S] = %20.13e", lambdaS);
            if (self.debugS) NSLog(@"Lagrange multiplier[*V] = %20.13e", lambdaV);
        }

        if        ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentP)) {
            for (NSUInteger i=0; i<numberOfColumnConstraints; i++) for (NSUInteger j=0; j<numberOfColumnConstraints; j++) eMatrix[i][j] -= lambdaS*eLambdaMatrix[i][j];
        } else if ((independentVariablesOfSystemPotential & independentT) && (independentVariablesOfSystemPotential & independentV)) {
            for (NSUInteger i=0; i<numberOfColumnConstraints; i++) for (NSUInteger j=0; j<numberOfColumnConstraints; j++) eMatrix[i][j] -= lambdaV*eLambdaMatrix[i][j];
        } else if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentV)) {
            for (NSUInteger i=0; i<numberOfColumnConstraints; i++) for (NSUInteger j=0; j<numberOfColumnConstraints; j++) eMatrix[i][j] -= lambdaS*eLambdaSMatrix[i][j] + lambdaV*eLambdaVMatrix[i][j];
        }

		//Compute e~ = (K^^T) E
		for (NSUInteger i=0; i<numberOfRowConstraints; i++) [MathSupport householderRowCol:HOUSEHOLDER_CALC_MODE_H2
																						 p:i
																						 l:i+1
																						 m:numberOfColumnConstraints-1
																						 v:cMatrix
																					  vRow:i
																						 h:&hVector[i]
																						 c:eMatrix
																				 cColStart:0
																				   cColEnd:numberOfColumnConstraints-1];

		// Form the last conCols - conRows of e^ = e~ K, i.e. compute e21^ = K2^^T E K1  and  e22^ = K2^^T E K2
		for (NSUInteger i=0; i<numberOfRowConstraints; i++) [MathSupport householderRowRow:HOUSEHOLDER_CALC_MODE_H2
																						 p:i
																						 l:i+1
																						 m:numberOfColumnConstraints-1
																						 v:cMatrix
																					  vRow:i
																						 h:&hVector[i]
																						 c:eMatrix
																				 cRowStart:numberOfRowConstraints
																				   cRowEnd:numberOfColumnConstraints-1];
		//Compute b ~ = (K^^T) B
		for (NSUInteger i=0; i<numberOfRowConstraints; i++) [MathSupport householderRowCol:HOUSEHOLDER_CALC_MODE_H2
																						 p:i
																						 l:i+1
																						 m:numberOfColumnConstraints-1
																						 v:cMatrix
																					  vRow:i
																						 h:&hVector[i]
																						 c:bMatrix
																				 cColStart:0
																				   cColEnd:0];
		//Compute b2^ = - b2~ - e21^ Y1^
		for (NSUInteger i=numberOfRowConstraints; i<numberOfColumnConstraints; i++)
			for (NSUInteger j=0; j<numberOfRowConstraints; j++) bMatrix[i][0] -= eMatrix[i][j]*yVector[j];
		// The least squares problem is now stored in the matrix and vector
		// e[numberOfRowConstraints:numberOfColumnConstraints-1][numberOfRowConstraints:numberOfColumnConstraints-1]
		// and b[numberOfRowConstraints:numberOfColumnConstraints-1][0]

		if (numberOfRowConstraints < numberOfColumnConstraints) {
            DoubleMatrix *aMatrixWrapper = [[DoubleMatrix alloc] initWithRowSize:numberOfColumnConstraints-numberOfRowConstraints
                                                               andWithColumnSize:numberOfColumnConstraints-numberOfRowConstraints];
			double **aMatrix = [aMatrixWrapper pointerToPointerToDouble];
            DoubleVector *hVectorWrapper = [[DoubleVector alloc] initWithSize:numberOfColumnConstraints-numberOfRowConstraints];
			double *hVector = [hVectorWrapper pointerToDouble];
            DoubleVector *gVectorWrapper = [[DoubleVector alloc] initWithSize:numberOfColumnConstraints-numberOfRowConstraints];
			double *gVector = [gVectorWrapper pointerToDouble];
            IntegerVector *pVectorWrapper = [[IntegerVector alloc] initWithSize:numberOfColumnConstraints-numberOfRowConstraints andInitialValue:0];
			NSInteger *pVector = [pVectorWrapper pointerToInteger];

			double fNorm = 0.0;
			for (NSUInteger i=0; i<(numberOfColumnConstraints-numberOfRowConstraints); i++) {
				for (NSUInteger j=0; j<(numberOfColumnConstraints-numberOfRowConstraints); j++) {
					aMatrix[i][j] = eMatrix[i+numberOfRowConstraints][j+numberOfRowConstraints];
					fNorm += aMatrix[i][j]*aMatrix[i][j];
				}
			}
			fNorm = sqrt(fNorm);
			double tolerance = self.PPPARAMETERS_RANK.doubleValue*fNorm, rNorm;
			NSInteger pseudoRank;
			[MathSupport hfti:aMatrix
							m:numberOfColumnConstraints-numberOfRowConstraints
							n:numberOfColumnConstraints-numberOfRowConstraints
							b:&bMatrix[numberOfRowConstraints]
						   nb:1
						  tau:tolerance
							k:&pseudoRank
						rnorm:&rNorm
							h:hVector
							g:gVector
							p:pVector];
			if (self.debugS) {
				NSLog(@"HFTI completed. Dimension = %lu, pseudo rank = %ld, fNorm = %g, rNorm = %g", numberOfColumnConstraints-numberOfRowConstraints,
					  pseudoRank, fNorm, rNorm);
				for (NSUInteger i=0; i<(numberOfColumnConstraints-numberOfRowConstraints); i++) NSLog(@"... diagonal element[%2.2lu][%2.2lu] = %13.6g",
																									  i, i, aMatrix[i][i]);
				NSLog(@"... ratio min/max diagonal elements = %13.6g, DBL_EPSILON = %13.6g",
					  aMatrix[numberOfColumnConstraints-numberOfRowConstraints-1][numberOfColumnConstraints-numberOfRowConstraints-1]/aMatrix[0][0], DBL_EPSILON);
			}
			if (self.PPOPTIONS_RANK.boolValue) {
				if ((numberOfColumnConstraints-numberOfRowConstraints) > pseudoRank) {
                    [self computeAndFillUsefulOutputQuantitiesForEquilibrateStateClass];
					exitStatus = @"failure, The hessian is rank deficient.";
					NSDictionary *dictionary = [NSDictionary dictionaryWithObjectsAndKeys:exitStatus, @"status", equilibrateState, @"results", nil];
					return dictionary;
				}
			}
		}

		for (NSUInteger i=0; i<numberOfRowConstraints; i++) bMatrix[i][0] = yVector[i];
		for (NSUInteger i=0; i<numberOfRowConstraints; i++) [MathSupport householderRowCol:HOUSEHOLDER_CALC_MODE_H2
																						 p:numberOfRowConstraints-1-i
																						 l:numberOfRowConstraints-i
																						 m:numberOfColumnConstraints-1
																						 v:cMatrix
																					  vRow:numberOfRowConstraints-1-i
																						 h:&hVector[numberOfRowConstraints-1-i]
																						 c:bMatrix
																				 cColStart:0
																				   cColEnd:0];
        DoubleVector *deltaMolesWrapper = [[DoubleVector alloc] initWithSize:numberOfColumnConstraints];
		double *deltaMoles = [deltaMolesWrapper pointerToDouble];
		double rNorm = 0.0, sNorm = 0.0, deltaT = 0.0, deltaP = 0.0;
        if ((independentVariablesOfSystemPotential & independentT) && (independentVariablesOfSystemPotential & independentP)) {
            for (NSUInteger i=0; i<numberOfColumnConstraints; i++) {
                deltaMoles[i] = bMatrix[i][0] - referenceMoles[i];
                sNorm += bMatrix[i][0]*bMatrix[i][0];
                rNorm += deltaMoles[i]*deltaMoles[i];
            }
        } else if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentP)) {
            for (NSUInteger i=0; i<(numberOfColumnConstraints-1); i++) {
                deltaMoles[i] = bMatrix[i][0];
                sNorm += pow(bMatrix[i][0]+referenceMoles[i], 2.0);
                rNorm += pow(deltaMoles[i], 2.0);
            }
            deltaT = bMatrix[numberOfColumnConstraints-1][0];
            sNorm += pow([equilibrateState T]+deltaT, 2.0);
            double scale = 1000.0;
            if ((scale = [equilibrateState T]) < 500.0) scale = 1000.0;
            rNorm += pow(deltaT/scale, 2.0);
        } else if ((independentVariablesOfSystemPotential & independentT) && (independentVariablesOfSystemPotential & independentV)) {
            for (NSUInteger i=0; i<(numberOfColumnConstraints-1); i++) {
                deltaMoles[i] = bMatrix[i][0];
                sNorm += pow(bMatrix[i][0]+referenceMoles[i], 2.0);
                rNorm += pow(deltaMoles[i], 2.0);
            }
            deltaP = bMatrix[numberOfColumnConstraints-1][0];
            sNorm += pow([equilibrateState P]+deltaP, 2.0);
            double scale = 1000.0;
            if ((scale = [equilibrateState P]) < 100.0) scale = 1000.0;
            rNorm += pow(deltaP/scale, 2.0);
        } else if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentV)) {
            for (NSUInteger i=0; i<(numberOfColumnConstraints-2); i++) {
                deltaMoles[i] = bMatrix[i][0];
                sNorm += pow(bMatrix[i][0]+referenceMoles[i], 2.0);
                rNorm += pow(deltaMoles[i], 2.0);
            }
            deltaT = bMatrix[numberOfColumnConstraints-2][0];
            sNorm += pow([equilibrateState T]+deltaT, 2.0);
            double scale = 1000.0;
            if ((scale = [equilibrateState T]) < 500.0) scale = 1000.0;
            rNorm += pow(deltaT/scale, 2.0);
            deltaP = bMatrix[numberOfColumnConstraints-1][0];
            sNorm += pow([equilibrateState P]+deltaP, 2.0);
            if ((scale = [equilibrateState P]) < 100.0) scale = 1000.0;
            rNorm += pow(deltaP/scale, 2.0);
        }

		rNorm = sqrt(rNorm);
		sNorm = sqrt(sNorm);
		quadraticIterations++;

		if (self.debugS) {
			NSLog(@"HFTI solution with sNorm = %g and rNorm = %g:", sNorm, rNorm);
			NSUInteger columnIndex = 0;
			for (NSString *key in [phaseList allKeys]) {
				NSUInteger count = [[phaseList objectForKey:key] count];
				EquilibrateStatePhase *phaseWrapper = [[equilibrateState phasesInSystem] objectForKey:key];
				id phase = [phaseWrapper phaseClassInstance];
				NSArray *indices = [phaseWrapper indicesOfNonZeroMolesOfEndmemberComponents];

				for (NSUInteger i=0; i<count; i++, columnIndex++) {
					NSString *component;
					if ([phase conformsToProtocol:@protocol(SolutionPhaseProtocol)])
						component = [[phase componentAtIndex:[[indices objectAtIndex:i] unsignedIntValue]] phaseName];
					else component = [NSString stringWithString:key];
					NSLog(@"... soln[%@][%@] = %13.6g, reference moles = %13.6g, delta = %13.6g",
						  [key stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
						  [component stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
                          ((independentVariablesOfSystemPotential & independentS) || (independentVariablesOfSystemPotential & independentV)) ? bMatrix[columnIndex][0]+referenceMoles[columnIndex] : bMatrix[columnIndex][0],
						  referenceMoles[columnIndex], deltaMoles[columnIndex]);
				}
			}
            if (independentVariablesOfSystemPotential & independentS) {
                NSLog(@"... soln[%@][%@] = %13.6g, reference T     = %13.6g, delta = %13.6g", @"T (°C)         ", @"               ", [equilibrateState T]+deltaT-273.15, [equilibrateState T]-273.15, deltaT);
            }
            if (independentVariablesOfSystemPotential & independentV) {
                NSLog(@"... soln[%@][%@] = %13.6g, reference P     = %13.6g, delta = %13.6g", @"P (MPa)        ", @"               ", ([equilibrateState P]+deltaP)/10.0, [equilibrateState P]/10.0, deltaP/10.0);
            }
		}
		/*
		if (quadraticIterations == 1) {
			NSLog(@" ");
			NSLog(@"Quad iter = %3.3d, rNorm = %20.13e", quadraticIterations, rNorm);
		} else {
			NSLog(@"Quad iter = %3.3d, rNorm = %20.13e", quadraticIterations, rNorm);
		}
		*/
		if (rNorm < self.PPPARAMETERS_QUAD_OPTIMAL.doubleValue*sNorm) {
			converged = YES;
			exitStatus = @"success, Optimal residual norm.";
			if (self.debugS) {
				NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><>");
				NSLog(@"Quadratic loop convergence with optimal rNorm.");
				NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><>");
			}
		} else if ((rNorm < self.PPPARAMETERS_QUAD_SUBOPTIMAL.doubleValue)
				   && (quadraticIterations > self.PPPARAMETERS_QUAD_MAX_ITERS.unsignedIntValue)) {
			converged = YES;
			exitStatus = @"success, Residual norm less than optimal.";
			if (self.debugS) {
				NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><><>");
				NSLog(@"Quadratic loop convergence with sub-optimal rNorm.");
				NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><><>");
			}
		} else if (quadraticIterations > self.PPPARAMETERS_QUAD_MAX_ITERS.unsignedIntValue) {
			converged = YES;
			exitStatus = [NSString localizedStringWithFormat:@"failure, Quadratic algorithm did not converge in %d iterations (rNorm: %g).",
						  self.PPPARAMETERS_QUAD_MAX_ITERS.unsignedIntValue,
						  rNorm];
			//exitStatus = [NSString stringWithString:@"failure, Quadratic algorithm did not converge in 100 iterations."];
			if (self.debugS) {
				NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><><>");
				NSLog(@"Quadratic loop did not converge in 100 iterations.");
				NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><><>");
			}
		} else if (numberOfRowConstraints == numberOfColumnConstraints) {
			converged = YES;
			exitStatus = @"success, No search direction. Row and column constraints are equal.";
			if (self.debugS) {
				NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><>");
				NSLog(@"Quadratic loop convergence. No search direction.");
				NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><>");
			}
		} else {
			// Do a linear search step
			double lambda = 1.0, stepSize = 0.10, pTotal = 0.0;
			NSInteger iter, status;

			do {
				double reltest = self.PPPARAMETERS_LINEAR_RELTEST.doubleValue;
				iter = self.PPPARAMETERS_LINEAR_MAX_ITERS.unsignedIntValue;
				status = [MathSupport min1d:&lambda st:&stepSize reltest:reltest ifn:&iter fnminval:&pTotal fn1d:^(double scaleLength, NSInteger *infeasible) {
					if ((scaleLength < -2.0) || (scaleLength > 2.0)) {
						*infeasible = TRUE;
						return 0.0;
					}

                    if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentV)) {
                        double tOld = [equilibrateState T];
                        double pOld = [equilibrateState P];
                        double tNew = tOld + scaleLength*deltaT;
                        double pNew = pOld + scaleLength*deltaP;
                        double referenceEntropy = [equilibrateState referenceEntropyOfSystem] + [equilibrateState correctionToReferenceEntropyOfSystem];
                        double referenceVolume  = [equilibrateState referenceVolumeOfSystem]  + [equilibrateState correctionToReferenceVolumeOfSystem];
                        double sumOfComputedEntropies = 0.0, sumOfComputedVolumes = 0.0;
                        NSUInteger iterations = 0;

                        if (self.debugS) NSLog(@"*.*.*.*.      Entry: S = %g V = %g, T = %g P = %g", referenceEntropy, referenceVolume, tNew-273.15, pNew);
                        do {
                            sumOfComputedEntropies = 0.0;
                            sumOfComputedVolumes   = 0.0;
                            double sumOfComputedHeatCapacity = 0.0, sumOfComputedDvdt = 0.0, sumOfComputedDvdp = 0.0;
                            NSUInteger columnIndex = 0;
                            *infeasible = FALSE;

                            for (NSString *key in [phaseList allKeys]) {
                                EquilibrateStatePhase *object = [phasesAvailable valueForKey:key];
                                double entropy = [object entropyFunctionAtT:tNew
                                                                       andP:pNew
                                           andDeltaMolesOfNonZeroComponents:&deltaMoles[columnIndex]
                                           andScalerMultiplierForCorrection:scaleLength];
                                if (entropy == compositionIsUnacceptable) { *infeasible = TRUE; return 0.0; }
                                double volume = [object volumeFunctionAtT:tNew
                                                                     andP:pNew
                                         andDeltaMolesOfNonZeroComponents:&deltaMoles[columnIndex]
                                         andScalerMultiplierForCorrection:scaleLength];
                                if (volume == compositionIsUnacceptable) { *infeasible = TRUE; return 0.0; }
                                double heatCapacity = [object heatCapacityFunctionAtT:tNew
                                                                                 andP:pNew
                                                     andDeltaMolesOfNonZeroComponents:&deltaMoles[columnIndex]
                                                     andScalerMultiplierForCorrection:scaleLength];
                                if (heatCapacity == compositionIsUnacceptable) { *infeasible = TRUE; return 0.0; }
                                double dvdp = [object dvdpFunctionAtT:tNew
                                                                 andP:pNew
                                     andDeltaMolesOfNonZeroComponents:&deltaMoles[columnIndex]
                                     andScalerMultiplierForCorrection:scaleLength];
                                if (dvdp == compositionIsUnacceptable) { *infeasible = TRUE; return 0.0; }
                                double dvdt = [object dvdtFunctionAtT:tNew
                                                                 andP:pNew
                                     andDeltaMolesOfNonZeroComponents:&deltaMoles[columnIndex]
                                     andScalerMultiplierForCorrection:scaleLength];
                                if (dvdt == compositionIsUnacceptable) { *infeasible = TRUE; return 0.0; }
                                if (self.debugV) NSLog(@"*.*.*.*.      Phase:%@ S = %g V = %g Cp = %g dvdp = %g dvdt = %g", [[object phaseClassInstance] phaseName], entropy, volume, heatCapacity, dvdp, dvdt);
                                sumOfComputedEntropies    += entropy;
                                sumOfComputedVolumes      += volume;
                                sumOfComputedHeatCapacity += heatCapacity;
                                sumOfComputedDvdp         += dvdp;
                                sumOfComputedDvdt         += dvdt;
                                columnIndex += [[phaseList objectForKey:key] count];
                            }
                            if (self.debugV) NSLog(@"*.*.*.*.      Iter = %lu S %20.13e V %20.13e Cp %20.13e dvdp %20.13e dvdt %20.13e",
                                              iterations, sumOfComputedEntropies, sumOfComputedVolumes, sumOfComputedHeatCapacity, sumOfComputedDvdp, sumOfComputedDvdt);

                            double d11 = sumOfComputedHeatCapacity/tNew;
                            double d12 = -sumOfComputedDvdt;
                            double d21 = sumOfComputedDvdt;
                            double d22 = sumOfComputedDvdp;
                            double determ = d11*d22 - d12*d21;
                            tNew -= (  (sumOfComputedEntropies-referenceEntropy)*d22 - (sumOfComputedVolumes-referenceVolume)*d12 )/determ;
                            pNew -= ( -(sumOfComputedEntropies-referenceEntropy)*d21 + (sumOfComputedVolumes-referenceVolume)*d11 )/determ;
                            if (self.debugV) NSLog(@"*.*.*.*.      T = %g P = %g", tNew-273.15, pNew);
                            if      (tNew <    300.0) tNew =    300.0;
                            else if (tNew >   5000.0) tNew =   5000.0;
                            if      (pNew <      1.0) pNew =      1.0; // was -10000.0
                            else if (pNew > 100000.0) pNew = 100000.0;

                            iterations++;
                        } while ((iterations < 100) &&
                                 (fabs(sumOfComputedVolumes-referenceVolume) > 1.0e-9*referenceVolume) &&
                                 (fabs(sumOfComputedEntropies-referenceEntropy) > 1.0e-9*referenceEntropy));
                        if (self.debugS) {
                            NSLog(@"...Iterations %lu", iterations);
                            NSLog(@"...Vref = %g V error = %g %% = %g", referenceVolume, sumOfComputedVolumes-referenceVolume, 100.0*(sumOfComputedVolumes-referenceVolume)/referenceVolume);
                            NSLog(@"...Sref = %g S error = %g %% = %g", referenceEntropy, sumOfComputedEntropies-referenceEntropy, 100.0*(sumOfComputedEntropies-referenceEntropy)/referenceEntropy);
                        }
                        if (iterations == 100) {
                            if (self.debugS) NSLog(@"*.*.*.*.      Exit: Infeasible S,V solution");
                            *infeasible = TRUE;
                            return 0.0;
                        }

                        [equilibrateState setT:tNew];
                        [equilibrateState setP:pNew];
                        if (self.debugS) NSLog(@"*.*.*.*.      Exit: T = %g and P = %g", tNew-273.15, pNew);

                    } else if (independentVariablesOfSystemPotential & independentS) {
                        double tOld = [equilibrateState T];
                        double pOld = [equilibrateState P];
                        double tNew = tOld + scaleLength*deltaT;
                        double referenceEntropy = [equilibrateState referenceEntropyOfSystem] + [equilibrateState correctionToReferenceEntropyOfSystem];
                        NSUInteger iterations = 0;
                        double sumOfComputedEntropies = 0.0;

                        if (self.debugS) NSLog(@"*.*.*.*.      Entry: S = %g, T = %g", referenceEntropy, tNew-273.15);
                        do {
                            sumOfComputedEntropies = 0.0;
                            double sumOfComputedHeatCapacities = 0.0;
                            NSUInteger columnIndex = 0;
                            *infeasible = FALSE;

                            for (NSString *key in [phaseList allKeys]) {
                                EquilibrateStatePhase *object = [phasesAvailable valueForKey:key];
                                double entropy      = [object entropyFunctionAtT:tNew
                                                                            andP:pOld
                                                andDeltaMolesOfNonZeroComponents:&deltaMoles[columnIndex]
                                                andScalerMultiplierForCorrection:scaleLength];
                                if (entropy == compositionIsUnacceptable) { *infeasible = TRUE; return 0.0; }
                                double heatCapacity = [object heatCapacityFunctionAtT:tNew
                                                                                 andP:pOld
                                                     andDeltaMolesOfNonZeroComponents:&deltaMoles[columnIndex]
                                                     andScalerMultiplierForCorrection:scaleLength];
                                if (heatCapacity == compositionIsUnacceptable) { *infeasible = TRUE; return 0.0; }
                                if (self.debugV) NSLog(@"*.*.*.*.      Phase:%@ S = %g Cp = %g", [[object phaseClassInstance] phaseName], entropy, heatCapacity);

                                sumOfComputedEntropies      += entropy;
                                sumOfComputedHeatCapacities += heatCapacity;
                                columnIndex += [[phaseList objectForKey:key] count];
                            }
                            if (self.debugV) NSLog(@"*.*.*.*.      Iter = %lu Entropy %20.13e, Cp = %20.13e", iterations, sumOfComputedEntropies, sumOfComputedHeatCapacities);

                            tNew -= (sumOfComputedEntropies-referenceEntropy)*tNew/sumOfComputedHeatCapacities;
                            if (self.debugV) NSLog(@"*.*.*.*.      T = %g", tNew-273.15);
                            if      (tNew <  300.0) tNew =  300.0;
                            else if (tNew > 5000.0) tNew = 5000.0;

                            iterations++;
                        } while ((iterations < 100) && (fabs(sumOfComputedEntropies-referenceEntropy) > 1.0e-9*referenceEntropy));

                        [equilibrateState setT:tNew];
                        if (self.debugS) NSLog(@"*.*.*.*.      Exit: T = %g", tNew-273.15);

                    } else if (independentVariablesOfSystemPotential & independentV) {
                        double tOld = [equilibrateState T];
                        double pOld = [equilibrateState P];
                        double pNew = pOld + scaleLength*deltaP;
                        double referenceVolume = [equilibrateState referenceVolumeOfSystem] + [equilibrateState correctionToReferenceVolumeOfSystem];
                        double sumOfComputedVolumes = 0.0;
                        NSUInteger iterations = 0;

                        if (self.debugS) NSLog(@"*.*.*.*.      Entry: V = %g, P = %g", referenceVolume, pNew);
                        do {
                            sumOfComputedVolumes = 0.0;
                            double sumOfComputedDvdp = 0.0;
                            NSUInteger columnIndex = 0;
                            *infeasible = FALSE;

                            for (NSString *key in [phaseList allKeys]) {
                                EquilibrateStatePhase *object = [phasesAvailable valueForKey:key];
                                double volume = [object volumeFunctionAtT:tOld
                                                                     andP:pNew
                                         andDeltaMolesOfNonZeroComponents:&deltaMoles[columnIndex]
                                         andScalerMultiplierForCorrection:scaleLength];
                                if (volume == compositionIsUnacceptable) { *infeasible = TRUE; return 0.0; }
                                double dvdp = [object dvdpFunctionAtT:tOld
                                                                 andP:pNew
                                     andDeltaMolesOfNonZeroComponents:&deltaMoles[columnIndex]
                                     andScalerMultiplierForCorrection:scaleLength];
                                if (dvdp   == compositionIsUnacceptable) { *infeasible = TRUE; return 0.0; }
                                if (self.debugV) NSLog(@"*.*.*.*.      Phase:%@ V = %g dvdp = %g", [[object phaseClassInstance] phaseName], volume, dvdp);

                                sumOfComputedVolumes += volume;
                                sumOfComputedDvdp    += dvdp;
                                columnIndex += [[phaseList objectForKey:key] count];
                            }
                            if (self.debugV) NSLog(@"*.*.*.*.      Iter = %lu Volume %20.13e, dvdp = %20.13e", iterations, sumOfComputedVolumes, sumOfComputedDvdp);

                            pNew -= (sumOfComputedVolumes-referenceVolume)/sumOfComputedDvdp;
                            if (self.debugV) NSLog(@"*.*.*.*.      P = %g", pNew);
                            if      (pNew <      1.0) pNew =      1.0; // was -10000.0
                            else if (pNew > 100000.0) pNew = 100000.0;

                            iterations++;
                        } while ((iterations < 100) && (fabs(sumOfComputedVolumes-referenceVolume) > 1.0e-9*referenceVolume));

                        [equilibrateState setP:pNew];
                        if (self.debugS) NSLog(@"*.*.*.*.      Exit: P = %g", pNew);
                    }

					*infeasible = FALSE;
					NSUInteger columnIndex = 0;
					double sumOfPotentials = 0.0;

					for (NSString *key in [phaseList allKeys]) {
						EquilibrateStatePhase *object = [phasesAvailable valueForKey:key];
						double potential = [object potentialFunctionFor:independentVariablesOfSystemPotential
																   andT:[equilibrateState T]
																   andP:[equilibrateState P]
									   andDeltaMolesOfNonZeroComponents:&deltaMoles[columnIndex]
									   andScalerMultiplierForCorrection:scaleLength];
						if (potential == compositionIsUnacceptable) {
							*infeasible = TRUE;
							return 0.0;
						}
						sumOfPotentials += potential;
						columnIndex += [[phaseList objectForKey:key] count];
					}
                    if (self.debugS) NSLog(@"*.*.*.*. Potential %20.13e scale %20.13e", sumOfPotentials, scaleLength);
					return sumOfPotentials;
				}];
				if (self.debugS) NSLog(@"Return from linear search with status %ld, lambda = %g, step = %g, iter = %ld, potential = %g",
								  status, lambda, stepSize, iter, pTotal);
				if (status == MIN1D_BAD_INITIAL) {
					lambda /= 2.0;
					if (fabs(lambda) < self.PPPARAMETERS_LINEAR_MIN_STEPLENGTH.doubleValue) {
						status = MIN1D_SUCCESS;
						if (!self.PPOPTIONS_ZERO_LINEAR.boolValue) {
                            [self computeAndFillUsefulOutputQuantitiesForEquilibrateStateClass];
							exitStatus = @"failure, The linear search direction has *zero* length.";
							NSDictionary *dictionary = [NSDictionary dictionaryWithObjectsAndKeys:exitStatus, @"status", equilibrateState, @"results", nil];
							return dictionary;
						}
					}
				} else if (status == MIN1D_ITERS_EXCEEDED) {
					status = MIN1D_SUCCESS;
				}

			} while (status != MIN1D_SUCCESS);

			// reassemble solution
			if (self.debugV) {
                if ((independentVariablesOfSystemPotential & independentT) && (independentVariablesOfSystemPotential & independentP)) {
                    for (NSUInteger i=0; i<numberOfColumnConstraints; i++) NSLog(@"... corr soln[%2.2lu] = %13.6g, reference moles = %13.6g",
                                                                                 i, referenceMoles[i]+lambda*deltaMoles[i], referenceMoles[i]);
                } else if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentP)) {
                    for (NSUInteger i=0; i<(numberOfColumnConstraints-1); i++) NSLog(@"... corr soln[%2.2lu] = %13.6g, reference moles = %13.6g",
                                                                                 i, referenceMoles[i]+lambda*deltaMoles[i], referenceMoles[i]);
                } else if ((independentVariablesOfSystemPotential & independentT) && (independentVariablesOfSystemPotential & independentV)) {
                    for (NSUInteger i=0; i<(numberOfColumnConstraints-1); i++) NSLog(@"... corr soln[%2.2lu] = %13.6g, reference moles = %13.6g",
                                                                                     i, referenceMoles[i]+lambda*deltaMoles[i], referenceMoles[i]);
                } else if ((independentVariablesOfSystemPotential & independentS) && (independentVariablesOfSystemPotential & independentV)) {
                    for (NSUInteger i=0; i<(numberOfColumnConstraints-2); i++) NSLog(@"... corr soln[%2.2lu] = %13.6g, reference moles = %13.6g",
                                                                                     i, referenceMoles[i]+lambda*deltaMoles[i], referenceMoles[i]);
                }
            }
			NSUInteger columnIndex = 0;
            NSString *removeMeToo = @"";
			for (NSString *key in [phaseList allKeys]) {
				EquilibrateStatePhase *object = [phasesAvailable valueForKey:key];
				[object addDeltaMolesToNonZeroMolesOfEndmemberComponents:&deltaMoles[columnIndex] andScalerMultiplierForCorrection:lambda];
				columnIndex += [[phaseList objectForKey:key] count];

				// remove the phase if the mass falls below tolerance
				if (![key isEqualToString:removeMeToo] && ([object mass] < self.PPPARAMETERS_MASSOUT.doubleValue)) {

                    if ([[[object phaseClassInstance] phaseName] hasSuffix:@"PseudoPhase"]) {
                        NSString *nameOfPseudoPhaseToBeDeleted = [[object phaseClassInstance] phaseName];
                        if (self.debugS) NSLog(@"A PseudoPhase has been selected for removal by the system: %@", nameOfPseudoPhaseToBeDeleted);
                        NSMutableSet *setOfPseudoPhasesRemainingInSystem = [NSMutableSet setWithCapacity:0];
                        for (NSString *nameOfPhaseInSystem in [[equilibrateState phasesInSystem] allKeys])
                            if (![nameOfPhaseInSystem isEqualToString:nameOfPseudoPhaseToBeDeleted] && [nameOfPhaseInSystem hasSuffix:@"PseudoPhase"]) [setOfPseudoPhasesRemainingInSystem addObject:nameOfPhaseInSystem];
                        if (self.debugS) NSLog(@"... there are %lu PseudoPhases remaining in the system.", [setOfPseudoPhasesRemainingInSystem count]);
                        if ([setOfPseudoPhasesRemainingInSystem count] == 1) {
                            removeMeToo = [setOfPseudoPhasesRemainingInSystem anyObject];
                            if (self.debugS) NSLog(@"... the additional PseudoPhase %@ will be removed too after the primary.",removeMeToo);
                        }
                    }

					[[equilibrateState phasesInSystem] removeObjectForKey:key];
					if ([[object phaseClassInstance] respondsToSelector:@selector(decrementInstanceCountOfPhase)]) [[object phaseClassInstance] decrementInstanceCountOfPhase];
					if (![[[object phaseClassInstance] phaseName] hasSuffix:@"PseudoPhase"]) [[equilibrateState potentialPhaseList] setObject:object forKey:key];
					lastPhaseRemovedFromSystem = [NSString stringWithString:key];
					[setOfPhasesRemovedFromSystem addObject:[NSString stringWithString:lastPhaseRemovedFromSystem]];

					double *phaseCompositionInElements  = [[object bulkCompositionInElements] pointerToDouble];
					if (![key isEqualToString:@"Liquid"]) {
						EquilibrateStatePhase *liquidWrapper = [[equilibrateState phasesInSystem] valueForKey:@"Liquid"];
						if (liquidWrapper) { // liquid is in the assemblage
							double *liquidCompositionInElements = [[liquidWrapper bulkCompositionInElements] pointerToDouble];
							[EquilibrateState accumulateIntoVectorOfElements:liquidCompositionInElements
															   aScaler:1.0
										   timesElementReferenceVector:phaseCompositionInElements];
							[liquidWrapper setMass:[[liquidWrapper phaseClassInstance] convertElementsToTotalMass:liquidCompositionInElements]];
						} else { // liquid is not in the assemblage, so the discarded phase is subtracted from the remaining system phases
                            if (self.debugV) {
                                NSLog(@"Before removal of phase %@ the composition of all other phases in the system is:", key);
                                for (NSString *tempPhaseName in [[equilibrateState phasesInSystem] allKeys]) {
                                    EquilibrateStatePhase *tempPhaseWrapper = [[equilibrateState phasesInSystem] objectForKey:tempPhaseName];
                                    double *bc = [[tempPhaseWrapper bulkCompositionInElements] pointerToDouble];
                                    NSLog(@"... Phase: %@", tempPhaseName);
                                    for (NSUInteger i=1; i<107; i++) if (bc[i] != 0.0) NSLog(@"... %@ moles are: %20.13e",
                                                                                             [[PhaseBase elementNameFromAtomicNumber:i] stringByPaddingToLength:2 withString:@" " startingAtIndex:0], bc[i]);
                                }
                            }
							if (![self distributeMolesOfElementsAmongstSystemPhases:phaseCompositionInElements
																	   withMultiplier:1.0
																 withExcludedPhases:[NSSet setWithObject:@"none"]
                                                            withNonNegativeSolution:NO]) {
								if (!self.PPOPTIONS_DISTRIBUTION.boolValue) {
                                    [self computeAndFillUsefulOutputQuantitiesForEquilibrateStateClass];
									exitStatus = @"failure, Molar redistribution was not successful (solid).";
									NSDictionary *dictionary = [NSDictionary dictionaryWithObjectsAndKeys:exitStatus, @"status", equilibrateState, @"results", nil];
									return dictionary;
								}
							}
                            if (self.debugV) {
                                NSLog(@"After removal of phase %@ the composition of all other phases in the system is:", key);
                                for (NSString *tempPhaseName in [[equilibrateState phasesInSystem] allKeys]) {
                                    EquilibrateStatePhase *tempPhaseWrapper = [[equilibrateState phasesInSystem] objectForKey:tempPhaseName];
                                    double *bc = [[tempPhaseWrapper bulkCompositionInElements] pointerToDouble];
                                    NSLog(@"... Phase: %@", tempPhaseName);
                                    for (NSUInteger i=1; i<107; i++) if (bc[i] != 0.0) NSLog(@"... %@ moles are: %20.13e",
                                                                                             [[PhaseBase elementNameFromAtomicNumber:i] stringByPaddingToLength:2 withString:@" " startingAtIndex:0], bc[i]);
                                }
                            }
						}


					} else { // special case for discarding the liquid phase
						if (![self distributeMolesOfElementsAmongstSystemPhases:phaseCompositionInElements
																   withMultiplier:1.0
															 withExcludedPhases:[NSSet setWithObject:@"none"]
                                                        withNonNegativeSolution:NO]) {
							if (!self.PPOPTIONS_DISTRIBUTION.boolValue) {
                                [self computeAndFillUsefulOutputQuantitiesForEquilibrateStateClass];
								exitStatus = @"failure, Molar redistribution was not successful (liquid).";
								NSDictionary *dictionary = [NSDictionary dictionaryWithObjectsAndKeys:exitStatus, @"status", equilibrateState, @"results", nil];
								return dictionary;
							}
						}
					}

					[EquilibrateState scaleVectorOfElements:phaseCompositionInElements aScaler:0.0];
					[object setMass:0.0];

					converged = NO;
					if (self.debugS) NSLog(@"The phase %@ has been removed from the system and placed on the potential phase list.", key);

                    if (![removeMeToo isEqualToString:@""]) { // The additional phase is a dependent PseudoPhase that must be removed because its abundance has just gone to zero
                        EquilibrateStatePhase *additionalPhaseObject = [phasesAvailable objectForKey:removeMeToo];
                        [[equilibrateState phasesInSystem] removeObjectForKey:removeMeToo];
                        //[[equilibrateState potentialPhaseList] setObject:additionalPhaseObject forKey:key];
                        [setOfPhasesRemovedFromSystem addObject:[NSString stringWithString:removeMeToo]];
                        [EquilibrateState scaleVectorOfElements:[[additionalPhaseObject bulkCompositionInElements] pointerToDouble] aScaler:0.0];
                        [additionalPhaseObject setMass:0.0];
                        if (self.debugS) NSLog(@"The phase %@ has been removed from the system and placed on the potential phase list.", removeMeToo);
                    }
				}

			}

			[historyPotentialCalculatedInLinearSearch addObject:[NSNumber numberWithDouble:pTotal]];
			if (!self.PPOPTIONS_DETECT_MINIMAL.boolValue &&
				[historyPotentialCalculatedInLinearSearch count] >= self.PPPARAMETERS_LINEAR_MINIMAL.unsignedIntValue ) {
				NSUInteger interval = self.PPPARAMETERS_LINEAR_MINIMAL.unsignedIntValue;
				double pTotalAverage = 0.0;
				for (NSUInteger i=0; i<interval; i++) pTotalAverage += [[historyPotentialCalculatedInLinearSearch objectAtIndex:i] doubleValue];
				pTotalAverage /= (double) interval;
				[historyPotentialCalculatedInLinearSearch removeObjectAtIndex:0];
				if (fabs(pTotalAverage-averageLastPotentialCalculatedInLinearSearch) < 10.0*DBL_EPSILON) {
					converged = YES;
					exitStatus = @"success, Minimal energy computed.";
                    convergenceByMinimalEnergy = YES;
					if (self.debugS) {
						NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><");
						NSLog(@"Quadratic loop convergence with minimal energy.");
						NSLog(@"<><><><><><><><><><><><><><><><><><><><><><><><");
					}
				}
				averageLastPotentialCalculatedInLinearSearch = pTotalAverage;
			}
			if (self.debugS) {
				if (quadraticIterations == 1) NSLog(@"Start:");
				NSLog(@"Quad iteration %3.3lu potential = %20.13e, %20.1f", quadraticIterations, pTotal, (pTotal-lastPotentialCalculatedInLinearSearch)/(pTotal*DBL_EPSILON));
			}
			lastPotentialCalculatedInLinearSearch = pTotal;
		}

        if (converged) {
			quadraticIterations = 0;

			// Test if phases are named correctly
			if (converged) {
				NSMutableDictionary *result = [NSMutableDictionary dictionaryWithCapacity:1];
				for (NSString *key in [equilibrateState phasesInSystem]) {
					EquilibrateStatePhase *phaseWrapper = [[equilibrateState phasesInSystem] valueForKey:key];
					id phase = [phaseWrapper phaseClassInstance];
					if ([phase respondsToSelector:@selector(nameOfPhaseWithComposition:)]) {
						NSString *nameOfPhase = [phase nameOfPhaseWithComposition:[[phase convertElementsToMoles:[[phaseWrapper bulkCompositionInElements] pointerToDouble]] pointerToDouble]];
						if (![nameOfPhase isEqualToString:key]) [result setObject:nameOfPhase forKey:key];
					}
				}
				if ([result count] > 0) {
					for (NSString *key in result) {
						if (self.debugS) NSLog(@"Found misnamed phase. %@ should be named %@", key, [result objectForKey:key]);
						EquilibrateStatePhase *misnamedPhaseWrapper = [[equilibrateState phasesInSystem] valueForKey:key];
						EquilibrateStatePhase *properlyNamedPhaseWrapper = [[equilibrateState potentialPhaseList] valueForKey:[result objectForKey:key]];
						if (properlyNamedPhaseWrapper) {
							id misnamedPhase = [misnamedPhaseWrapper phaseClassInstance];
							id properlyNamedPhase = [properlyNamedPhaseWrapper phaseClassInstance];

							[properlyNamedPhaseWrapper setMass:[misnamedPhaseWrapper mass]];
							[misnamedPhaseWrapper setMass:0.0];

							double *misnamedPhaseComposition = [[misnamedPhaseWrapper bulkCompositionInElements] pointerToDouble];
							double *properlyNamedPhaseComposition = [[properlyNamedPhaseWrapper bulkCompositionInElements] pointerToDouble];
							for (NSUInteger i=0; i<107; i++) {
								properlyNamedPhaseComposition[i] = misnamedPhaseComposition[i];
								misnamedPhaseComposition[i] = 0.0;
							}

							if ([misnamedPhase respondsToSelector:@selector(decrementInstanceCountOfPhase)]) [misnamedPhase decrementInstanceCountOfPhase];
							[[equilibrateState phasesInSystem] removeObjectForKey:key];
							[[equilibrateState potentialPhaseList] setObject:misnamedPhaseWrapper forKey:key];

							if ([properlyNamedPhase respondsToSelector:@selector(incrementInstanceCountOfPhase)]) [properlyNamedPhase incrementInstanceCountOfPhase];
							[[equilibrateState potentialPhaseList] removeObjectForKey:[result objectForKey:key]];
							[[equilibrateState phasesInSystem] setObject:properlyNamedPhaseWrapper forKey:[result objectForKey:key]];
						} else {
							properlyNamedPhaseWrapper = [[equilibrateState phasesInSystem] valueForKey:[result objectForKey:key]];

							NSLog(@"%@ %@ %f %f", key, [result objectForKey:key], [equilibrateState T]-273.15, [equilibrateState P]);
							NSLog(@"%@ %@", [result objectForKey:key],
								  [[properlyNamedPhaseWrapper phaseClassInstance] nameOfPhaseWithComposition:[[properlyNamedPhaseWrapper bulkCompositionInElements] pointerToDouble]]);
							/*
							 // switch only if the other phase is also misnamed
							 if (![key isEqualToString:[[properlyNamedPhaseWrapper phaseClassInstance] nameOfPhaseWithComposition:[properlyNamedPhaseWrapper bulkCompositionInElements]]]) {
							 double mass = [misnamedPhaseWrapper mass];
							 [misnamedPhaseWrapper setMass:[properlyNamedPhaseWrapper mass]];
							 [properlyNamedPhaseWrapper setMass:mass];

							 double *misnamedPhaseComposition = [misnamedPhaseWrapper bulkCompositionInElements];
							 double *properlyNamedPhaseComposition = [properlyNamedPhaseWrapper bulkCompositionInElements];

							 for (NSUInteger i=0; i<107; i++) {
							 double moles = misnamedPhaseComposition[i];
							 misnamedPhaseComposition[i] = properlyNamedPhaseComposition[i];
							 properlyNamedPhaseComposition[i] = moles;
							 }
							 }
							 */
						}

						if (self.debugS) NSLog(@"... phases have been switched.");
					}
				}
			}

			if (self.PPOPTIONS_QUAD_ITERS.boolValue) {
				// If there has occurred a quadratic failure, then terminate without checking on saturation state condition
				if ([exitStatus hasPrefix:@"success"]) converged = converged && ![self evaluateSaturationState:lastPhaseRemovedFromSystem];
			} else converged = converged && ![self evaluateSaturationState:lastPhaseRemovedFromSystem];

			// Test for phase unmixing and add additional phase if required
			if (converged) {
				NSDictionary *result = nil;
				for (NSString *key in [equilibrateState phasesInSystem]) {
					EquilibrateStatePhase *phaseWrapper = [[equilibrateState phasesInSystem] valueForKey:key];
					id phase = [phaseWrapper phaseClassInstance];
					if ([phase respondsToSelector:@selector(checkForAndDetermineCompositionOfCoexistingImmisciblePhase:andT:andP:)]) {
						result = [phase checkForAndDetermineCompositionOfCoexistingImmisciblePhase:[[phase convertElementsToMoles:[[phaseWrapper bulkCompositionInElements] pointerToDouble]] pointerToDouble]
																							  andT:[equilibrateState T]
																							  andP:[equilibrateState P]];
						if ([[result objectForKey:@"additionalPhaseDetected"] boolValue]) break;
					}
				}

				if (result) {
					if ([[result objectForKey:@"additionalPhaseDetected"] boolValue] && ![[result objectForKey:@"nameOfCoexistingPhase"] isEqualToString:lastPhaseRemovedFromSystem] ) {
						NSString *nameOfCoexistingPhase = [result objectForKey:@"nameOfCoexistingPhase"];
						double coexistingPhaseAffinity = [[result objectForKey:@"affinity"] doubleValue];
						NSArray *coexistingPhaseComposition = [result objectForKey:@"composition"];

						if (self.debugS) NSLog(@"Phase %@ will be added to the assemblage (scaled affinity = %f)", nameOfCoexistingPhase, coexistingPhaseAffinity);

						EquilibrateStatePhase *phaseWrapper = [[equilibrateState potentialPhaseList] objectForKey:nameOfCoexistingPhase];
						id phase = [phaseWrapper phaseClassInstance];
						[[equilibrateState potentialPhaseList] removeObjectForKey:nameOfCoexistingPhase];

						NSUInteger na = [coexistingPhaseComposition count];
                        DoubleVector *xWrapper = [[DoubleVector alloc] initWithSize:na];
						double *x = [xWrapper pointerToDouble];
						double molesIn = self.PPPARAMETERS_MOLESIN.doubleValue;
						for (NSUInteger i=0; i<na; i++) x[i] = [[coexistingPhaseComposition objectAtIndex:i] doubleValue]*molesIn;
						[phaseWrapper setBulkCompositionInElements:[phase convertMolesToElements:x]];

						double *phaseCompositionInElements  = [[phaseWrapper bulkCompositionInElements] pointerToDouble];
						double massOfPhase = [[phaseWrapper phaseClassInstance] convertElementsToTotalMass:phaseCompositionInElements];
						[phaseWrapper setMass:massOfPhase];

						EquilibrateStatePhase *liquidWrapper;
						if ((liquidWrapper = [[equilibrateState phasesInSystem] objectForKey:@"Liquid"])) {
							double *liquidCompositionInElements = [[liquidWrapper bulkCompositionInElements] pointerToDouble];
							[EquilibrateState accumulateIntoVectorOfElements:liquidCompositionInElements
															   aScaler:-1.0
										   timesElementReferenceVector:phaseCompositionInElements];
							if (![[liquidWrapper phaseClassInstance] testPermissibleValuesOfComponents:[[[liquidWrapper phaseClassInstance] convertElementsToMoles:liquidCompositionInElements] pointerToDouble]]) {
								[EquilibrateState accumulateIntoVectorOfElements:liquidCompositionInElements
																   aScaler:1.0
											   timesElementReferenceVector:phaseCompositionInElements];
								if (![self distributeMolesOfElementsAmongstSystemPhases:phaseCompositionInElements
																		 withMultiplier:-1.0
																	 withExcludedPhases:[NSSet setWithObject:@"Liquid"]
                                                                withNonNegativeSolution:NO]) {
									if (!self.PPOPTIONS_DISTRIBUTION.boolValue) {
                                        [self computeAndFillUsefulOutputQuantitiesForEquilibrateStateClass];
                                        exitStatus = @"failure, Molar redistribution was not successful (solids, 1).";
                                        NSDictionary *dictionary = [NSDictionary dictionaryWithObjectsAndKeys:exitStatus, @"status", equilibrateState, @"results", nil];
                                        return dictionary;
                                    }
								}
							} else [liquidWrapper setMass:[liquidWrapper mass]-massOfPhase];
						} else {
							if (![self distributeMolesOfElementsAmongstSystemPhases:phaseCompositionInElements
																	 withMultiplier:-1.0
																 withExcludedPhases:[NSSet setWithObject:@"none"]
                                                            withNonNegativeSolution:NO]) {
								if (!self.PPOPTIONS_DISTRIBUTION.boolValue) {
                                    [self computeAndFillUsefulOutputQuantitiesForEquilibrateStateClass];
                                    exitStatus = @"failure, Molar redistribution was not successful (solids, 2).";
                                    NSDictionary *dictionary = [NSDictionary dictionaryWithObjectsAndKeys:exitStatus, @"status", equilibrateState, @"results", nil];
                                    return dictionary;
                                }
							}
						}

						[[equilibrateState phasesInSystem] setObject:phaseWrapper forKey:nameOfCoexistingPhase];
						if ([[phaseWrapper phaseClassInstance] respondsToSelector:@selector(incrementInstanceCountOfPhase)]) [[phaseWrapper phaseClassInstance] incrementInstanceCountOfPhase];

						converged = NO;
					}
				}
			}

			for (NSString *removedPhase in setOfPhasesRemovedFromSystem) {
				if (self.debugS) NSLog(@"Removed phase %@ %lu times", [removedPhase stringByPaddingToLength:15 withString:@" " startingAtIndex:0],
								  [setOfPhasesRemovedFromSystem countForObject:removedPhase]);
				if ([setOfPhasesRemovedFromSystem countForObject:removedPhase] >= self.PPPARAMETERS_CYCLES_MAX.unsignedIntValue) {
					converged = YES;
					exitStatus = [NSString localizedStringWithFormat:@"failure, Addition/removal cycling for phase %@ detected.", removedPhase];
				}
			}
			lastPhaseRemovedFromSystem = @"";

		}

	} // end while loop on !converged

    if (([equilibrateState correctionToReferenceVolumeOfSystem] != 0.0) || ([equilibrateState correctionToReferenceEntropyOfSystem] != 0.0))
        exitStatus = [[exitStatus stringByReplacingOccurrencesOfString:@"success" withString:@"failure"] stringByAppendingString:@" (S or V mismatch)"];

	[self computeAndFillUsefulOutputQuantitiesForEquilibrateStateClass];
	NSDictionary *dictionary = [NSDictionary dictionaryWithObjectsAndKeys:exitStatus, @"status", equilibrateState, @"results", nil];
	return dictionary;
}

@end
