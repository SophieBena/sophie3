//
//  PhaseBase.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/14/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>

@class DoubleVector;

@interface PhaseBase : NSObject <NSSecureCoding> {
	NSString *phaseFormula;
	NSString *phaseName;
    DoubleVector *formulaAsElementArray;
	double mw;
    double entropyFromRobieEtAl1979;
}

/**
 Chemical formula of the phase.  A custom setter routine is defined to assign the property.
 This routine automatically initializes the formula conversion array and molecular weight.
 */
@property (nonatomic, readwrite, copy) NSString *phaseFormula;
/**
 Name of the phase.
 */
@property (readwrite, copy) NSString *phaseName;
/**
 Molecular weight of the phase.  Computed when formula is set.  Reeadonly.
 */
@property (readwrite, assign) double mw;
/**
 Number of moles of each elememnt calculated from the phase formula. Computed when formula is set.  Reeadonly.
 */
@property (readwrite) DoubleVector *formulaAsElementArray;
/**
 Entropy of the elements calculated from phase stoichiometry using values tabulated by
 Robie, Hemingway and Fisher (1979) USGS Bull 1452
 */
@property (readwrite, assign) double entropyFromRobieEtAl1979;

/**
 Name of element from element index
 */
+(NSString *)elementNameFromAtomicNumber:(NSUInteger)atomicNumber;

/**
 Moles of elements (standard order) => Moles of phase
 */
-(double)convertElementsToMolesOfPhase:(double *)e;
/**
 Moles of elements (standard order) => Total mass of the phase (g)
 */
-(double)convertElementsToMassOfPhase:(double *)e;

@end
