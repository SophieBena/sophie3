//
//  StixrudeSolutionPhase.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 3/10/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "SolutionPhaseProtocol.h"
#import "PhaseBase.h"

#define NA 6
#define NR 5
#define NS 16

@interface StixrudeSolutionPhase : PhaseBase <SolutionPhaseProtocol> {
@protected
	BOOL computeMixingQuantities;
	NSArray *endmembers;
	double dArray[NA];
	double adjustmentCoefficient[NA];
	double wArray[NA][NA];
	double siteMultiplicity[NS];
	double siteSpeciesCoefficients[NS][NA];
	NSUInteger na, nr, ns;
	double nAtoms;
	void (^convertMolesOfElementsToMolesOfSpecies)(double *, double *);
	NSString *(^convertMolesOfSpeciesToDisplayFormula)(double *);
    void (^correctActivityCoefficients)(double *, double *);

    NSUInteger numberOfIterationsAllowedInSaturationMethod;
    double minimumAffinityDiffereneceToAllowConvergenceInSaturationMethod;
}

@property NSUInteger numberOfIterationsAllowedInSaturationMethod;
@property double minimumAffinityDiffereneceToAllowConvergenceInSaturationMethod;

/**
 Initializes an instance of a StixrudeSolutionPhase class.

 @param name
 Name of phase.

 @param speciesProperties
 Length na (<=NA) of array containing StixrudeProperties objects.

 @param speciesWeights
 Length na (<=NA) of array containing NSNumber objects. Each object holds a double that represents the van Laar weighting coefficient for the species.

 @param speciesXfactor
 Length na (<=NA) of array containing NSNumber objects. Each object holds a double that represents the van Laar weighting coefficient for the species.

 @param wIn
 Length na*(na-1)/2 [<=NA*(NA-1)] of array containing NSNumber objects.  Each object holds a double that represents site redundancy in the saturation/affinity method.

 @param sites
 Length ns (<=NS) of array of objects.  Each object is an NSArray that contains:
 - at index 0: NSNumber object holding a double which represents the site multiplicity.
 - at index 1: An NSArray of length na of NSNumber objects. Each holds a double that represents the species coefficient for the site.

 @param nAtomIn
 Number of atoms in the solution formula unit.

 @param convertor
 Code block that implements conversion of moles of elements vector (1st argument) to moles of species vector (2nd argument).

 @param display
 Code block that implements formula display of solution.  Input vector of moles of species.  Return NSString with formula.

 @param adjuster
 Code block called from affinityAndCompositionFromLiquidChemicalPotentialSum to adjust activity coefficients in teh course of estimatng saturation compositions.  Used to avoid regions of composition space.

 @return
 Pointer to the instantiated object.
 */
-(id)initWithPhaseName:(NSString *)name
		   withSpecies:(NSArray *)speciesProperties
	withSpeciesWeights:(NSArray *)speciesWeights
	withSpeciesXfactor:(NSArray *)speciesXfactor
			withWarray:(NSArray *)wIn
			 withSites:(NSArray *)sites
			WithNAtoms:(NSUInteger)nAtomIn
  withElementConvertor:(void (^)(double *, double *))convertor
 andWithFormulaDisplay:(NSString *(^)(double *))display
andWithSaturationStateAdjuster:(void (^)(double *, double *))adjuster;

#undef NA
#undef NR
#undef NS

@end
