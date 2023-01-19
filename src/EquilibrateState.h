//
//  EquilibrateState.h
//  derived from:
//  MeltsState.h
//  originally part of the PhasePlot package
//
//  Created by Mark Ghiorso on 6/20/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>

@class DoubleVector;

static const NSUInteger independentT = 001LU;
static const NSUInteger independentS = 002LU;
static const NSUInteger independentP = 004LU;
static const NSUInteger independentV = 010LU;
static const NSUInteger independentH = 020LU;

static const double compositionIsUnacceptable = 123456789.0;

@interface EquilibrateStatePhase : NSObject {
	id     phaseClassInstance;
	double mass;
    DoubleVector *bulkCompositionInElements;
	double affinity;
}

@property (readwrite, strong) id phaseClassInstance;
@property (readwrite, assign) double mass;
@property (readwrite, assign) DoubleVector *bulkCompositionInElements;
@property (readwrite, assign) double affinity;

/**
 Returns the molar sum of the bulkCompositionInElements array
 */
-(double)sumOfMolesOfElements;

/**
 Determines the indices of non-zero entries for moles of components in the phase.

 @return
 An array of indices of non-zero components.
 */
-(NSArray *)indicesOfNonZeroMolesOfEndmemberComponents;

/**
 Determines non-zero entries for moles of components and composition/potential energy arrays for the phase.

 Proceedure:
 - First determine if the phase is a solution or a stoichiometric phase.
 - For each component of the solution or for the stoichiometric phase report the:
   - Total moles.
   - A pointer to a vector that represents the stoichiometry in terms of elements.
   - Chemical potential.
   - Hessian (row/column), if the phase is a solution, otherwise the entry is missing:
     - An array that holds the derivatives of the chemical potential with repect to all solution components whose molar abundances are non-zero.
 - Repeat the above as required for all the endmembers of the solution (see return section for format).

 @param constraint
   Specifies constraints for calculation of potential energy functions:
   - independentT & independentP =>> derivatives of the Gibbs free energy
   - independentT & independentV =>> derivatives of the Helmholtz energy
   - independentS & independentP =>> derivatives of the enthalpy
   - independentS & independentV =>> derinatives of the internal energy
   - independentH & independentP =>> derivatives of the negative entropy

 @param t
   Temperature in K.

 @param p
   Pressure in bars.

 @return
   An NSArray object with the following structure:
   - NSArray object: one entry for each non-zero component or one entry for the stoichiometric phase:
     - NSNumber object: Contains a double value.  This the number of moles of the component/phase.
     - NSValue object: Contains a pointer to double value. This is an array of length 107 containing stoiciometric abundances of the elements in one mole of the component/phase.
     - NSNumber object: Contains a double value.  This the first derivative of the potential function with respect to the component/phase.
     - (only for components) NSArray object: Contains NSNumber objects that hold a double value. These values are the elements of the Hessian (second derivative matrix of the potential function) row for the non-zero components in the solution.

     - (additional entries are included to handle constraint vectors, see comments in code).

 The number of objects in the array is one for a stoichiometric phase and equal to the number of components with non-zero molar abundances for a solution.
 */
-(NSArray *)nonZeroMolesOfEndmemberComponents:(NSUInteger)constraint andT:(double)t andP:(double)p;
/**
 Computes the potential function for the phase.

 @param constraint
 Specifies constraints for calculation of potential energy functions:
 - independentT & independentP =>> derivatives of the Gibbs free energy
 - independentT & independentV =>> derivatives of the Helmholtz energy
 - independentS & independentP =>> derivatives of the enthalpy
 - independentS & independentV =>> derinatives of the internal energy
 - independentH & independentP =>> derivatives of the negative entropy

 @param t
 Temperature in K.

 @param p
 Pressure in bars.

 @param deltaMoles
 Corrections to non-zero molar abundances of components (or the moles of the phase). Usually determined by a quadratic minimization solution returned from HFTI.

 @param scaleFractor
 A proposed scaler multiplier for computation of the molar abundances of each component or the phase, as: mNew = mOld + scaleFactor * deltaMoles.

 @return
 The value of the potential or the number compositionIsUnacceptable if the computed composition is infeasible.
 */
-(double)potentialFunctionFor:(NSUInteger)constraint andT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor;
-(double)entropyFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor;
-(double)volumeFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor;
-(double)enthalpyFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor;
-(double)heatCapacityFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor;
-(double)dvdpFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor;
-(double)dvdtFunctionAtT:(double)t andP:(double)p andDeltaMolesOfNonZeroComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor;

/**
 Resets molar abundances of phases and components in phases by updating the bulkCompositionInElements array

 @param deltaMoles
 Corrections to non-zero molar abundances of components (or the moles of the phase). Usually determined by a quadratic minimization solution returned from HFTI.

 @param scaleFractor
 A proposed scaler multiplier for computation of the molar abundances of each component or the phase, as: mNew = mOld + scaleFactor * deltaMoles.
 */
-(void)addDeltaMolesToNonZeroMolesOfEndmemberComponents:(double *)deltaMoles andScalerMultiplierForCorrection:(double)scaleFractor;

@end

static const NSUInteger FO2_NONE = 001LU;
static const NSUInteger FO2_HM   = 002LU;
static const NSUInteger FO2_NNO  = 004LU;
static const NSUInteger FO2_QFM  = 010LU;
static const NSUInteger FO2_COH  = 020LU;
static const NSUInteger FO2_IW   = 040LU;

@interface EquilibrateState : NSObject {
	double mass;                                    // Mass of the system
	DoubleVector *bulkCompositionInOxides;          // current bulk composition (moles of oxides)
    NSUInteger nOxides;                             // number of oxides in system
    DoubleVector *bulkCompositionInElements;        // current bulk composition (moles of elements)

	NSMutableDictionary * phasesInSystem;           // A collection of EquilibrateStatePhase objects
	NSMutableDictionary * potentialPhaseList;       // A collection of EquilibrateStatePhase objects

	double T;                                       // current temperature (K)
	double P;                                       // current pressure (bars)
	double referenceEnthalpyOfSystem;               // Total reference enthalpy of the system (J)
	double referenceEntropyOfSystem;                // Total reference entropy of the system (J/K)
	double referenceVolumeOfSystem;                 // Total reference volume of the system (J/bar)

    double correctionToReferenceEntropyOfSystem;    // feasibility correction to total reference entropy of the system (J/K)
    double correctionToReferenceVolumeOfSystem;     // feasibility correction to total reference volume of the system (J/bar)

	double     fo2;                                 // current value of fo2 (numeric, base 10 log)
	NSUInteger fo2Path;                             // current value of fo2 path (i.e. FO2_NONE, etc
	double     fo2Delta;                            // offset from fo2Path
	double     referenceOxygen;                     // reference value of O2 content in the system

	BOOL   isIsenthalpic;                           // isenthalpic mode (YES/NO)
	BOOL   isIsentropic;                            // isentropic mode (YES/NO)
	BOOL   isIsochoric;                             // isochoric mode (YES/NO)
}

@property (readwrite, assign) double mass;
@property (readwrite, assign) DoubleVector *bulkCompositionInOxides;
@property (readwrite, assign) DoubleVector *bulkCompositionInElements;
@property (readwrite, strong) NSMutableDictionary *phasesInSystem;
@property (readwrite, strong) NSMutableDictionary *potentialPhaseList;
@property (readwrite, assign) double T;
@property (readwrite, assign) double P;
@property (readwrite, assign) double referenceEntropyOfSystem;
@property (readwrite, assign) double referenceVolumeOfSystem;
@property (readwrite, assign) double referenceEnthalpyOfSystem;
@property (readwrite, assign) double correctionToReferenceEntropyOfSystem;
@property (readwrite, assign) double correctionToReferenceVolumeOfSystem;
@property (readwrite, assign) double fo2;
@property (readwrite, assign) NSUInteger fo2Path;
@property (readwrite, assign) double fo2Delta;
@property (readwrite, assign) double referenceOxygen;

/**
 elements[] += scaler * reference[].

 Takes an input vector of lenth 107 containig mole numbers of elements indexed on atomic number and
 adds to that vector a scaler multple of another vector of length 107 which represents a phase stoichiometry.

 @param elements
   double vector of length 107. Input and output.

 @param scaler
   input only

 @param reference
   double vector of length 107. Input only.
 */
+(void)accumulateIntoVectorOfElements:(double *)elements aScaler:(double)scaler timesElementReferenceVector:(double *)reference;

/**
 elements[] *= scaler.

 Takes an input vector of lenth 107 containig mole numbers of elements indexed on atomic number and
 multiplies each element of that vector by a scaler.

 @param elements
 input/output double vector of length 107. Input and output.

 @param scaler
 input only
 */
+(void)scaleVectorOfElements:(double *)elements aScaler:(double)scaler;

/**
 Initialize with the number of oxides (system components)

 @param na
   Number of components in teh system (in the usual case the number of oxides)
 */
-(id)initWithComponentNumber:(NSUInteger)na;

/**
 Determines the number of non-zero entries in the element array instance variable bulkCompositionInElements

 @return
   number of non-zero entries
 */
-(NSUInteger)numberOfNonZeroElements;

/**
 Generates and return a table of element indices for non-zero entries in the element array instance variable bulkCompositionInElements

 @return
   An NSArray of NSNumber objects each containing a NSUInteger index
 */
-(NSArray *)hashTableOfNonZeroElementEntries;

/**
 Formats an instance of the EquilibrateState class as an XML document

 @return
 Pointer to an NSXMLDocument with root element EquilibrateState hat contains the current instance
 */
-(NSXMLElement *)getEquilibrateStateAsXMLElement:(id)delegate selector:(SEL)convertMolesOfElementsToGramsOfOxides;

@end
